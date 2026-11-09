#!/usr/bin/env python3
"""Independent leakage-controlled feature-selection audit.

IMPORTANT: This is a NEW nested-CV experiment, not a reconstruction of the
historical VIF and forward-wrapper selection in the manuscript. Only outer-fold
holdout metrics are reported. Independent published test sets are never loaded.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_selection import SelectKBest, f_classif
from sklearn.metrics import accuracy_score, matthews_corrcoef, roc_auc_score, make_scorer
from sklearn.model_selection import GridSearchCV, StratifiedKFold
from sklearn.pipeline import Pipeline


def evaluate(data: Path, folds: int, trees: int, k_values: list[int]):
    matrix = pd.read_csv(data / 'scavenpred_train_88.csv')
    labels = pd.read_csv(data / 'scavenpred_train_labels.csv')
    lookup = dict(zip(labels.seq.astype(str), labels['class'].astype(int)))
    assert matrix.seq.is_unique and labels.seq.is_unique
    assert matrix.seq.map(lookup).notna().all(), 'Training labels missing'
    assert (matrix['class'].astype(int) == matrix.seq.map(lookup).astype(int)).all()
    features = [c for c in matrix.columns if c not in ('seq','class')]
    X = matrix[features].apply(pd.to_numeric, errors='raise')
    y = matrix['class'].astype(int)
    assert X.notna().all().all(), 'Missing descriptor data'
    outer = StratifiedKFold(n_splits=folds,shuffle=True,random_state=42)
    reports=[]
    for fold,(i_train,i_valid) in enumerate(outer.split(X,y),1):
        # The selector is inside the Pipeline: each GridSearchCV inner fit learns
        # its selection only from that inner TRAIN split.
        model=Pipeline([
            ('selector', SelectKBest(score_func=f_classif)),
            ('rf',RandomForestClassifier(n_estimators=trees,random_state=42,n_jobs=-1)),
        ])
        inner=StratifiedKFold(n_splits=3,shuffle=True,random_state=42+fold)
        grid=GridSearchCV(model,{'selector__k':k_values,'rf__max_depth':[7,10],
                                 'rf__min_samples_leaf':[1,4]},
                          scoring=make_scorer(matthews_corrcoef),cv=inner,
                          n_jobs=1,refit=True)
        grid.fit(X.iloc[i_train],y.iloc[i_train])
        pred=grid.predict(X.iloc[i_valid])
        proba=grid.predict_proba(X.iloc[i_valid])[:,1]
        selected=[x for x,chosen in zip(features,grid.best_estimator_.named_steps['selector'].get_support()) if chosen]
        reports.append({'outer_fold':fold,'n_train':len(i_train),'n_validation':len(i_valid),
                        'accuracy':float(accuracy_score(y.iloc[i_valid],pred)),
                        'mcc':float(matthews_corrcoef(y.iloc[i_valid],pred)),
                        'roc_auc':float(roc_auc_score(y.iloc[i_valid],proba)),
                        'selected_features':selected,'best_inner_params':grid.best_params_})
        print('Fold',fold,'MCC',round(reports[-1]['mcc'],4),'selected',selected)
    return {'method':'Nested, train-fold-only SelectKBest + RandomForest; NOT original historical wrapper',
            'independent_test_used':False,'outer_folds':reports,
            'mean_outer_mcc':float(np.mean([r['mcc'] for r in reports])),
            'mean_outer_accuracy':float(np.mean([r['accuracy'] for r in reports]))}


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--data',type=Path,default=Path(__file__).resolve().parents[1]/'data'/'restored')
    parser.add_argument('--folds',type=int,default=5)
    parser.add_argument('--trees',type=int,default=200)
    parser.add_argument('--k-values',type=int,nargs='+',default=[5,10,14])
    parser.add_argument('--output',type=Path,default=Path('nested_cv_results.json'))
    args=parser.parse_args()
    assert all(k>0 and k<=88 for k in args.k_values)
    out=evaluate(args.data,args.folds,args.trees,args.k_values)
    args.output.write_text(json.dumps(out,indent=2))
    print('Mean MCC',out['mean_outer_mcc'])
