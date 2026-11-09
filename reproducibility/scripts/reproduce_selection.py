#!/usr/bin/env python3
"""Executable reproduction of the corrected ScavenPred training-only feature pipeline.

Welch unpaired t-test (p<0.05), class-signal-ordered Pearson redundancy
filter (|r|>=0.7), iterative VIF (<10) with class-correlation tie-break,
and forward Random Forest 5-fold CV wrapper optimizing MCC.

Independent test records are NEVER loaded during feature selection.
Only finalized, disjoint, training-peptide feature data are used.
The forward-wrapper CV results are selection scores, not unbiased external
performance. Validate the final frozen model on the external test separately.
"""
from __future__ import annotations
import argparse
import json
import warnings
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import ttest_ind
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import matthews_corrcoef, make_scorer
from sklearn.model_selection import StratifiedKFold, cross_val_score
from statsmodels.stats.outliers_influence import variance_inflation_factor


def reconstruct(data: Path, max_features: int=12, trees: int=100, folds: int=5):
    # No test-data read is permitted in this function.
    train = pd.read_csv(data/'scavenpred_train_88.csv')
    X = train.drop(columns=['seq','class']).apply(pd.to_numeric,errors='raise')
    y = train['class'].astype(int)
    assert len(train)==1648 and train['seq'].is_unique
    assert X.shape[1]==88 and X.notna().all().all()

    statistics=[]
    for feature in X.columns:
        pos = X.loc[y==1,feature].values
        neg = X.loc[y==0,feature].values
        _,p = ttest_ind(pos,neg,equal_var=False)
        corr = float(np.corrcoef(X[feature],y)[0,1])
        statistics.append({'feature':feature,'welch_p':float(p),'correlation_with_class':corr})
    significant = sorted([s for s in statistics if s['welch_p']<0.05],
                         key=lambda s:-abs(s['correlation_with_class']))

    # Retain more class-discriminating feature when a pair is redundant.
    retained=[]
    for row in significant:
        f=row['feature']
        if all(abs(X[f].corr(X[kept]))<0.70 for kept in retained):
            retained.append(f)

    corr_map={r['feature']:abs(r['correlation_with_class']) for r in statistics}
    candidates=retained.copy()
    removals=[]
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore',category=RuntimeWarning)
        while len(candidates)>1:
            matrix=np.column_stack([np.ones(len(X)),X[candidates].to_numpy(dtype=float)])
            vif=np.asarray([variance_inflation_factor(matrix,i+1)
                            for i in range(len(candidates))])
            largest=float(np.max(vif))
            if largest<10:
                break
            # The historical procedure removes the weakest class-correlated
            # feature among effectively tied high-VIF features.
            tied=np.where(np.isclose(vif,largest,rtol=1e-6,atol=1e-8))[0]
            idx=min(tied,key=lambda i:corr_map[candidates[i]])
            removals.append({'feature':candidates[idx],'vif':float(vif[idx])})
            candidates.pop(idx)

    splitter=StratifiedKFold(n_splits=folds,shuffle=True,random_state=42)
    scorer=make_scorer(matthews_corrcoef)
    selected=[]
    wrapper=[]
    for step in range(max_features):
        trials=[]
        for f in candidates:
            if f in selected:
                continue
            fields=selected+[f]
            rf=RandomForestClassifier(n_estimators=trees,max_depth=7,random_state=42,n_jobs=1)
            scores=cross_val_score(rf,X[fields],y,cv=splitter,scoring=scorer,n_jobs=1)
            trials.append({'added':f,'subset':fields,'mean_mcc':float(scores.mean()),
                           'sd_mcc':float(scores.std())})
        if not trials:break
        best=max(trials,key=lambda r:r['mean_mcc'])
        improvement=(best['mean_mcc']-wrapper[-1]['mean_mcc']) if wrapper else float('inf')
        if wrapper and improvement<=0:
            wrapper.append({'step':step+1,'decision':'STOP',**best})
            break
        selected.append(best['added'])
        wrapper.append({'step':step+1,'decision':'ADD',**best})
        print(f"Step {step+1}: add {best['added']}; mean CV MCC={best['mean_mcc']:.6f}",flush=True)
    return {'method':'Training-only: Welch p<0.05, correlation |r|<0.70, VIF<10, RF100 forward MCC',
            'independent_test_used':False,
            'training_peptides':len(train),'initial_descriptors':len(X.columns),
            'significant_descriptors':len(significant),'correlation_filtered':len(retained),
            'vif_filtered':len(candidates),'vif_candidates':candidates,
            'vif_removal_log':removals,'wrapper':wrapper,'final_features':selected}


if __name__ == '__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--data',type=Path,default=Path(__file__).resolve().parents[1]/'data'/'restored')
    p.add_argument('--out',type=Path,default=Path('reproduced_feature_selection.json'))
    p.add_argument('--max-features',type=int,default=12)
    p.add_argument('--trees',type=int,default=100)
    p.add_argument('--folds',type=int,default=5)
    a=p.parse_args()
    result=reconstruct(a.data,max_features=a.max_features,trees=a.trees,folds=a.folds)
    a.out.write_text(json.dumps(result,indent=2,allow_nan=True))
    print('Summary',json.dumps({k:v for k,v in result.items() if k in ['initial_descriptors','significant_descriptors','correlation_filtered','vif_filtered','final_features']},indent=2))
