#!/usr/bin/env python3
"""Strict nested CV: supervised feature selection and RF tuning inside INNER folds.
Independent-test files are not read. All decisions fitted anew in training folds.
"""
from pathlib import Path
import argparse, json, hashlib, platform, sys
import numpy as np
import pandas as pd
from sklearn.pipeline import Pipeline
from sklearn.feature_selection import SelectKBest, f_classif
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedKFold, GridSearchCV
from sklearn.metrics import accuracy_score, matthews_corrcoef, roc_auc_score, average_precision_score, confusion_matrix,make_scorer
import sklearn, scipy

CONFIG = {'scavenpred': ('scavenpred_train_88.csv',None), 'ace-ipred':('ace_train_88.csv','ace_train_labels.csv')}

def load(root,name):
    feature_file,label_file=CONFIG[name]
    df=pd.read_csv(root/feature_file)
    seq=df.seq.astype(str).str.upper().str.strip()
    assert seq.is_unique
    if label_file:
        labels=pd.read_csv(root/label_file)
        labels.seq=labels.seq.astype(str).str.upper().str.strip()
        assert labels.seq.is_unique
        y=seq.map(labels.set_index('seq')['class'])
        assert y.notna().all() and len(labels)==len(df)
    else: y=df['class']
    X=df.drop(columns=['Unnamed: 0','seq','class','index'],errors='ignore').astype(float)
    assert X.shape[1]==88 and X.notna().all().all() and set(y.unique())=={0,1}
    return seq,X,y.astype(int)

def run(root,name,out,trees=100,seed=42):
    seq,X,y=load(root,name)
    outer=StratifiedKFold(n_splits=5,shuffle=True,random_state=seed)
    oof=[]; folds=[]
    for i,(a,b) in enumerate(outer.split(X,y),1):
        inner=StratifiedKFold(n_splits=3,shuffle=True,random_state=seed+i)
        estimator=Pipeline([('select',SelectKBest(score_func=f_classif)),('rf',RandomForestClassifier(n_estimators=trees,random_state=seed,n_jobs=1))])
        grid=GridSearchCV(estimator,{'select__k':[5,10], 'rf__max_depth':[7,10],'rf__min_samples_leaf':[1,4]},cv=inner,scoring=make_scorer(matthews_corrcoef),n_jobs=1,refit=True,error_score='raise')
        grid.fit(X.iloc[a],y.iloc[a]);pred=grid.predict(X.iloc[b]);score=grid.predict_proba(X.iloc[b])[:,list(grid.classes_).index(1)]
        choices=list(X.columns[grid.best_estimator_.named_steps['select'].get_support()])
        fold={'fold':i,'train_n':len(a),'validation_n':len(b),'outer_train_sha256':hashlib.sha256('\n'.join(seq.iloc[a]).encode()).hexdigest(),'outer_validation_sha256':hashlib.sha256('\n'.join(seq.iloc[b]).encode()).hexdigest(),'best_params':grid.best_params_,'chosen_features':choices,'accuracy':float(accuracy_score(y.iloc[b],pred)),'mcc':float(matthews_corrcoef(y.iloc[b],pred)),'roc_auc':float(roc_auc_score(y.iloc[b],score)),'average_precision':float(average_precision_score(y.iloc[b],score)),'confusion_matrix_0_1':confusion_matrix(y.iloc[b],pred,labels=[0,1]).tolist()}
        folds.append(fold)
        for seqid,truth,prediction,prob in zip(seq.iloc[b],y.iloc[b],pred,score): oof.append({'sequence':seqid,'fold':i,'truth':int(truth),'prediction':int(prediction),'positive_probability':float(prob)})
        print(name,'fold',i,'mcc',round(fold['mcc'],3),'features',','.join(choices),flush=True)
    table=pd.DataFrame(oof)
    assert len(table)==len(X) and table.sequence.is_unique and sorted(table.sequence)==sorted(seq)
    summary={'model':name,'design':'strict 5x3 nested CV with select-k-best and RF hyperparameter tuning in each inner fold','n':len(table),'outer_accuracy':float(accuracy_score(table.truth,table.prediction)),'outer_mcc':float(matthews_corrcoef(table.truth,table.prediction)),'outer_roc_auc':float(roc_auc_score(table.truth,table.positive_probability)),'outer_pr_auc':float(average_precision_score(table.truth,table.positive_probability)),'confusion_matrix_0_1':confusion_matrix(table.truth,table.prediction,labels=[0,1]).tolist(),'independent_test_loaded':False,'train_only_input_files':[CONFIG[name][0]]+([CONFIG[name][1]] if CONFIG[name][1] else []),'settings':{'outer_folds':5,'inner_folds':3,'selection_k':[5,10],'max_depth':[7,10],'min_samples_leaf':[1,4],'rf_trees':trees,'seed':seed},'versions':{'python':platform.python_version(),'sklearn':sklearn.__version__,'numpy':np.__version__,'pandas':pd.__version__,'scipy':scipy.__version__},'folds':folds}
    out.mkdir(parents=True,exist_ok=True)
    table.to_csv(out/(name+'_strict_nested_oof.csv'),index=False)
    (out/(name+'_strict_nested_report.json')).write_text(json.dumps(summary,indent=2))
    print('DONE',name,{k:v for k,v in summary.items() if k in ['n','outer_accuracy','outer_mcc','outer_roc_auc','outer_pr_auc']},flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--data',type=Path,required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--model',choices=['scavenpred','ace-ipred','both'],default='both');p.add_argument('--trees',type=int,default=100)
    a=p.parse_args()
    for name in CONFIG if a.model=='both' else [a.model]:run(a.data,name,a.out,trees=a.trees)
