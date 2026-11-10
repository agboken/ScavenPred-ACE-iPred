#!/usr/bin/env python3
"""Grouped 5x3 nested CV retaining Welch, Pearson, VIF, RF forward wrapper.

Every supervised prefilter and wrapper operates only on the corresponding
outer training partition. The 3-fold grouped wrapper scores select features,
but must not be presented as unbiased generalisation estimates. Only grouped,
held-out outer predictions contribute to the final reported CV metrics.
External test data are never loaded.
"""
from pathlib import Path
import argparse,json,time,warnings,platform
import numpy as np,pandas as pd
from scipy.stats import ttest_ind
from statsmodels.stats.outliers_influence import variance_inflation_factor
from sklearn.model_selection import StratifiedGroupKFold
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score,matthews_corrcoef,roc_auc_score,average_precision_score,confusion_matrix
from strict_group_nested_validation import dataset,make_groups
import sklearn

def screen(X,y,candidate_limit=12):
    y=np.asarray(y); statistics=[]
    for f in X.columns:
        p=float(ttest_ind(X.loc[y==1,f].to_numpy(), X.loc[y==0,f].to_numpy(), equal_var=False).pvalue)
        corr=float(np.corrcoef(X[f].to_numpy(),y)[0,1])
        if np.isfinite(p) and p<.05 and np.isfinite(corr):
            statistics.append((f,abs(corr),p))
    statistics.sort(key=lambda row:(-row[1],row[0]))
    retained=[]
    for f,_,_ in statistics:
        if all(abs(float(X[f].corr(X[g])))<.70 for g in retained):retained.append(f)
    candidates=retained.copy()
    signal={f:c for f,c,_ in statistics}
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore',category=RuntimeWarning)
        while len(candidates)>1:
            matrix=np.column_stack([np.ones(len(X)),X[candidates].to_numpy(float)])
            vifs=np.asarray([variance_inflation_factor(matrix,i+1) for i in range(len(candidates))])
            largest=np.nanmax(vifs)
            if largest<10:break
            tied=np.where(np.isclose(vifs,largest,rtol=1e-6,atol=1e-8))[0]
            drop=min(tied,key=lambda i:(signal[candidates[i]],candidates[i]))
            candidates.pop(drop)
    return candidates[:candidate_limit],{'welch':len(statistics),'pearson':len(retained),'vif':len(candidates)}

def rf_parameters(name,trees):
    return dict(n_estimators=trees,max_depth=7 if name=='scavenpred' else 10,
                min_samples_leaf=1 if name=='scavenpred' else 4,
                max_features='sqrt' if name=='scavenpred' else 1.0,
                class_weight=None if name=='scavenpred' else 'balanced',
                random_state=42,n_jobs=1)

def evaluate(data,out,name,max_features,trees,candidate_limit):
    seq,X,y=dataset(data,name)
    y=y.to_numpy();groups=make_groups(seq)
    outer=StratifiedGroupKFold(n_splits=5,shuffle=True,random_state=42)
    oof=[];details=[]
    for fold,(tr,va) in enumerate(outer.split(X,y,groups),1):
        assert set(groups[tr]).isdisjoint(set(groups[va]))
        Xt=X.iloc[tr].reset_index(drop=True); yt=y[tr]; gt=groups[tr]
        candidates,counts=screen(Xt,yt,candidate_limit)
        if not candidates:raise RuntimeError('No candidates in outer fold '+str(fold))
        inner=StratifiedGroupKFold(n_splits=3,shuffle=True,random_state=42+fold)
        inner_splits=list(inner.split(Xt,yt,gt))
        for it,iv in inner_splits:assert set(gt[it]).isdisjoint(set(gt[iv]))
        selected=[];steps=[];cache={};previous=-np.inf
        def score(features):
            key=tuple(features)
            if key not in cache:
                scores=[]
                for it,iv in inner_splits:
                    model=RandomForestClassifier(**rf_parameters(name,trees))
                    model.fit(Xt.iloc[it][features],yt[it])
                    pred=model.predict(Xt.iloc[iv][features])
                    scores.append(matthews_corrcoef(yt[iv],pred))
                cache[key]=float(np.mean(scores))
            return cache[key]
        for step in range(min(max_features,len(candidates))):
            options=[(score(selected+[f]),f) for f in candidates if f not in selected]
            best,f=max(options,key=lambda z:(z[0],z[1]))
            if selected and best<=previous:break
            selected.append(f);previous=best
            steps.append({'feature_added':f,'inner_group_cv_mcc_selection_score':best})
        model=RandomForestClassifier(**rf_parameters(name,1000 if name=='scavenpred' else 200))
        model.fit(Xt[selected],yt)
        p=model.predict(X.iloc[va][selected])
        proba=model.predict_proba(X.iloc[va][selected])[:,list(model.classes_).index(1)]
        for idx,pred,prob in zip(va,p,proba):
            oof.append({'sequence':seq.iloc[idx],'group':int(groups[idx]),'fold':fold,
                        'actual':int(y[idx]),'prediction':int(pred),'positive_probability':float(prob)})
        details.append({'fold':fold,'training_n':len(tr),'validation_n':len(va),
                        'validation_group_disjoint':True,'screen_counts':counts,
                        'candidate_features':candidates,'selected_features':selected,
                        'selection_steps':steps,
                        'outer_accuracy':float(accuracy_score(y[va],p)),
                        'outer_mcc':float(matthews_corrcoef(y[va],p)),
                        'inner_models_fitted':len(cache)*3,
                        'outer_confusion_matrix':confusion_matrix(y[va],p,labels=[0,1]).tolist()})
        print(name,'fold',fold,'MCC',details[-1]['outer_mcc'],flush=True)
    df=pd.DataFrame(oof)
    assert len(df)==len(X) and df.sequence.is_unique and df.fold.nunique()==5
    report={'model':name,'training_n':len(X),'group_n':len(np.unique(groups)),
            'independent_test_read':False,'candidate_limit':candidate_limit,
            'max_features':max_features,'inner_trees':trees,
            'oof_accuracy':float(accuracy_score(df.actual,df.prediction)),
            'oof_mcc':float(matthews_corrcoef(df.actual,df.prediction)),
            'oof_roc_auc':float(roc_auc_score(df.actual,df.positive_probability)),
            'oof_pr_auc':float(average_precision_score(df.actual,df.positive_probability)),
            'confusion_matrix':confusion_matrix(df.actual,df.prediction,labels=[0,1]).tolist(),
            'folds':details,
            'versions':{'python':platform.python_version(),'sklearn':sklearn.__version__,
                        'numpy':np.__version__,'pandas':pd.__version__},
            'note':'Welch/Pearson/VIF fitted to OUTER training partition. Inner wrapper scores select features and are not independent generalization estimates. Final OOF scores use untouched OUTER validation groups.'}
    out.mkdir(parents=True,exist_ok=True)
    df.to_csv(out/(name+'_wrapper_group_oof.csv'),index=False)
    (out/(name+'_wrapper_group_report.json')).write_text(json.dumps(report,indent=2))
    print(json.dumps({k:report[k] for k in ['model','oof_accuracy','oof_mcc','oof_roc_auc','oof_pr_auc','confusion_matrix']},indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--data',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--model',choices=['scavenpred','ace-ipred','both'],default='both')
    p.add_argument('--inner-trees',type=int,default=40)
    p.add_argument('--candidate-limit',type=int,default=12)
    args=p.parse_args()
    for name in ['scavenpred','ace-ipred'] if args.model=='both' else [args.model]:
        evaluate(args.data,args.out,name,5 if name=='scavenpred' else 10,
                 args.inner_trees,args.candidate_limit)
