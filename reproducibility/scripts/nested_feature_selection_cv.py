#!/usr/bin/env python3
"""Nested five-fold evaluation, independently screening features within each outer training fold.
This is a NEW audit, not a recreation of undocumented historical model development.
Inner three-fold scores are used only to choose features, not to estimate performance.
"""
import argparse,json
from pathlib import Path
import numpy as np,pandas as pd
from scipy.stats import ttest_ind
from statsmodels.stats.outliers_influence import variance_inflation_factor
from sklearn.model_selection import StratifiedKFold,cross_val_score
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score,matthews_corrcoef,roc_auc_score,confusion_matrix,make_scorer

def dataset(folder,name):
    if name=='scavenpred':
        df=pd.read_csv(folder/'scavenpred_train_88.csv')
        y=df['class'].astype(int).to_numpy()
    else:
        df=pd.read_csv(folder/'ace_train_88.csv')
        labs=pd.read_csv(folder/'ace_train_labels.csv')
        labs['seq']=labs.seq.astype(str).str.upper().str.strip()
        assert labs.seq.is_unique
        y=df.seq.astype(str).str.upper().str.strip().map(labs.set_index('seq')['class'])
        assert y.notna().all()
        y=y.to_numpy(dtype=int)
    seq=df.seq.astype(str).str.upper().str.strip()
    assert seq.is_unique and set(np.unique(y))=={0,1}
    X=df.drop(columns=['seq','class','Unnamed: 0'],errors='ignore').astype(float)
    assert X.shape[1]==88 and X.notna().all().all()
    return X,y,seq

def screen(X,y,limit):
    stats=[]
    for f in X.columns:
        p=ttest_ind(X.loc[y==1,f],X.loc[y==0,f],equal_var=False).pvalue
        corr=np.corrcoef(X[f].to_numpy(),y)[0,1]
        if np.isfinite(p) and p<.05 and np.isfinite(corr):stats.append((f,abs(corr),p))
    stats.sort(key=lambda t:(-t[1],t[0]))
    selected=[]
    for f,c,p in stats:
        if all(abs(X[f].corr(X[g]))<.7 for g in selected):selected.append(f)
    strength={f:c for f,c,_ in stats}
    candidates=selected.copy()
    while len(candidates)>1:
        matrix=np.column_stack([np.ones(len(X)),X[candidates].values])
        with np.errstate(divide='ignore',invalid='ignore'):
            vif=np.asarray([variance_inflation_factor(matrix,i+1) for i in range(len(candidates))])
        worst=np.nanmax(vif)
        if worst<10:break
        tied=[i for i,v in enumerate(vif) if np.isclose(v,worst,rtol=1e-6,atol=1e-8)]
        candidates.pop(min(tied,key=lambda i:strength[candidates[i]]))
    return candidates[:limit],dict(welch_pass=len(stats),correlation_retained=len(selected),vif_retained=len(candidates))

def run(folder,name,out,max_features,trees,candidate_limit,seed):
    X,y,seq=dataset(folder,name)
    outer=StratifiedKFold(n_splits=5,shuffle=True,random_state=seed)
    folds=[];preds=[]
    for fold,(tr,te) in enumerate(outer.split(X,y),1):
        XT=X.iloc[tr].reset_index(drop=True);yt=y[tr]
        candidates,counts=screen(XT,yt,candidate_limit)
        assert candidates
        inner=StratifiedKFold(n_splits=3,shuffle=True,random_state=seed+fold)
        selected=[];inner_steps=[];best=-float('inf')
        for step in range(max_features):
            trials=[]
            for f in candidates:
                if f in selected:continue
                fields=selected+[f]
                clf=RandomForestClassifier(n_estimators=trees,max_depth=7 if name=='scavenpred' else 10,
                    min_samples_leaf=1 if name=='scavenpred' else 4,
                    class_weight=None if name=='scavenpred' else 'balanced',
                    random_state=seed,n_jobs=1)
                scores=cross_val_score(clf,XT[fields],yt,cv=inner,scoring=make_scorer(matthews_corrcoef),n_jobs=1)
                trials.append((float(scores.mean()),f))
            if not trials:break
            score,f=max(trials,key=lambda x:(x[0],x[1]))
            if selected and score<=best:break
            best=score;selected.append(f);inner_steps.append(dict(added=f,inner_mcc=score))
        final=RandomForestClassifier(n_estimators=1000 if name=='scavenpred' else 200,
            max_depth=7 if name=='scavenpred' else 10,
            min_samples_leaf=1 if name=='scavenpred' else 4,
            max_features='sqrt' if name=='scavenpred' else 1.0,
            class_weight=None if name=='scavenpred' else 'balanced',
            n_jobs=1,random_state=42)
        final.fit(XT[selected],yt)
        probability=final.predict_proba(X.iloc[te][selected])[:,list(final.classes_).index(1)]
        predicted=final.predict(X.iloc[te][selected])
        folds.append(dict(fold=fold,n_train=len(tr),n_outer_test=len(te),features=selected,
            outer_accuracy=float(accuracy_score(y[te],predicted)),
            outer_mcc=float(matthews_corrcoef(y[te],predicted)),
            outer_roc_auc=float(roc_auc_score(y[te],probability)),
            outer_confusion_matrix=confusion_matrix(y[te],predicted,labels=[0,1]).tolist(),
            screening_counts=counts,inner_selection=inner_steps))
        for j,idx in enumerate(te):
            preds.append(dict(fold=fold,sequence=str(seq.iloc[idx]),actual=int(y[idx]),
                prediction=int(predicted[j]),positive_probability=float(probability[j])))
        print(name,'fold',fold,'selected',selected,'outer MCC',folds[-1]['outer_mcc'],flush=True)
    frame=pd.DataFrame(preds);out.mkdir(parents=True,exist_ok=True)
    frame.to_csv(out/f'{name}_nested_oof.csv',index=False)
    summary=dict(model=name,evaluation='5 outer / 3 inner folds',screening_scope='outer training only',
        independent_test_used=False,n_training_peptides=len(X),
        settings=dict(seed=seed,inner_trees=trees,candidate_limit=candidate_limit,
            max_selected_features=max_features,outer_final_trees=1000 if name=='scavenpred' else 200),
        out_of_fold_accuracy=float(accuracy_score(frame.actual,frame.prediction)),
        out_of_fold_mcc=float(matthews_corrcoef(frame.actual,frame.prediction)),
        out_of_fold_roc_auc=float(roc_auc_score(frame.actual,frame.positive_probability)),
        mean_fold_mcc=float(np.mean([f['outer_mcc'] for f in folds])),folds=folds,
        scope='New nested evaluation; does not establish historical test-set selection provenance.')
    (out/f'{name}_nested_summary.json').write_text(json.dumps(summary,indent=2))
    print(json.dumps({k:v for k,v in summary.items() if k!='folds'},indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--data',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--model',choices=['scavenpred','ace-ipred','both'],default='both')
    p.add_argument('--max-features',type=int,default=4)
    p.add_argument('--trees',type=int,default=30)
    p.add_argument('--candidate-limit',type=int,default=10)
    p.add_argument('--seed',type=int,default=42)
    a=p.parse_args()
    names=['scavenpred','ace-ipred'] if a.model=='both' else [a.model]
    for name in names:run(a.data,name,a.out,a.max_features,a.trees,a.candidate_limit,a.seed)
