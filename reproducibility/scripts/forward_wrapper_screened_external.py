#!/usr/bin/env python3
"""External evaluation after training-only Welch/Pearson/VIF and grouped RF forward wrapper.

The same wrapper as grouped_nested_forward_wrapper_cv.py is fitted on the
complete training data. The external cohort is sequence-screened before
evaluation; its labels do not enter feature selection or model fitting.
"""
import argparse,json
from pathlib import Path
import numpy as np,pandas as pd
from sklearn.model_selection import StratifiedGroupKFold
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score,matthews_corrcoef,roc_auc_score,average_precision_score,confusion_matrix
from strict_group_nested_validation import dataset,make_groups
from grouped_nested_forward_wrapper_cv import screen,rf_parameters
import sklearn

def related(a,b):
    if len(a)==len(b) and len(a)>=5:
        if sum(x!=y for x,y in zip(a,b))<=int(np.floor(.2*len(a)+1e-9)):return True
    return min(len(a),len(b))>=5 and (a in b or b in a)

def run(data,out,name,inner_trees=40,candidate_limit=12):
    seq,X,y=dataset(data,name);groups=make_groups(seq);y=y.to_numpy()
    filename='scavenpred_test_88.csv' if name=='scavenpred' else 'ace_test_88.csv'
    test=pd.read_csv(data/filename)
    test.seq=test.seq.astype(str).str.strip().str.upper()
    assert test.seq.is_unique and set(test.seq).isdisjoint(set(seq))
    candidates,counts=screen(X,y,candidate_limit)
    inner=StratifiedGroupKFold(n_splits=3,shuffle=True,random_state=42)
    splits=list(inner.split(X,y,groups))
    cache={};selected=[];steps=[];previous=-np.inf
    def score(features):
        key=tuple(features)
        if key not in cache:
            scores=[]
            for it,iv in splits:
                model=RandomForestClassifier(**rf_parameters(name,inner_trees))
                model.fit(X.iloc[it][features],y[it])
                scores.append(matthews_corrcoef(y[iv],model.predict(X.iloc[iv][features])))
            cache[key]=float(np.mean(scores))
        return cache[key]
    for _ in range(min(5 if name=='scavenpred' else 10,len(candidates))):
        options=[(score(selected+[f]),f) for f in candidates if f not in selected]
        best,f=max(options,key=lambda v:(v[0],v[1]))
        if selected and best<=previous:break
        selected.append(f);previous=best;steps.append({'feature':f,'inner_group_mcc':best})
    # Reconstruct a selected mean-amino-acid descriptor only from training
    # sequences/descriptors, never from external labels.
    for f in selected:
        if f not in test.columns:
            aa='ACDEFGHIKLMNPQRSTVWY'
            A=np.asarray([[s.count(c)/len(s) for c in aa] for s in seq])
            weights,_,rank,_=np.linalg.lstsq(A,X[f].to_numpy(),rcond=None)
            residual=float(np.max(np.abs(A@weights-X[f].to_numpy())))
            assert rank==20 and residual<1e-9,(f,residual)
            B=np.asarray([[s.count(c)/len(s) for c in aa] for s in test.seq])
            test[f]=B@weights
    overlap=np.asarray([any(related(q,t) for t in seq) for q in test.seq])
    kept=~overlap
    model=RandomForestClassifier(**rf_parameters(name,1000 if name=='scavenpred' else 200))
    model.fit(X[selected],y)
    pred=model.predict(test.loc[kept,selected])
    prob=model.predict_proba(test.loc[kept,selected])[:,list(model.classes_).index(1)]
    truth=test.loc[kept,'class'].astype(int).to_numpy()
    out.mkdir(parents=True,exist_ok=True)
    full=pd.DataFrame({'sequence':test.seq,'excluded_related_to_train':overlap,
                       'actual':test['class'].astype(int),'prediction':np.nan,'positive_probability':np.nan})
    full.loc[kept,'prediction']=pred
    full.loc[kept,'positive_probability']=prob
    full.to_csv(out/(name+'_wrapper_external_predictions.csv'),index=False)
    report={'model':name,'training_n':len(X),'external_before_screening':len(test),
            'excluded_related':int(overlap.sum()),'external_evaluated':int(kept.sum()),
            'selected_features':selected,'selection_steps':steps,'screening_counts':counts,
            'accuracy':float(accuracy_score(truth,pred)),
            'mcc':float(matthews_corrcoef(truth,pred)),
            'roc_auc':float(roc_auc_score(truth,prob)),
            'pr_auc':float(average_precision_score(truth,prob)),
            'confusion_matrix':confusion_matrix(truth,pred,labels=[0,1]).tolist(),
            'test_labels_used_in_fitting':False,'sklearn':sklearn.__version__}
    (out/(name+'_wrapper_external_report.json')).write_text(json.dumps(report,indent=2))
    print(name,json.dumps(report,indent=2),flush=True)

if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--data',type=Path,required=True)
    parser.add_argument('--out',type=Path,required=True)
    parser.add_argument('--model',choices=['both','scavenpred','ace-ipred'],default='both')
    parser.add_argument('--inner-trees',type=int,default=40)
    parser.add_argument('--candidate-limit',type=int,default=12)
    args=parser.parse_args()
    for name in ['scavenpred','ace-ipred'] if args.model=='both' else [args.model]:
        run(args.data,args.out,name,args.inner_trees,args.candidate_limit)
