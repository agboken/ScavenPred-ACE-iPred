#!/usr/bin/env python3
"""Strict nested CV with sequence-related groups disjoint between folds.
Unsupervised groups: >=80% same-length identity or exact substring of length>=5.
All supervised selection and tuning performed inside inner training partitions.
"""
from pathlib import Path
import argparse, json, hashlib, platform
from collections import defaultdict,Counter
import numpy as np,pandas as pd
from sklearn.pipeline import Pipeline
from sklearn.feature_selection import SelectKBest,f_classif
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedGroupKFold,GridSearchCV
from sklearn.metrics import accuracy_score,matthews_corrcoef,roc_auc_score,average_precision_score,confusion_matrix,make_scorer
import sklearn,scipy

CONFIG={'scavenpred':('scavenpred_train_88.csv',None),'ace-ipred':('ace_train_88.csv','ace_train_labels.csv')}
def dataset(root,name):
 f,l=CONFIG[name];df=pd.read_csv(root/f); seq=df.seq.astype(str).str.upper().str.strip()
 assert seq.is_unique
 if l:
  labels=pd.read_csv(root/l);labels.seq=labels.seq.astype(str).str.upper().str.strip()
  assert labels.seq.is_unique and len(labels)==len(df)
  y=seq.map(labels.set_index('seq')['class']);assert y.notna().all()
 else:y=df['class']
 X=df.drop(columns=['seq','class','Unnamed: 0','index'],errors='ignore').astype(float)
 assert X.shape[1]==88 and X.notna().all().all() and set(y.unique())=={0,1}
 return seq,X,y.astype(int)
def make_groups(seq):
 ss=seq.tolist();n=len(ss);parent=list(range(n));size=[1]*n
 def root(i):
  while parent[i]!=i:parent[i]=parent[parent[i]];i=parent[i]
  return i
 def link(a,b):
  a,b=root(a),root(b)
  if a!=b:
   if size[a]<size[b]:a,b=b,a
   parent[b]=a;size[a]+=size[b]
 same_len=defaultdict(list)
 for i,s in enumerate(ss):same_len[len(s)].append(i)
 for length,ids in same_len.items():
  if length<5:continue
  tolerance=int(np.floor(length*.2+1e-9))
  for j,i in enumerate(ids):
   for k in ids[j+1:]:
    if sum(c!=d for c,d in zip(ss[i],ss[k]))<=tolerance:link(i,k)
 ix=defaultdict(list)
 for i,s in enumerate(ss):
  if len(s)>=5:ix[s].append(i)
 for j,s in enumerate(ss):
  for length in range(5,len(s)):
   for a in range(len(s)-length+1):
    for i in ix.get(s[a:a+length],[]):link(i,j)
 roots=[root(i) for i in range(n)];mapping={r:j for j,r in enumerate(sorted(set(roots)))}
 return np.array([mapping[r] for r in roots])
def run(root,name,out,trees=100,seed=42):
 seq,X,y=dataset(root,name);groups=make_groups(seq)
 outer=StratifiedGroupKFold(n_splits=5,shuffle=True,random_state=seed)
 oof=[];details=[]
 for fold,(tr,va) in enumerate(outer.split(X,y,groups),1):
  gtr=groups[tr];gva=groups[va];assert set(gtr).isdisjoint(set(gva))
  model=Pipeline([('select',SelectKBest(f_classif)),('rf',RandomForestClassifier(n_estimators=trees,random_state=seed,n_jobs=1))])
  inner=StratifiedGroupKFold(n_splits=3,shuffle=True,random_state=seed+fold)
  grid=GridSearchCV(model,{'select__k':[5,10],'rf__max_depth':[7,10],'rf__min_samples_leaf':[1,4]},scoring=make_scorer(matthews_corrcoef),cv=inner,n_jobs=1,error_score='raise',refit=True)
  grid.fit(X.iloc[tr],y.iloc[tr],groups=gtr)
  pred=grid.predict(X.iloc[va]);prob=grid.predict_proba(X.iloc[va])[:,list(grid.classes_).index(1)]
  chosen=list(X.columns[grid.best_estimator_.named_steps['select'].get_support()])
  item={'fold':fold,'n_train':len(tr),'n_validation':len(va),'n_validation_groups':len(set(gva)),'group_disjoint':True,'features':chosen,'best_params':grid.best_params_,'accuracy':float(accuracy_score(y.iloc[va],pred)),'mcc':float(matthews_corrcoef(y.iloc[va],pred)),'roc_auc':float(roc_auc_score(y.iloc[va],prob)),'average_precision':float(average_precision_score(y.iloc[va],prob))}
  details.append(item)
  for s,truth,p,pb in zip(seq.iloc[va],y.iloc[va],pred,prob):oof.append({'sequence':s,'fold':fold,'truth':int(truth),'prediction':int(p),'positive_probability':float(pb)})
  print(name,'fold',fold,'n',len(va),'mcc',round(item['mcc'],3),flush=True)
 oof=pd.DataFrame(oof)
 assert len(oof)==len(seq) and oof.sequence.is_unique and set(oof.sequence)==set(seq)
 report={'model':name,'design':'5x3 group-disjoint nested CV; unsupervised near-identity and shared-substring groups','training_n':len(seq),'group_n':len(set(groups)),'oof_accuracy':float(accuracy_score(oof.truth,oof.prediction)),'oof_mcc':float(matthews_corrcoef(oof.truth,oof.prediction)),'oof_roc_auc':float(roc_auc_score(oof.truth,oof.positive_probability)),'oof_pr_auc':float(average_precision_score(oof.truth,oof.positive_probability)),'confusion_matrix':confusion_matrix(oof.truth,oof.prediction,labels=[0,1]).tolist(),'independent_test_read':False,'folds':details,'versions':{'python':platform.python_version(),'sklearn':sklearn.__version__,'pandas':pd.__version__}}
 out.mkdir(parents=True,exist_ok=True);oof.to_csv(out/(name+'_group_nested_oof.csv'),index=False);(out/(name+'_group_nested_report.json')).write_text(json.dumps(report,indent=2))
 print('DONE',name,{k:report[k] for k in ['group_n','oof_accuracy','oof_mcc','oof_roc_auc','oof_pr_auc']},flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--data',type=Path,required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--model',choices=['both','scavenpred','ace-ipred'],default='both');a=p.parse_args()
 for name in CONFIG if a.model=='both' else [a.model]:run(a.data,name,a.out)
