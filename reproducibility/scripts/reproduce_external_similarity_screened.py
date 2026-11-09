"""Train exclusively on original training data; screen external by sequence relatedness before scoring.
Selected hyperparameters and features fitted in CV on the training partition only.
"""
import json
import argparse
from pathlib import Path
import numpy as np,pandas as pd
from sklearn.pipeline import Pipeline
from sklearn.feature_selection import SelectKBest,f_classif
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedGroupKFold,GridSearchCV
from sklearn.metrics import accuracy_score,matthews_corrcoef,roc_auc_score,average_precision_score,confusion_matrix,make_scorer
from strict_group_nested_validation import dataset,make_groups

parser=argparse.ArgumentParser()
parser.add_argument('--data',type=Path,required=True)
parser.add_argument('--scavenpred-test',type=Path,required=True)
parser.add_argument('--out',type=Path,required=True)
args=parser.parse_args()
folder=args.data;out=args.out;out.mkdir(parents=True,exist_ok=True)
def related(a,b):
 if len(a)==len(b) and len(a)>=5 and sum(x!=y for x,y in zip(a,b))<=int(np.floor(.2*len(a)+1e-9)):return True
 if min(len(a),len(b))>=5 and (a in b or b in a):return True
 return False
for name,test_path in [('scavenpred',str(args.scavenpred_test)),('ace-ipred',str(folder/'ace_test_88.csv'))]:
 seq,X,y=dataset(folder,name)
 test=pd.read_csv(test_path);test.seq=test.seq.astype(str).str.upper().str.strip()
 assert test.seq.is_unique and not(set(seq)&set(test.seq))
 # Recover any missing global-mean descriptor using TRAINING sequence features only.
 for missing in [c for c in X.columns if c not in test.columns]:
  alphabet='ACDEFGHIKLMNPQRSTVWY'
  matrix=np.array([[s.count(c)/len(s) for c in alphabet] for s in seq])
  weights,_,rank,_=np.linalg.lstsq(matrix,X[missing].to_numpy(float),rcond=None)
  residual=np.max(np.abs(matrix@weights-X[missing].to_numpy(float)))
  assert rank==len(alphabet) and residual<1e-10,(missing,residual)
  matrix_test=np.array([[s.count(c)/len(s) for c in alphabet] for s in test.seq])
  test[missing]=matrix_test@weights
  print(name,'reconstructed',missing,'training-only residue averaging residual',residual,flush=True)
 Xtest=test[X.columns].astype(float);truth=test['class'].astype(int)
 overlap=np.array([any(related(q,t) for t in seq) for q in test.seq]);keep=~overlap
 print(name,'independent total',len(test),'similarity_excluded',overlap.sum(),'remaining',keep.sum(),flush=True)
 groups=make_groups(seq)
 pipeline=Pipeline([('select',SelectKBest(f_classif)),('rf',RandomForestClassifier(n_estimators=100,random_state=42,n_jobs=1))])
 inner=StratifiedGroupKFold(n_splits=3,shuffle=True,random_state=42)
 grid=GridSearchCV(pipeline,{'select__k':[5,10],'rf__max_depth':[7,10],'rf__min_samples_leaf':[1,4]},cv=inner,scoring=make_scorer(matthews_corrcoef),n_jobs=1,refit=True,error_score='raise')
 grid.fit(X,y,groups=groups)
 pred=grid.predict(Xtest.loc[keep]);proba=grid.predict_proba(Xtest.loc[keep])[:,1];target=truth.loc[keep]
 report={'model':name,'n_train':len(X),'n_external_total':len(test),'n_removed_related':int(overlap.sum()),'n_external_screened':int(keep.sum()),'selected_features':list(X.columns[grid.best_estimator_.named_steps['select'].get_support()]),'params':grid.best_params_,'accuracy':float(accuracy_score(target,pred)),'mcc':float(matthews_corrcoef(target,pred)),'roc_auc':float(roc_auc_score(target,proba)),'pr_auc':float(average_precision_score(target,proba)),'confusion_matrix_0_1':confusion_matrix(target,pred,labels=[0,1]).tolist(),'test_labels_not_used_for_training':True}
 result=pd.DataFrame({'sequence':test.seq,'excluded_related_to_train':overlap,'actual':truth,'prediction':np.nan,'positive_probability':np.nan})
 result.loc[keep,'prediction']=pred;result.loc[keep,'positive_probability']=proba
 result.to_csv(out/(name+'_external_full_audit.csv'),index=False)
 (out/(name+'_external_report.json')).write_text(json.dumps(report,indent=2))
 print('DONE',name,report,flush=True)
