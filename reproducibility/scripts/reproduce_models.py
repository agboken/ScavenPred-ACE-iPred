#!/usr/bin/env python3
"""Reproduce published frozen RF models from the public, independently held-out data.

This script does not reselect features; see nested_feature_selection.py for a
separate train-only, exploratory leakage-controlled selection audit.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score, matthews_corrcoef, roc_auc_score, confusion_matrix

CONFIG = {
 'scavenpred': {
  'train': 'scavenpred_train_88.csv', 'test': 'scavenpred_test_88.csv',
  'features': ['E5','ProtFP3','BLOSUM6','BLOSUM3','Z3'],
  'params': dict(n_estimators=1000,max_depth=7,random_state=42,n_jobs=-1),
  'expected_test': 256,
 },
 'ace-ipred': {
  'train': 'ace_train_88.csv', 'test': 'ace_test_88.csv',
  'labels': 'ace_train_labels.csv',
  'features': ['BLOSUM4','BLOSUM1','KF5','VHSE2','BLOSUM10','SVGER2','SVGER3','ProtFP7','SVGER4','MSWHIM3'],
  'params': dict(n_estimators=200,max_depth=10,min_samples_leaf=4,max_features=1.0,class_weight='balanced',random_state=42,n_jobs=-1),
  'expected_test': 253,
 },
}


def normalized(series):
    return series.astype(str).str.strip().str.upper()


def reproduce(data_dir: Path, name: str):
    cfg = CONFIG[name]
    train = pd.read_csv(data_dir / cfg['train'])
    test = pd.read_csv(data_dir / cfg['test'])
    train['seq'] = normalized(train['seq'])
    test['seq'] = normalized(test['seq'])
    assert len(train) == len(set(train.seq)), 'Duplicate training sequences detected'
    assert len(test) == len(set(test.seq)), 'Duplicate independent-test sequences detected'
    assert not (set(train.seq) & set(test.seq)), 'Exact train/test sequence overlap detected'
    assert len(test) == cfg['expected_test'], 'Unexpected test-set size'
    if name == 'ace-ipred':
        labels = pd.read_csv(data_dir / cfg['labels'])
        label_map = dict(zip(normalized(labels.seq), labels['class']))
        y_train = train.seq.map(label_map)
        assert y_train.notna().all(), 'Missing ACE training labels'
    else:
        y_train = train['class']
    assert test['class'].notna().all(), 'Missing independent-test labels'
    columns = cfg['features']
    assert train[columns].notna().all().all() and test[columns].notna().all().all()
    model = RandomForestClassifier(**cfg['params'])
    model.fit(train[columns], y_train.astype(int))
    predictions = model.predict(test[columns]).astype(int)
    probability = model.predict_proba(test[columns])[:, list(model.classes_).index(1)]
    truth = test['class'].astype(int)
    result = dict(
        model=name,
        n_train=len(train),n_test=len(test),
        features=columns,
        accuracy=float(accuracy_score(truth,predictions)),
        mcc=float(matthews_corrcoef(truth,predictions)),
        roc_auc=float(roc_auc_score(truth,probability)),
        confusion_matrix=confusion_matrix(truth,predictions,labels=[0,1]).tolist(),
        random_seed=42,
        note='Frozen preselected-feature reproduction; does not rerun historical feature selection.',
    )
    return result, pd.DataFrame(dict(sequence=test.seq, actual=truth,prediction=predictions,positive_probability=probability))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--data',type=Path,default=Path(__file__).resolve().parents[1]/'data'/'restored')
    parser.add_argument('--output',type=Path,default=Path('reproduced_results'))
    parser.add_argument('--model', choices=['scavenpred','ace-ipred','both'], default='both')
    args=parser.parse_args()
    args.output.mkdir(parents=True,exist_ok=True)
    names=CONFIG.keys() if args.model=='both' else [args.model]
    all_results=[]
    for name in names:
        result,predictions=reproduce(args.data,name)
        all_results.append(result)
        predictions.to_csv(args.output/f'{name}_independent_predictions.csv',index=False)
        print(json.dumps(result,indent=2))
    (args.output/'metrics.json').write_text(json.dumps(all_results,indent=2))
