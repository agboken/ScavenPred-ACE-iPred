# Final feature-selection and validation audit — 10 November 2026

## Locked protocol

The primary corrected evaluation is **five-fold nested, sequence-group-disjoint cross-validation**. Only the 88-descriptor curated training matrices are read. Peptide groups are formed without reference to class labels: (a) sequences of equal length with at least 80% positional identity for lengths >=5; (b) exact full-sequence substring/superstring relations of length >=5. Transitive connected sequences remain in the same group. Both outer five-fold and inner three-fold splits use `StratifiedGroupKFold` and are group-disjoint. A `Pipeline(SelectKBest(f_classif), RandomForestClassifier)` performs all supervised feature selection and fitting separately within each inner-training fold; grid parameters `k=5/10`, `max_depth=7/10`, `min_samples_leaf=1/4`; 100 trees; selection scoring by MCC; seed 42. Only held-out outer predictions contribute to the reported metrics.

Executable source: [strict_group_nested_validation.py](../scripts/strict_group_nested_validation.py)

## Verified results

| Metric | ScavenPred | ACE-iPred |
|---|---:|---:|
| Peptides | 1,648 | 1,933 |
| Sequence groups | 1,580 | 1,718 |
| Out-of-fold accuracy | **70.33%** | **75.06%** |
| MCC | **0.403345** | **0.498451** |
| ROC-AUC | **0.781950** | **0.819098** |
| PR-AUC (average precision) | **0.774145** | **0.820272** |
| TN / FP / FN / TP | 672 / 199 / 290 / 487 | 652 / 247 / 235 / 799 |

### Group-disjoint outer-fold results

| Model | Fold | Validation n | Accuracy | MCC |
|---|---:|---:|---:|---:|
| ScavenPred | 1 | 330 | 0.696970 | 0.390761 |
| ScavenPred | 2 | 331 | 0.709970 | 0.416875 |
| ScavenPred | 3 | 329 | 0.696049 | 0.389319 |
| ScavenPred | 4 | 329 | 0.711246 | 0.420854 |
| ScavenPred | 5 | 329 | 0.702128 | 0.401918 |
| ACE-iPred | 1 | 388 | 0.770619 | 0.538383 |
| ACE-iPred | 2 | 386 | 0.735751 | 0.470137 |
| ACE-iPred | 3 | 386 | 0.735751 | 0.469491 |
| ACE-iPred | 4 | 386 | 0.738342 | 0.473321 |
| ACE-iPred | 5 | 387 | 0.772610 | 0.542070 |

### Independent verification

The generated out-of-fold files were separately reloaded and verified: all 1,648 / 1,933 peptide sequences appear exactly once, all five outer-fold results report group-disjoint partitions, total fold counts equal dataset counts, and recomputed accuracy/MCC/ROC-AUC/PR-AUC match the saved reports to numerical precision. ACE labels are aligned by sequence.

### Reproduce

```bash
python reproducibility/scripts/unpack_data.py
python reproducibility/scripts/strict_group_nested_validation.py --data reproducibility/data/csv --out group_nested_results --model both
```

The corrected fold-contained evaluation supersedes performance claims based on the superseded feature-selection CV method. The original frozen predictor weights and legacy external-test calculations are left unchanged; this is an independently reproducible *new evaluation protocol*, not a claim about earlier experimental chronology.
