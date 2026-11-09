# Final corrected validation — 10 November 2026

## Scope
This is the corrected, executable validation using published sequence descriptor matrices. It supersedes the *performance interpretation* of the exploratory wrapper CV. Frozen historical test results remain available as archival results.

## Leakage controls
1. Five outer folds are generated using sequence-related groups computed without outcome labels. Groups join same-length peptides at >=80% identity (for length >=5) and substring/superstring pairs sharing an exact substring of length >=5.
2. Within each outer training set, three-fold grouped inner CV selects features and Random Forest parameters. **Every supervised feature selection step is fitted exclusively inside the respective inner training partition**, via an sklearn `Pipeline` and `GridSearchCV`.
3. Outer validation folds do not contribute to model tuning or selection. Every sequence has exactly one outer out-of-fold prediction.
4. The external cohorts are not loaded in the nested CV code. External scoring fits a final training model and filters close train/test sequence relations before reporting metrics. Independent-test labels are not fitted.

## Group-disjoint nested CV (out-of-fold results)

| Measure | ScavenPred | ACE-iPred |
|---|---:|---:|
| Training n | 1,648 | 1,933 |
| Sequence-related groups | 1,580 | 1,718 |
| Accuracy | 0.7032766990 | 0.7506466632 |
| MCC | 0.4033449280 | 0.4984507203 |
| ROC-AUC | 0.7819500655 | 0.8190978371 |
| PR-AUC | 0.7741452753 | 0.8202718191 |
| Confusion matrix [TN, FP; FN, TP] | [672, 199; 290, 487] | [652, 247; 235, 799] |

## Similarity-screened external evaluation

| Measure | ScavenPred | ACE-iPred |
|---|---:|---:|
| Original external n | 256 | 253 |
| Excluded as close to training | 6 | 5 |
| Evaluated n | 250 | 248 |
| Accuracy | 0.812 | 0.750 |
| MCC | 0.52603799995 | 0.46136366140 |
| ROC-AUC | 0.86035525321 | 0.84250489237 |
| PR-AUC | 0.95473403834 | 0.92948605429 |
| Confusion matrix [TN, FP; FN, TP] | [41, 13; 34, 162] | [54, 19; 43, 132] |

## Reproduction
- Group-disjoint nested CV: [strict_group_nested_validation.py](../scripts/strict_group_nested_validation.py)
- Similarity-screened external evaluation: [reproduce_external_similarity_screened.py](../scripts/reproduce_external_similarity_screened.py)
- Dataset restoration and dependencies: [reproducibility README](../README.md)

## Scientific interpretation
The out-of-fold results use a new, strictly fold-contained, grouped nested-CV model-development protocol. The external results use a sequence-screened existing external cohort, **not a prospectively collected new test set**. Thus the robust claim is that the corrected evaluation is leakage-controlled by construction and reproducible from published code and matrices. The old exploratory notebook does not constitute evidence for this corrected protocol.

The selector in this strict verification is `SelectKBest(f_classif)` in a nested sklearn pipeline; it is **not identical** to the earlier Welch/correlation/VIF/forward-wrapper research procedure. Manuscript Methods must accurately identify the method used for each reported result. The external evaluation must not be described as a historical untouched test if earlier exploratory results were consulted.

Do not substitute these metrics for a different experimental protocol without updating Methods, sample counts, model settings, and all corresponding tables and figures.
