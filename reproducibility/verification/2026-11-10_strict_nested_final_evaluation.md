# Final strict leakage-controlled nested evaluation — 10 November 2026

## Protocol and verified outcome

This is a **new, end-to-end fold-contained evaluation** of ScavenPred and ACE-iPred on the curated training datasets. All supervised feature selection, model-parameter search and fitting happen within training partitions using a scikit-learn `Pipeline` embedded inside `GridSearchCV`. The outer five validation folds are used only for performance measurement. No independent-test dataset is read by the script.

The full runnable code is [strict_nested_validation.py](../scripts/strict_nested_validation.py). The algorithm uses **univariate ANOVA SelectKBest**, with `k ∈ {5,10}`; Random Forest `max_depth ∈ {7,10}` and `min_samples_leaf ∈ {1,4}`; 100 trees, three inner stratified folds and five outer stratified folds, random seed 42. Model selection uses inner-fold MCC. Each outer held-out peptide is predicted exactly once; summary metrics are calculated from those out-of-fold predictions.

| Metric | ScavenPred | ACE-iPred |
|---|---:|---:|
| Curated training samples | 1,648 | 1,933 |
| Unique out-of-fold predictions | 1,648 | 1,933 |
| Accuracy | **0.712379** | **0.737196** |
| MCC | **0.421646** | **0.470990** |
| ROC-AUC | **0.786466** | **0.815765** |
| Average precision (PR-AUC) | **0.775106** | **0.820776** |
| Confusion matrix (TN, FP, FN, TP) | 672, 199, 275, 502 | 632, 267, 241, 793 |

### Five outer folds

| Model | Fold | Accuracy | MCC | ROC-AUC | PR-AUC |
|---|---:|---:|---:|---:|---:|
| ScavenPred | 1 | 0.693939 | 0.386693 | 0.780276 | 0.766673 |
| ScavenPred | 2 | 0.739394 | 0.476992 | 0.813182 | 0.810865 |
| ScavenPred | 3 | 0.696970 | 0.390761 | 0.761494 | 0.747030 |
| ScavenPred | 4 | 0.729483 | 0.456206 | 0.800260 | 0.801358 |
| ScavenPred | 5 | 0.702128 | 0.404420 | 0.783760 | 0.766538 |
| ACE-iPred | 1 | 0.746770 | 0.490193 | 0.816009 | 0.818506 |
| ACE-iPred | 2 | 0.720930 | 0.438122 | 0.805314 | 0.811622 |
| ACE-iPred | 3 | 0.710594 | 0.418357 | 0.790338 | 0.795947 |
| ACE-iPred | 4 | 0.743523 | 0.484535 | 0.827225 | 0.853172 |
| ACE-iPred | 5 | 0.764249 | 0.524880 | 0.840958 | 0.836780 |

### Reproduction

Restore the versioned data archive first:

```bash
python reproducibility/scripts/unpack_data.py
python reproducibility/scripts/strict_nested_validation.py --data reproducibility/data/restored --out strict_results --model both --trees 100
```

### Verified controls

1. Only versioned training feature matrices and ACE training labels are loaded; independent-test files are not used.
2. Every fold fits its feature selector and hyperparameter selection solely on the corresponding outer-training data; the selector is inside the inner-CV Pipeline.
3. Both training datasets have 88 descriptor columns and unique peptide sequences. ACE labels are matched by sequence, not row number.
4. All 1,648 and 1,933 peptides respectively receive exactly one held-out prediction. Summary metrics were independently recomputed from the out-of-fold CSVs and matched the generated reports.

### Implementation identity

This is the **final corrected nested evaluation protocol** with ANOVA feature screening and fold-contained RF parameter tuning. The public web predictors continue to use their pre-existing frozen models, so these new nested metrics must be attributed to this **new evaluated pipeline**, not represented as a rerun of the old global Welch/VIF/wrapper procedure or proof of a prospectively unopened historical external test.
