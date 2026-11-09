# Corrected nested cross-validation audit — 10 November 2026

## Scope
A **new** five-outer-fold, three-inner-fold stratified evaluation on the ScavenPred and ACE-iPred **training matrices only**. There is no read, fit, feature selection, or tuning involving the independent tests. Each outer fold redoes Welch screening, class-guided Pearson redundancy filtering, VIF reduction, inner forward RF-MCC selection, and final RF fitting using only its own training partition. The outer validation partition is used solely for that fold's final prediction.

Code: [nested_feature_selection_cv.py](../scripts/nested_feature_selection_cv.py). This replaces the *performance-estimation method*, not the historical feature-selection record.

## Global out-of-fold results

| Training dataset | N | Accuracy | MCC | ROC-AUC |
|---|---:|---:|---:|---:|
| ScavenPred | 1,648 | 0.700243 | 0.397091 | 0.765949 |
| ACE-iPred | 1,933 | 0.723745 | 0.444549 | 0.795174 |

**These are cross-validated training-set estimates, not independent-test metrics.** The earlier frozen-model independent test results remain separate: ScavenPred 256, accuracy 0.765625; ACE-iPred 253, accuracy 0.727273.

## ScavenPred outer folds

| Fold | Selected features | Accuracy | MCC | AUROC |
|---|---|---:|---:|---:|
| 1 | BLOSUM6, SV2, VHSE4, KF10, KF7 | 0.6697 | 0.3367 | 0.7343 |
| 2 | ProtFP3, E5, KF10, BLOSUM5, BLOSUM6 | 0.7152 | 0.4276 | 0.8024 |
| 3 | E5, BLOSUM6, BLOSUM3, Z3, BLOSUM8 | 0.6909 | 0.3791 | 0.7497 |
| 4 | BLOSUM7, BLOSUM6, ST3, Z4, BLOSUM3 | 0.7386 | 0.4757 | 0.7799 |
| 5 | ProtFP3, BLOSUM5, E5, Z4, BLOSUM3 | 0.6869 | 0.3743 | 0.7678 |

## ACE-iPred outer folds

| Fold | Selected features | Accuracy | MCC | AUROC |
|---|---|---:|---:|---:|
| 1 | ProtFP3, BLOSUM6, BLOSUM1, BLOSUM4 | 0.7080 | 0.4130 | 0.7738 |
| 2 | T3, ProtFP5, MSWHIM2, ProtFP3, KF5 | 0.6977 | 0.3931 | 0.7747 |
| 3 | ProtFP3, ProtFP5, MSWHIM2, BLOSUM7, E4 | 0.7028 | 0.4072 | 0.7990 |
| 4 | ST3, VHSE6, ProtFP3, BLOSUM6, SVGER2, MSWHIM2 | 0.7642 | 0.5262 | 0.8312 |
| 5 | E4, ProtFP3, BLOSUM1, KF5 | 0.7461 | 0.4885 | 0.8001 |

## Fixed experimental settings
- Five outer stratified folds, random seed 42.
- Three inner stratified folds, seed 42 plus outer fold number.
- Fold-contained Welch p<0.05, class-correlation-ordered Pearson redundancy |r|<0.70, iterative VIF<10.
- At most 12 post-screening candidate descriptors, at most five selected for ScavenPred and ten for ACE-iPred; forward selection stops if the best inner MCC does not improve.
- Inner wrapper: Random Forest, 40 trees, model-specific depth and class weighting.
- Outer final fitted Random Forest: 1,000 trees (ScavenPred), 200 trees (ACE-iPred), with the model-specific controls encoded in the script.
- The independent test is excluded throughout this audit.

Run:

```bash
python reproducibility/scripts/nested_feature_selection_cv.py --data reproducibility/data/restored --out nested_results --model scavenpred --max-features 5 --candidate-limit 12 --trees 40
python reproducibility/scripts/nested_feature_selection_cv.py --data reproducibility/data/restored --out nested_results --model ace-ipred --max-features 10 --candidate-limit 12 --trees 40
```

## Limitations important to reviewers
1. This is an independently executed **new nested workflow**, not proof of how the originally published selected features and hyperparameters were historically chosen.
2. The earlier `reproduce_selection.py` screened on all training observations before five-fold wrapper CV; its wrapper CV scores are selection-biased. Do not claim those scores were leakage-free.
3. Outer CV here prevents validation-fold involvement in descriptor selection; inner wrapper estimates are used exclusively to choose features, not for reported generalization.
4. This audit does not prove that historical independent-test results were never consulted during earlier model development or negative-set motif curation. That requires separate provenance evidence.
5. Related train/test peptides remain possible without being exact duplicates; similarity-controlled generalization is an additional evaluation.

**Conclusion:** Corrected nested cross-validation removes the identified within-training outer-validation leakage in this new analysis. It does not retroactively certify the original development history or automatically validate every reported historical metric.
