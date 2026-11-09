# Reconciliation of ScavenPred and ACE-iPred evaluation results

**Verified 10 November 2026.** Multiple correctly archived evaluations were previously displayed without enough distinction. The two **canonical corrected evaluation** protocols are the sequence-group-disjoint nested CV and similarity-screened external evaluation below.

## Evaluation comparison

| Protocol | ScavenPred n | ScavenPred accuracy | MCC | ROC-AUC | ACE-iPred n | ACE-iPred accuracy | MCC | ROC-AUC |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Legacy frozen-model independent test | 256 | 76.5625% | 0.377006 | 0.823609 | 253 | 72.7273% | 0.416194 | 0.785218 |
| Earlier training-only nested CV (not group-disjoint) | 1,648 | 70.0243% | 0.397091 | 0.765949 | 1,933 | 72.3745% | 0.444549 | 0.795174 |
| **Canonical corrected group-disjoint nested CV** | **1,648** | **70.3277%** | **0.403345** | **0.781950** | **1,933** | **75.0647%** | **0.498451** | **0.819098** |
| **Canonical corrected similarity-screened external evaluation** | **250** | **81.20%** | **0.526038** | **0.860355** | **248** | **75.00%** | **0.461364** | **0.842505** |

The corresponding canonical PR-AUC values are:
- Group-disjoint nested CV: ScavenPred **0.774145**, ACE-iPred **0.820272**
- Similarity-screened external: ScavenPred **0.954734**, ACE-iPred **0.929486**

## Exactly which results to use

1. **Cross-validation:** [Final grouped nested CV](2026-11-10_final_group_disjoint_nested_cv.md). This evaluates every training sample once as an outer-fold holdout and groups closely related peptide sequences. It uses inner grouped CV, `SelectKBest(f_classif)`, and RF grid tuning. It is **not** the historical Welch/VIF/forward-wrapper pipeline.
2. **External:** [Final similarity-screened external](2026-11-10_final_external_similarity_screened.md). It excludes 6 ScavenPred and 5 ACE-iPred close relatives of training sequences before external scoring, trains on the original training cohort using grouped inner CV, and evaluates a *newly selected RF model* on 250 / 248 sequences.

## Where older numbers came from

- [Frozen-model rerun](2026-11-10_independent_test_rerun.json): preselected features of an earlier frozen RF evaluated on 256 / 253. Not the newly selected model and not similarity-screened.
- [Earlier corrected nested-CV audit](2026-11-10_corrected_nested_cv.md): different feature-selection method/split protocol, no sequence-group isolation. Kept only for audit history.
- The web-deployed ScavenPred and ACE-iPred predictors are not automatically updated by these experimental reruns. **Do not claim the externally validated revised model is the live web model** unless separately deployed and checked.

## Arithmetic verification

- Grouped ScavenPred: (672 + 487) / 1,648 = 0.703276699
- Grouped ACE-iPred: (652 + 799) / 1,933 = 0.750646663
- Screened ScavenPred: (41 + 162) / 250 = 0.812
- Screened ACE-iPred: (54 + 132) / 248 = 0.75

Each canonical result matches its recorded confusion matrix and saved local prediction-level evaluation. All canonical values are traceable to the corresponding published scripts and final GitHub reports.

**Manuscript rule:** use a single consistent model/protocol per table and ensure its sample counts, methods, features, and reported metrics are paired. Do not mix the old frozen RF or earlier nested-CV numbers with the canonical corrected validation.
