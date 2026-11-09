# Leakage re-audit — 10 November 2026

Scope: public frozen-model scripts, published training/test feature matrices, and published ScavenPred feature-selection implementation.

## Confirmed

| Item | ScavenPred | ACE-iPred |
|---|---:|---:|
| Training sequences | 1,648 | 1,933 |
| Independent-test sequences | 256 | 253 |
| Exact sequence overlap between partitions | 0 | 0 |
| Within-partition exact duplicates | 0 in both | 0 in both |

- The frozen model reproduction script fits its Random Forest on training descriptors/labels and evaluates on a distinct independent-test matrix. The test matrix is not passed to `fit`.
- The ACE training labels are joined to the training feature matrix by peptide sequence.
- Published matrices contain no missing values.
- The frozen features and RF hyperparameters are predefined in the reproduction script. This **does not prove** they were chosen without previously consulting the independent test set during historical model development.

## Leakage issue found: ScavenPred feature-selection CV

The public `reproducibility/scripts/reproduce_selection.py` performs univariate Welch tests, class-correlation ranking, Pearson redundancy filtering and VIF elimination on **all 1,648 training peptides before** passing subsets to five-fold cross-validation inside the forward wrapper.

Thus the validation folds contribute their labels to the feature-screening process that determines the wrapper's candidates, a form of *within-training cross-validation leakage*. These five-fold wrapper scores are **not unbiased held-out performance estimates**. Repeatedly using the same folds to select additions also biases the highest wrapper score upward.

This issue does not, by itself, imply leakage into the distinct 256-peptide independent test. It does mean the present code must **not** be described as fully leakage-free five-fold feature-selection evaluation.

## Relatedness checks

- ScavenPred: at least 8 training/test same-length peptide pairs differ by only one substitution; at least 6 independent-test sequences participate in a length >=5 exact substring/superstring relation across partitions.
- ACE-iPred: 171 same-length train/test peptide pairs differ by only one substitution (many are short peptides); at least 5 test sequences have a length >=5 substring/superstring relation to a training sequence.

These findings are **not exact duplicate leakage**. They are a reason to separately report similarity-controlled evaluation, especially for short peptides.

## Remaining provenance questions

- Was the test partition fixed before motif-based negative curation, descriptor screening and hyperparameter selection?
- Were motif exclusion rules for putative negatives derived exclusively from the training partition?
- Were independent-test labels or scores ever used to change features, thresholds, or hyperparameters?
- Can the originally reported five-fold CV pipeline be rerun with **all** supervised screening and feature selection fitted separately inside each outer training fold?

## Required correction

1. For a valid five-fold performance estimate, place **Welch screening, class-correlation ranking, redundancy filtering, VIF, wrapper selection and RF fitting** inside each **outer training fold**. Use inner CV within each outer training fold to select features/model settings.
2. Evaluate once on each untouched outer validation fold and aggregate the five held-out scores. Never use outer validation outcomes to revise the pipeline.
3. Use the independent test only once after selecting model settings from training data; disclose if it was historically consulted more often.
4. Distinguish independent-test frozen-model reproduction from a full leakage-free reconstruction of model-development history.

**Current conclusion:** exact train/test separation and fitting of the frozen RFs are verified. Full elimination of information leakage in feature-selection CV and in the original historical development workflow is **not yet verified**. No revised claims should be made until corrected nested evaluation is executed and documented.
