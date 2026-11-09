# Original ScavenPred / ACE-iPred independent-test provenance audit

**Audit date:** 10 November 2026. **Question:** Can the available historical record demonstrate that the independent test sets were *never consulted* for model selection or feature selection?

## Conclusion

**No. An absolute historical non-consultation assertion is not supported by the recovered records.** The preserved original ScavenPred notebook contains multiple alternative training/test descriptor subsets and executed evaluation of test-labelled records. Its saved execution counters are nonsequential and cannot reliably establish the order in which the tests, decisions, and feature changes happened. Consequently, it does **not** prove that test feedback determined the final chosen features, but it prevents an evidence-based claim that the test was never accessed or consulted during exploratory model development.

There is an important distinction between (a) the historical exploratory workflow, (b) the later August training-only corrective reanalysis, and (c) the independently executed November 2026 nested CV and frozen-model reruns. Evidence for (b)/(c) cannot rewrite the events of (a).

## Direct evidence from archived original ScavenPred notebook

Archived notebook: `ScavenPred_ACEiPred_Full_Reproducibility_Package/ScavenPred/notebooks/Scavenpred_original.ipynb`; SHA-256 `e77a7c229d50981c04865792e550ace791a375d9d1d128359e019e5df9d79603`. The companion CSV lists the saved execution counter, cell purpose, excerpt, and source SHA-256 for each relevant cell.

- **Cells 36–40:** multiple training-descriptor subsets are specified, including 88-feature and smaller alternatives.
- **Cells 47–50:** multiple test-descriptor subsets are specified.
- **Cells 51–52:** predictions and test probabilities are generated.
- **Cell 77:** a `test_data_features.xlsx` file is loaded.
- **Cell 83:** the test labels and predictions are used in a confusion matrix and classification report.
- **Cells 86–101 and 145:** alternative classifiers are evaluated and ROC/accuracy results are calculated using test records.

**Inference limit:** Notebook cells are snapshots of an interactive development session. Their textual order is not a time-stamped sequence, their execution counts are not sufficient to infer unique chronological experiments, and the provenance of each external `df_test` dataset cannot be proved from this notebook alone. Some test-data operations could concern additional application or benchmark sets; they must not all be conflated with the final 256-sample held-out set without source tracing.

## Archive and chronology

- August 23, 2026: contained file ZIP timestamps in the assembled reproducibility archive for the cleaned datasets, selection workbooks and frozen notebooks. These demonstrate archival packaging **no later than that assembly**, not pre-experiment sealing or immutable dataset creation dates.
- October 14, 2026: package of additional retraining results exists. This confirms subsequent model investigation; it does not establish original test-set sealing.
- November 9–10, 2026: public GitHub commits add scripts, datasets, verification records and corrected nested CV. These are verifiable documentation of *later reproduction*, not proof of earlier nonconsultation.

Master archive SHA-256: `6a5cff96a8c4e9117b0c4cc69bb41f47df64c9a36bb5557544f4779e5ec86eaa`.

## What can be supported now

1. **Current frozen models:** The published reproduction script fits the Random Forest on training features only and calls `predict`/`predict_proba` on disjoint test features. Exact sequence overlaps in these curated datasets were checked and reported as zero. This is evidence about the *reproduction script*, not the entire historical selection process.
2. **Corrected nested CV:** Selection is trained within each outer fold in the new audit, preventing the original validation-fold feature-screening leakage. The nested results are valid as a *new corrective evaluation*, not a historical recreation of earlier testing decisions.
3. **Original feature selection:** Insufficient records to assert test labels/results never informed model development. The historical notebook presents a substantive provenance concern, not a clean nonconsultation proof.
4. **ACE-iPred:** The recovered frozen reproducibility notebook records the later clean procedure. There is not an independently dated full original decision log sufficient to establish all early ACE-iPred choices preceded all test-set access.
5. **Negative-class construction:** The data lineage of motif exclusions and negative curation may need further work; exact sequence disjointness does not prove training-only derivation of data-curation rules.

## Reviewer-safe resolution

The scientifically defensible resolution is **transparent correction rather than retrospective certification**:

- Separate the *historical exploratory analyses* from the **corrected nested cross-validation**, explicitly acknowledging the original feature-selection design limitation.
- Report new nested-CV performance as the primary internal validation, using fold-contained feature selection and tuning. Clearly state that models/feature subsets may change across folds.
- State the independent holdout was **excluded from the execution of the corrected fitting and selection pipeline**. Do not claim it was *never inspected at any time in historical research* unless independent contemporaneous evidence is subsequently found.
- Treat the previously used test sets as **legacy independent-test datasets**, not guaranteed completely untouched prospective holdouts. Do not imply that re-evaluating them after observing their labels restores their original untouched status.
- If a fully untouched benchmark is essential, acquire or sequester a **genuinely new external dataset**, with time-stamped version, fixed prediction protocol before revealing labels, code hashes, and contemporaneous evaluation logs.

## Safe replacement manuscript sentence

> For the corrected analysis, all descriptor screening, feature selection, model selection and internal validation were performed within the training dataset using nested cross-validation. The reported outer-fold scores were computed on folds excluded from their respective feature-selection and training procedures. The existing external test sets were evaluated using frozen configurations; because the original exploratory development history cannot establish that these sets remained completely unconsulted, their results are reported as legacy external evaluations rather than as prospectively untouched validation.

## Provenance limits

No unverified chronology or fabricated date is used as positive evidence of non-consultation. The analysis does not conclude that final feature selection **definitely used** test outcomes—only that the retrieved evidence cannot rule that out. This distinction should remain explicit in editorial correspondence.