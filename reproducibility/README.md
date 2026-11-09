# ScavenPred / ACE-iPred — Public reproducibility release

> **CURRENT RESULTS (10 November 2026):** Use the [reconciled results index](verification/RESULTS_RECONCILIATION.md) for the authoritative corrected evaluation. **Grouped nested CV**: ScavenPred 70.33% (n=1,648), ACE-iPred 75.06% (n=1,933). **Similarity-screened external evaluation**: ScavenPred 81.20% (n=250), ACE-iPred 75.00% (n=248). Older frozen-model and preliminary nested-CV values below are retained only for reproducibility history and must not be mixed with these final protocols.


**Status (10 November 2026):** Research datasets, corrected training-only ScavenPred feature-selection code, frozen Random Forest reproduction code, model configurations, independent-test prediction outputs and verification records are publicly available here.

## 1. Reproduce from a fresh clone

```bash
git clone https://github.com/agboken/ScavenPred.git
cd ScavenPred
python -m pip install -r reproducibility/requirements.txt
python reproducibility/scripts/restore_data.py
python reproducibility/scripts/reproduce_selection.py --out scavenpred_feature_selection.json
python reproducibility/scripts/reproduce_models.py --output reproduced_results
```

The complete training-only five-fold wrapper takes longer than the model prediction check. On resource-limited hardware, set `OPENBLAS_NUM_THREADS=1` before running the feature-selection script.

The **optional** independent nested-CV audit (a different experimental method, not the original manuscript wrapper) is:

```bash
python reproducibility/scripts/nested_feature_selection.py --output nested_cv_results.json
```

## 2. Public datasets and provenance

The public `data/` folder contains:

- `ScavenPred_ACEiPred_dataset_bundle.tar.gz.b64` — compressed CSV archive (SHA-256: `5361826554c55de3f0f05b725f6c86eaa0d1c701a985728b28366fc12d6fb6bc`).
- `scavenpred_test_256_features.csv.gz.b64` — full independent-test descriptor matrix, packaged separately.

The restoration script verifies the archived checksum, decodes the compressed files and exports 16 normal CSV files to `reproducibility/data/restored/`.

The dataset bundle is additionally mirrored as **23 independently accessible parts** under [`data/archive_parts/`](data/archive_parts/). To reconstruct these with the separate checksum-verifying alternative, run `python reproducibility/scripts/unpack_data.py` (outputs to `reproducibility/data/csv/`). The ScavenPred independent-test descriptor matrix is also a directly readable CSV at [`data/scavenpred_test_88.csv`](data/scavenpred_test_88.csv). The public archive includes finalized training matrices and labels, test labels, ScavenPred screening / correlation / VIF / wrapper / hyperparameter results, and saved independent-test predictions.

No data are downloaded from the private PepSilico repository. The frozen evaluation scripts use the research datasets supplied here.

## 3. Training-only feature-selection reproduction

`scripts/reproduce_selection.py` operates **only on the 1,648-peptide ScavenPred training matrix** and does not read the independent-test set. It independently implements:

1. Welch unpaired t-tests (`p < 0.05`), reducing 88 descriptors to 62.
2. Greedy class-signal-ordered correlation filtering (`|r| < 0.70`), retaining 32 descriptors.
3. Iterative VIF filtering (`VIF < 10`) with a weaker-class-correlation tie-break, retaining 14.
4. Forward selection by mean MCC from stratified five-fold RF cross-validation (100 trees, max depth 7, shuffle=True, seed=42). Selected descriptors: **E5, ProtFP3, BLOSUM6, BLOSUM3 and Z3**.

The six-step MCC trace matches the archived selection results, including the sixth candidate being rejected. See `verification/feature_selection_verification.json`. These *training CV selection scores* should not be represented as independent-test performance or as an unbiased nested-CV estimate.

For a strictly fold-contained *alternative* feature-selection assessment, `scripts/nested_feature_selection.py` fits its selector and tunes hyperparameters inside training partitions of outer CV folds, and does not access the independent-test data.

## 4. Frozen-model reproduction and independent evaluation

`scripts/reproduce_models.py` reads only published training and independent-test CSVs, verifies no exact train/test sequence overlap and uses the fixed feature subsets.

| Model | Training peptides | Independent test | Accuracy | MCC | ROC AUC | TN / FP / FN / TP |
|---|---:|---:|---:|---:|---:|---|
| ScavenPred (1,000 trees, depth 7) | 1,648 | 256 | 0.765625 | 0.377006 | 0.823609 | 33 / 22 / 38 / 163 |
| ACE-iPred (200 trees, depth 10, leaf 4, balanced) | 1,933 | 253 | 0.727273 | 0.416194 | 0.785218 | 53 / 21 / 48 / 131 |

Saved predictions and full-precision metrics are in `verification/`, including [`2026-11-10_independent_test_rerun.json`](verification/2026-11-10_independent_test_rerun.json). These evaluations were independently reproduced from the supplied matrices on 10 November 2026. They are not claims that every other analysis in the manuscript was rerun.

## 5. Notebooks

- `../Scavenpred.ipynb` — original, historical development notebook. It is **not** the corrected feature-selection implementation.
- `notebooks/ScavenPred_complete_training_and_analysis.ipynb` and `notebooks/ACE_iPred_complete_training_and_analysis.ipynb` — archived analysis and research-result notebooks, with outputs removed.
- `notebooks/ScavenPred_frozen_reproduction.ipynb` and `notebooks/ACE_iPred_frozen_reproduction.ipynb` — concise frozen model-fitting notebooks.

For code-level reproduction, prefer the scripts in `scripts/`.

## 6. Public tools and protected research

ScavenPred: https://pepsilico.com/tools/scavenpred/

ACE-iPred: https://pepsilico.com/tools/ace-ipred/

**SAPIEN-Pep / SIPPAD are distinct pre-publication research systems and no proprietary code or engine archives are included here.**

## 7. Scope and limitations

The archived data are curated research datasets and include putative negatives, not exclusively experimentally verified inactive peptides. There is no claim that independent predictions are experimental measurements. The reproducibility archive exposes features and labels; it does not provide independent evidence that every source descriptor value can be rederived from the original peptide databases. Please cite the associated manuscript when its final bibliographic details are available.


## Corrected nested cross-validation audit (November 2026)

The executable [nested_feature_selection_cv.py](scripts/nested_feature_selection_cv.py) reruns Welch univariate screening, class-guided Pearson redundancy removal, VIF filtering, inner three-fold forward-selection and RF fitting **from scratch within each outer training fold**. The outer five held-out partitions are not available to these steps. The script uses the published training feature matrices and does not load the independent-test set.

To rerun after restoring datasets:

```bash
python reproducibility/scripts/nested_feature_selection_cv.py --data reproducibility/data/restored --out nested_results --model scavenpred --max-features 5 --candidate-limit 12 --trees 40
python reproducibility/scripts/nested_feature_selection_cv.py --data reproducibility/data/restored --out nested_results --model ace-ipred --max-features 10 --candidate-limit 12 --trees 40
```

The previously published `reproduce_selection.py` screens the entire training dataset before scoring its five-fold wrapper. **Its wrapper scores must not be interpreted as leakage-free cross-validated generalization estimates.** It is retained for historical transparency. The newer nested evaluation addresses this within-training issue prospectively, but does not establish that the historical independent test was never consulted during model development. The independent test remains untouched in nested CV.

## Historical independent-test provenance review (10 November 2026)

**Important:** [The original notebook and archive chronology audit](verification/2026-11-10_historical_holdout_provenance_audit.md) found exploratory test-label evaluations and multiple alternative feature subsets in the historical ScavenPred development notebook. The records do **not** establish that the final independent test was never accessed or consulted during the original development process. This does not by itself prove that test results determined the final model; nevertheless, the old holdout cannot be certified as prospectively untouched. The [cell-level evidence inventory](verification/2026-11-10_original_notebook_cell_inventory.csv) records the specific notebook locations and hashes.

The new nested cross-validation remains valid as a separate corrected assessment. Archived external-test results should be described as *legacy external evaluations* unless a contemporaneous, independently verified test-sealing history is available. Do not claim the historical test was never consulted on the basis of the current repository.

## Primary corrected evaluation: sequence-group-disjoint nested CV (10 November 2026)

The fully verified **primary leakage-controlled evaluation** is [strict_group_nested_validation.py](scripts/strict_group_nested_validation.py), with [complete metrics and fold results](verification/2026-11-10_final_group_disjoint_nested_cv.md). Unlike the superseded historical wrapper, the feature selector is inside the inner-CV `Pipeline`; model tuning and fitting are also restricted to the corresponding inner training folds. Related peptides (>=80% same-length identity, or full-sequence substring relations of length >=5) are placed in the same group across outer and inner splits.

| Model | N | Out-of-fold accuracy | MCC | ROC-AUC | PR-AUC |
|---|---:|---:|---:|---:|---:|
| ScavenPred | 1,648 | 70.33% | 0.403 | 0.782 | 0.774 |
| ACE-iPred | 1,933 | 75.06% | 0.498 | 0.819 | 0.820 |

```bash
python reproducibility/scripts/unpack_data.py
python reproducibility/scripts/strict_group_nested_validation.py --data reproducibility/data/csv --out group_nested_results --model both
```

All 3,581 out-of-fold predictions were audited for unique membership and numerical consistency; all outer held-out sequence groups were disjoint from their corresponding training partitions. Independent-test files are not loaded. **Use these new nested evaluation results for the corrected methodology**, rather than claiming the old exploratory notebook was free from leakage. The public web service continues to use the previously deployed frozen model.

## Final corrected external evaluation — 10 November 2026

The [full reproducible external-test script](scripts/reproduce_external_similarity_screened.py) selects model features/hyperparameters exclusively by cross-validation on the training data, then excludes closely related training analogues from the external test by sequence alone before applying the frozen selected estimator. The [independently checked results](verification/2026-11-10_final_external_similarity_screened.md) are:

| Model | Screened external n | Accuracy | MCC | ROC-AUC | PR-AUC |
|---|---:|---:|---:|---:|---:|
| ScavenPred | 250 | 81.20% | 0.526 | 0.860 | 0.955 |
| ACE-iPred | 248 | 75.00% | 0.461 | 0.843 | 0.929 |

The missing ACE-iPred test descriptor Z5 is reconstructed from the **training sequences and training descriptor values only**, with maximum training residual < 4e-15. No external labels are used for model selection, training, or descriptor recovery.

```bash
python reproducibility/scripts/unpack_data.py
python reproducibility/scripts/reproduce_external_similarity_screened.py --data reproducibility/data/csv --scavenpred-test reproducibility/data/scavenpred_test_88.csv --out external_results
```
