# ScavenPred and ACE-iPred — Reproducibility

**Primary corrected method (10 November 2026): grouped nested five-fold CV retaining the original Welch / Pearson / VIF / forward-wrapper selection.** The corresponding external evaluation uses the same selection method, fitted exclusively on training data, with close training relatives excluded from the external test by sequence criteria.

**Automated verification passed:** [GitHub Actions run 42731088361](https://github.com/agboken/ScavenPred-ACE-iPred/actions/runs/42731088361) independently regenerated all five grouped outer folds for both models, both similarity-screened external predictions, the four reported MCC values, and the saved confusion matrices. The run also provides a downloadable prediction-level artifact.

## Complete sequential research notebooks

Both notebooks contain **inline executable Python code**, not merely calls to older notebooks: checksum-controlled dataset restoration; 88-descriptor loading and activity-label checks; class/length/distribution plots; sequence-relatedness grouping and disjoint-fold audits; train-only Welch, Pearson, VIF screening; 12-candidate, three-fold inner MCC forward-wrapper search; five-fold outer model evaluation; held-out confusion matrices, ROC/PR curves and descriptor-selection frequency; full-training wrapper selection; similarity-screened external evaluation; per-sequence MCC checks; Random Forest descriptor importance and optional SHAP interpretation.

- [ScavenPred: full sequential computational notebook](notebooks/ScavenPred_complete_training_and_analysis.ipynb)
- [ACE-iPred: full sequential computational notebook](notebooks/ACE_iPred_complete_training_and_analysis.ipynb)

**Upstream provenance boundary:** The curated 88-descriptor matrices are provided, but the complete original UniProt protein collection, fragmentation procedure and experimentally sourced activity records are not reconstructed from these files. Thus the notebooks reproduce all available **current modelling, validation, visualisation and analysis** stages, not unprovided raw source collection.

The verified scripts in `scripts/` remain the canonical reference implementation; the notebooks embed their functions to document and execute each stage in sequence. The [lightweight notebook preflight](../../.github/workflows/current-notebooks-preflight.yml) checks syntax, data restoration, grouping and preliminary feature screening; the full model runs remain covered by the [grouped-wrapper verification workflow](../../.github/workflows/grouped-forward-wrapper.yml).

## Primary corrected results

| Evaluation | ScavenPred | ACE-iPred |
|---|---:|---:|
| Training peptides | 1,648 | 1,933 |
| Grouped nested CV accuracy | **70.15%** | **73.98%** |
| Grouped nested CV MCC | **0.399791** | **0.477642** |
| Grouped nested CV ROC-AUC | **0.770657** | **0.811789** |
| Grouped nested CV PR-AUC | **0.763040** | **0.818212** |
| Screened external peptides | **250** | **248** |
| Screened external accuracy | **78.40%** | **76.21%** |
| Screened external MCC | **0.438924** | **0.491399** |
| Screened external ROC-AUC | **0.832389** | **0.843366** |
| Screened external PR-AUC | **0.940418** | **0.921615** |

**Confusion matrices**, in TN / FP / FN / TP order: nested CV ScavenPred **675 / 196 / 296 / 481**, ACE-iPred **655 / 244 / 259 / 775**; external ScavenPred **36 / 18 / 36 / 160**, ACE-iPred **56 / 17 / 42 / 133**.

- [Grouped nested CV report and fold-selected features](verification/2026-11-10_grouped_nested_forward_wrapper.md)
- [Matching forward-wrapper external test report](verification/2026-11-10_grouped_forward_wrapper_external.md)
- [Group-disjoint wrapper verification workflow](../../actions/workflows/grouped-forward-wrapper.yml)
- [Reproducible grouped nested forward-wrapper script](scripts/grouped_nested_forward_wrapper_cv.py)
- [Reproducible forward-wrapper external script](scripts/forward_wrapper_screened_external.py)

## Reproduce from a fresh clone

```bash
git clone https://github.com/agboken/ScavenPred-ACE-iPred.git
cd ScavenPred-ACE-iPred
python -m pip install -r reproducibility/requirements.txt
python reproducibility/scripts/unpack_data.py

OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
python reproducibility/scripts/grouped_nested_forward_wrapper_cv.py \
  --data reproducibility/data/csv \
  --out grouped_wrapper_results --model both \
  --inner-trees 40 --candidate-limit 12

OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
python reproducibility/scripts/forward_wrapper_screened_external.py \
  --data reproducibility/data/csv \
  --out grouped_wrapper_results --model both \
  --inner-trees 40 --candidate-limit 12
```

All training and test descriptor CSVs are reconstructed and checksum-verified from the public archive. **scikit-learn 1.8.0** is pinned in `requirements.txt`. The public [GitHub Actions workflow](../../actions/workflows/grouped-forward-wrapper.yml) reruns and checks the metrics and provides full prediction files as downloadable run artifacts.

## Method and leakage controls

- Training data are grouped from peptide sequences without using activity labels: for length ≥5, equal-length sequences sharing ≥80% positional identity, or substring/superstring relations of length ≥5, are kept in one group. Relatedness is transitive.
- Each of five outer groups of folds is held out from **all** descriptor screening and model selection. Only outer training data are supplied to the Welch (`p<0.05`), Pearson (`|r|<0.70`), VIF (`<10`), and forward-wrapper stages.
- The forward wrapper adds positions/descriptors according to **inner grouped three-fold mean MCC**, with at most 12 candidate features, at most five selected ScavenPred features or ten ACE-iPred features. Inner wrapper scores are selection criteria; they are not unbiased generalization estimates.
- The final model in each outer iteration is trained on that outer training fold, then predicts the untouched outer validation group. Reported cross-validation scores come only from concatenated outer out-of-fold predictions.
- For the separate external assessment, the same forward-wrapper method selects features on the **full training dataset**, then fits the final RF; test activity labels never enter feature selection, tuning, or fitting. Exactly 6 ScavenPred and 5 ACE-iPred external peptides are excluded by sequence-based similarity screening.
- The CV and external evaluations correspond to the same feature-selection **procedure**, although selected feature subsets vary between outer folds. The full-training final models are different from previously deployed web predictors.

## Final fitted external models

ScavenPred's full-training wrapper selected **ProtFP3, BLOSUM6, BLOSUM3, E5, ProtFP2**.

ACE-iPred's full-training wrapper selected **ProtFP3, BLOSUM1, BLOSUM6, BLOSUM4, BLOSUM5**.

Both final models use the original fixed model-specific Random Forest settings, with forward selection using 40-tree inner evaluation and final fitting using 1,000 or 200 trees, respectively.

## Datasets

The 23 base64 parts under [`data/archive_parts`](data/archive_parts/) reconstruct the checksum-verified public dataset bundle (SHA-256 `5361826554c55de3f0f05b725f6c86eaa0d1c701a985728b28366fc12d6fb6bc`). Its CSV files include ScavenPred training `scavenpred_train_88.csv` (1,648 sequences), ACE-iPred training `ace_train_88.csv` (1,933 sequences), ACE training labels, test descriptor matrices and feature screening/provenance tables. All training and evaluation matrices contain 88 descriptor families before selection (the ACE external matrix originally had one missing descriptor, which is not selected in the forward-wrapper models).

The `data/` directory retains archived supplementary records and a direct ScavenPred test-feature CSV. No proprietary SAPIEN-Pep or SIPPAD engine files are distributed.

## Historical analyses retained for transparency

The following are **not** the primary wrapper-based manuscript results:

- [Previously published grouped `SelectKBest` CV](verification/2026-11-10_final_group_disjoint_nested_cv.md): ScavenPred 70.33%, ACE-iPred 75.06%.
- [Previously published `SelectKBest` external evaluation](verification/2026-11-10_final_external_similarity_screened.md): ScavenPred 81.20%, ACE-iPred 75.00%.
- [Original frozen-model independent tests](verification/2026-11-10_independent_test_rerun.json): ScavenPred 76.56% (256), ACE-iPred 72.73% (253). These belong to different frozen models and the historical [frozen-model CI](../../actions/workflows/reproducibility.yml).
- [Exploratory notebook](notebooks/archive/Scavenpred_original_archived.ipynb) and the [historical holdout provenance audit](verification/2026-11-10_historical_holdout_provenance_audit.md).

Do **not** mix metrics, feature subsets or cohort sizes between evaluation protocols. In particular, the live [ScavenPred](https://pepsilico.com/tools/scavenpred/) and [ACE-iPred](https://pepsilico.com/tools/ace-ipred/) predictor endpoints have not been replaced with the new wrapper-selected models by this research analysis.

## Scientific scope

The negative cohorts include *putative* negatives; absence of a motif does not establish experimental inactivity. The models predict free-radical-scavenging or ACE-inhibitory peptide activities, not experimentally observed values. Cite the associated manuscript when final publication information is available.
