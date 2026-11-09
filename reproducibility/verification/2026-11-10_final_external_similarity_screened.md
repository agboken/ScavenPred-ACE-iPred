# Final external evaluation: training-only feature selection and similarity screening

## Procedure

The [reproducible evaluation script](../scripts/reproduce_external_similarity_screened.py) loads the corrected training matrices, groups related **training** peptides using an unsupervised sequence-based rule, then uses three-fold StratifiedGroupKFold with a SelectKBest + Random Forest Pipeline to select features and hyperparameters. The external-test labels are not supplied to model fitting, feature screening or hyperparameter search.

Before scoring, external peptides are excluded when (i) they share >=80% residue identity with an equal-length training sequence of length >=5, or (ii) either sequence is an exact substring of the other of length >=5. The exclusion uses peptide sequences, not activity labels. No training peptides are added or deleted as a consequence of external performance.

The ACE test matrix contained 87 of the 88 training descriptors (missing Z5). Z5 is a mean amino-acid descriptor. It was reconstructed using a linear amino-acid frequency model fitted **only to ACE training sequences and their Z5 values**, with maximal training reconstruction residual **3.8e-15**. External test activity labels were not used for the reconstruction.

## Reproduced results

| Metric | ScavenPred | ACE-iPred |
|---|---:|---:|
| Training samples | 1,648 | 1,933 |
| External initial sample size | 256 | 253 |
| Excluded by sequence relatedness | 6 | 5 |
| Evaluated external peptides | **250** | **248** |
| Accuracy | **81.20%** | **75.00%** |
| MCC | **0.526038** | **0.461364** |
| ROC-AUC | **0.860355** | **0.842505** |
| PR-AUC (average precision) | **0.954734** | **0.929486** |
| TN / FP / FN / TP | 41 / 13 / 34 / 162 | 54 / 19 / 43 / 132 |

**ScavenPred fitted descriptors:** BLOSUM6, F2, E5, ProtFP3, ProtFP5, SV1, ST5, SVGER2, SVGER8, T2. Depth 10, leaf 4, 100 trees.

**ACE-iPred fitted descriptors:** BLOSUM7, KF4, E4, ProtFP3, ProtFP5, ST3, ST5, SVGER8, VHSE1, Z1. Depth 7, leaf 1, 100 trees.

### Reproduce

```bash
python reproducibility/scripts/unpack_data.py
python reproducibility/scripts/reproduce_external_similarity_screened.py --data reproducibility/data/csv --scavenpred-test reproducibility/data/scavenpred_test_88.csv --out external_results
```

The scripts regenerate per-peptide external results and a machine-readable metric report. This evaluation uses the revised, fold-contained feature-selection pipeline; the previously deployed web models are not modified by this verification.
