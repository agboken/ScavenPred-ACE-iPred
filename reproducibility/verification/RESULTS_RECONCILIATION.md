# Reconciliation of ScavenPred and ACE-iPred evaluations

**Primary corrected experiments: grouped nested CV with the original Welch/Pearson/VIF/forward-wrapper feature-selection method and its matching similarity-screened external evaluation.** These results supersede the previous `SelectKBest` evaluation for the revised wrapper-based manuscript. Always pair metrics with the correct model and cohort.

| Evaluation protocol | ScavenPred n | Accuracy | MCC | ROC-AUC | PR-AUC | ACE-iPred n | Accuracy | MCC | ROC-AUC | PR-AUC |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| **Primary wrapper grouped nested CV** | **1,648** | **70.15%** | **0.399791** | **0.770657** | **0.763040** | **1,933** | **73.98%** | **0.477642** | **0.811789** | **0.818212** |
| **Primary wrapper similarity-screened external** | **250** | **78.40%** | **0.438924** | **0.832389** | **0.940418** | **248** | **76.21%** | **0.491399** | **0.843366** | **0.921615** |
| Previous grouped nested CV with `SelectKBest` | 1,648 | 70.33% | 0.403345 | 0.781950 | 0.774145 | 1,933 | 75.06% | 0.498451 | 0.819098 | 0.820272 |
| Previous screened external with `SelectKBest` | 250 | 81.20% | 0.526038 | 0.860355 | 0.954734 | 248 | 75.00% | 0.461364 | 0.842505 | 0.929486 |
| Historical frozen RF independent test | 256 | 76.56% | 0.377006 | 0.823609 | — | 253 | 72.73% | 0.416194 | 0.785218 | — |
| Early non-grouped nested CV audit | 1,648 | 70.02% | 0.397091 | 0.765949 | — | 1,933 | 72.37% | 0.444549 | 0.795174 | — |

## MCC verification from primary confusion matrices

| Evaluation | ScavenPred [TN,FP;FN,TP] | ACE-iPred [TN,FP;FN,TP] |
|---|---|---|
| Wrapper grouped CV | [675,196;296,481] | [655,244;259,775] |
| Wrapper screened external | [36,18;36,160] | [56,17;42,133] |

- Grouped ScavenPred accuracy = (675+481)/1648 = 0.70145631
- Grouped ACE-iPred accuracy = (655+775)/1933 = 0.73978272
- External ScavenPred accuracy = (36+160)/250 = 0.784
- External ACE-iPred accuracy = (56+133)/248 = 0.76209677

All primary MCCs, ROC-AUCs, and PR-AUCs were recomputed from the corresponding prediction-level outputs; the outer MCC is calculated from **aggregated outer-fold confusion counts**, not an arithmetic mean of fold MCCs.

## Authoritative reproducibility references

- [Grouped nested forward-wrapper CV](2026-11-10_grouped_nested_forward_wrapper.md)
- [Matching forward-wrapper external evaluation](2026-11-10_grouped_forward_wrapper_external.md)
- [Nested wrapper source code](../scripts/grouped_nested_forward_wrapper_cv.py)
- [External wrapper source code](../scripts/forward_wrapper_screened_external.py)
- [Automated verification](../../.github/workflows/grouped-forward-wrapper.yml)

## Why historical numbers differ

The old frozen RF models, the corrected `SelectKBest` models and the new forward-wrapper models have **different selected descriptors and/or model parameters and datasets**. The `SelectKBest` models are preserved as prior alternative research analyses, not as manuscript results for a claim about the forward wrapper.

The new forward-wrapper selected final external descriptor sets are:
- **ScavenPred:** ProtFP3, BLOSUM6, BLOSUM3, E5, ProtFP2
- **ACE-iPred:** ProtFP3, BLOSUM1, BLOSUM6, BLOSUM4, BLOSUM5

The public PepSilico predictor interfaces were not automatically redeployed with these research models. No historical evaluation is evidence that the corrected wrapper-based protocol was used at an earlier date.

**For revised manuscript tables, always use the first two rows and the documented wrapper pipeline, not the superseded model configurations.**
