# Similarity-screened external validation using the forward-wrapper-selected models

**Date:** 10 November 2026. This evaluation is paired with the [grouped nested forward-wrapper CV](2026-11-10_grouped_nested_forward_wrapper.md), and **not** with the separate `SelectKBest` experiments.

The complete curated training data were subjected to Welch screening (p<0.05), class-signal-ordered Pearson correlation filtering (|r|<0.70), VIF<10, and forward RF wrapper based on three-fold sequence-group-disjoint training CV. Candidate limit 12, RF 40 trees during wrapper selection, and a fixed final RF with 1,000 trees (ScavenPred) or 200 trees (ACE-iPred). The test sets were filtered solely from peptide sequences to exclude same-length >=80% identity for length >=5 and full-sequence substring/superstring matches length >=5. No test activity labels entered descriptor selection or model fitting.

| Metric | ScavenPred | ACE-iPred |
|---|---:|---:|
| Training peptides | 1,648 | 1,933 |
| Starting external cohort | 256 | 253 |
| Excluded as close to training | 6 | 5 |
| External peptides evaluated | **250** | **248** |
| Accuracy | **0.784000** | **0.762096774194** |
| MCC | **0.438924355585** | **0.491398685358** |
| ROC-AUC | **0.832388510960** | **0.843365949119** |
| PR-AUC | **0.940417591016** | **0.921615205486** |
| Confusion matrix [TN,FP;FN,TP] | **[36,18;36,160]** | **[56,17;42,133]** |

The aggregate MCC values were independently recomputed from the saved prediction-level outputs.

## Final full-training selected descriptors

**ScavenPred (5):** ProtFP3, BLOSUM6, BLOSUM3, E5, ProtFP2.

**ACE-iPred (5):** ProtFP3, BLOSUM1, BLOSUM6, BLOSUM4, BLOSUM5.

The selected features differ from those in the earlier frozen deployed RF models and from the newer `SelectKBest` external evaluation. The selection method—not the predictive application—has been consistently matched between this external run and the new grouped nested wrapper CV.

[Executable external code](../scripts/forward_wrapper_screened_external.py) · [Grouped forward-wrapper CV](2026-11-10_grouped_nested_forward_wrapper.md) · [Combined GitHub Actions verification](../../.github/workflows/grouped-forward-wrapper.yml)

```bash
python reproducibility/scripts/unpack_data.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
python reproducibility/scripts/forward_wrapper_screened_external.py \
  --data reproducibility/data/csv --out grouped_wrapper_results \
  --model both --inner-trees 40 --candidate-limit 12
```

These external cohorts are the existing curated tests, not prospectively gathered new cohorts. The earlier `SelectKBest` external results remain archival and must not be paired with this wrapper-based CV in any manuscript table.
