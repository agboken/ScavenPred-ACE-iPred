# Grouped nested CV retaining Welch, Pearson, VIF and RF forward-wrapper selection

**Evaluation date:** 10 November 2026. **Protocol:** five outer, three inner `StratifiedGroupKFold` folds, fixed seed 42, unsupervised grouping of peptide sequences using the existing sequence-relationship grouping function. Feature selection runs anew inside every **outer training partition**: Welch p<0.05; Pearson |r|<0.70; VIF<10; at most 12 ranked candidate descriptors; forward selection by inner grouped-CV mean MCC. Maximum selected descriptors: 5 for ScavenPred, 10 for ACE-iPred. Wrapper scoring uses 40-tree RF, and final outer training fits use 1,000 trees (ScavenPred) or 200 trees (ACE-iPred) with original fixed model-specific parameters. All outer validation groups are untouched until scoring. External test data are never read.

**Clarification:** The preliminary screening is fitted to all records in an **outer training** partition before the inner wrapper search. Therefore inner wrapper MCC is a model-selection score, not an unbiased inner generalization estimate. The externally reported *outer out-of-fold metrics* remain leakage-controlled because none of the outer validation examples participate in any selection stage. This is a direct nested assessment of the original selection sequence, not a claim that screening was separately re-fitted within each inner wrapper fit.

| Metric | ScavenPred | ACE-iPred |
|---|---:|---:|
| Training peptides | 1,648 | 1,933 |
| Related-sequence groups | 1,580 | 1,718 |
| Accuracy (outer OOF) | **0.701456310680** | **0.739782721159** |
| MCC (outer OOF) | **0.399790720859** | **0.477641534351** |
| ROC-AUC | **0.770656666179** | **0.811789049944** |
| PR-AUC / average precision | **0.763040290543** | **0.818211575943** |
| Confusion matrix [TN,FP;FN,TP] | **[675,196;296,481]** | **[655,244;259,775]** |

The MCC values above were separately recomputed from the combined outer predictions, **not averaged from the five fold MCCs**. Every training peptide appears once in the saved prediction records, and all five outer validation groups are disjoint from their corresponding outer training groups.

### Outer-fold selected features

| Fold | ScavenPred features | ACE-iPred features |
|---|---|---|
| 1 | BLOSUM6, BLOSUM3, E5, VHSE1, VHSE4 | ProtFP3, MSWHIM2, ProtFP5, BLOSUM6, SVGER5, BLOSUM4, VHSE6, T3 |
| 2 | SVGER2, BLOSUM6, BLOSUM3, E5, Z4 | ProtFP3, E4, BLOSUM1, BLOSUM6, BLOSUM4, SVGER2, T5 |
| 3 | E5, BLOSUM6, BLOSUM3, ST5, SVGER2 | ST3, BLOSUM6, VHSE6, MSWHIM2, ProtFP3, ProtFP5, T3, SVGER5, BLOSUM4 |
| 4 | BLOSUM6, E5, SVGER2, Z4, Z3 | E4, SVGER2, BLOSUM6, T5, BLOSUM4 |
| 5 | E5, BLOSUM6, BLOSUM7, VHSE4, Z4 | ProtFP3, ProtFP5, MSWHIM2, BLOSUM5, KF5 |

## Comparison to prior grouped nested SelectKBest evaluation

| Model | New forward-wrapper accuracy | Previous SelectKBest accuracy | New forward-wrapper MCC | Previous SelectKBest MCC |
|---|---:|---:|---:|---:|
| ScavenPred | **70.15%** | 70.33% | **0.400** | 0.403 |
| ACE-iPred | **73.98%** | 75.06% | **0.478** | 0.498 |

The two procedures are different algorithms. Do **not** combine these new grouped-CV values with external evaluation of a model whose features were chosen by `SelectKBest`. A wrapper-consistent external evaluation must be performed and independently verified before the manuscript is switched to the wrapper-based grouped protocol.

## Reproduce

```bash
python reproducibility/scripts/unpack_data.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
python reproducibility/scripts/grouped_nested_forward_wrapper_cv.py \
  --data reproducibility/data/csv --out grouped_wrapper_results \
  --model both --inner-trees 40 --candidate-limit 12
```

[Executable grouped wrapper](../scripts/grouped_nested_forward_wrapper_cv.py) · [GitHub Actions verification](../../.github/workflows/grouped-forward-wrapper.yml) (workflow publishes full out-of-fold predictions and per-fold reports as an artifact).

All model selection occurred inside outer training data; this evaluation is not a retrospective statement about the old notebook.
