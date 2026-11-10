# ScavenPred and ACE-iPred

Open research code, curated peptide datasets, and reproducibility workflows for free-radical-scavenging (ScavenPred) and ACE-inhibitory (ACE-iPred) peptide activity prediction.

## Directly accessible research datasets

Both prediction tasks now have their **own visible CSV files** under [reproducibility/data](reproducibility/data). These are restored from the repository's checksum-verified source archive.

| Tool | Training descriptors | Training labels | Independent test descriptors |
|---|---|---|---|
| ScavenPred | [scavenpred_train_88.csv](reproducibility/data/scavenpred_train_88.csv) | [scavenpred_train_labels.csv](reproducibility/data/scavenpred_train_labels.csv) | [scavenpred_test_88.csv](reproducibility/data/scavenpred_test_88.csv) |
| ACE-iPred | [ace_train_88.csv](reproducibility/data/ace_train_88.csv) | [ace_train_labels.csv](reproducibility/data/ace_train_labels.csv) | [ace_test_88.csv](reproducibility/data/ace_test_88.csv) |

These are **curated numerical descriptor matrices**, not full raw source-protein retrieval histories. ACE-iPred's independent descriptor file has 253 rows before sequence-relatedness screening; the corrected external evaluation includes 248 peptides. ScavenPred's test file has 256 before screening and 250 after.

## Reproducibility

Start with [reproducibility/README.md](reproducibility/README.md). Consult the [results reconciliation](reproducibility/verification/RESULTS_RECONCILIATION.md) for the corrected grouped nested cross-validation and similarity-screened external evaluation, and use the linked scripts for executable reproduction. Historical experiments remain clearly identified and are not interchangeable with the corrected protocols.

## Public prediction tools

- ScavenPred: https://pepsilico.com/tools/scavenpred/
- ACE-iPred: https://pepsilico.com/tools/ace-ipred/

The live web predictors and corrected experimental validation scripts are distinct implementations. Model predictions are not experimental measurements.

## Research scope

SAPIEN-Pep and SIPPAD are separate research systems, and their protected code is not included.

## Citation

Please cite the associated manuscript once final publication details are available.
