## Evaluating Sequence and Structural Similarity Metrics for Predicting Shared Paralog Functions

This repository contains the data processing, dataset construction, analysis code, and feature evaluations for the paper titled: **Evaluating Sequence and Structural Similarity Metrics for Predicting Shared Paralog Functions**, currently *in preparation*.

The main datasets, along with all computed paralog pair features, can be found in the `./data` directory.



### Sequence similarity features computing overview:

Scripts in `./compute_features` are used to compute all sequence similarity features presented in the paper, for human and yeast paralog pairs, including metrics such as sequence identity, alignment scores, and other relevant features.

> [!WARNING]
> They require external softfware (e.g., Foldseek) and databases (e.g., AF-DB).

| Script                                 | Brief description                                        |
|:---------------------------------------|:---------------------------------------------------------|
| {species}_ProtT5_distances.ipynb       | Compute distances between ProtT5 embeddings (*per-protein*) extracted from UniProt |
| {species}_esm2_distances.py            | Compute ESM2 embeddings and distances between them |
| {species}_compare_pairs_seq.py         | Compute sequence similarity searches using MMseqs2 |
| {species}_compare_pairs_struct.py      | Compute predicted structure similarity searches using Foldseek |



### Main data processing and analysis notebooks overview:

Main notebooks used to build the four datasets (PPI-SL X Yeast-Human), analyze their contents, explore the relationships between the features, and evaluate the predictive power of individual and combined features on these datasets.

| Notebook                               | Figures                 | Brief description                                        |
|:---------------------------------------|:------------------------|:---------------------------------------------------------|
| 1a_human_prepare_datasets.ipynb        | Fig. 3B, 4B, 5B         | Analyse features relationships in all human paralog pairs & build PPI and SL dataset for human paralog pairs |
| 1b_yeast_prepare_datasets.ipynb        | Fig. S1                 | Analyse features relationships in all yeast paralog pairs & build PPI and SL dataset for yeast paralog pairs |
| 2_datasets_overlap.ipynb               | Fig. 2                  | Analyse overlaps between PPI-SL labels   |
| 3_evaluate_features.ipynb              | Fig. 3, 4, 5, 6, S2     | Evaluate individual and combined features on the four datasets   |

