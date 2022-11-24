# reanalysis_scRNA_seq_benchmark 
[![DOI](https://zenodo.org/badge/451781506.svg)](https://zenodo.org/badge/latestdoi/451781506)

A reanalysis the results of the paper [A practical solution to pseudoreplication bias in single-cell studies](https://www.nature.com/articles/s41467-021-21038-1)

Recently Zimmerman et al. proposed the use of mixed models over pseudobulk aggregation approaches, reporting improved performance on a novel simulation approach of hierarchical single-cell expression data. However, their reported results could not prove the superiority of mixed models as they are based on separate calculations of type 1 (performance of the models on non-differentially expressed genes) and type 2 error (performance on differentially expressed genes). To correctly benchmark the models, a reanalysis using a balanced measure of performance, considering both the type 1 and type 2 errors (both the differentially and non-differentially expressed genes), is necessary.

This analysis was conducted using a modified version of hierarchicell which returns the Matthews correlation coefficient performance metric as well as the type 1 error rates, uses the same simulated data across approaches and has checkpointing capabilities (so runs can continue from where they left off if aborted or crashed) is available at: https://github.com/neurogenomics/hierarchicell.

## Results
Click [here](https://al-murphy.github.io/reanalysis_scRNA_seq_benchmark/analysis.html) 
to view the results of the benchmarking analysis.

Our results demonstrate that pseudobulk approaches are far from being "too conservative" and are, in fact, the best performing models based on this simulated dataset for the analysis of single-cell expression data.
