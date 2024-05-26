# MSc-thesis on "Bias correction and differential motif activity in ATAC-seq data"

## setup
- R version and library have to be specified in the `config.yaml` file  
  (e.g., `R: "R_LIBS_USER=/path/to/library /path/to/R/executable"`)

## Content 
### bulkATAC-seq: Evaluation of the weight and/or insertion model performance on the 10 experimental bulk ATAC-seq datasets and Cusanovich dataset
**Weight model**
- `code/`:
  - `00-get_dat.R`: read the ATAC-seq fragments in BED/BAM format and save it as a list of fragments for samples
  - `01-dat_ttl/wgt.R`: get the original (`ttl`) or weighted (`wgt`, either only fragment weighted or fragment+peak weighted) peak-level accessibility counts in `SummarizedExperiment` object
  - `02-dif.R`: differential accessible motif activity analysis using either `chromVAR-limma` or `limma-voom` 
  - `03-pkd.R`: differential accessiblity test on peak levels using `limma-voom`
  
  **Insertion model**
  - `04-mms.R`: Calculate the motif scores and background scores for each dataset
  - `05-mmd.R`: Using the motif scores and background scores to calculate the z-score and apply `limma` or directly apply `limma-voom` on the motif scores.

### scATAC-seq: Evaluation of the fragment weight correction on the 6 scATAC-seq datasets
- `code/`:
  - `00-frg.R`: read the ATAC-seq fragments in BED/BAM format and save it as a list of fragments for barcodes/cells
  - `01-get_ttl/wgt.R`: get the original (`ttl`) or weighted (`wgt`, only fragment weighted) peak-level accessibility counts in `SummarizedExperiment` object
  - `02-fil.R`: standard quality control and dimensionality reduction in scATAC-seq
  - `03-dif.R`: differential accessible motif analysis using `limma-voom`
  - `04-pkd.R`: differential accessible analysis on peak levels using `limma-voom`
  - `05-sta.R`: evaluation statistics on the peak-level accessibility counts, e.g. coefficient of variation (CV), difference ratio
 
### loess-comparison: Investigate the subsampling and loess parameters and the weight model on semi-simulated CTCF/MAZ-FL/GC datasets
- `code/`:
  - `00-sim_frg.R`:  
  - `01-weight.R`: 

