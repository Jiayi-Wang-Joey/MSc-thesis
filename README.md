# MSc-thesis on "Bias correction and differential motif activity in ATAC-seq data"

## setup
- R version and library have to be specified in the `config.yaml` file  
  (e.g., `R: "R_LIBS_USER=/path/to/library /path/to/R/executable"`)

## Content 
### bulkATAC-seq: Evaluation of the weight and insertion model performance on the 10 experimental bulk ATAC-seq datasets
**Weight model**
- `00-get_dat.R`: read the ATAC-seq fragments in BED/BAM format and save it as a list of fragments for each dataset
- `01-dat_ttl.R` and `get_wgt.R`: get the original (`ttl`) or weighted (`wgt`, either only fragment weighted or fragment+peak weighted) peak-level accessibility counts in `SummarizedExperiment` object
- `02-dif.R`:
  - Differential accessible motif activity analysis using either `chromVAR-limma` or `limma-voom` 
  - **in:** `02-dif-<dif>.R`<br>: `dif`=`chromVAR-limma` or `limma-voom`
- `03-pkd.R`: Peak-level accessiblity test using `limma-voom`
**Insertion model**
  
