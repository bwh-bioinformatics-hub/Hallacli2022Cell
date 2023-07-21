# Hallacli2022Cell
Bioinformatics code repository for Hallacli et al. Cell 2022

## ./src/main.sh
The primary bash script includes multiple steps of the pipeline:
- process each sample (e.g. QC, mapping, and gene expression quantification), see _RNAseq.lite.sh
- merge all samples to generate gene expression matrix,
- perform sample level QC / outlier detection, see _normQC.R
- DE analysis, see _DE.R and _DE2.R

## ./src/main.R
The R script to concatenate multi-timepoint consecutive analysis result together. 
