# TF binding motifs in Shank3 region

#### input files:
`shank3_allfeatures.txt`, and `shank3_up2K.txt` are BED like files that provide the coordinates of SHANK3 gene. These files are taken from GENCODE V31 HG38 annotations with custom scripts. Alternatively, you can manually produce these annotation files.

#### output file:
`shank3_region_TF_motifs.txt` all matching motifs with motifs' coordinates in the SHANK3 region. Subsequent columns indicate the binding motif's tf names, symbols(if available), descriptions(if available), and species.

#### R script:
`shank3.R` is the script that takes previously mentioned inputs and produce the output as noted above. Note a number of Bioconductor packages are required. Rscript was tested under R/3.6.1.

Run `Rscript shank3.R` in bash or load the R script file in Rstudio. 