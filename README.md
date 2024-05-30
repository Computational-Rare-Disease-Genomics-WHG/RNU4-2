# RNU4-2
 Code relating to the RNU4-2 paper

## SNV_depletion
Contains code to calculate regional depletion of SNVs within snRNA genes (as per Figure 4). 
- calculate_depletion.R contains code to tally SNVs per window and calculate the proportion of observed SNVs.
- run_rscript.sh does some processing of the input data (rearranging columns etc) for input to the calculate_depletion.R script. This includes retrieving variants from UK Biobank within snRNA regions using bcftools and formatting the files for input into the R script.
- plotting_ouptut.Rmd contains code used to create plots in the paper, and to create supplementary table 8, using output of the run_rscript.sh script. 