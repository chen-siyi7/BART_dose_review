args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=2L)stop('Usage: Rscript replay.R REPRODUCIBILITY_DIRECTORY NEW_OUTPUT_DIRECTORY')
library(bartDoseReview)
replay_results(args[1L],args[2L])
