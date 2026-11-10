
args<-commandArgs(trailingOnly=TRUE);stopifnot(length(args)%in%c(1L,2L))
script<-sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE));root<-normalizePath(dirname(script))
source(file.path(root,'analysis_helpers.R'))
out<-if(length(args)==2L)args[2]else file.path(root,'inputs')
ast_prepare_inputs(args[1],out)
cat('Participant joins, endpoints, missing-outcome order and prediction weights checked; no sampling\n')
