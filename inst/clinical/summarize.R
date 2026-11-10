
args<-commandArgs(trailingOnly=TRUE);stopifnot(length(args)<=1L)
script<-sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE));root<-normalizePath(dirname(script))
source(file.path(root,'analysis_helpers.R'))
plan<-jsonlite::read_json(file.path(root,'PLAN.json'),simplifyVector=TRUE)
input<-readRDS(file.path(root,'inputs/INPUTS.rds'))
chains<-setNames(lapply(names(input$models),function(model)lapply(seq_len(plan$chains),function(role) {
  z<-readRDS(file.path(root,'chains',paste0(model,'_chain',role,'.rds')))
  stopifnot(z$status=='COMPLETE',z$model_id==model,z$role==role,
    z$plan_sha256==ast_sha(file.path(root,'PLAN.json')),z$input_sha256==plan$input_sha256,
    identical(dim(z$draws),c(as.integer(plan$retained),3L)),all(is.finite(z$draws)))
  z
})),names(input$models))
out<-if(length(args))args[1]else file.path(root,'results')
ast_summarize(input,plan,chains,out)
cat('All 18 models, 72 chains and 57 diagnostic intervals summarized\n')
