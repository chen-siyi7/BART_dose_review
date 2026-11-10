

a<-commandArgs(trailingOnly=TRUE);stopifnot(length(a)%in%c(2L,3L))
script<-sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE));root<-normalizePath(file.path(dirname(script),'..'))
for(n in c('reference.R','audits.R'))source(file.path(root,'R',n))
stem<-a[2];stopifnot(!file.exists(paste0(stem,'.rds')),!file.exists(paste0(stem,'.json')))
input<-read_json(file.path(root,'data/base/compact_reference/bank_v2/INPUTS.json'))
u<-Filter(function(u) u$id==a[1],input$units);stopifnot(length(u)==1L);u<-u[[1L]]
c<-Filter(function(c) c$id==u$catalog_id,input$catalogs)[[1L]]
n<-unlist(c$pattern_counts);X<-do.call(rbind,lapply(c$X_unique,unlist));law<-catalog_law(X,c$cutpoints,n)
proposals<-if(length(a)==3L) as.integer(a[3]) else 1000000L;stopifnot(proposals>0L)
set.seed(2061109000L)
if(u$status=='IDENTITY_SKIPPED') q<-rep(u$identity,proposals) else {
  v<-sample_reference(proposals,2048L,50L,prepare_reference(law,nrow(X)),unlist(u$successes),n,n/sum(n),u$offset,'aggregated');q<-v$q
}
saveRDS(q,paste0(stem,'.rds'))
delta<-.01/(2*sum(vapply(input$units,function(v)v$bank==u$bank,TRUE)))
bands<-probability_intervals(q,u$endpoint,delta)
if(u$status=='IDENTITY_SKIPPED') {
  gate<-as.integer(if(u$endpoint=='T')u$identity<.45 else u$identity>.3)
  bands$mean_lower<-bands$mean_upper<-u$identity;bands$gate_lower<-bands$gate_upper<-gate
}
save_json(c(list(unit=u$id,language='R',proposals=proposals,seed=2061109000L,scalar_error=delta),bands),paste0(stem,'.json'))
