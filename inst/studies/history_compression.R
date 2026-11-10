args <- commandArgs(trailingOnly = TRUE)
script_arg <- grep('^--file=',commandArgs(),value=TRUE)
script_dir <- if(length(script_arg)) dirname(normalizePath(sub('^--file=','',script_arg[1]))) else getwd()
root <- if(length(args)) normalizePath(args[1]) else dirname(script_dir)
out <- if(length(args)>1L) args[2] else file.path(root,'outputs','history_compression')
dir.create(out,recursive=TRUE,showWarnings=FALSE)
for(f in c('reference.R','decisions.R','certificates.R')) source(file.path(root,'R',f))
x <- jsonlite::fromJSON(file.path(root,'data/base/compact_certificate/INPUTS.json'),simplifyVector=FALSE)$records
rows <- lapply(x,function(z) {
  p <- prepare_certificate(z); a <- minimum_certificate(p)
  data.frame(id=z$id, bank=z$bank, case_id=z$case_id, history=paste(z$bank,z$case_id,sep=':'),
    policy=z$policy,stage=z$stage,mode=z$mode,certified=a$status=='CERTIFIED',
    full=length(p$available),minimum=a$cost,
    stringsAsFactors=FALSE,row.names=NULL)
})
d <- do.call(rbind,rows)
d$saved <- d$full-d$minimum
d$fraction_saved <- d$saved/d$full
stopifnot(nrow(d)==292L, length(unique(d$history))==32L,
          sum(d$certified)==223L, sum(d$saved,na.rm=TRUE)==743L,
          sum(d$saved>0,na.rm=TRUE)==218L,
          all(d$minimum[d$certified]<=d$full[d$certified]))
original <- jsonlite::fromJSON(file.path(root,'data/base/compact_certificate/RECORDS.json'),simplifyVector=FALSE)
stopifnot(length(original)==nrow(d))
for(i in seq_along(original)) {
  stopifnot(identical(d$certified[i],original[[i]]$status=='CERTIFIED'))
  if(d$certified[i])stopifnot(d$minimum[i]==original[[i]]$cost,
                             d$minimum[i]==original[[i]]$best_greedy_cost)
}
write.csv(d,file.path(out,'REPORT_COMPRESSION.csv'),row.names=FALSE)
hist_summary <- function(a,scope) {
  r <- lapply(split(a,a$history),function(g) {
    z <- g[g$certified,]
    data.frame(scope=scope,history=g$history[1],bank=g$bank[1],reports=nrow(g),
      certified=nrow(z),certified_fraction=mean(g$certified),
      full_median=if(nrow(z))median(z$full)else NA_real_,
      minimum_median=if(nrow(z))median(z$minimum)else NA_real_,
      mean_saved=if(nrow(z))mean(z$saved)else NA_real_,
      median_saved=if(nrow(z))median(z$saved)else NA_real_,
      mean_fraction_saved=if(nrow(z))mean(z$fraction_saved)else NA_real_,
      reduced_fraction=if(nrow(z))mean(z$saved>0)else NA_real_)
  })
  do.call(rbind,r)
}
h <- rbind(hist_summary(d,'all'),hist_summary(d[d$mode=='reference',],'reference'),
           hist_summary(d[d$mode!='identity',],'excluding_identity'))
write.csv(h,file.path(out,'HISTORY_COMPRESSION.csv'),row.names=FALSE)
overall <- function(a,group) {
  z <- a[a$certified,]; h <- hist_summary(a,group); hv <- h[!is.na(h$mean_saved),]
  data.frame(group=group,histories=length(unique(a$history)),reports=nrow(a),certified=nrow(z),
    unresolved=sum(!a$certified),histories_with_certificates=nrow(hv),
    minimum_median=median(z$minimum),full_median=median(z$full),
    minimum_min=min(z$minimum),minimum_max=max(z$minimum),
    saved_median=median(z$saved),fraction_saved_median=median(z$fraction_saved),
    histories_ever_reduced=sum(hv$reduced_fraction>0),
    equal_history_mean_saved=mean(hv$mean_saved),
    history_mean_saved_median=median(hv$mean_saved),
    history_mean_saved_min=min(hv$mean_saved),history_mean_saved_max=max(hv$mean_saved),
    equal_history_fraction_saved=mean(hv$mean_fraction_saved),
    history_fraction_saved_median=median(hv$mean_fraction_saved),
    history_fraction_saved_min=min(hv$mean_fraction_saved),history_fraction_saved_max=max(hv$mean_fraction_saved),
    total_saved=sum(z$saved),total_before=sum(z$full))
}
ss <- list(all=d,reference=d[d$mode=='reference',],excluding_identity=d[d$mode!='identity',],
  IA1=d[d$stage=='IA1',],Final=d[d$stage=='Final',])
for(b in unique(d$bank))ss[[b]]<-d[d$bank==b,]
for(b in unique(d$bank))for(s in unique(d$stage))ss[[paste(b,s,sep=':')]]<-d[d$bank==b & d$stage==s,]
summary <- do.call(rbind,lapply(names(ss),function(n)overall(ss[[n]],n)))
stopifnot(summary$histories_ever_reduced[summary$group=='all']==32L,
  summary$histories_with_certificates[summary$group=='reference']==25L,
  summary$histories_ever_reduced[summary$group=='reference']==25L,
  summary$histories_ever_reduced[summary$group=='excluding_identity']==32L)
write.csv(summary,file.path(out,'COMPRESSION_SUMMARY.csv'),row.names=FALSE)
print(summary,digits=6,row.names=FALSE)
cat('\nHistogram of minimum lengths by stage:\n');print(with(d[d$certified,],table(stage,minimum)))
cat('\nHistories with no reference certificate:\n');print(h[h$scope=='reference' & h$certified==0,c('history','reports')],row.names=FALSE)
cat('\nSavings by history, all reports:\n');print(h[h$scope=='all',c('history','reports','certified','mean_saved','mean_fraction_saved','reduced_fraction')],row.names=FALSE)
jsonlite::write_json(list(scope='Post-result descriptive summaries of fixed reports, with equal weighting of histories for history means; no population inference or confidence intervals.',
  summaries=summary,history=h,report_counts=as.data.frame(table(d$bank,d$mode))),file.path(out,'HISTORY_COMPRESSION.json'),pretty=TRUE,auto_unbox=TRUE,digits=17,na='null')
