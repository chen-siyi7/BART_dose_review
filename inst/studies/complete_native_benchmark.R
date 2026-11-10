

a<-commandArgs(trailingOnly=TRUE);stopifnot(length(a)==1L)
script<-sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE));root<-normalizePath(file.path(dirname(script),'..'))
out<-normalizePath(a[1]);stopifnot(!file.exists(file.path(out,'RESULTS.json')))
for(n in c('reference.R','decisions.R','audits.R'))source(file.path(root,'R',n))
plan<-read_json(file.path(out,'PLAN.json'));work<-read_json(file.path(root,'data/worked_review/INPUTS.json'))
prep_times<-unlist(read_json(file.path(out,'PREPROCESSING.json')));rows<-list()
for(u in work$units) for(method in c('direct','aggregated')) rows[[length(rows)+1L]]<-read_json(file.path(out,paste0(method,'_',u$endpoint,u$dose,'.json')))
for(r in rows)stopifnot(file.exists(file.path(out,paste0(r$method,'_',r$endpoint,r$dose,'.bin'))))
methods<-lapply(c('direct','aggregated'),function(method) {
  rr<-Filter(function(x)x$method==method,rows)
  t0<-proc.time()[3];actions<-setNames(lapply(c('IA1','Final'),function(s)work_actions(rr,s)),c('IA1','Final'));at<-unname(proc.time()[3]-t0)
  setup<-sum(prep_times[unique(vapply(work$units,function(u)u$catalog_id,''))])
  list(method=method,seconds=setup+sum(vapply(rr,function(r)r$seconds,0))+at,catalog_seconds=setup,actions=actions,
    accepted_min=min(vapply(rr,function(r)r$accepted,0)),accepted_max=max(vapply(rr,function(r)r$accepted,0)))
})
lib<-Sys.getenv('BART_LIBRARY',unset='');if(nzchar(lib)).libPaths(c(normalizePath(lib),.libPaths()))
stopifnot(as.character(packageVersion('BART'))=='2.9.9')
save_json(list(status='FROZEN_BEFORE_NATIVE_MEASUREMENT',original_plan_sha256=digest::digest(file=file.path(out,'PLAN.json'),algo='sha256'),
  reason='Initial native process stopped at version preflight: system BART 2.9.10. Use existing BART 2.9.9 library for the specified fits.',
  reference_measurements='Complete and retained without rerun',BART=as.character(packageVersion('BART')),
  script_sha256=digest::digest(file=script,algo='sha256')),file.path(out,'NATIVE_CONTINUATION_PLAN.json'))
native<-file.path(out,'native_2_9_9');stopifnot(!dir.exists(native));dir.create(native)
rs<-file.path(R.home('bin'),'Rscript');t0<-proc.time()[3];ct<-numeric(4)
for(role in 1:4) {
  t<-proc.time()[3];status<-system2(rs,c(shQuote(file.path(root,'scripts/fit_chain.R')),shQuote(file.path(root,'data/worked_review/case03.rds')),
    shQuote(file.path(native,paste0('chain',role,'.rds'))),role),stdout=file.path(native,paste0('chain',role,'.log')),stderr=file.path(native,paste0('chain',role,'.log')))
  stopifnot(status==0L);ct[role]<-unname(proc.time()[3]-t);cat('Chain',role,ct[role],'seconds\n');flush.console()
}
status<-system2(rs,c(shQuote(file.path(root,'scripts/diagnose.R')),shQuote(native),shQuote(file.path(native,'DIAGNOSTICS.csv'))),
  stdout=file.path(native,'diagnostics.log'),stderr=file.path(native,'diagnostics.log'));stopifnot(status==0L)
dr<-read.csv(file.path(native,'DIAGNOSTICS.csv'));nr<-list()
for(j in 1:5)for(ep in c('T','E')) {
  mn<-dr[dr$dose==j&dr$series==paste0(ep,'_mean'),];g<-dr[dr$dose==j&dr$series==paste0(ep,'_gate'),]
  nr[[length(nr)+1L]]<-list(dose=j,endpoint=ep,mean=mn$mean,mean_lower=mn$lower,mean_upper=mn$upper,gate_lower=g$lower,gate_upper=g$upper)
}
na<-setNames(lapply(c('IA1','Final'),function(s)work_actions(nr,s)),c('IA1','Final'))
methods[[3L]]<-list(method='native_four_chains',seconds=unname(proc.time()[3]-t0),chain_process_seconds=ct,intervals=nr,actions=na,qualified_quantities=sum(dr$qualified))
save_json(list(status='COMPLETE',language='R',plan_sha256=digest::digest(file=file.path(out,'PLAN.json'),algo='sha256'),methods=methods,endpoint_results=rows,
  scope='R implementation measurements; fixed prior and likelihood, new reference streams and specified four-chain BART 2.9.9 fits'),file.path(out,'RESULTS.json'))
save_json(list(status='PASS',session=capture.output(sessionInfo())),file.path(out,'COMPLETE.json'))
print(methods)
