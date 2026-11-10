
argv<-commandArgs(trailingOnly=TRUE);stopifnot(length(argv)==1L,!file.exists(argv[1]))
script<-sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE));root<-normalizePath(file.path(dirname(script),'..'))
out<-argv[1];dir.create(out,recursive=TRUE);out<-normalizePath(out)
for(f in c('reference.R','decisions.R','feasibility.R'))source(file.path(root,'R',f))
save_json<-function(x,p)jsonlite::write_json(x,p,auto_unbox=TRUE,pretty=TRUE,digits=17,na='null',null='null')
sha<-function(p)digest::digest(file=p,algo='sha256');elapsed<-function()unname(proc.time()[3])
plan<-list(status='FROZEN_BEFORE_GENERATION_AND_MEASUREMENT',time_UTC=format(Sys.time(),tz='UTC'),
 purpose='Descriptive feasibility and three new fixed reviews; no population-level improvement criterion',
 designs=list(repeated_profiles=c(6L,8L,10L,12L),continuous_profiles=c(12L,20L,30L),
              counts_repeated=5L,counts_continuous=1L,cuts=c('coarse: at most three observed gaps per variable','fine: every observed gap')),
 limits=list(seconds=30,states=100000L,partitions=10000L,products=2000000L),
 guard_pilot='An initial catalog-only pilot used limits of 6000 states and 4000 partitions; all three planned review catalogs stopped before patient-outcome generation. Limits were increased before any review outcome was generated or sampled. All pilot results are retained.',
 reviews=list(list(id=1L,H=6L,cuts='fine',scenario='separated',seed=611080001L),
              list(id=2L,H=8L,cuts='coarse',scenario='separated',seed=611080002L),
              list(id=3L,H=8L,cuts='coarse',scenario='near-boundary',seed=611080003L)),
 truth=list(T_separated=c(.08,.13,.22,.32,.62),E_separated=c(.12,.70,.32,.22,.08),
            T_boundary=c(.44,.47,.48,.49,.52),E_boundary=c(.38,.40,.41,.40,.39),probit_factor_coefficient=.35),
 proposals=1000000L,trees=50L,chunk=2048L,reference_bank_error=.01,scalar_error=.01/60,
 native=list(chains=4L,burnin=8000L,retained_per_chain=32000L,version='2.9.9',offset=0),
 order='All catalogs; each review references in alternating endpoint order, followed by four sequential BART chains',
 sources=setNames(lapply(c(script,file.path(root,'R/feasibility.R'),file.path(root,'R/reference.R'),
                          file.path(root,'R/decisions.R'),file.path(root,'scripts/fit_review.R')),sha),
                  c('feasibility_study.R','feasibility.R','reference.R','decisions.R','fit_review.R')),
 scope='Fixed seeds and resource limits; failures retained; object memory measures stored catalog state, not process peak RSS')
save_json(plan,file.path(out,'PLAN.json'));cat('Study plan frozen\n');flush.console()
make_X<-function(H) {
  set.seed(611081000L+H)
  x<-matrix(runif(H*5L),H,5L);x[,1]<-sort(x[,1]);x[,4]<-rep(0:1,length.out=H)
  colnames(x)<-paste0('x',1:5);x
}
make_cuts<-function(X,density)lapply(1:ncol(X),function(v) {
  xs<-sort(unique(X[,v]));gaps<-(xs[-1]+xs[-length(xs)])/2
  if(density=='coarse'&&length(gaps)>3L)gaps<-gaps[unique(round(seq(1,length(gaps),length.out=3L)))]
  gaps
})
catalogs<-list();records<-list();k<-0L
for(kind in c('repeated','continuous'))for(H in if(kind=='repeated')c(6L,8L,10L,12L)else c(12L,20L,30L))for(density in c('coarse','fine')) {
  id<-paste(kind,H,density,sep='_');X<-make_X(H);cuts<-make_cuts(X,density);n<-rep(if(kind=='repeated')5L else 1L,H)
  start<-elapsed();z<-catalog_law_limited(X,cuts,n)
  if(z$status=='COMPLETE') {
    a<-prepare_reference(z$law,H)

    if(H<=8L) {
      ordinary<-catalog_law(X,cuts,n);aa<-vapply(z$law$probabilities,as.numeric,0);bb<-vapply(ordinary$probabilities,as.numeric,0)
      stopifnot(setequal(names(aa),names(bb)),max(abs(aa-bb[names(aa)]))<1e-12)
    }
    catalogs[[id]]<-list(X=X,cuts=cuts,counts=n,a=a,setup_seconds=z$seconds,law=z$law)
    prepared_bytes<-as.numeric(object.size(a));partitions<-length(a$parts)
  } else {prepared_bytes<-NA_real_;partitions<-NA_integer_}
  k<-k+1L;records[[k]]<-data.frame(id=id,kind=kind,H=H,patients_per_dose=sum(n),density=density,
    cuts=sum(lengths(cuts)),status=z$status,reason=if(is.null(z$reason))''else z$reason,seconds=z$seconds,
    states=z$states,products=z$products,partitions=partitions,memo_bytes=z$memo_bytes,prepared_bytes=prepared_bytes)
  saveRDS(catalogs,file.path(out,'CATALOGS.rds'));write.csv(do.call(rbind,records),file.path(out,'FEASIBILITY.csv'),row.names=FALSE,na='')
  cat(id,z$status,round(z$seconds,2),'seconds;',partitions,'partitions\n');flush.console()
}
results<-list();all_endpoint_rows<-list();summaries<-list();si<-0L
for(spec in plan$reviews) {
  name<-sprintf('review%02d',spec$id);folder<-file.path(out,name);dir.create(folder)
  key<-paste('repeated',spec$H,spec$cuts,sep='_')
  if(is.null(catalogs[[key]])){results[[name]]<-list(status='CATALOG_UNAVAILABLE');next}
  catz<-catalogs[[key]];X<-catz$X;H<-nrow(X);n<-catz$counts
  set.seed(spec$seed);dat<-list();successes<-list();units<-list();k<-0L
  riskT<-if(spec$scenario=='separated')plan$truth$T_separated else plan$truth$T_boundary
  riskE<-if(spec$scenario=='separated')plan$truth$E_separated else plan$truth$E_boundary
  XT<-as.numeric(scale(X[,1]));XE<-as.numeric(scale(X[,2]))
  for(j in 1:5) {
    index<-rep(seq_len(H),n);d<-data.frame(dose=j,profile=index,X[index,,drop=FALSE])
    d$T<-rbinom(length(index),1,pnorm(qnorm(riskT[j])+.35*XT[index]))
    d$E<-rbinom(length(index),1,pnorm(qnorm(riskE[j])+.35*XE[index]))
    dat[[j]]<-d
    for(ep in c('T','E')){k<-k+1L;units[[k]]<-list(dose=j,endpoint=ep,successes=as.numeric(tapply(d[[ep]],d$profile,sum)))}
  }
  case<-list(id=spec$id,X=X,cuts=catz$cuts,counts=n,dat=do.call(rbind,dat),units=units,fit_seed=611082000L+100L*spec$id)
  saveRDS(case,file.path(folder,'CASE.rds'));write.csv(case$dat,file.path(folder,'PATIENTS.csv'),row.names=FALSE)
  rows<-list();k<-0L
  for(i in seq_along(units)) {
    u<-units[[i]]
    for(method in if(i%%2)c('direct','aggregated')else c('aggregated','direct')) {
      set.seed(611083000L+spec$id*100L+i*2L+as.integer(method=='aggregated'))
      start<-elapsed();v<-sample_reference(plan$proposals,plan$chunk,50L,catz$a,u$successes,n,n/sum(n),0,method)
      b<-probability_intervals(v$q,u$endpoint,plan$scalar_error);sec<-elapsed()-start;k<-k+1L
      rows[[k]]<-c(list(method=method,dose=u$dose,endpoint=u$endpoint,seconds=sec),b)
      writeBin(v$q,file.path(folder,paste0(method,'_',u$endpoint,u$dose,'.bin')),size=8,endian='little')
      save_json(rows,file.path(folder,'ENDPOINTS.json'))
      cat(name,method,u$endpoint,u$dose,round(sec,2),'seconds;',length(v$q),'accepted\n');flush.console()
    }
  }
  methods<-list()
  for(method in c('direct','aggregated')) {
    r<-Filter(function(v)v$method==method,rows);start<-elapsed()
    actions<-setNames(lapply(c('IA1','Final'),function(s)work_actions(r,s)),c('IA1','Final'))
    methods[[method]]<-list(seconds=catz$setup_seconds+sum(vapply(r,function(v)v$seconds,0))+elapsed()-start,
      accepted_min=min(vapply(r,function(v)v$accepted,0)),accepted_max=max(vapply(r,function(v)v$accepted,0)),actions=actions)
  }

  for(stage in c('IA1','Final'))stopifnot(length(intersect(methods$direct$actions[[stage]],methods$aggregated$actions[[stage]]))>0L)
  rs<-file.path(R.home('bin'),'Rscript');native<-file.path(folder,'native');dir.create(native);start<-elapsed()
  for(role in 1:4) {
    status<-system2(rs,c(shQuote(file.path(root,'scripts/fit_review.R')),shQuote(file.path(folder,'CASE.rds')),
      shQuote(file.path(native,paste0('chain',role,'.rds'))),role),
      stdout=file.path(native,paste0('chain',role,'.log')),stderr=file.path(native,paste0('chain',role,'.log')))
    stopifnot(status==0L);cat(name,'native chain',role,'complete\n');flush.console()
  }
  chain<-lapply(1:4,function(role)readRDS(file.path(native,paste0('chain',role,'.rds'))))
  d<-study_diagnostics(chain,plan$scalar_error);write.csv(d,file.path(native,'DIAGNOSTICS.csv'),row.names=FALSE,na='')
  nr<-rows_from_diagnostics(d);na<-setNames(lapply(c('IA1','Final'),function(s)work_actions(nr,s)),c('IA1','Final'))
  methods$native<-list(seconds=elapsed()-start,actions=na,qualified=sum(d$qualified))

  d95<-study_diagnostics(chain,.05/20);write.csv(d95,file.path(native,'DIAGNOSTICS_95.csv'),row.names=FALSE,na='')
  zz<-diagnostic_intervals(d95);fitted95<-setNames(lapply(c('IA1','Final'),function(s)interval_actions(s,zz$gates,zz$scores)),c('IA1','Final'))
  for(stage in c('IA1','Final')) {
    ref<-methods$aggregated$actions[[stage]];fit<-fitted95[[stage]];si<-si+1L
    summaries[[si]]<-data.frame(review=name,H=H,scenario=spec$scenario,stage=stage,
      reference_size=length(ref),fitted_size=length(fit),contained=all(ref%in%fit),
      established_miss=!any(ref%in%fit),unresolved=any(ref%in%fit)&&!all(ref%in%fit))
  }
  results[[name]]<-list(status='COMPLETE',spec=spec,catalog_id=key,partitions=length(catz$a$parts),
    setup_seconds=catz$setup_seconds,methods=methods,endpoint_results=rows,fitted95=fitted95)
  save_json(results[[name]],file.path(folder,'RESULTS.json'));save_json(results,file.path(out,'RESULTS.json'))
  write.csv(do.call(rbind,summaries),file.path(out,'ACCURACY.csv'),row.names=FALSE)
}
save_json(list(status='PASS',plan_sha256=sha(file.path(out,'PLAN.json')),reviews=length(results),
  limits_stop_without_partial_law=TRUE,session=capture.output(sessionInfo())),file.path(out,'COMPLETE.json'))
cat('Three-review feasibility study complete\n')
