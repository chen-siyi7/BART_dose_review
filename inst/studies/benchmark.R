
argv<-commandArgs(trailingOnly=TRUE);stopifnot(length(argv)%in%c(1L,2L))
script<-sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE))
root<-normalizePath(file.path(dirname(script),'..'));out<-argv[1];mode<-if(length(argv)==2L) argv[2] else 'all'
stopifnot(mode%in%c('short','workload','all'),!file.exists(out))
dir.create(out,recursive=TRUE);out<-normalizePath(out)
for(n in c('reference.R','decisions.R')) source(file.path(root,'R',n))
read_json<-function(p) jsonlite::fromJSON(p,simplifyVector=FALSE)
save_json<-function(x,p) jsonlite::write_json(x,p,pretty=TRUE,auto_unbox=TRUE,digits=17,na='null',null='null')
sha<-function(p) digest::digest(file=p,algo='sha256')
elapsed<-function() unname(proc.time()[3])
plan<-list(status='FROZEN_BEFORE_R_MEASUREMENTS',language='R',R=R.version.string,platform=R.version$platform,
  selection='All original 42 model-history pairs; previously chosen repeated-grid case 3 for the full workload',
  methods=c('direct','aggregated'),trees=50L,chunk=2048L,short_proposals=1024L,short_repetitions=5L,
  full_proposals=1000000L,full_repetitions=1L,delta=.01/20,seed_base=2061109000L,
  order='Short batches then full workload; reference methods alternate by endpoint; four sequential native chains follow',
  aggregation='Independent Gaussian sums for each unique full-point leaf mask, conditional on multinomial partition counts',
  source_hashes=setNames(lapply(c(file.path(root,'R/reference.R'),file.path(root,'R/decisions.R'),script,
    file.path(root,'scripts/fit_chain.R'),file.path(root,'scripts/diagnose.R')),sha),
    c('reference.R','decisions.R','benchmark.R','fit_chain.R','diagnose.R')),
  inference='Single-machine descriptive timing. Models, cases and proposal counts fixed before this R measurement.')
save_json(plan,file.path(out,'PLAN.json'))
cat('Frozen R benchmark plan\n');flush.console()
input<-read_json(file.path(root,'data/base/compact_reference/bank_v2/INPUTS.json'))
prep<-new.env(parent=emptyenv());prep_times<-numeric()
for(cat in input$catalogs) {
  t0<-elapsed();X<-do.call(rbind,lapply(cat$X_unique,unlist))
  law<-catalog_law(X,cat$cutpoints,unlist(cat$pattern_counts));a<-prepare_reference(law,nrow(X))
  expected<-setNames(vapply(cat$probability_strings,function(w) as.numeric(gmp::as.bigq(w)),0),
     vapply(cat$partitions,function(x) paste(sort(unlist(x)),collapse=','),''))
  actual<-vapply(law$probabilities,as.numeric,0)
  stopifnot(setequal(names(actual),names(expected)),max(abs(actual-expected[names(actual)]))<1e-12)
  assign(cat$id,list(a=a,cat=cat),prep);prep_times[cat$id]<-elapsed()-t0
}
save_json(as.list(prep_times),file.path(out,'PREPROCESSING.json'))
cat('Prepared',length(ls(prep)),'catalogs\n');flush.console()
if(mode%in%c('short','all')) {
  times<-list();k<-0L
  for(i in seq_along(input$units)) {
    u<-input$units[[i]];if(u$status!='TIMED') next
    z<-get(u$catalog_id,prep);counts<-unlist(z$cat$pattern_counts)
    for(rep in 1:5) for(method in (if((i+rep)%%2) c('direct','aggregated') else c('aggregated','direct'))) {
      set.seed(plan$seed_base+i*20L+rep*2L+as.integer(method=='aggregated'))
      t0<-elapsed();v<-sample_reference(1024L,2048L,50L,z$a,unlist(u$successes),counts,counts/sum(counts),u$offset,method)
      k<-k+1L;times[[k]]<-data.frame(id=u$id,bank=u$bank,case_id=u$case_id,policy=u$policy,method=method,
        repetition=rep,seconds=elapsed()-t0,accepted=length(v$q))
    }
    if(i%%20L==0L) {cat('Short unit',i,'/',length(input$units),'\n');flush.console()}
  }
  tt<-do.call(rbind,times);write.csv(tt,file.path(out,'SHORT_TIMING.csv'),row.names=FALSE)
  keys<-unique(tt[c('bank','case_id','policy')]);pair<-list()
  for(i in seq_len(nrow(keys))) {
    x<-keys[i,];rr<-tt[tt$bank==x$bank&tt$case_id==x$case_id&tt$policy==x$policy,,drop=FALSE]
    units<-Filter(function(u) u$bank==x$bank&&u$case_id==x$case_id&&u$policy==x$policy,input$units)
    setup<-sum(prep_times[unique(vapply(units,function(u) u$catalog_id,''))])
    total<-vapply(c('direct','aggregated'),function(method) {
      r<-rr[rr$method==method,,drop=FALSE]; setup+sum(tapply(r$seconds,r$id,median))
    },0)
    pair[[i]]<-data.frame(x,direct_seconds=total[1],aggregated_seconds=total[2],ratio=total[2]/total[1],row.names=NULL)
  }
  pp<-do.call(rbind,pair);write.csv(pp,file.path(out,'SHORT_PAIR_TOTALS.csv'),row.names=FALSE)
  summaries<-lapply(c('small','deeper','crossed','all'),function(bank) {
    p<-if(bank=='all') pp else pp[pp$bank==bank,,drop=FALSE]
    data.frame(bank=bank,pairs=nrow(p),direct_ms=1000*median(p$direct_seconds),aggregated_ms=1000*median(p$aggregated_seconds),ratio=median(p$ratio),aggregation_faster=sum(p$ratio<1))
  })
  ss<-do.call(rbind,summaries);write.csv(ss,file.path(out,'SHORT_SUMMARY.csv'),row.names=FALSE);print(ss);flush.console()
}
if(mode%in%c('workload','all')) {
  work<-read_json(file.path(root,'data/worked_review/INPUTS.json'));rows<-list();k<-0L
  for(i in seq_along(work$units)) {
    u<-work$units[[i]];z<-get(u$catalog_id,prep);counts<-unlist(z$cat$pattern_counts)
    for(method in (if(i%%2) c('direct','aggregated') else c('aggregated','direct'))) {
      seed<-plan$seed_base+10000L+i*10L+as.integer(method=='aggregated');set.seed(seed)
      t0<-elapsed();v<-sample_reference(1000000L,2048L,50L,z$a,unlist(u$successes),counts,counts/sum(counts),u$offset,method)
      b<-probability_intervals(v$q,u$endpoint,plan$delta);sec<-elapsed()-t0
      k<-k+1L;rows[[k]]<-c(list(method=method,unit_id=u$id,dose=u$dose,endpoint=u$endpoint,seed=seed,seconds=sec),b,
        list(log_likelihood_bound=v$log_likelihood_bound,unit_envelope_accepts_same_uniform=v$unit_envelope_accepts_same_uniform))
      stem<-paste0(method,'_',u$endpoint,u$dose);writeBin(v$q,file.path(out,paste0(stem,'.bin')),size=8,endian='little')
      save_json(rows[[k]],file.path(out,paste0(stem,'.json')))
      cat(stem,round(sec,3),'seconds;',length(v$q),'accepted\n');flush.console()
    }
  }
  methods<-lapply(c('direct','aggregated'),function(method) {
    r<-Filter(function(x) x$method==method,rows)
    t0<-elapsed();actions<-setNames(lapply(c('IA1','Final'),function(s) work_actions(r,s)),c('IA1','Final'));action_sec<-elapsed()-t0
    setup<-sum(prep_times[unique(vapply(work$units,function(u) u$catalog_id,''))])
    list(method=method,seconds=setup+sum(vapply(r,function(x)x$seconds,0))+action_sec,catalog_seconds=setup,actions=actions,
      accepted_min=min(vapply(r,function(x)x$accepted,0)),accepted_max=max(vapply(r,function(x)x$accepted,0)))
  })
  native<-file.path(out,'native');dir.create(native);rs<-file.path(R.home('bin'),'Rscript')
  t0<-elapsed();chain_seconds<-numeric(4)
  for(role in 1:4) {
    start<-elapsed();status<-system2(rs,c(shQuote(file.path(root,'scripts/fit_chain.R')),
      shQuote(file.path(root,'data/worked_review/case03.rds')),shQuote(file.path(native,paste0('chain',role,'.rds'))),role),
      stdout=file.path(native,paste0('chain',role,'.log')),stderr=file.path(native,paste0('chain',role,'.log')))
    stopifnot(status==0L);chain_seconds[role]<-elapsed()-start;cat('Native chain',role,chain_seconds[role],'seconds\n');flush.console()
  }
  status<-system2(rs,c(shQuote(file.path(root,'scripts/diagnose.R')),shQuote(native),shQuote(file.path(native,'DIAGNOSTICS.csv'))),
      stdout=file.path(native,'diagnostics.log'),stderr=file.path(native,'diagnostics.log'));stopifnot(status==0L)
  dr<-read.csv(file.path(native,'DIAGNOSTICS.csv'));nr<-list();k<-0L
  for(j in 1:5) for(ep in c('T','E')) {
    mn<-dr[dr$dose==j&dr$series==paste0(ep,'_mean'),];g<-dr[dr$dose==j&dr$series==paste0(ep,'_gate'),]
    k<-k+1L;nr[[k]]<-list(dose=j,endpoint=ep,mean=mn$mean,mean_lower=mn$lower,mean_upper=mn$upper,gate_lower=g$lower,gate_upper=g$upper)
  }
  na<-setNames(lapply(c('IA1','Final'),function(s) work_actions(nr,s)),c('IA1','Final'))
  methods[[3]]<-list(method='native_four_chains',seconds=elapsed()-t0,chain_process_seconds=chain_seconds,
    intervals=nr,actions=na,qualified_quantities=sum(dr$qualified))
  save_json(list(status='COMPLETE',language='R',plan_sha256=sha(file.path(out,'PLAN.json')),methods=methods,endpoint_results=rows,
    scope='R implementation measurements; fixed prior and likelihood, new random streams, no selection by measured speed'),file.path(out,'RESULTS.json'))
  print(methods)
}
save_json(list(status='PASS',mode=mode,session=capture.output(sessionInfo())),file.path(out,'COMPLETE.json'))
