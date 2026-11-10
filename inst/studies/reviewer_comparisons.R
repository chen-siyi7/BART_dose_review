
argv<-commandArgs(trailingOnly=TRUE);stopifnot(length(argv)==1L,!file.exists(argv[1]))
script<-sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE))
root<-normalizePath(file.path(dirname(script),'..'));out<-argv[1];dir.create(out,recursive=TRUE);out<-normalizePath(out)
for(f in c('reference.R','decisions.R','feasibility.R','prior_trees.R','prior_tree_audit.R'))source(file.path(root,'R',f))
savej<-function(x,p)jsonlite::write_json(x,p,auto_unbox=TRUE,pretty=TRUE,digits=17,na='null',null='null')
sha<-function(p)digest::digest(file=p,algo='sha256');elapsed<-function()unname(proc.time()[3])
cases<-lapply(1:3,function(i)readRDS(file.path(root,'data/R_feasibility',sprintf('review%02d/CASE.rds',i))))
plan<-list(status='FROZEN_BEFORE_NEW_VALIDATION_AND_MEASUREMENT',time_UTC=format(Sys.time(),tz='UTC'),
  purpose='Reviewer-requested descriptive comparisons; existing fixed reviews selected before new measurements',
  cases='All three previously reported reviews; review 3 is reused for its previously unresolved near-threshold decisions',
  methods=c('tree_prior','direct','aggregated'),trees=50L,benchmark_proposals=50000L,
  chunk=256L,reference_error=.01,benchmark_scalar_error=.01/180,
  budget_review=3L,budget_proposals=c(100000L,1000000L,2000000L),budget_method='direct',budget_scalar_error=.01/60,
  order='Reviews in numerical order; methods rotate across dose-endpoints; the budget study follows the benchmark',
  timing='One sequential R run per route and review; all design preparation, prior generation, rejection, intervals and action scan included; descriptive machine-specific timings',
  seed_base=611091000L,existing_patient_outcomes=TRUE,new_native_fits=FALSE,
  earlier_run='An initial partial run is retained as R_reviewer_comparisons_pilot. The endpoint timer was moved past the final interval construction. Cases, seeds, methods and budgets are unchanged; final times come from this complete rerun.',
  sources=setNames(lapply(c(script,file.path(root,'R/prior_trees.R'),file.path(root,'R/prior_tree_audit.R'),
    file.path(root,'R/reference.R'),file.path(root,'R/decisions.R')),sha),
    c('reviewer_comparisons.R','prior_trees.R','prior_tree_audit.R','reference.R','decisions.R')),
  session=capture.output(sessionInfo()),threads=as.list(Sys.getenv(c('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','VECLIB_MAXIMUM_THREADS'))))
savej(plan,file.path(out,'PLAN.json'));cat('Reviewer-comparison plan frozen\n');flush.console()
checks<-audit_prior_trees(root,out);cat('Complete-tree prior and posterior checks passed\n');flush.console()
all_rows<-list();summaries<-list();k<-0L;h<-0L
for(i in 1:3) {
  case<-cases[[i]];folder<-file.path(out,sprintf('review%02d',i));dir.create(folder)
  start<-elapsed();tree<-prepare_tree_prior(case$X,case$cuts,case$counts);tree_setup<-elapsed()-start
  start<-elapsed();law<-catalog_law_limited(case$X,case$cuts,case$counts);stopifnot(law$status=='COMPLETE')
  catalog<-prepare_reference(law$law,nrow(case$X));catalog_setup<-elapsed()-start
  saveRDS(list(tree=tree,catalog=catalog,tree_setup=tree_setup,catalog_setup=catalog_setup),file.path(folder,'PREPARATION.rds'))
  rows<-list()
  for(u in seq_along(case$units)) {
    unit<-case$units[[u]];order<-plan$methods[((seq_len(3L)+u+i-3L)%%3L)+1L]
    for(method in order) {
      set.seed(plan$seed_base+i*100L+u*3L+match(method,plan$methods))
      a<-if(method=='tree_prior')tree else catalog
      draw<-sample_reference_checkpoints(plan$benchmark_proposals,plan$chunk,50L,a,unit$successes,
        case$counts,case$counts/sum(case$counts),0,method,plan$benchmark_scalar_error,unit$endpoint)
      band<-draw$reports[[1]];r<-c(list(review=i,method=method,dose=unit$dose,endpoint=unit$endpoint),band,
        list(tree_attempts=draw$tree_attempts,valid_tree_attempts=draw$valid_tree_attempts,visited_nodes=draw$visited_nodes))
      stem<-paste(method,unit$endpoint,unit$dose,sep='_');saveRDS(draw,file.path(folder,paste0(stem,'.rds')))
      k<-k+1L;all_rows[[k]]<-r;rows[[length(rows)+1L]]<-r
      savej(rows,file.path(folder,'ENDPOINTS.json'))
      cat('Review',i,method,paste0(unit$endpoint,unit$dose),round(band$seconds,2),'s;',band$accepted,'accepted\n');flush.console()
    }
  }
  for(method in plan$methods) {
    r<-Filter(function(x)x$method==method,rows);start<-elapsed()
    actions<-setNames(lapply(c('IA1','Final'),function(s)work_actions(r,s)),c('IA1','Final'))
    h<-h+1L;summaries[[h]]<-list(review=i,method=method,setup_seconds=if(method=='tree_prior')tree_setup else catalog_setup,
      seconds=(if(method=='tree_prior')tree_setup else catalog_setup)+sum(vapply(r,function(x)x$seconds,0))+elapsed()-start,
      accepted_min=min(vapply(r,function(x)x$accepted,0)),accepted_max=max(vapply(r,function(x)x$accepted,0)),actions=actions)
  }
  savej(summaries,file.path(out,'BENCHMARK_SUMMARY.json'))
}
savej(all_rows,file.path(out,'BENCHMARK_ENDPOINTS.json'))
case<-cases[[3]];start<-elapsed();law<-catalog_law_limited(case$X,case$cuts,case$counts)
stopifnot(law$status=='COMPLETE');catalog<-prepare_reference(law$law,nrow(case$X));setup<-elapsed()-start
budget<-file.path(out,'budget');dir.create(budget);rows<-list()
for(u in seq_along(case$units)) {
  unit<-case$units[[u]];set.seed(plan$seed_base+1000L+u)
  draw<-sample_reference_checkpoints(plan$budget_proposals,2048L,50L,catalog,unit$successes,case$counts,
    case$counts/sum(case$counts),0,'direct',plan$budget_scalar_error,unit$endpoint)
  saveRDS(draw,file.path(budget,paste0(unit$endpoint,unit$dose,'.rds')))
  for(j in seq_along(plan$budget_proposals))rows[[length(rows)+1L]]<-c(list(dose=unit$dose,endpoint=unit$endpoint),draw$reports[[j]])
  savej(rows,file.path(budget,'ENDPOINTS.json'))
  cat('Budget review',paste0(unit$endpoint,unit$dose),round(draw$reports[[3]]$seconds,2),'s;',length(draw$q),'accepted\n');flush.console()
}
bs<-lapply(plan$budget_proposals,function(N) {
  r<-Filter(function(x)x$proposals==N,rows);start<-elapsed()
  actions<-setNames(lapply(c('IA1','Final'),function(s)work_actions(r,s)),c('IA1','Final'))
  list(proposals=N,seconds=setup+sum(vapply(r,function(x)x$seconds,0))+elapsed()-start,
    accepted_min=min(vapply(r,function(x)x$accepted,0)),accepted_max=max(vapply(r,function(x)x$accepted,0)),
    largest_gate_width=max(vapply(r,function(x)x$gate_upper-x$gate_lower,0)),
    largest_efficacy_mean_width=max(vapply(Filter(function(x)x$endpoint=='E',r),function(x)x$mean_upper-x$mean_lower,0)),
    actions=actions)
})
savej(bs,file.path(out,'BUDGET_SUMMARY.json'))
savej(list(status='COMPLETE',plan_sha256=sha(file.path(out,'PLAN.json')),benchmark=summaries,budget=bs,
  validation=checks,scope='Fixed-history descriptive comparisons; simultaneous error allocated separately to benchmark and budget families'),file.path(out,'RESULTS.json'))
cat('All reviewer comparisons complete\n')
