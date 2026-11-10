

audit_reviewer_comparisons <- function(root,out) {
  folder<-file.path(root,'data/R_reviewer_comparisons')
  plan<-read_json(file.path(folder,'PLAN.json'));result<-read_json(file.path(folder,'RESULTS.json'))
  sha<-function(p)digest::digest(file=p,algo='sha256')
  stopifnot(result$status=='COMPLETE',sha(file.path(folder,'PLAN.json'))==result$plan_sha256,
            plan$benchmark_proposals==50000L,plan$benchmark_scalar_error==.01/180,
            plan$budget_scalar_error==.01/60)
  for(n in names(plan$sources)) {
    p<-if(n=='reviewer_comparisons.R')file.path(root,'scripts',n)else file.path(root,'R',n)
    stopifnot(source_hash_matches(p,plan$sources[[n]],root),
      source_hash_matches(file.path(folder,'frozen_sources',n),plan$sources[[n]],root))
  }
  same<-function(actual,expected,fields)for(f in fields) {
    b<-expected[[f]];if(is.null(b))b<-NA_real_
    stopifnot(isTRUE(all.equal(as.numeric(actual[[f]]),as.numeric(b),tolerance=1e-11)))
  }
  fields<-c('accepted','gate_count','mean','mean_lower','mean_upper','gate_lower','gate_upper')
  endpoint_checks<-0L;action_checks<-0L
  read_draw<-function(p,endpoint,delta) {
    d<-readRDS(p);stopifnot(length(d$q)==length(d$accepted_at),all(is.finite(d$q)),
      all(d$q>=0&d$q<=1),all(diff(d$accepted_at)>0),all(d$accepted_at>=1),
      all(d$accepted_at==floor(d$accepted_at)))
    for(b in d$reports) {
      stopifnot(all(d$accepted_at<=max(vapply(d$reports,function(x)x$proposals,0))))
      same(probability_intervals(d$q[d$accepted_at<=b$proposals],endpoint,delta),b,fields)
      endpoint_checks<<-endpoint_checks+1L
    };d
  }
  totals<-read_json(file.path(folder,'BENCHMARK_SUMMARY.json'))
  stopifnot(length(totals)==9L,length(read_json(file.path(folder,'BENCHMARK_ENDPOINTS.json')))==90L)
  for(i in 1:3) {
    cf<-file.path(folder,sprintf('review%02d',i));rows<-read_json(file.path(cf,'ENDPOINTS.json'))
    stopifnot(length(rows)==30L)
    case<-readRDS(file.path(root,'data/R_feasibility',sprintf('review%02d/CASE.rds',i)))
    prep<-readRDS(file.path(cf,'PREPARATION.rds'))
    stopifnot(identical(prep$tree$counts,as.numeric(case$counts)),prep$tree$H==nrow(case$X))
    for(e in rows) {
      stem<-paste(e$method,e$endpoint,e$dose,sep='_')
      d<-read_draw(file.path(cf,paste0(stem,'.rds')),e$endpoint,plan$benchmark_scalar_error)
      stopifnot(length(d$reports)==1L,d$reports[[1]]$proposals==plan$benchmark_proposals)
      same(d$reports[[1]],e,c(fields,'proposals','seconds'))
      if(e$method=='tree_prior')stopifnot(d$tree_attempts>=d$valid_tree_attempts,
                                         d$valid_tree_attempts>=50000*50)
    }
    for(method in unlist(plan$methods)) {
      rr<-Filter(function(x)x$method==method,rows)
      s<-Filter(function(x)x$review==i&&x$method==method,totals);stopifnot(length(s)==1L);s<-s[[1]]
      stopifnot(s$accepted_min==min(vapply(rr,function(x)x$accepted,0)),
                s$accepted_max==max(vapply(rr,function(x)x$accepted,0)),
                s$seconds>=s$setup_seconds+sum(vapply(rr,function(x)x$seconds,0))-1e-7)
      for(stage in c('IA1','Final')) {
        stopifnot(setequal(work_actions(rr,stage),unlist(s$actions[[stage]])));action_checks<-action_checks+1L
      }
    }
  }
  rows<-read_json(file.path(folder,'budget/ENDPOINTS.json'));stopifnot(length(rows)==30L)
  for(endpoint in c('T','E'))for(dose in 1:5) {
    d<-read_draw(file.path(folder,'budget',paste0(endpoint,dose,'.rds')),endpoint,plan$budget_scalar_error)
    stopifnot(identical(as.numeric(names(d$reports)),as.numeric(unlist(plan$budget_proposals))))
    for(b in d$reports) {
      r<-Filter(function(x)x$endpoint==endpoint&&x$dose==dose&&x$proposals==b$proposals,rows)
      stopifnot(length(r)==1L);same(b,r[[1]],c(fields,'seconds'))
    }
  }
  budget<-read_json(file.path(folder,'BUDGET_SUMMARY.json'));stopifnot(length(budget)==3L)
  for(s in budget) {
    rr<-Filter(function(x)x$proposals==s$proposals,rows)
    stopifnot(s$accepted_min==min(vapply(rr,function(x)x$accepted,0)),
              s$accepted_max==max(vapply(rr,function(x)x$accepted,0)))
    for(stage in c('IA1','Final')) {
      stopifnot(setequal(work_actions(rr,stage),unlist(s$actions[[stage]])));action_checks<-action_checks+1L
    }
    stopifnot(abs(s$largest_gate_width-max(vapply(rr,function(x)x$gate_upper-x$gate_lower,0)))<1e-12,
      abs(s$largest_efficacy_mean_width-max(vapply(Filter(function(x)x$endpoint=='E',rr),
        function(x)x$mean_upper-x$mean_lower,0)))<1e-12)
  }
  stopifnot(result$validation$status=='PASS')
  value<-list(status='PASS',source_hashes=TRUE,
    source_verification='Original hashes or recorded comment-free source transformations',
    plan_hash=TRUE,matched_reviews=3L,
    routes=3L,reconstructed_endpoint_intervals=endpoint_checks,reconstructed_action_sets=action_checks,
    fixed_budget_prefixes=3L,new_timing=FALSE,pilot_used_as_final=FALSE)
  save_json(value,file.path(out,'REVIEWER_COMPARISON_CHECKS.json'));value
}
