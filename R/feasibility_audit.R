
audit_feasibility_study <- function(data,out) {
  folder<-file.path(data,'R_feasibility');plan<-read_json(file.path(folder,'PLAN.json'))
  complete<-read_json(file.path(folder,'COMPLETE.json'));stopifnot(complete$status=='PASS',complete$reviews==3L)
  stopifnot(digest::digest(file=file.path(folder,'PLAN.json'),algo='sha256')==complete$plan_sha256)
  ff<-read.csv(file.path(folder,'FEASIBILITY.csv'));stopifnot(nrow(ff)==14L)
  cats<-readRDS(file.path(folder,'CATALOGS.rds'))
  for(c in cats)stopifnot(abs(sum(vapply(c$law$probabilities,as.numeric,0))-1)<1e-12)
  c<-cats[['repeated_6_coarse']];stopped<-catalog_law_limited(c$X,c$cuts,c$counts,max_states=0L)
  stopifnot(stopped$status=='RESOURCE_LIMIT',is.null(stopped$law))
  baseline<-catalog_law(c$X,c$cuts,c$counts)
  a<-vapply(baseline$probabilities,as.numeric,0);b<-vapply(c$law$probabilities,as.numeric,0)
  stopifnot(setequal(names(a),names(b)),max(abs(a-b[names(a)]))<1e-12)
  result<-read_json(file.path(folder,'RESULTS.json'));checks<-0L;accuracy<-list();k<-0L
  for(name in names(result)) {
    r<-result[[name]];stopifnot(r$status=='COMPLETE');cf<-file.path(folder,name)
    case<-readRDS(file.path(cf,'CASE.rds'))
    for(u in case$units) {
      d<-case$dat[case$dat$dose==u$dose,,drop=FALSE]
      stopifnot(identical(as.numeric(tapply(d[[u$endpoint]],d$profile,sum)),as.numeric(u$successes)))
    }
    for(e in r$endpoint_results) {
      p<-file.path(cf,paste0(e$method,'_',e$endpoint,e$dose,'.bin'))
      q<-readBin(p,'double',n=file.info(p)$size/8L,size=8,endian='little')
      actual<-probability_intervals(q,e$endpoint,plan$scalar_error)
      for(field in c('accepted','gate_count','mean','mean_lower','mean_upper','gate_lower','gate_upper')) {
        expected<-e[[field]];if(is.null(expected))expected<-NA_real_
        stopifnot(isTRUE(all.equal(as.numeric(actual[[field]]),as.numeric(expected),tolerance=1e-11)))
      };checks<-checks+1L
    }
    for(method in c('direct','aggregated')) {
      rows<-Filter(function(v)v$method==method,r$endpoint_results)
      for(stage in c('IA1','Final'))stopifnot(setequal(work_actions(rows,stage),unlist(r$methods[[method]]$actions[[stage]])))
    }
    chains<-lapply(1:4,function(role)readRDS(file.path(cf,'native',paste0('chain',role,'.rds'))))
    dr<-study_diagnostics(chains,plan$scalar_error);saved<-read.csv(file.path(cf,'native/DIAGNOSTICS.csv'))
    for(field in names(dr))stopifnot(isTRUE(all.equal(dr[[field]],saved[[field]],tolerance=1e-10)))
    rows<-rows_from_diagnostics(dr)
    for(stage in c('IA1','Final'))stopifnot(setequal(work_actions(rows,stage),unlist(r$methods$native$actions[[stage]])))
    d95<-study_diagnostics(chains,.05/20);z<-diagnostic_intervals(d95)
    for(stage in c('IA1','Final')) {
      fit<-interval_actions(stage,z$gates,z$scores);ref<-unlist(r$methods$aggregated$actions[[stage]])
      stopifnot(setequal(fit,unlist(r$fitted95[[stage]])));k<-k+1L
      accuracy[[k]]<-data.frame(review=name,H=r$spec$H,scenario=r$spec$scenario,stage=stage,
        reference_size=length(ref),fitted_size=length(fit),contained=all(ref%in%fit),
        established_miss=!any(ref%in%fit),unresolved=any(ref%in%fit)&&!all(ref%in%fit))
    }
  }
  actual<-do.call(rbind,accuracy);expected<-read.csv(file.path(folder,'ACCURACY.csv'))
  stopifnot(isTRUE(all.equal(actual,expected,check.attributes=FALSE)))
  write.csv(actual,file.path(out,'NEW_REVIEW_ACCURACY.csv'),row.names=FALSE)
  list(status='PASS',design_configurations=14L,completed_catalogs=sum(ff$status=='COMPLETE'),
       limited_catalogs=sum(ff$status=='RESOURCE_LIMIT'),new_reviews=3L,endpoint_reports=checks,
       native_chains=12L,no_partial_law_on_limit=TRUE,baseline_recursion_agrees=TRUE)
}

audit_weighted_example <- function(data,out) {
  rr<-read_json(file.path(data,'base/compact_certificate/INPUTS.json'))$records
  w<-Filter(function(x)x$bank=='deeper'&&x$case_id==3&&x$stage=='Final'&&x$mode=='reference',rr)
  stopifnot(length(w)==1L);z<-prepare_certificate(w[[1]])
  costs<-setNames(rep(1,length(z$atoms)),z$atoms);costs['Tgate:5']<-4
  minimum<-minimum_certificate(z,costs);oracle<-exhaustive_certificate(z,costs)
  displayed<-c('Emean:1','Tgate:2','Emean:2','Emean:3','Emean:4','Tgate:5')
  stopifnot(minimum$cost==6,oracle$cost==6,sum(costs[displayed])==9,
            identical(certificate_actions(z,displayed),'dose:2'),
            identical(certificate_actions(z,minimum$atoms),'dose:2'))
  value<-list(status='PASS',id=w[[1]]$id,weights=as.list(costs),minimum=minimum,
    displayed_cost=sum(costs[displayed]),subset_oracle_cost=oracle$cost,
    scope='Illustrative reporting preference; weights are not measured review effort')
  save_json(value,file.path(out,'WEIGHTED_REPORT_EXAMPLE.json'));value
}
