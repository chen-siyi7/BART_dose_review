read_json <- function(p) jsonlite::fromJSON(p,simplifyVector=FALSE)
save_json <- function(x,p) jsonlite::write_json(x,p,pretty=TRUE,auto_unbox=TRUE,digits=17,na='null',null='null')
fraction <- function(x) as.numeric(gmp::as.bigq(as.character(x)))
source_hash_matches <- function(path, expected, root) {
  actual <- digest::digest(file=path,algo='sha256')
  if(identical(actual,expected)) return(TRUE)
  root <- normalizePath(root)
  path <- normalizePath(path)
  record <- file.path(root,'SOURCE_CODE_PROVENANCE.json')
  if(!startsWith(path,paste0(root,'/')) || !file.exists(record)) return(FALSE)
  entries <- jsonlite::read_json(record,simplifyVector=FALSE)$files
  entry <- entries[[substring(path,nchar(root)+2L)]]
  revision<-!is.null(entry)&&identical(entry$change,'interval_and_exact_score_handling')&&
    identical(entry$validation$status,'PASS')&&
    identical(entry$validation$file,'REVISION_VALIDATION.json')&&
    file.exists(file.path(root,entry$validation$file))&&
    identical(digest::digest(file=file.path(root,entry$validation$file),algo='sha256'),entry$validation$sha256)
  !is.null(entry) && isTRUE(entry$comment_free) &&
    (entry$change %in% c('comments_only','comments_and_source_verification')||revision) &&
    identical(actual,entry$distributed_sha256) &&
    identical(expected,entry$original_sha256)
}
verify_manifest <- function(root) {
  manifest<-read_json(file.path(root,'MANIFEST.json'))
  actual<-list.files(root,recursive=TRUE,all.files=TRUE,no..=TRUE)
  actual<-actual[!file.info(file.path(root,actual))$isdir]
  stopifnot(setequal(actual,c(names(manifest),'MANIFEST.json')))
  for(n in names(manifest)) {
    p<-file.path(root,n);m<-manifest[[n]]
    stopifnot(!nzchar(Sys.readlink(p)),file.info(p)$size==m$bytes,digest::digest(file=p,algo='sha256')==m$sha256)
  }
  list(status='PASS',members=length(manifest),complete_inventory=TRUE)
}
audit_prior <- function(data,out) {
  bank<-read_json(file.path(data,'base/compact_reference/bank_v2/INPUTS.json'))


  exact_ids<-vapply(c('small','deeper','crossed'),function(b) Filter(function(c) c$bank==b,bank$catalogs)[[1L]]$id,'')
  rows<-list()
  for(i in seq_along(bank$catalogs)) {
    c<-bank$catalogs[[i]];X<-do.call(rbind,lapply(c$X_unique,unlist))
    exact<-c$id%in%exact_ids;law<-catalog_law(X,c$cutpoints,unlist(c$pattern_counts),exact=exact)
    expected<-setNames(c$probability_strings,vapply(c$partitions,function(p) paste(sort(unlist(p)),collapse=','),''))
    stopifnot(setequal(names(law$probabilities),names(expected)))
    if(exact) for(k in names(expected)) stopifnot(law$probabilities[[k]]==gmp::as.bigq(expected[[k]])) else {
      actual<-vapply(law$probabilities,as.numeric,0);target<-vapply(expected,function(w) fraction(w),0)
      stopifnot(max(abs(actual-target[names(actual)]))<1e-12)
    }
    rows[[i]]<-data.frame(id=c$id,partitions=length(law$probabilities),exact_rational=exact)
  }
  write.csv(do.call(rbind,rows),file.path(out,'PRIOR_CATALOGS.csv'),row.names=FALSE)
  list(status='PASS',catalogs=length(rows),rational_catalogs=length(exact_ids))
}
audit_prediction_points <- function(data) {
  input<-read_json(file.path(data,'worked_review/INPUTS.json'));c<-input$catalogs[[1L]]
  X<-do.call(rbind,lapply(c$X_unique,unlist));n<-unlist(c$pattern_counts);H<-nrow(X)
  extra<-rbind(vapply(c$cutpoints,function(v) unlist(v)[1L],0),apply(X,2L,max)+10)
  original<-catalog_law(X,c$cutpoints,n,exact=TRUE)
  full<-catalog_law(rbind(X,extra),c$cutpoints,c(n,0,0),exact=TRUE)
  marginal<-list();mask<-2^H-1L
  for(k in names(full$probabilities)) {
    leaves<-as.integer(strsplit(k,',',fixed=TRUE)[[1L]]);leaves<-bitwAnd(leaves,mask);leaves<-sort(leaves[leaves>0])
    key<-paste(leaves,collapse=',');if(is.null(marginal[[key]])) marginal[[key]]<-gmp::as.bigq(0)
    marginal[[key]]<-marginal[[key]]+full$probabilities[[k]]
  }
  stopifnot(setequal(names(marginal),names(original$probabilities)))
  for(k in names(marginal)) stopifnot(marginal[[k]]==original$probabilities[[k]])
  list(status='PASS',training_profiles=H,prediction_only_profiles=2L,exact_rational_marginal_preserved=TRUE)
}
audit_certificates <- function(data,out) {
  inputs<-read_json(file.path(data,'base/compact_certificate/INPUTS.json'))$records
  expected<-read_json(file.path(data,'base/compact_certificate/RECORDS.json'))
  rows<-list();subsets<-0L
  for(i in seq_along(inputs)) {
    z<-prepare_certificate(inputs[[i]]);e<-expected[[i]]
    stopifnot(setequal(certificate_actions(z),unlist(inputs[[i]]$full_actions)))
    a<-minimum_certificate(z);stopifnot(a$status==e$status)
    if(a$status=='CERTIFIED') {
      b<-exhaustive_certificate(z);subsets<-subsets+b$subsets
      d<-vapply(c('dose-forward','dose-reverse','target-mean-first','descending-cost'),function(o) deletion_certificate(z,o)$cost,0)
      stopifnot(a$cost==e$cost,b$cost==a$cost,min(d)==a$cost,
        identical(certificate_actions(z,a$atoms),z$target))

      for(atom in a$atoms) stopifnot(!identical(certificate_actions(z,setdiff(a$atoms,atom)),z$target))
    }
    rows[[i]]<-data.frame(id=inputs[[i]]$id,status=a$status,cost=a$cost,full=length(z$available))
    if(i%%40L==0L) {cat('Certificates',i,'/',length(inputs),'\n');flush.console()}
  }
  rr<-do.call(rbind,rows);write.csv(rr,file.path(out,'CERTIFICATES.csv'),row.names=FALSE)
  list(status='PASS',reports=nrow(rr),certifiable=sum(rr$status=='CERTIFIED'),unresolved=sum(rr$status!='CERTIFIED'),
       median_minimum=median(rr$cost,na.rm=TRUE),median_full=median(rr$full[rr$status=='CERTIFIED']),subset_rectangles=subsets)
}
audit_frontiers <- function(data,out) {
  source<-file.path(data,'base/explanation/source/independent_refit_20261028')
  cells<-read.csv(file.path(source,'CELLS.csv'),stringsAsFactors=FALSE)
  pairs<-read.csv(file.path(data,'base/explanation/report/PAIR_ROWS.csv'),stringsAsFactors=FALSE)
  targets<-read.csv(file.path(data,'base/explanation/report/TARGET_ROWS.csv'),stringsAsFactors=FALSE)
  key<-function(x) paste(x$case_id,x$refit,sep='|');by<-split(cells,key(cells));tr<-split(targets,key(targets))
  checks<-0L;changed<-0L;ranking<-0L
  for(i in seq_len(nrow(pairs))) {
    p<-pairs[i,];if(p$status!='valid') next
    cc<-by[[key(p)]];cc<-cc[order(cc$position),];eta<-max(abs(cc$q_original-cc$q_refit))
    actual<-frontier_costs('Final',cc$dose,matrix(cc$k_original,ncol=1),cc$q_original,eta)
    t<-tr[[key(p)]];wanted<-actual[paste0('dose:',t$target)]
    stopifnot(all(wanted==vapply(t$frontier_cost,fraction,0)),sum(abs(cc$k_original-cc$k_refit))==p$observed_count_L1)

    keep<-cc$k_original>=300;rho<-ifelse(keep,cc$k_original-299,pmax(0,300-cc$k_original));oracle<-setNames(rep(Inf,length(cc$dose)+1L),c('dose:0',paste0('dose:',cc$dose)))
    bits<-all_masks(length(cc$dose));lo<-pmax(0,cc$q_original-eta);hi<-pmin(1,cc$q_original+eta)
    for(b in seq_len(nrow(bits))) {
      member<-which(bits[b,]);cost<-sum(rho[bits[b,]!=keep])
      if(!length(member)) oracle['dose:0']<-min(oracle['dose:0'],cost)
      for(r in member) {
        q<-lo;q[r]<-hi[r];winner<-member[which.max(q[member])]
        if(winner==r) {name<-paste0('dose:',cc$dose[r]);oracle[name]<-min(oracle[name],cost)}
      }
    }
    stopifnot(all(actual==oracle[names(actual)]))
    checks<-checks+nrow(t);changed<-changed+as.integer(p$changed==1L)
    ranking<-ranking+as.integer(p$changed==1L&&actual[paste0('dose:',p$refit_action)]==0)
  }
  stopifnot(checks==6893L,changed==140L,ranking==136L)
  list(status='PASS',planned_pairs=nrow(pairs),valid_pairs=sum(pairs$status=='valid'),target_costs=checks,changed=changed,ranking_only=ranking)
}
audit_intercept <- function(path,out) {
  truth<-read.csv(file.path(path,'EXACT_POSTERIORS.csv'));jobs<-read.csv(file.path(path,'INPUTS.csv'))
  counts<-read.csv(file.path(path,'COUNTS.csv'));diagnostics<-read.csv(file.path(path,'DIAGNOSTICS.csv'))
  for(i in seq_len(nrow(truth))) {
    r<-truth[i,];logf<-function(t) dnorm(t,sd=1.5,log=TRUE)+(if(r$s>0)r$s*pnorm(t,log.p=TRUE)else 0)+(if(r$s<r$n)(r$n-r$s)*pnorm(-t,log.p=TRUE)else 0)
    shift<-optimize(logf,c(-10,10),maximum=TRUE)$objective;f<-function(t) exp(logf(t)-shift)
    integral<-function(fun,l=-Inf,u=Inf) integrate(fun,l,u,rel.tol=1e-11,abs.tol=1e-12,subdivisions=1000L)$value
    norm<-integral(f);v<-c(integral(function(t)pnorm(t)*f(t))/norm,integral(f,u=qnorm(.45))/norm,integral(f,l=qnorm(.3))/norm)
    stopifnot(max(abs(v-as.numeric(r[c('mean','T_gate','E_gate')])) )<1e-10)
  }
  rows<-list();k<-0L
  for(i in seq_len(nrow(jobs))) {
    job<-jobs[i,];cc<-counts[counts$case_id==job$case_id,];dr<-diagnostics[diagnostics$case_id==job$case_id,]
    by<-function(j,ep) {x<-cc[cc$dose==j&cc$endpoint==ep,];truth[truth$n==x$n&truth$s==x$s,]}
    gate<-vapply(1:5,function(j) by(j,'T')$T_gate>=.3&&(job$stage=='Final'||by(j,'E')$E_gate>=.3),TRUE)
    q<-vapply(1:5,function(j) by(j,'E')$mean,0)
    action<-if(job$stage=='Final') paste0('dose:',if(any(gate)) which(gate)[which.max(q[gate])] else 0L) else mask_label(gate,1:5)
    for(mode in c('uniform','decision')) {
      z<-diagnostic_intervals(dr,mode);labels<-interval_actions(job$stage,z$gates,z$scores)
      k<-k+1L;rows[[k]]<-data.frame(case_id=job$case_id,method=mode,exact_action=action,actions=paste(sort(labels),collapse=';'),
        miss=as.integer(!action%in%labels),singleton=as.integer(length(labels)==1L),wrong_singleton=as.integer(length(labels)==1L&&!action%in%labels))
    }
  }
  rr<-do.call(rbind,rows);e<-read.csv(file.path(path,'ACTION_SETS.csv'),stringsAsFactors=FALSE)
  for(n in names(rr)) stopifnot(identical(as.character(rr[[n]]),as.character(e[[n]])))
  write.csv(rr,file.path(out,'INTERCEPT_ACTIONS.csv'),row.names=FALSE)
  list(status='PASS',histories=nrow(jobs),quadrature_values=nrow(truth),decision_singletons=sum(rr$singleton[rr$method=='decision']),uniform_singletons=sum(rr$singleton[rr$method=='uniform']),misses=sum(rr$miss))
}
audit_forest_counts <- function(data,out) {
  root<-file.path(data,'base/compact_certificate/original_reports');dirs<-list.dirs(root,recursive=FALSE)
  rows<-list();k<-0L
  for(path in dirs) {
    records<-read_json(file.path(path,'ACTION_SETS.json'))
    for(r in records) {
      ref<-unlist(r$reference);native<-unlist(r$native)
      bounds<-c(as.integer(!any(ref%in%native)),as.integer(!all(ref%in%native)))
      wrong<-if(length(native)==1L) bounds else c(0L,0L)
      stopifnot(identical(as.numeric(bounds),as.numeric(c(r$miss_lower,r$miss_upper))),
        identical(as.numeric(wrong),as.numeric(c(r$wrong_singleton_lower,r$wrong_singleton_upper))))
      k<-k+1L;rows[[k]]<-data.frame(bank=basename(path),case_id=r$case_id,policy=r$policy,stage=r$stage,mode=r$mode,
         reference_singleton=as.integer(length(ref)==1L),native_singleton=as.integer(length(native)==1L),miss_lower=bounds[1],miss_upper=bounds[2],wrong_lower=wrong[1],wrong_upper=wrong[2])
    }
  }
  rr<-do.call(rbind,rows);write.csv(rr,file.path(out,'FOREST_ACCURACY.csv'),row.names=FALSE)
  list(status='PASS',stage_method_rows=nrow(rr),banks=length(dirs))
}
audit_workload <- function(path) {
  result<-read_json(file.path(path,'RESULTS.json'));rr<-result$endpoint_results
  for(r in rr) {
    q<-readBin(file.path(path,paste0(r$method,'_',r$endpoint,r$dose,'.bin')),'double',n=r$accepted,size=8,endian='little')
    b<-probability_intervals(q,r$endpoint,.01/20)
    for(n in names(b)) stopifnot(abs(b[[n]]-r[[n]])<1e-10)
  }
  for(m in result$methods) {
    rows<-if(m$method=='native_four_chains') m$intervals else Filter(function(r) r$method==m$method,rr)
    for(s in c('IA1','Final')) stopifnot(setequal(work_actions(rows,s),unlist(m$actions[[s]])))
  }
  list(status='PASS',sample_series=length(rr),methods=length(result$methods))
}
audit_sampler <- function() {

  law<-catalog_law(matrix(0,1L,1L),list(.5),30,exact=TRUE);a<-prepare_reference(law,1L)
  logf<-function(x) dnorm(x,sd=1.5,log=TRUE)+12*pnorm(x,log.p=TRUE)+18*pnorm(-x,log.p=TRUE)
  shift<-optimize(logf,c(-10,10),maximum=TRUE)$objective;f<-function(x)exp(logf(x)-shift)
  norm<-integrate(f,-Inf,Inf,rel.tol=1e-10)$value
  truth<-integrate(function(x) pnorm(x)*f(x),-Inf,Inf,rel.tol=1e-10)$value/norm
  means<-numeric()
  for(method in c('direct','aggregated')) {
    set.seed(171101L+as.integer(method=='aggregated'))
    z<-sample_reference(50000L,2048L,50L,a,12,30,1,0,method);q<-z$q
    stopifnot(length(q)>100L,abs(mean(q)-truth)<6*sd(q)/sqrt(length(q)))
    means[method]<-mean(q)
  }

  gates<-lapply(1:3,function(j)rbind(c(.7,.9),c(.1,.2)))
  z<-prepare_certificate(list(stage='Final',target='dose:2',gates=gates,scores=list(c(.1,.2),c(.7,.8),c(.4,.5))))
  costs<-setNames(c(1,5,2,3,7,1),z$atoms)
  m<-minimum_certificate(z,costs);e<-exhaustive_certificate(z,costs);stopifnot(m$cost==e$cost)

  z<-prepare_certificate(list(stage='Final',target='dose:1',gates=gates,scores=list(c(1,1),c(0,1),c(0,1))))
  stopifnot(minimum_certificate(z)$cost==exhaustive_certificate(z)$cost)
  regions<-path_regions(c(900,900),c(.3,.2),c(900,900),c(.2,.3));p<-earliest_path(regions[['2']])
  stopifnot(p$t==gmp::as.bigq(1,2),!p$closed,earliest_path(regions[['2']],2L)$t==1)
  list(status='PASS',stump_quadrature_mean=truth,posterior_means=as.list(means),weighted_subset_oracle=TRUE,endpoint_ties=TRUE,affine_path_boundary=TRUE)
}
