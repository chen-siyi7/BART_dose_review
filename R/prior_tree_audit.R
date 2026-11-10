audit_prior_trees <- function(root,out) {
  cats<-readRDS(file.path(root,'data/R_feasibility/CATALOGS.rds'))
  records<-list();i<-0L
  for(id in c('repeated_6_fine','repeated_8_coarse','continuous_20_fine')) {
    c<-cats[[id]];a<-prepare_tree_prior(c$X,c$cuts,c$counts)
    set.seed(611090100L+i);z<-draw_tree_attempts(40000L,a)
    keys<-vapply(split(z$leaf_mask,z$leaf_root),function(v)paste(sort(v),collapse=','),'')
    w<-vapply(c$law$probabilities,as.numeric,0);obs<-setNames(numeric(length(w)),names(w))
    tab<-table(keys)/length(z$valid);stopifnot(all(names(tab)%in%names(w)))
    obs[names(tab)]<-as.numeric(tab)
    tolerance<-sqrt(log(2*length(w)/1e-8)/(2*length(z$valid)))
    stopifnot(max(abs(obs-w))<=tolerance,
      abs(length(z$valid)/40000-as.numeric(c$law$normalizer))<=sqrt(log(2/1e-8)/80000))
    for(r in z$valid[seq_len(min(200L,length(z$valid)))]) {
      masks<-z$leaf_mask[z$leaf_root==r]
      stopifnot(all(tree_training_counts(masks,a)>=5),sum(masks)==a$root_mask)
    }
    i<-i+1L;records[[i]]<-list(design=id,attempts=40000L,valid=length(z$valid),
      partitions=length(w),maximum_probability_error=max(abs(obs-w)),simultaneous_tolerance=tolerance)
  }
  c<-cats[['repeated_6_fine']];a<-prepare_tree_prior(c$X,c$cuts,c$counts)
  set.seed(611090200L);f<-prior_tree_fields(20000L,50L,a)
  d<-ncol(f);alpha<-1e-6;z<-qnorm(1-alpha/(4*d));v<-diag(var(f))
  stopifnot(all(abs(colMeans(f))<=z*1.5/sqrt(nrow(f))),
    all((nrow(f)-1)*v/2.25>=qchisq(alpha/(4*d),nrow(f)-1)),
    all((nrow(f)-1)*v/2.25<=qchisq(1-alpha/(4*d),nrow(f)-1)))


  a<-prepare_tree_prior(matrix(0,1,1),list(.5),12L)
  set.seed(611090201L)
  draw<-sample_reference_checkpoints(20000L,256L,2L,a,7L,12L,1,0,'tree_prior',.001,'E')
  logf<-function(x)dnorm(x,sd=1.5,log=TRUE)+7*pnorm(x,log.p=TRUE)+5*pnorm(-x,log.p=TRUE)
  shift<-optimize(logf,c(-10,10),maximum=TRUE)$objective
  den<-integrate(function(x)exp(logf(x)-shift),-12,12,rel.tol=1e-11)$value
  mu<-integrate(function(x)pnorm(x)*exp(logf(x)-shift),-12,12,rel.tol=1e-11)$value/den
  gate<-integrate(function(x)exp(logf(x)-shift),qnorm(.3),12,rel.tol=1e-11)$value/den
  band<-draw$reports[[1]]
  stopifnot(mu>=band$mean_lower,mu<=band$mean_upper,gate>=band$gate_lower,gate<=band$gate_upper)
  result<-list(status='PASS',partition_laws=records,Gaussian_marginals=TRUE,
    intercept_posterior=list(mean=mu,gate=gate,intervals=band),
    scope='Simulation checks supplement the complete-tree conditioning argument; no catalog is used by the sampler')
  jsonlite::write_json(result,file.path(out,'PRIOR_TREE_CHECKS.json'),pretty=TRUE,auto_unbox=TRUE,digits=17)
  result
}
