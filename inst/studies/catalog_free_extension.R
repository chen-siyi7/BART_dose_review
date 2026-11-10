

argv <- commandArgs(trailingOnly=TRUE)
stopifnot(length(argv)==2L,dir.exists(argv[1]),!file.exists(argv[2]))
root <- normalizePath(argv[1]); out <- argv[2]
dir.create(out,recursive=TRUE); out <- normalizePath(out)
for(f in c('reference.R','decisions.R','prior_trees.R')) source(file.path(root,'R',f))
savej <- function(x,p) jsonlite::write_json(x,p,auto_unbox=TRUE,pretty=TRUE,digits=17,na='null',null='null')
sha <- function(p) digest::digest(file=p,algo='sha256')
elapsed <- function() unname(proc.time()[3])
old <- read.csv(file.path(root,'data/R_feasibility/FEASIBILITY.csv'),stringsAsFactors=FALSE)
failed <- old$id[old$status=='RESOURCE_LIMIT']
plan <- list(status='FROZEN_BEFORE_NEW_GENERATION_AND_MEASUREMENT',
  time_UTC=format(Sys.time(),tz='UTC'),
  purpose='Descriptive test of complete-tree references beyond completed partition catalogs',
  prior_designs=failed,prior_tree_attempts=40000L,
  review_designs=c('repeated_12_fine','continuous_30_fine'),
  review_scenarios=c('separated','near-boundary'),
  selection='Largest repeated-profile and distinct-profile designs in the previous fourteen-design grid; both previously hit catalog resource limits',
  patient_generation='Same fixed profile/cut recipes and response templates as the previous feasibility study; newly generated independent outcomes',
  trees=50L,proposals_per_endpoint=50000L,chunk=256L,
  review_family_error=.01,scalar_error=.01/80,
  minimum_leaf_count=5L,offset=0,
  seed_base=611092000L,
  timing='One sequential R run per fixed review, including tree setup, prior generation, likelihood rejection, interval construction and decision scan; no native fitting and no partition enumeration',
  scope='Descriptive fixed-history extension, not clinical validation or a population-level trial assessment; memory is prepared object size, not process peak RSS',
  sources=setNames(lapply(c(sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE)),file.path(root,'R',c('reference.R','decisions.R','prior_trees.R'))),sha),
                   c('catalog_free_extension.R','reference.R','decisions.R','prior_trees.R')),
  session=capture.output(sessionInfo()))
savej(plan,file.path(out,'PLAN.json'))
cat('Catalog-free extension plan frozen\n');flush.console()
make_X <- function(H) {
  set.seed(611081000L+H)
  x <- matrix(runif(H*5L),H,5L); x[,1] <- sort(x[,1]); x[,4] <- rep(0:1,length.out=H)
  colnames(x) <- paste0('x',1:5);x
}
make_cuts <- function(X,density) lapply(1:ncol(X),function(v) {
  xs <- sort(unique(X[,v])); gaps <- (xs[-1]+xs[-length(xs)])/2
  if(density=='coarse'&&length(gaps)>3L) gaps <- gaps[unique(round(seq(1,length(gaps),length.out=3L)))]
  gaps
})
make_design <- function(id) {
  r <- old[old$id==id,,drop=FALSE];stopifnot(nrow(r)==1L)
  X <- make_X(r$H);cuts <- make_cuts(X,r$density)
  counts <- rep(if(r$kind=='repeated')5L else 1L,r$H)
  list(id=id,X=X,cuts=cuts,counts=counts)
}
designs <- setNames(lapply(failed,make_design),failed)
saveRDS(designs,file.path(out,'DESIGNS.rds'))
prior <- list()
for(i in seq_along(designs)) {
  d <- designs[[i]];start <- elapsed();a <- prepare_tree_prior(d$X,d$cuts,d$counts)
  setup <- elapsed()-start;set.seed(plan$seed_base+i)
  start <- elapsed();z <- draw_tree_attempts(plan$prior_tree_attempts,a);seconds <- elapsed()-start
  stopifnot(length(z$valid)>0L)
  for(r in z$valid[seq_len(min(200L,length(z$valid)))]) {
    masks <- z$leaf_mask[z$leaf_root==r]
    stopifnot(all(tree_training_counts(masks,a)>=5),sum(masks)==a$root_mask)
  }
  prior[[i]] <- data.frame(design=d$id,H=a$H,patients_per_dose=sum(d$counts),
     tree_attempts=plan$prior_tree_attempts,valid_trees=length(z$valid),
     valid_fraction=length(z$valid)/plan$prior_tree_attempts,
     setup_seconds=setup,sampling_seconds=seconds,prepared_bytes=as.numeric(object.size(a)),
     checked_valid_trees=min(200L,length(z$valid)))
  write.csv(do.call(rbind,prior),file.path(out,'PRIOR_FEASIBILITY.csv'),row.names=FALSE)
  cat(d$id,'valid trees',length(z$valid),'/',plan$prior_tree_attempts,'\n');flush.console()
}
truth <- list(separated=list(T=c(.08,.13,.22,.32,.62),E=c(.12,.70,.32,.22,.08)),
              `near-boundary`=list(T=c(.44,.47,.48,.49,.52),E=c(.38,.40,.41,.40,.39)))
rows_all <- summaries <- list();review <- 0L
for(id in plan$review_designs) for(scenario in plan$review_scenarios) {
  review <- review+1L;d <- designs[[id]];X <- d$X;n <- d$counts;H <- nrow(X)
  set.seed(plan$seed_base+100L+review)
  z <- list(T=as.numeric(scale(X[,1])),E=as.numeric(scale(X[,2])))
  dat <- units <- list();k <- 0L
  for(j in 1:5) {
    index <- rep(seq_len(H),n)
    dd <- data.frame(dose=j,profile=index,X[index,,drop=FALSE])
    for(ep in c('T','E')) dd[[ep]] <- rbinom(length(index),1,pnorm(qnorm(truth[[scenario]][[ep]][j])+.35*z[[ep]][index]))
    dat[[j]] <- dd
    for(ep in c('T','E')) {k <- k+1L;units[[k]] <- list(dose=j,endpoint=ep,successes=as.numeric(tapply(dd[[ep]],dd$profile,sum)))}
  }
  folder <- file.path(out,sprintf('review%02d',review));dir.create(folder)
  case <- c(d,list(review=review,scenario=scenario,dat=do.call(rbind,dat),units=units))
  saveRDS(case,file.path(folder,'CASE.rds'));write.csv(case$dat,file.path(folder,'PATIENTS.csv'),row.names=FALSE)
  start <- elapsed();a <- prepare_tree_prior(X,d$cuts,n);setup <- elapsed()-start
  rows <- list()
  for(u in seq_along(units)) {
    unit <- units[[u]];set.seed(plan$seed_base+1000L+review*100L+u)
    draw <- sample_reference_checkpoints(plan$proposals_per_endpoint,plan$chunk,50L,a,
      unit$successes,n,n/sum(n),0,'tree_prior',plan$scalar_error,unit$endpoint)
    band <- draw$reports[[1L]]
    rows[[u]] <- c(list(review=review,design=id,scenario=scenario,dose=unit$dose,endpoint=unit$endpoint),band,
      list(tree_attempts=draw$tree_attempts,valid_tree_attempts=draw$valid_tree_attempts,visited_nodes=draw$visited_nodes))
    saveRDS(draw,file.path(folder,paste0(unit$endpoint,unit$dose,'.rds')))
    savej(rows,file.path(folder,'ENDPOINTS.json'))
    cat('Review',review,id,scenario,paste0(unit$endpoint,unit$dose),band$accepted,'accepted;',round(band$seconds,2),'s\n');flush.console()
  }
  start <- elapsed();actions <- setNames(lapply(c('IA1','Final'),function(s)work_actions(rows,s)),c('IA1','Final'))
  seconds <- setup+sum(vapply(rows,function(r)r$seconds,0))+elapsed()-start
  summaries[[review]] <- list(review=review,design=id,scenario=scenario,H=H,patients_per_dose=sum(n),
    setup_seconds=setup,seconds=seconds,prepared_bytes=as.numeric(object.size(a)),
    accepted_min=min(vapply(rows,function(r)r$accepted,0)),accepted_max=max(vapply(rows,function(r)r$accepted,0)),
    actions=actions)
  rows_all <- c(rows_all,rows)
  savej(summaries,file.path(out,'REVIEW_SUMMARY.json'))
  savej(rows_all,file.path(out,'ENDPOINTS.json'))
}
savej(list(status='COMPLETE',plan_sha256=sha(file.path(out,'PLAN.json')),
  prior_feasibility=do.call(rbind,prior),reviews=summaries,
  scope=plan$scope),file.path(out,'RESULTS.json'))
cat('Catalog-free extension complete\n')
