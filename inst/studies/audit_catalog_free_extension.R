
args <- commandArgs(trailingOnly=TRUE);stopifnot(length(args)%in%c(2L,3L))
root <- normalizePath(args[1]);out <- normalizePath(args[2])
report <- if(length(args)==3L)args[3]else out
dir.create(report,recursive=TRUE,showWarnings=FALSE);report <- normalizePath(report)
for(f in c('reference.R','decisions.R','prior_trees.R'))source(file.path(root,'R',f))
source_audit <- new.env(parent=baseenv())
sys.source(file.path(root,'R/audits.R'),envir=source_audit)
readj <- function(p)jsonlite::read_json(p,simplifyVector=FALSE)
plan <- readj(file.path(out,'PLAN.json'));result <- readj(file.path(out,'RESULTS.json'))
stopifnot(result$status=='COMPLETE',identical(result$plan_sha256,digest::digest(file=file.path(out,'PLAN.json'),algo='sha256')))
for(f in c('reference.R','decisions.R','prior_trees.R'))stopifnot(source_audit$source_hash_matches(file.path(root,'R',f),plan$sources[[f]],root))
old <- read.csv(file.path(root,'data/R_feasibility/FEASIBILITY.csv'),stringsAsFactors=FALSE)
designs <- readRDS(file.path(out,'DESIGNS.rds'))
stopifnot(setequal(names(designs),old$id[old$status=='RESOURCE_LIMIT']))
for(d in designs) {
  H <- nrow(d$X);set.seed(611081000L+H)
  X <- matrix(runif(H*5L),H,5L);X[,1] <- sort(X[,1]);X[,4] <- rep(0:1,length.out=H)
  colnames(X) <- paste0('x',1:5)
  stopifnot(identical(X,d$X))
  rr <- old[old$id==d$id,,drop=FALSE]
  cuts <- lapply(seq_len(ncol(X)),function(v) {
    xs <- sort(unique(X[,v]));gaps <- (xs[-1]+xs[-length(xs)])/2
    if(rr$density=='coarse'&&length(gaps)>3L)gaps <- gaps[unique(round(seq(1,length(gaps),length.out=3L)))]
    gaps
  })
  stopifnot(identical(cuts,d$cuts),identical(d$counts,rep(if(rr$kind=='repeated')5L else 1L,H)))
}
pr <- read.csv(file.path(out,'PRIOR_FEASIBILITY.csv'),stringsAsFactors=FALSE)
for(i in seq_along(designs)) {
  d <- designs[[i]];a <- prepare_tree_prior(d$X,d$cuts,d$counts)
  set.seed(plan$seed_base+i);z <- draw_tree_attempts(plan$prior_tree_attempts,a)
  r <- pr[pr$design==d$id,,drop=FALSE]
  stopifnot(nrow(r)==1L,r$valid_trees==length(z$valid),r$tree_attempts==z$attempts,
    abs(r$valid_fraction-length(z$valid)/z$attempts)<1e-12)
  for(k in z$valid[seq_len(min(200L,length(z$valid)))]) {
    masks <- z$leaf_mask[z$leaf_root==k]
    stopifnot(all(tree_training_counts(masks,a)>=5),sum(masks)==a$root_mask)
  }
}
summaries <- readj(file.path(out,'REVIEW_SUMMARY.json'));checks <- 0L
truth <- list(separated=list(T=c(.08,.13,.22,.32,.62),E=c(.12,.70,.32,.22,.08)),
              `near-boundary`=list(T=c(.44,.47,.48,.49,.52),E=c(.38,.40,.41,.40,.39)))
for(i in seq_along(summaries)) {
  summary <- summaries[[i]];folder <- file.path(out,sprintf('review%02d',i))
  case <- readRDS(file.path(folder,'CASE.rds'));rows <- readj(file.path(folder,'ENDPOINTS.json'))
  stopifnot(case$id==summary$design,case$scenario==summary$scenario,
    identical(case$X,designs[[case$id]]$X),identical(case$cuts,designs[[case$id]]$cuts))
  set.seed(plan$seed_base+100L+i)
  zprofile <- list(T=as.numeric(scale(case$X[,1])),E=as.numeric(scale(case$X[,2])))
  index <- rep(seq_len(nrow(case$X)),case$counts)
  for(j in 1:5)for(ep in c('T','E')) {
    generated <- rbinom(length(index),1,pnorm(qnorm(truth[[case$scenario]][[ep]][j])+.35*zprofile[[ep]][index]))
    stored <- case$dat[case$dat$dose==j,ep]
    stopifnot(identical(as.numeric(generated),as.numeric(stored)))
  }
  for(u in case$units) {
    z <- readRDS(file.path(folder,paste0(u$endpoint,u$dose,'.rds')))
    stopifnot(length(z$q)==length(z$accepted_at),all(z$q>=0&z$q<=1),
      !is.unsorted(z$accepted_at,strictly=TRUE),all(z$accepted_at<=plan$proposals_per_endpoint),
      all(z$accepted_at>=1L))
    expected <- probability_intervals(z$q,u$endpoint,plan$scalar_error)
    stored <- rows[[which(vapply(rows,function(r)r$dose==u$dose&&r$endpoint==u$endpoint,TRUE))]]
    for(n in names(expected)) {
      if(length(expected[[n]])==1L&&is.na(expected[[n]]))stopifnot(is.null(stored[[n]]))
      else stopifnot(isTRUE(all.equal(expected[[n]],stored[[n]],tolerance=1e-12,check.attributes=FALSE)))
    }
    checks <- checks+1L
  }
  for(stage in c('IA1','Final'))stopifnot(identical(work_actions(rows,stage),unlist(summary$actions[[stage]],use.names=FALSE)))
}
frozen <- file.path(report,'frozen_sources');dir.create(frozen,showWarnings=FALSE)
generation <- file.path(dirname(sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE))),'catalog_free_extension.R')
stopifnot(source_audit$source_hash_matches(generation,plan$sources[['catalog_free_extension.R']],root))
file.copy(c(generation,file.path(root,'R',c('reference.R','decisions.R','prior_trees.R'))),frozen,overwrite=TRUE)
answer <- list(status='PASS',prior_designs=nrow(pr),replayed_tree_attempts=sum(pr$tree_attempts),
  source_verification='Original hashes or recorded comment-free source transformations',
  reviews=length(summaries),endpoint_reports=checks,
  scalar_error=plan$scalar_error,total_quantity_error=checks*2*plan$scalar_error,
  saved_draw_intervals_reproduced=TRUE,action_sets_reproduced=TRUE,
  scope='Audits fixed-design feasibility, saved posterior draws and reported reference sets; no claim about clinical performance')
jsonlite::write_json(answer,file.path(report,'AUDIT.json'),auto_unbox=TRUE,pretty=TRUE,digits=17)
cat('Catalog-free extension audit PASS\n')
