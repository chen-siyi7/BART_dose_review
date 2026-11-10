

a<-commandArgs(trailingOnly=TRUE);stopifnot(length(a)==3L)
script<-sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE));root<-normalizePath(file.path(dirname(script),'..'))
case<-as.integer(a[1]);role<-as.integer(a[2]);stem<-a[3]
stopifnot(case%in%1:900,role%in%1:2,!file.exists(paste0(stem,'.json')),!file.exists(paste0(stem,'.rds')))
lib<-Sys.getenv('BART_LIBRARY',unset='');if(nzchar(lib)) .libPaths(c(normalizePath(lib),.libPaths()))
suppressPackageStartupMessages(library(BART));stopifnot(as.character(packageVersion('BART'))=='2.9.9')
path<-file.path(root,'data/extended/comparative_reference_20261103');p<-jsonlite::fromJSON(file.path(path,'PLAN.json'))
jobs<-read.csv(file.path(path,'INPUTS.csv'));job<-jobs[jobs$case_id==case,]
input<-file.path(path,'inputs',basename(job$input_path));stopifnot(digest::digest(file=input,algo='sha256')==job$input_sha256)
z<-readRDS(input);seed<-as.integer(p$fit_seed_base+case*8L+role)
RNGkind('Mersenne-Twister','Inversion','Rejection');set.seed(seed);rng0<-.Random.seed
warnings<-character();start<-proc.time()[3]
fit<-tryCatch(withCallingHandlers({
  dat<-z$dat;X<-as.matrix(dat[paste0('x',1:5)]);T<-E<-matrix(NA_real_,1000L,5L)
  one<-function(x,y) {
    invisible(capture.output(f<-BART::pbart(x.train=x,y.train=y)))
    stopifnot(nrow(f$prob.train)==1000L,all(is.finite(f$prob.train)));f$prob.train
  }
  if(z$layout=='separate_dose_X') for(j in 1:5) {
    k<-which(dat$dose==j);T[,j]<-rowMeans(one(X[k,,drop=FALSE],dat$T[k]));E[,j]<-rowMeans(one(X[k,,drop=FALSE],dat$E[k]))
  } else {
    if(z$layout=='pooled_logdose_X') X<-cbind(logdose=log(c(10,20,35,55,70))[dat$dose],X)
    tp<-one(X,dat$T);ep<-one(X,dat$E)
    for(j in 1:5) {k<-which(dat$dose==j);T[,j]<-rowMeans(tp[,k,drop=FALSE]);E[,j]<-rowMeans(ep[,k,drop=FALSE])}
  };list(T=T,E=E)
},warning=function(w) {warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')}),error=function(e)list(error=conditionMessage(e)))
out<-list(case_id=case,role=role,seed=seed,status=if(is.null(fit$error))'valid' else 'failed_fit',error=fit$error,
  stage=z$stage,layout=z$layout,candidates=z$candidates,warnings=warnings,initial_rng=rng0,final_rng=.Random.seed,elapsed=unname(proc.time()[3]-start))
if(is.null(fit$error)) {
  out$kT<-colSums(fit$T<.45);out$kE<-colSums(fit$E>.3);out$q<-colMeans(fit$E);out$pT<-colMeans(fit$T)
  keep<-out$kT>=300;if(z$stage!='Final')keep<-keep&out$kE>=300
  out$action<-if(z$stage=='Final') paste0('dose:',if(any(keep)) which(keep)[which.max(out$q[keep])] else 0L) else paste0('mask:',paste(which(keep),collapse=':'))
}
saveRDS(c(out,list(fit=fit)),paste0(stem,'.rds'))
jsonlite::write_json(out,paste0(stem,'.json'),auto_unbox=TRUE,pretty=TRUE,digits=17,na='null',null='null')
cat(out$status,'\n')
