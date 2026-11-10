
argv<-commandArgs(trailingOnly=TRUE);stopifnot(length(argv)==3L)
lib<-Sys.getenv('BART_LIBRARY',unset='');if(nzchar(lib)).libPaths(c(normalizePath(lib),.libPaths()))
suppressPackageStartupMessages(library(BART));stopifnot(as.character(packageVersion('BART'))=='2.9.9')
z<-readRDS(argv[1]);role<-as.integer(argv[3]);stopifnot(role%in%1:4)
RNGkind('Mersenne-Twister','Inversion','Rejection');seed<-z$fit_seed+role;set.seed(seed)
xi<-matrix(NA_real_,nrow=5L,ncol=max(lengths(z$cuts)))
for(v in 1:5)xi[v,seq_along(z$cuts[[v]])]<-z$cuts[[v]]
T<-E<-matrix(NA_real_,32000L,5L);start<-unname(proc.time()[3]);warnings<-character()
withCallingHandlers({for(j in 1:5)for(ep in c('T','E')) {
  d<-z$dat[z$dat$dose==j,,drop=FALSE]
  invisible(capture.output(f<-BART::pbart(x.train=as.matrix(d[paste0('x',1:5)]),y.train=d[[ep]],
    xinfo=xi,rm.const=FALSE,ntree=50L,k=2,base=.95,power=2,sparse=FALSE,augment=FALSE,binaryOffset=0,
    nskip=8000L,ndpost=32000L,nkeeptrain=32000L,nkeeptest=0L,nkeeptreedraws=0L)))
  stopifnot(nrow(f$prob.train)==32000L,ncol(f$prob.train)==nrow(d),all(is.finite(f$prob.train)))
  if(ep=='T')T[,j]<-rowMeans(f$prob.train)else E[,j]<-rowMeans(f$prob.train)
}},warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')})
saveRDS(list(T=T,E=E,role=role,seed=seed,elapsed=unname(proc.time()[3])-start,warnings=warnings),argv[2])
