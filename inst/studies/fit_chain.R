options(warn=1)
a<-commandArgs(trailingOnly=TRUE);stopifnot(length(a)==3L)
input<-a[1];out<-a[2];role<-as.integer(a[3]);stopifnot(role%in%1:4)
lib<-Sys.getenv('BART_LIBRARY',unset='');if(nzchar(lib)).libPaths(c(normalizePath(lib),.libPaths()))
suppressPackageStartupMessages(library(BART));stopifnot(as.character(packageVersion('BART'))=='2.9.9')
z<-readRDS(input);RNGkind('Mersenne-Twister','Inversion','Rejection');set.seed(2107111000L+role)
T<-E<-matrix(NA_real_,32000L,5L);started<-proc.time()[3];warnings<-character()
withCallingHandlers({
for(j in 1:5)for(ep in c('T','E')){
 d<-z$dat[z$dat$dose==j,,drop=FALSE]
 invisible(capture.output(f<-BART::pbart(x.train=as.matrix(d[paste0('x',1:5)]),y.train=d[[ep]],
  ntree=50L,k=2,base=.95,power=2,sparse=FALSE,augment=FALSE,binaryOffset=0,
  nskip=8000L,ndpost=32000L,nkeeptrain=32000L,nkeeptest=0L,nkeeptreedraws=0L)))
 stopifnot(identical(dim(f$prob.train),c(32000L,30L)),all(is.finite(f$prob.train)))
 if(ep=='T')T[,j]<-rowMeans(f$prob.train)else E[,j]<-rowMeans(f$prob.train)
}
},warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')})
saveRDS(list(role=role,seed=2107111000L+role,T=T,E=E,elapsed=unname(proc.time()[3]-started),warnings=warnings),out)
