

options(warn=1)
a<-commandArgs(trailingOnly=TRUE);stopifnot(length(a)==1L,!file.exists(a[1]))
script<-sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE))
package<-normalizePath(file.path(dirname(script),'..'))
root<-a[1];dir.create(root,recursive=TRUE);dir.create(file.path(root,'raw'))
p<-jsonlite::fromJSON(file.path(package,'data/extended/exact_reference_benchmark_20261103/PLAN.json'))
targets<-jsonlite::fromJSON(file.path(package,'data/PUBLISHED_TARGETS.json'))
RNGkind('Mersenne-Twister','Inversion','Rejection');set.seed(as.integer(p$data_seed))
jobs<-list();case<-0L;counts<-list()
for(profile in 1:5)for(stage in c('IA1','Final'))for(rep in 1:30) {
 case<-case+1L;n<-if(stage=='IA1')12L else 18L
 st<-rbinom(5,n,targets$toxicity[profile,]);se<-rbinom(5,n,targets$efficacy[profile,])
 jobs[[case]]<-data.frame(case_id=case,profile=profile,stage=stage,n=n)
 counts[[case]]<-data.frame(case_id=case,endpoint=rep(c('T','E'),each=5),dose=rep(1:5,2),n=n,s=c(st,se))
}
jobs<-do.call(rbind,jobs);cc<-do.call(rbind,counts)
write.csv(jobs,file.path(root,'INPUTS.csv'),row.names=FALSE);write.csv(cc,file.path(root,'COUNTS.csv'),row.names=FALSE)
saveRDS(.Random.seed,file.path(root,'DATA_RNG_FINAL.rds'))
prior_var<-2.25
exact<-function(n,s) {
 logf<-function(t)dnorm(t,0,sqrt(prior_var),log=TRUE)+(if(s>0)s*pnorm(t,log.p=TRUE)else 0)+(if(s<n)(n-s)*pnorm(t,lower.tail=FALSE,log.p=TRUE)else 0)
 shift<-optimize(logf,c(-10,10),maximum=TRUE)$objective
 f<-function(t)exp(logf(t)-shift)
 integral<-function(fun,l=-Inf,u=Inf)integrate(fun,l,u,rel.tol=1e-11,abs.tol=1e-12,subdivisions=1000L)$value
 norm<-integral(f)
 c(mean=integral(function(t)pnorm(t)*f(t))/norm,T_gate=integral(f,u=qnorm(.45))/norm,E_gate=integral(f,l=qnorm(.3))/norm,normalizer_scaled=norm,log_shift=shift)
}
truth<-list();ix<-0L
for(n in c(12L,18L))for(s in 0:n){ix<-ix+1L;truth[[ix]]<-as.data.frame(as.list(c(n=n,s=s,exact(n,s))))}
write.csv(do.call(rbind,truth),file.path(root,'EXACT_POSTERIORS.csv'),row.names=FALSE)
N<-nrow(cc);draws<-2000L;burn<-1000L;states<-vector('list',4L);warnings<-character()
for(chain in 1:4) {
 set.seed(as.integer(p$chain_seed_base+chain));f<-rep(0,N);keep<-matrix(NA_real_,draws,N);mask<-outer(cc$n,1:18,`>=`);ones<-outer(cc$s,1:18,`>=`)
 t0<-proc.time()
 for(iter in seq_len(burn+draws)) {
  u<-matrix(runif(N*18L),N,18L)

  prob<-u*ifelse(ones,pnorm(f),pnorm(-f));stopifnot(all(prob>0&prob<1))
  z<-matrix(qnorm(as.vector(prob),lower.tail=TRUE),N,18L)
  z[ones]<--z[ones];z<-z+f
  mean_f<-rowSums(z*mask)/(cc$n+1/prior_var)
  f<-rnorm(N,mean_f,1/sqrt(cc$n+1/prior_var));stopifnot(all(is.finite(f)))
  if(iter>burn)keep[iter-burn,]<-pnorm(f)
 }
 path<-file.path(root,'raw',paste0('chain',chain,'.rds'));saveRDS(list(chain=chain,probabilities=keep,counts=cc,initial_state=rep(0,N),final_rng=.Random.seed,cpu=unname(sum((proc.time()-t0)[1:2])),elapsed=unname((proc.time()-t0)[3])),path)
 states[[chain]]<-keep;cat('chain',chain,'complete\n');flush.console()
}
allrows<-list();ix<-0L;metadata<-list();stat<-function(fun,v)tryCatch(suppressWarnings(fun(v)),error=function(e)NA_real_)
for(case in jobs$case_id) {
 ser<-array(NA_real_,c(draws,4L,20L));si<-0L
 for(type in c('T_mean','E_mean','T_gate','E_gate'))for(dose in 1:5) {
  endpoint<-substr(type,1,1);col<-which(cc$case_id==case&cc$endpoint==endpoint&cc$dose==dose)
  x<-sapply(states,function(z)z[,col]);if(type=='T_gate')x<-1*(x<.45);if(type=='E_gate')x<-1*(x>.3)
  si<-si+1L;ser[,,si]<-x;rh<-stat(posterior::rhat,x);ess<-stat(posterior::ess_bulk,x);mc<-stat(posterior::mcse_mean,x)
  se<-sqrt(stat(posterior::mcse_mean,x[,1:2])^2+stat(posterior::mcse_mean,x[,3:4])^2)
  paired<-if(is.finite(se)&&se>0)abs(mean(x[,1:2])-mean(x[,3:4]))/se else NA_real_
  constant<-min(x)==max(x);qualified<-!constant&&all(is.finite(c(rh,ess,mc,paired)))&&rh<=1.05&&ess>=100&&mc>0&&paired<=6
  zcrit<-qnorm(1-.05/40)
  ix<-ix+1L;allrows[[ix]]<-data.frame(case_id=case,series=type,dose=dose,mean=mean(x),rhat=rh,ess=ess,mcse=mc,constant=constant,paired_z=paired,qualified=qualified,precision_filtered=qualified&&mc<=.002,lower=if(qualified)max(0,mean(x)-zcrit*mc)else 0,upper=if(qualified)min(1,mean(x)+zcrit*mc)else 1)
 }
 path<-file.path(root,'raw',sprintf('case%04d_series.bin',case));writeBin(as.double(ser),path,size=8,endian='little')
 metadata[[as.character(case)]]<-list(path=path,sha256=digest::digest(file=path,algo='sha256'),dimensions=dim(ser))
}
write.csv(do.call(rbind,allrows),file.path(root,'DIAGNOSTICS.csv'),row.names=FALSE,na='')
jsonlite::write_json(metadata,file.path(root,'SERIES_META.json'),auto_unbox=TRUE,pretty=TRUE)
jsonlite::write_json(list(status='COMPLETE',histories=300L,independent_chains=4L,burn=1000L,retained=2000L,series=ix,scope='Proper finite-offset stump-only probit model; not the native forest or an adaptive trial'),file.path(root,'R_COMPLETE.json'),auto_unbox=TRUE,pretty=TRUE)
