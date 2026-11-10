a<-commandArgs(trailingOnly=TRUE);stopifnot(length(a)==2L);inputdir<-a[1];out<-a[2]
obj<-lapply(1:4,function(i)readRDS(file.path(inputdir,paste0('chain',i,'.rds'))))
safe<-function(f,x)tryCatch(suppressWarnings(f(x)),error=function(e)NA_real_)
rows<-list();k<-0L;delta<-.01/20
for(type in c('T_mean','E_mean','T_gate','E_gate'))for(j in 1:5){
 x<-sapply(obj,function(v)switch(type,T_mean=v$T[,j],E_mean=v$E[,j],T_gate=as.numeric(v$T[,j]<.45),E_gate=as.numeric(v$E[,j]>.30)))
 mn<-mean(x);rh<-safe(posterior::rhat,x);ess<-safe(posterior::ess_bulk,x);mc<-safe(posterior::mcse_mean,x)
 l<-safe(posterior::mcse_mean,x[,1:2]);r<-safe(posterior::mcse_mean,x[,3:4]);se<-sqrt(l*l+r*r)
 paired<-if(is.finite(se)&&se>0)abs(mean(x[,1:2])-mean(x[,3:4]))/se else NA_real_
 ok<-min(x)!=max(x)&&all(is.finite(c(rh,ess,mc,paired)))&&rh<=1.05&&ess>=100&&mc>0&&paired<=6
 mult<-qnorm(1-delta/2);k<-k+1L
 rows[[k]]<-data.frame(series=type,dose=j,mean=mn,rhat=rh,ess=ess,mcse=mc,paired_z=paired,qualified=ok,
                      lower=if(ok)max(0,mn-mult*mc)else 0,upper=if(ok)min(1,mn+mult*mc)else 1)
}
write.csv(do.call(rbind,rows),out,row.names=FALSE,na='')
