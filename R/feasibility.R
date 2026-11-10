

catalog_law_limited <- function(X,cuts,counts,minimum=5L,max_seconds=30,
                                max_states=100000L,max_partitions=10000L,max_products=2000000L) {
  X<-as.matrix(X);counts<-as.numeric(counts);H<-nrow(X);p<-ncol(X)
  cuts<-lapply(cuts,as.numeric)
  stopifnot(H>=1L,H<=30L,minimum>=1L,all(counts>=0),all(counts==floor(counts)),
            sum(counts)>=minimum,length(counts)==H,length(cuts)==p,all(lengths(cuts)>0),
            all(is.finite(X)),all(vapply(cuts,function(v)all(is.finite(v))&&!is.unsorted(v,strictly=TRUE),TRUE)))
  start<-unname(proc.time()[3]);states<-0L;products<-0;bits<-bitwShiftL(1L,seq_len(H)-1L)
  memo<-new.env(hash=TRUE,parent=emptyenv())
  limit<-function(reason)stop(structure(list(message=reason,call=NULL),class=c('catalog_limit','error','condition')))
  check<-function()if(unname(proc.time()[3])-start>max_seconds)limit('elapsed-time limit')
  visit<-function(mask,lo,hi,depth) {
    tag<-paste(c(mask,lo,hi,depth),collapse=':')
    if(exists(tag,memo,inherits=FALSE))return(get(tag,memo,inherits=FALSE))
    check();states<<-states+1L;if(states>max_states)limit('node-state limit')
    members<-which(bitwAnd(mask,bits)!=0L)
    if(sum(counts[members])<minimum)return(list())
    good<-which(lo<=hi);grow<-if(!length(good))0 else .95/(depth+1)^2
    out<-setNames(list(1-grow),as.character(mask))
    for(v in good)for(cc in seq.int(lo[v],hi[v])) {
      check();left<-members[X[members,v]<cuts[[v]][cc]];right<-setdiff(members,left)
      if(min(sum(counts[left]),sum(counts[right]))<minimum)next
      lm<-sum(bits[left]);rm<-sum(bits[right]);lh<-hi;lh[v]<-cc-1L;rl<-lo;rl[v]<-cc+1L
      a<-visit(lm,lo,lh,depth+1L);b<-visit(rm,rl,hi,depth+1L)
      products<<-products+length(a)*length(b);if(products>max_products)limit('partition-product limit')
      mult<-grow/(p*(hi[v]-lo[v]+1L))
      for(ka in names(a))for(kb in names(b)) {
        key<-paste(sort(as.integer(strsplit(paste(ka,kb,sep=','),',',fixed=TRUE)[[1L]])),collapse=',')
        if(is.null(out[[key]]))out[[key]]<-0
        out[[key]]<-out[[key]]+mult*a[[ka]]*b[[kb]]
        if(length(out)>max_partitions)limit('partition-count limit')
      }
    }
    assign(tag,out,memo);out
  }
  ans<-tryCatch({
    raw<-visit(2^H-1L,rep(1L,p),lengths(cuts),0L);normalizer<-Reduce(`+`,raw,init=0)
    stopifnot(normalizer>0)
    list(status='COMPLETE',law=list(probabilities=lapply(raw,function(v)v/normalizer),normalizer=normalizer,
                                  regions=length(ls(memo))))
  },catalog_limit=function(e)list(status='RESOURCE_LIMIT',reason=conditionMessage(e)))
  ans$seconds<-unname(proc.time()[3])-start;ans$states<-states;ans$products<-products
  ans$memo_bytes<-as.numeric(object.size(as.list(memo)))
  ans
}

study_diagnostics <- function(chains,delta) {
  safe<-function(f,x)tryCatch(suppressWarnings(f(x)),error=function(e)NA_real_)
  rows<-list();k<-0L
  for(type in c('T_mean','E_mean','T_gate','E_gate'))for(j in 1:5) {
    x<-sapply(chains,function(v)switch(type,T_mean=v$T[,j],E_mean=v$E[,j],
                                      T_gate=as.numeric(v$T[,j]<.45),E_gate=as.numeric(v$E[,j]>.30)))
    mn<-mean(x);rh<-safe(posterior::rhat,x);ess<-safe(posterior::ess_bulk,x);mc<-safe(posterior::mcse_mean,x)
    l<-safe(posterior::mcse_mean,x[,1:2]);r<-safe(posterior::mcse_mean,x[,3:4]);se<-sqrt(l*l+r*r)
    agree<-if(is.finite(se)&&se>0)abs(mean(x[,1:2])-mean(x[,3:4]))/se else NA_real_
    ok<-min(x)!=max(x)&&all(is.finite(c(rh,ess,mc,agree)))&&rh<=1.05&&ess>=100&&mc>0&&agree<=6
    mult<-qnorm(1-delta/2);k<-k+1L
    rows[[k]]<-data.frame(series=type,dose=j,mean=mn,rhat=rh,ess=ess,mcse=mc,paired_z=agree,qualified=ok,
      precision_filtered=ok&&mc<=.002,lower=if(ok)max(0,mn-mult*mc)else 0,upper=if(ok)min(1,mn+mult*mc)else 1)
  }
  do.call(rbind,rows)
}

rows_from_diagnostics <- function(dr) {
  rows<-list();k<-0L
  for(j in 1:5)for(ep in c('T','E')) {
    mn<-dr[dr$dose==j&dr$series==paste0(ep,'_mean'),];g<-dr[dr$dose==j&dr$series==paste0(ep,'_gate'),]
    k<-k+1L;rows[[k]]<-list(dose=j,endpoint=ep,mean=mn$mean,mean_lower=mn$lower,mean_upper=mn$upper,
                           gate_lower=g$lower,gate_upper=g$upper)
  }
  rows
}
