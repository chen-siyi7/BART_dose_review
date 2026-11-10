

catalog_law <- function(X, cuts, counts, minimum=5L, exact=FALSE) {
  X <- as.matrix(X); counts <- as.numeric(counts); H <- nrow(X); p <- ncol(X)
  cuts <- lapply(cuts, as.numeric)
  stopifnot(H>=1L,H<=12L,length(counts)==H,all(is.finite(X)),
            all(counts>=0),all(counts==floor(counts)),sum(counts)>=minimum,
            length(cuts)==p,all(lengths(cuts)>0L),
            all(vapply(cuts,function(x) all(is.finite(x))&&!is.unsorted(x,strictly=TRUE),TRUE)))
  zero <- if(exact) gmp::as.bigq(0) else 0
  one <- if(exact) gmp::as.bigq(1) else 1
  bits <- bitwShiftL(1L,seq_len(H)-1L)
  memo <- new.env(hash=TRUE,parent=emptyenv())
  visit <- function(mask,lo,hi,depth) {
    tag <- paste(c(mask,lo,hi,depth),collapse=':')
    if(exists(tag,memo,inherits=FALSE)) return(get(tag,memo,inherits=FALSE))
    members <- which(bitwAnd(mask,bits)!=0L)
    if(sum(counts[members])<minimum) return(list())
    good <- which(lo<=hi)
    grow <- if(!length(good)) zero else if(exact) gmp::as.bigq(19,20*(depth+1)^2) else .95/(depth+1)^2
    out <- setNames(list(one-grow),as.character(mask))
    for(v in good) for(cc in seq.int(lo[v],hi[v])) {
      left <- members[X[members,v]<cuts[[v]][cc]]; right <- setdiff(members,left)
      if(min(sum(counts[left]),sum(counts[right]))<minimum) next
      lm <- sum(bits[left]); rm <- sum(bits[right])
      lh <- hi; lh[v] <- cc-1L; rl <- lo; rl[v] <- cc+1L
      a <- visit(lm,lo,lh,depth+1L); b <- visit(rm,rl,hi,depth+1L)
      mult <- grow/(p*(hi[v]-lo[v]+1L))
      for(ka in names(a)) for(kb in names(b)) {
        key <- paste(sort(as.integer(strsplit(paste(ka,kb,sep=','),',',fixed=TRUE)[[1L]])),collapse=',')
        if(is.null(out[[key]])) out[[key]] <- zero
        out[[key]] <- out[[key]]+mult*a[[ka]]*b[[kb]]
      }
    }
    assign(tag,out,memo);out
  }
  raw <- visit(2^H-1L,rep(1L,p),lengths(cuts),0L)
  norm <- Reduce(`+`,raw,init=zero);stopifnot(norm>0)
  list(probabilities=lapply(raw,function(w) w/norm),normalizer=norm,regions=length(ls(memo)))
}

prepare_reference <- function(law,H) {
  parts <- lapply(names(law$probabilities),function(k) as.integer(strsplit(k,',',fixed=TRUE)[[1L]]))
  masks <- sort(unique(unlist(parts))); bits<-bitwShiftL(1L,seq_len(H)-1L)
  membership <- outer(masks,bits,function(a,b) as.numeric(bitwAnd(a,b)!=0L))
  incidence <- t(vapply(parts,function(x) as.numeric(masks%in%x),numeric(length(masks))))
  weights <- vapply(law$probabilities,as.numeric,0)
  stopifnot(abs(sum(weights)-1)<1e-12,all(weights>=0))
  list(parts=parts,masks=masks,membership=membership,incidence=incidence,
       weights=weights,H=H,lengths=lengths(parts))
}

prior_fields <- function(n,trees,a,method=c('aggregated','direct')) {
  method<-match.arg(method); tau<-3/(2*sqrt(trees))
  if(method=='aggregated') {
    counts <- t(rmultinom(n,trees,a$weights))


    leaf_counts <- counts %*% a$incidence
    z <- sqrt(leaf_counts)*matrix(rnorm(length(leaf_counts)),nrow=n)
    return(tau*(z %*% a$membership))
  }
  labels <- sample.int(length(a$parts),n*trees,replace=TRUE,prob=a$weights)
  masks <- unlist(a$parts[labels],use.names=FALSE)
  ids <- rep(rep(seq_len(n),each=trees),a$lengths[labels])
  z <- rnorm(length(masks))*tau
  rowsum(z*a$membership[match(masks,a$masks),,drop=FALSE],ids,reorder=FALSE)
}

probability_intervals <- function(q,endpoint,delta) {
  stopifnot(all(is.finite(q)),all(q>=0&q<=1),delta>0,delta<1,endpoint%in%c('T','E'))
  n<-length(q); k<-sum(if(endpoint=='T') q<.45 else q>.30)
  out<-list(accepted=n,gate_count=k,mean=if(n) mean(q) else NA_real_,
            mean_lower=0,mean_upper=1,gate_lower=0,gate_upper=1)
  if(n<2L) return(out)
  radius<-sqrt(2*var(q)*log(4/delta)/n)+7*log(4/delta)/(3*(n-1))
  out$mean_lower<-max(0,mean(q)-radius);out$mean_upper<-min(1,mean(q)+radius)
  out$gate_lower<-if(k==0) 0 else qbeta(delta/2,k,n-k+1)
  out$gate_upper<-if(k==n) 1 else qbeta(1-delta/2,k+1,n-k)
  out
}

sample_reference <- function(proposals,chunk,trees,a,successes,counts,prediction_weights,offset,method) {
  s<-as.numeric(successes); n<-as.numeric(counts); f<-n-s; w<-as.numeric(prediction_weights)
  stopifnot(length(s)==a$H,length(n)==a$H,length(w)==a$H,all(n>=0),all(s>=0&s<=n),
            all(s==floor(s)),all(n==floor(n)),all(w>=0),abs(sum(w)-1)<1e-12,is.finite(offset))
  log_bound<-sum(ifelse(s>0,s*log(ifelse(s>0,s/pmax(n,1),1)),0)+
                 ifelse(f>0,f*log(ifelse(f>0,f/pmax(n,1),1)),0))
  accepted<-vector('list',ceiling(proposals/chunk));unit_accepts<-0L
  for(i in seq_along(accepted)) {
    size<-min(chunk,proposals-(i-1L)*chunk)
    field<-prior_fields(size,trees,a,method)+offset
    loglik<-numeric(size)
    if(any(s>0)) loglik<-loglik+as.vector(pnorm(field[,s>0,drop=FALSE],log.p=TRUE)%*%s[s>0])
    if(any(f>0)) loglik<-loglik+as.vector(pnorm(-field[,f>0,drop=FALSE],log.p=TRUE)%*%f[f>0])
    stopifnot(max(loglik-log_bound)<1e-9)
    logu<-log(runif(size));hit<-logu<loglik-log_bound
    unit_accepts<-unit_accepts+sum(logu<loglik)
    accepted[[i]]<-as.vector(pnorm(field[hit,,drop=FALSE])%*%w)
  }
  list(q=unlist(accepted,use.names=FALSE),log_likelihood_bound=log_bound,
       unit_envelope_accepts_same_uniform=unit_accepts)
}
