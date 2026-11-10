
target_cost <- function(doses,gates,q,target,eta=0) {
  z<-membership_costs(gates)
  if(target==0L) return(sum(z$rho[z$keep]))
  r<-match(target,doses);stopifnot(!is.na(r),eta>=0,eta<=1)
  lo<-pmax(0,q-eta);hi<-min(1,q[r]+eta)
  blockers<-setdiff(which(z$keep&(lo>hi|(lo==hi&seq_along(doses)<r))),r)
  (if(z$keep[r])0 else z$rho[r])+sum(z$rho[blockers])
}
frontier_curve <- function(doses,gates,q) {
  breaks<-sort(unique(c(0,1,as.vector(outer(q,q,`-`))/2)))
  breaks<-breaks[breaks>=0&breaks<=1]
  rows<-list()
  for(i in seq_along(breaks)) {
    eta<-breaks[i];costs<-frontier_costs('Final',doses,gates,q,eta)
    rows[[length(rows)+1L]]<-data.frame(eta=eta,side='point',target=names(costs),cost=as.numeric(costs))
    if(i<length(breaks)) {
      mid<-(eta+breaks[i+1L])/2;costs<-frontier_costs('Final',doses,gates,q,mid)
      rows[[length(rows)+1L]]<-data.frame(eta=eta,side='open-right',target=names(costs),cost=as.numeric(costs))
    }
  };do.call(rbind,rows)
}
frontier_witness <- function(doses,gates,q,target,eta=0) {
  z<-membership_costs(gates);K<-gates;score<-pmax(0,q-eta)
  r<-if(target==0L)0L else match(target,doses);stopifnot(!is.na(r))
  if(r==0L) remove<-which(z$keep) else {
    score[r]<-min(1,q[r]+eta)
    if(!z$keep[r]) K[r,]<-pmax(K[r,],300)
    remove<-setdiff(which(z$keep&(score>score[r]|(score==score[r]&seq_along(doses)<r))),r)
  }
  for(j in remove) K[j,which.min(K[j,])]<-299L
  cost<-sum(abs(K-gates));stopifnot(cost==target_cost(doses,gates,q,target,eta))
  keep<-apply(K,1L,function(g)all(g>=300));winner<-if(any(keep))doses[which(keep)[which.max(score[keep])]] else 0L
  stopifnot(winner==target)
  list(target=target,cost=cost,counts=K,scores=score)
}


binary_rational <- function(text) {


  one<-function(s) {
    s<-tolower(trimws(s));neg<-startsWith(s,'-');s<-sub('^[+-]','',s)
    parts<-strsplit(s,'e',fixed=TRUE)[[1L]];exponent<-if(length(parts)>1L)as.integer(parts[2L])else 0L
    decimal<-strsplit(parts[1L],'.',fixed=TRUE)[[1L]];places<-if(length(decimal)>1L)nchar(decimal[2L])else 0L
    digits<-sub('^0+','',gsub('.','',parts[1L],fixed=TRUE));if(!nzchar(digits))digits<-'0'
    num<-gmp::as.bigz(digits);den<-gmp::as.bigz(10)^max(0L,places-exponent)
    if(exponent>places)num<-num*gmp::as.bigz(10)^(exponent-places)
    value<-gmp::as.bigq(num,den);if(value==0)return(value)
    e<-floor(log2(as.numeric(value)))
    power<-function(k)if(k>=0)gmp::as.bigq(gmp::as.bigz(2)^k)else gmp::as.bigq(1,gmp::as.bigz(2)^(-k))
    while(value<power(e))e<-e-1L
    while(value>=power(e+1L))e<-e+1L
    unit<-power(52L-e);scaled<-value*unit;nn<-gmp::numerator(scaled);dd<-gmp::denominator(scaled)
    m<-nn%/%dd;remainder<-nn%%dd
    if(2L*remainder>dd||(2L*remainder==dd&&m%%2L==1L))m<-m+1L
    answer<-gmp::as.bigq(m)/unit;if(neg)-answer else answer
  }
  do.call(base::c,lapply(as.character(text),one))
}
affine_restrict <- function(bounds,constant,slope,weak) {
  if(is.null(bounds))return(NULL)
  if(slope==0) return(if(constant>0||(constant==0&&weak))bounds else NULL)
  z<--constant/slope
  if(slope>0) {
    if(z>bounds$lo) {bounds$lo<-z;bounds$lc<-weak} else if(z==bounds$lo) bounds$lc<-bounds$lc&&weak
  } else {
    if(z<bounds$hi) {bounds$hi<-z;bounds$hc<-weak} else if(z==bounds$hi) bounds$hc<-bounds$hc&&weak
  }
  if(bounds$lo>bounds$hi||(bounds$lo==bounds$hi&&!(bounds$lc&&bounds$hc)))NULL else bounds
}
path_regions <- function(k0,q0,k1,q1) {
  J<-length(k0);q0<-gmp::as.bigq(q0);q1<-gmp::as.bigq(q1);regions<-vector('list',J+1L);names(regions)<-as.character(0:J)
  bits<-all_masks(J)
  for(b in seq_len(nrow(bits))) {
    bounds<-list(lo=gmp::as.bigq(0),lc=TRUE,hi=gmp::as.bigq(1),hc=TRUE)
    for(j in seq_len(J)) {
      sign<-if(bits[b,j])1 else -1
      bounds<-affine_restrict(bounds,gmp::as.bigq(sign*(k0[j]-300)),gmp::as.bigq(sign*(k1[j]-k0[j])),bits[b,j])
      if(is.null(bounds))break
    }
    if(is.null(bounds))next
    members<-which(bits[b,]);targets<-if(length(members))members else 0L
    for(r in targets) {
      rr<-bounds
      if(r>0L)for(j in setdiff(members,r)) {
        rr<-affine_restrict(rr,q0[r]-q0[j],q1[r]-q0[r]-q1[j]+q0[j],r<j)
        if(is.null(rr))break
      }
      if(!is.null(rr))regions[[as.character(r)]][[length(regions[[as.character(r)]])+1L]]<-rr
    }
  };regions
}
earliest_path <- function(regions,g=0L) {
  if(!length(regions))return(NULL)
  if(g==0L) {
    lo<-Reduce(min,lapply(regions,function(r)r$lo))
    return(list(t=lo,closed=any(vapply(regions,function(r)r$lo==lo&&r$lc,TRUE))))
  }
  times<-list()
  for(r in regions) {
    z<-r$lo*g;m<- -((-gmp::numerator(z))%/%gmp::denominator(z));t<-gmp::as.bigq(m,g)
    if(t==r$lo&&!r$lc) t<-t+gmp::as.bigq(1,g)
    if(t<r$hi||(t==r$hi&&r$hc))times[[length(times)+1L]]<-t
  }
  if(!length(times))NULL else list(t=Reduce(min,times),closed=TRUE)
}
