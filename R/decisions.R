
all_masks <- function(J) if(J==0L) matrix(logical(),1L,0L) else as.matrix(expand.grid(rep(list(c(FALSE,TRUE)),J)))
mask_label <- function(bits,doses) paste0('mask:',paste(doses[bits],collapse=':'))
action_alphabet <- function(stage,doses) {
  if(stage=='Final') c('dose:0',paste0('dose:',doses)) else
    apply(all_masks(length(doses)),1L,mask_label,doses=doses)
}
interval_actions <- function(stage,gates,scores,doses=seq_len(length(gates)),cutoff=.30) {
  gg<-lapply(gates,function(g) if(stage=='Final') g[1L,,drop=FALSE] else g)
  mandatory<-vapply(gg,function(g) all(g[,1L]>=cutoff),TRUE)
  possible<-vapply(gg,function(g) all(g[,2L]>=cutoff),TRUE)
  if(stage!='Final') {
    bits<-all_masks(length(doses));ok<-apply(bits,1L,function(b) all((!mandatory|b)&(!b|possible)))
    return(sort(apply(bits[ok,,drop=FALSE],1L,mask_label,doses=doses)))
  }
  out<-if(any(mandatory)) character() else 'dose:0'
  for(r in which(possible)) {
    j<-setdiff(which(mandatory),r)
    if(all(scores[r,2L]>scores[j,1L]|(scores[r,2L]==scores[j,1L]&r<j))) out<-c(out,paste0('dose:',doses[r]))
  }
  sort(out)
}
work_actions <- function(rows,stage) {
  by<-function(j,e) rows[[which(vapply(rows,function(r) r$dose==j&&r$endpoint==e,TRUE))]]
  gates<-lapply(1:5,function(j) do.call(rbind,lapply(c('T','E'),function(e) {r<-by(j,e);c(r$gate_lower,r$gate_upper)})))
  scores<-t(vapply(1:5,function(j) {r<-by(j,'E');c(r$mean_lower,r$mean_upper)},numeric(2)))
  interval_actions(stage,gates,scores)
}
return_values <- function(x) {
  doses<-as.integer(unlist(x$candidates));gates<-cbind(as.integer(unlist(x$kT)),as.integer(unlist(x$kE)))
  if(x$stage=='Final') gates<-gates[,1L,drop=FALSE]
  list(doses=doses,gates=gates,q=as.numeric(unlist(x$q)))
}
membership_costs <- function(gates,kappa=300L) {
  keep<-apply(gates,1L,function(g) all(g>=kappa))
  rho<-vapply(seq_len(nrow(gates)),function(j) if(keep[j]) min(gates[j,]-kappa+1L) else sum(pmax(0,kappa-gates[j,])),0)
  list(keep=keep,rho=rho)
}
frontier_costs <- function(stage,doses,gates,q,eta=0) {
  z<-membership_costs(gates);keep<-z$keep;rho<-z$rho
  if(stage!='Final') {
    bits<-all_masks(length(doses))
    cost<-apply(bits,1L,function(b) sum(rho[b!=keep]))
    return(setNames(cost,apply(bits,1L,mask_label,doses=doses)))
  }
  lo<-pmax(0,q-eta);hi<-pmin(1,q+eta)
  cost<-vapply(seq_along(doses),function(r) {
    blockers<-which(keep&(lo>hi[r]|(lo==hi[r]&seq_along(doses)<r)))
    blockers<-setdiff(blockers,r)
    (if(keep[r]) 0 else rho[r])+sum(rho[blockers])
  },0)
  setNames(c(sum(rho[keep]),cost),c('dose:0',paste0('dose:',doses)))
}
frontier <- function(stage,doses,gates,q,budget,eta) {
  v<-frontier_costs(stage,doses,gates,q,eta);sort(names(v)[v<=budget])
}
coordinate_actions <- function(stage,doses,gates,q,budget,eta) {
  gg<-lapply(seq_along(doses),function(j) cbind(as.numeric(gates[j,]-budget>=300),as.numeric(gates[j,]+budget>=300)))
  interval_actions(stage,gg,cbind(pmax(0,q-eta),pmin(1,q+eta)),doses,cutoff=.5)
}
pair_score <- function(a,b) {
  if(a$status!='valid'||b$status!='valid') return(Inf)
  aa<-return_values(a);bb<-return_values(b);stopifnot(identical(aa$doses,bb$doses))
  max(sum(abs(aa$gates-bb$gates))/1000,if(a$stage=='Final'&&length(aa$doses)>1L) max(abs(aa$q-bb$q)) else 0)
}
kl_interval <- function(mean,n,tail) {
  stopifnot(mean>=0,mean<=1,n>0,tail>0,tail<1)
  target<-log(1/tail)/n
  divergence<-function(p) {
    if(p<=0) return(if(mean==0) 0 else Inf)
    if(p>=1) return(if(mean==1) 0 else Inf)
    (if(mean>0) mean*log(mean/p) else 0)+(if(mean<1) (1-mean)*log((1-mean)/(1-p)) else 0)
  }
  lo<-0;hi<-1
  if(mean>0) {a<-0;b<-mean;for(i in 1:90) {p<-(a+b)/2;if(divergence(p)>target)a<-p else b<-p};lo<-max(0,a-1e-12)}
  if(mean<1) {a<-mean;b<-1;for(i in 1:90) {p<-(a+b)/2;if(divergence(p)>target)b<-p else a<-p};hi<-min(1,b+1e-12)}
  c(lo,hi)
}
diagnostic_intervals <- function(rows,mode='decision') {
  band<-function(kind,j) {
    r<-rows[rows$series==kind&rows$dose==j,,drop=FALSE];stopifnot(nrow(r)==1L)
    ok<-if(mode=='uniform') r$precision_filtered else r$qualified
    if(isTRUE(ok)) c(r$lower,r$upper) else c(0,1)
  }
  gates<-lapply(1:5,function(j) do.call(rbind,lapply(c('T','E'),function(e) {
    gate<-band(paste0(e,'_gate'),j);mu<-band(paste0(e,'_mean'),j)
    bounds<-if(e=='T') c(max(0,1-mu[2]/.45),min(1,(1-mu[1])/.55)) else c(max(0,(mu[1]-.3)/.7),min(1,mu[2]/.3))
    out<-c(max(gate[1],bounds[1]),min(gate[2],bounds[2]));stopifnot(out[1]<=out[2]+1e-12);out
  })))
  scores<-t(vapply(1:5,function(j) band('E_mean',j),numeric(2)))
  list(gates=gates,scores=scores)
}
