

prepare_tree_prior <- function(X,cuts,counts,minimum=5L) {
  X<-as.matrix(X);counts<-as.numeric(counts);H<-nrow(X);p<-ncol(X)
  cuts<-lapply(cuts,as.numeric)
  stopifnot(H>=1L,H<=30L,p>=1L,length(counts)==H,all(is.finite(X)),
    all(counts>=0),all(counts==floor(counts)),sum(counts)>=minimum,
    length(cuts)==p,all(lengths(cuts)>0L),
    all(vapply(cuts,function(x)all(is.finite(x))&&!is.unsorted(x,strictly=TRUE),TRUE)))
  bits<-bitwShiftL(1L,seq_len(H)-1L)
  left<-lapply(seq_len(p),function(v)vapply(cuts[[v]],function(c)as.integer(sum(bits[X[,v]<c])),0L))
  tab<-NULL
  if(H<=16L) {
    masks<-0:(2^H-1L);tab<-numeric(length(masks))
    for(h in seq_len(H))tab<-tab+counts[h]*as.numeric(bitwAnd(masks,bits[h])!=0L)
  }
  list(H=H,p=p,bits=bits,counts=counts,minimum=minimum,ncuts=lengths(cuts),
       left=left,count_table=tab,root_mask=as.integer(2^H-1L))
}
tree_training_counts <- function(mask,a) {
  if(!is.null(a$count_table))return(a$count_table[mask+1L])
  total<-numeric(length(mask))
  for(h in seq_len(a$H))total<-total+a$counts[h]*as.numeric(bitwAnd(mask,a$bits[h])!=0L)
  total
}
draw_tree_attempts <- function(n,a) {
  n<-as.integer(n);valid<-rep(TRUE,n);root<-seq_len(n);mask<-rep(a$root_mask,n)
  lo<-matrix(1L,n,a$p);hi<-matrix(rep(a$ncuts,each=n),n,a$p)
  leaf_root<-leaf_mask<-list();depth<-0L;k<-0L;nodes<-0L
  while(length(root)) {
    keep<-valid[root];root<-root[keep];mask<-mask[keep]
    lo<-lo[keep,,drop=FALSE];hi<-hi[keep,,drop=FALSE]
    if(!length(root))break
    nodes<-nodes+length(root)
    grows<-rowSums(lo<=hi)>0L & runif(length(root))<.95/(depth+1)^2
    stops<-which(!grows)
    if(length(stops)){k<-k+1L;leaf_root[[k]]<-root[stops];leaf_mask[[k]]<-mask[stops]}
    split<-which(grows)
    if(!length(split))break
    variable<-sample.int(a$p,length(split),replace=TRUE)
    lower<-lo[cbind(split,variable)];upper<-hi[cbind(split,variable)]
    bad<-upper<lower
    if(any(bad))valid[root[split[bad]]]<-FALSE
    keep<-!bad & valid[root[split]]
    split<-split[keep];variable<-variable[keep];lower<-lower[keep];upper<-upper[keep]
    if(!length(split))break
    cut<-lower+floor(runif(length(split))*(upper-lower+1L))
    full_left<-integer(length(split))
    for(v in seq_len(a$p)){ix<-which(variable==v);full_left[ix]<-a$left[[v]][cut[ix]]}
    left<-bitwAnd(mask[split],full_left);right<-bitwXor(mask[split],left)
    good<-tree_training_counts(left,a)>=a$minimum & tree_training_counts(right,a)>=a$minimum
    if(any(!good))valid[root[split[!good]]]<-FALSE
    keep<-good & valid[root[split]]
    split<-split[keep];variable<-variable[keep];cut<-cut[keep]
    left<-left[keep];right<-right[keep];nn<-length(split)
    if(!nn)break
    next_lo<-rbind(lo[split,,drop=FALSE],lo[split,,drop=FALSE])
    next_hi<-rbind(hi[split,,drop=FALSE],hi[split,,drop=FALSE])
    next_hi[cbind(seq_len(nn),variable)]<-cut-1L
    next_lo[cbind(nn+seq_len(nn),variable)]<-cut+1L
    root<-rep(root[split],2L);mask<-c(left,right);lo<-next_lo;hi<-next_hi;depth<-depth+1L
  }
  roots<-unlist(leaf_root,use.names=FALSE);masks<-unlist(leaf_mask,use.names=FALSE)
  keep<-valid[roots]
  list(valid=which(valid),leaf_root=roots[keep],leaf_mask=masks[keep],attempts=n,nodes=nodes)
}
prior_tree_fields <- function(n,trees,a,stats=NULL) {
  need<-as.integer(n*trees);used<-0L;field<-matrix(0,n,a$H)
  attempts_so_far<-valid_so_far<-0
  while(used<need) {

    rate<-if(attempts_so_far)max(.01,valid_so_far/attempts_so_far)else .8
    attempts<-min(65536L,max(64L,as.integer(ceiling(1.05*(need-used)/rate))))
    z<-draw_tree_attempts(attempts,a)
    attempts_so_far<-attempts_so_far+attempts;valid_so_far<-valid_so_far+length(z$valid)
    if(!is.null(stats)) {
      stats$attempts<-stats$attempts+z$attempts;stats$valid<-stats$valid+length(z$valid)
      stats$nodes<-stats$nodes+z$nodes
    }
    take<-min(length(z$valid),need-used)
    if(!take)next
    map<-integer(attempts);map[z$valid[seq_len(take)]]<-used+seq_len(take)
    id<-map[z$leaf_root];keep<-id>0L;id<-id[keep];masks<-z$leaf_mask[keep]
    membership<-outer(masks,a$bits,function(x,b)as.numeric(bitwAnd(x,b)!=0L))
    leaf_values<-rnorm(length(id))*3/(2*sqrt(trees))
    sums<-rowsum(membership*leaf_values,(id-1L)%/%trees+1L,reorder=TRUE)
    rows<-as.integer(rownames(sums));field[rows,]<-field[rows,,drop=FALSE]+sums
    used<-used+take
  }
  field
}



sample_reference_checkpoints <- function(checkpoints,chunk,trees,a,successes,counts,
                                        prediction_weights,offset,method,delta,endpoint) {
  stopifnot(all(diff(checkpoints)>0),all(checkpoints==floor(checkpoints)),chunk>=1L,
            method%in%c('tree_prior','direct','aggregated'))
  s<-as.numeric(successes);n<-as.numeric(counts);f<-n-s;w<-as.numeric(prediction_weights)
  stopifnot(length(s)==a$H,length(n)==a$H,length(w)==a$H,all(n>=0),all(s>=0&s<=n),
            all(s==floor(s)),all(n==floor(n)),
            all(w>=0),abs(sum(w)-1)<1e-12,is.finite(offset))
  bound<-sum(ifelse(s>0,s*log(ifelse(s>0,s/pmax(n,1),1)),0)+
             ifelse(f>0,f*log(ifelse(f>0,f/pmax(n,1),1)),0))
  stats<-new.env(parent=emptyenv());stats$attempts<-stats$valid<-stats$nodes<-0
  accepted<-list();hit_at<-list();k<-0L;done<-0L;reports<-list();start<-unname(proc.time()[3])
  for(target in checkpoints) {
    while(done<target) {
      size<-min(chunk,target-done)
      field<-(if(method=='tree_prior')prior_tree_fields(size,trees,a,stats)else prior_fields(size,trees,a,method))+offset
      ll<-numeric(size)
      if(any(s>0))ll<-ll+as.vector(pnorm(field[,s>0,drop=FALSE],log.p=TRUE)%*%s[s>0])
      if(any(f>0))ll<-ll+as.vector(pnorm(-field[,f>0,drop=FALSE],log.p=TRUE)%*%f[f>0])
      stopifnot(all(is.finite(ll)),max(ll-bound)<1e-9)
      hit<-which(log(runif(size))<ll-bound);k<-k+1L
      accepted[[k]]<-as.vector(pnorm(field[hit,,drop=FALSE])%*%w)
      hit_at[[k]]<-done+hit;done<-done+size
    }
    q<-unlist(accepted,use.names=FALSE)
    band<-probability_intervals(q,endpoint,delta)
    reports[[as.character(target)]]<-c(list(proposals=target,seconds=unname(proc.time()[3])-start),band)
  }
  list(q=unlist(accepted,use.names=FALSE),accepted_at=unlist(hit_at,use.names=FALSE),
       reports=reports,tree_attempts=stats$attempts,valid_tree_attempts=stats$valid,
       visited_nodes=stats$nodes,log_likelihood_bound=bound)
}
