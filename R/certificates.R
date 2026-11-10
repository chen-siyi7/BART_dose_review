
prepare_certificate <- function(record) {
  gates<-lapply(record$gates,function(row) if(is.matrix(row)) row else do.call(rbind,lapply(row,unlist)))
  scores<-do.call(rbind,lapply(record$scores,unlist));J<-length(gates)
  atoms<-unlist(lapply(seq_len(J),function(j) c(paste0('Tgate:',j),
    if(record$stage!='Final') paste0('Egate:',j) else if(J>1L) paste0('Emean:',j))))
  bands<-lapply(atoms,function(a) {
    v<-strsplit(a,':',fixed=TRUE)[[1L]];j<-as.integer(v[2])
    switch(v[1],Tgate=gates[[j]][1L,],Egate=gates[[j]][2L,],Emean=scores[j,])
  });names(bands)<-atoms
  stopifnot(all(vapply(bands,function(b) length(b)==2L&&all(is.finite(b))&&0<=b[1]&&b[1]<=b[2]&&b[2]<=1,TRUE)))
  list(stage=record$stage,gates=gates,scores=scores,target=record$target,atoms=atoms,bands=bands,
       available=atoms[!vapply(bands,function(x) identical(as.numeric(x),c(0,1)),TRUE)])
}
certificate_actions <- function(z,selected=z$atoms) {
  stopifnot(all(selected%in%z$atoms))
  g<-z$gates;q<-z$scores
  for(j in seq_along(g)) {
    if(!paste0('Tgate:',j)%in%selected) g[[j]][1L,]<-c(0,1)
    if(!paste0('Egate:',j)%in%selected) g[[j]][2L,]<-c(0,1)
    if(!paste0('Emean:',j)%in%selected) q[j,]<-c(0,1)
  }
  interval_actions(z$stage,g,q)
}
certifies <- function(z,selected=z$atoms) {
  if(is.null(z$target)) return(FALSE)
  band<-function(a) if(a%in%selected) z$bands[[a]] else c(0,1)
  J<-length(z$gates)
  if(z$stage!='Final') {
    members<-as.integer(strsplit(sub('^mask:','',z$target),':',fixed=TRUE)[[1L]])
    for(j in seq_len(J)) {
      g<-rbind(band(paste0('Tgate:',j)),band(paste0('Egate:',j)))
      if(j%in%members) {if(!all(g[,1]>=.3))return(FALSE)} else if(!any(g[,2]<.3)) return(FALSE)
    };return(TRUE)
  }
  r<-as.integer(sub('^dose:','',z$target))
  if(r==0L) return(all(vapply(seq_len(J),function(j) band(paste0('Tgate:',j))[2]<.3,TRUE)))
  if(band(paste0('Tgate:',r))[1]<.3) return(FALSE)
  lower<-band(paste0('Emean:',r))[1]
  for(j in setdiff(seq_len(J),r)) {
    if(band(paste0('Tgate:',j))[2]<.3) next
    upper<-band(paste0('Emean:',j))[2]
    if(!(lower>upper||(lower==upper&&r<j))) return(FALSE)
  };TRUE
}
minimum_certificate <- function(z,costs=setNames(rep(1,length(z$atoms)),z$atoms)) {
  stopifnot(all(z$atoms%in%names(costs)),all(is.finite(costs)),all(costs>0))
  fail<-list(status='UNCERTIFIABLE',cost=NA_real_,atoms=character())
  if(!certifies(z)) return(fail)
  least<-function(options) if(!length(options)) NULL else options[[which.min(vapply(options,function(a) sum(costs[a]),0))]]
  J<-length(z$gates);passes<-fails<-vector('list',J)
  for(j in seq_len(J)) {
    aa<-c(paste0('Tgate:',j),if(z$stage!='Final') paste0('Egate:',j))
    gg<-do.call(rbind,z$bands[aa])
    if(all(gg[,1]>=.3)) passes[[j]]<-aa
    options<-lapply(aa[gg[,2]<.3],function(a) a)
    if(length(options)) fails[[j]]<-least(options)
  }
  if(z$stage!='Final') {
    members<-as.integer(strsplit(sub('^mask:','',z$target),':',fixed=TRUE)[[1L]])
    selected<-unlist(lapply(seq_len(J),function(j) if(j%in%members) passes[[j]] else fails[[j]]),use.names=FALSE)
  } else {
    r<-as.integer(sub('^dose:','',z$target))
    if(r==0L) selected<-unlist(fails,use.names=FALSE) else {
      options<-list()
      for(retain in c(FALSE,TRUE)) {
        target_atom<-paste0('Emean:',r)
        if(retain&&!target_atom%in%z$atoms) next
        lower<-if(retain) z$scores[r,1L] else 0
        sel<-c(passes[[r]],if(retain) target_atom);possible<-TRUE
        for(j in setdiff(seq_len(J),r)) {
          alt<-list();if(!is.null(fails[[j]])) alt[[length(alt)+1L]]<-fails[[j]]
          if(lower==1&&r<j) alt[[length(alt)+1L]]<-character()
          upper<-z$scores[j,2L]
          if(lower>upper||(lower==upper&&r<j)) alt[[length(alt)+1L]]<-paste0('Emean:',j)
          if(!length(alt)) {possible<-FALSE;break}
          sel<-c(sel,least(alt))
        }
        if(possible) options[[length(options)+1L]]<-sel
      }
      if(!length(options)) return(fail)
      selected<-least(options)
    }
  }
  selected<-z$atoms[z$atoms%in%selected];stopifnot(certifies(z,selected))
  list(status='CERTIFIED',cost=sum(costs[selected]),atoms=selected)
}
deletion_certificate <- function(z,order='dose-forward') {
  if(!certifies(z)) return(list(cost=NA_real_,atoms=character()))
  aa<-z$available;sequence<-aa
  if(order=='dose-reverse') sequence<-aa[order(-as.integer(sub('.*:','',aa)),seq_along(aa)%%2L)]
  if(order=='target-mean-first'&&z$stage=='Final'&&z$target!='dose:0') {
    a<-sub('dose:','Emean:',z$target);sequence<-c(intersect(a,aa),setdiff(aa,a))
  }
  for(a in sequence) {trial<-setdiff(aa,a);if(certifies(z,trial)) aa<-trial}
  list(cost=length(aa),atoms=aa)
}
exhaustive_certificate <- function(z,costs=setNames(rep(1,length(z$atoms)),z$atoms)) {
  if(is.null(z$target)||!identical(certificate_actions(z),z$target)) return(list(cost=NA_real_,atoms=character(),subsets=0L))
  aa<-z$available;bits<-all_masks(length(aa));best<-Inf;sel<-NULL

  for(i in seq_len(nrow(bits))) {
    candidate<-aa[bits[i,]];value<-sum(costs[candidate]);if(value>best) next
    if(identical(certificate_actions(z,candidate),z$target)&&value<best) {best<-value;sel<-candidate}
  }
  list(cost=best,atoms=sel,subsets=nrow(bits))
}
