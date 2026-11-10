

a<-commandArgs(trailingOnly=TRUE);stopifnot(length(a)==1L,!file.exists(a[1]))
script<-sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE));root<-normalizePath(file.path(dirname(script),'..'))
for(n in c('decisions.R','sensitivity.R','audits.R'))source(file.path(root,'R',n))
data<-file.path(root,'data/base/explanation');cells<-read.csv(file.path(data,'source/independent_refit_20261028/CELLS.csv'),colClasses='character')
for(n in c('case_id','refit','position','dose','k_original','k_refit'))cells[[n]]<-as.integer(cells[[n]])
targets<-read.csv(file.path(data,'report/TARGET_ROWS.csv'),colClasses='character')
key<-function(x)paste(x$case_id,x$refit,sep='|');cc<-split(cells,key(cells));tt<-split(targets,key(targets));checks<-0L
gcd<-function(a,b)if(b==0L)a else Recall(b,a%%b)
for(k in names(tt)) {
  c<-cc[[k]];c<-c[order(c$position),];t<-tt[[k]]
  regions<-path_regions(c$k_original,binary_rational(c$q_original),c$k_refit,binary_rational(c$q_refit))
  movement<-sum(abs(c$k_original-c$k_refit));g<-Reduce(gcd,abs(c$k_original-c$k_refit),init=0L)
  for(i in seq_len(nrow(t))) {
    r<-t[i,];j<-if(r$target=='0')0L else match(as.integer(r$target),c$dose)
    for(kind in c('continuous','lattice')) {
      value<-earliest_path(regions[[as.character(j)]],if(kind=='lattice')g else 0L)
      reachable<-as.integer(r[[paste0(kind,'_reachable')]])
      stopifnot(reachable==as.integer(!is.null(value)))
      if(!is.null(value)) {
        cost_name<-if(kind=='lattice')'lattice_cost' else 'continuous_cost_infimum'
        if(value$t!=gmp::as.bigq(r[[paste0(kind,'_first_t')]]))stop('Path mismatch ',k,' target ',r$target,' ',kind,' R=',as.character(value$t),' expected=',r[[paste0(kind,'_first_t')]])
        stopifnot(value$t==gmp::as.bigq(r[[paste0(kind,'_first_t')]]),value$closed==as.logical(as.integer(r[[paste0(kind,'_first_t_attained')]])),
          movement*value$t==gmp::as.bigq(r[[cost_name]]))
      }
    };checks<-checks+1L
  }
}
save_json(list(status='PASS',target_paths=checks,exact_rational=TRUE),a[1]);cat(checks,'path rows passed\n')
