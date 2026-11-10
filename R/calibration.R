
comparison_sets <- function(a,scaling=1) {
  if(a$status!='valid') return(NULL)
  z<-return_values(a);full<-action_alphabet(a$stage,z$doses)
  grid<-c(0,.0005,.001,.002,.005,.01,.02,.05,.1,.2,.5,1,length(z$doses)*ncol(z$gates))
  lapply(grid,function(x) {
    budget<-round(1000*x);eta<-if(a$stage=='Final') min(1,scaling*x) else 0
    j<-frontier(a$stage,z$doses,z$gates,z$q,budget,eta)
    b<-coordinate_actions(a$stage,z$doses,z$gates,z$q,budget,eta)
    costs<-membership_costs(z$gates)
    warn<-min(costs$rho)<=budget||(a$stage=='Final'&&length(frontier(a$stage,z$doses,z$gates,z$q,0,eta))>1L)
    w<-if(warn) full else a$action
    stopifnot(all(j%in%b),all(j%in%w),a$action%in%j)
    list(joint=sort(j),coordinate=sort(b),warning=sort(w))
  })
}
replay_comparison <- function(path,out) {
  jobs<-read.csv(file.path(path,'INPUTS.csv'),stringsAsFactors=FALSE)
  pair<-lapply(jobs$case_id,function(i) lapply(1:2,function(role) read_json(file.path(path,'raw',sprintf('case%04d_role%d.json',i,role)))))
  names(pair)<-jobs$case_id
  scales<-c(.25,1,4);scale_names<-c('1/4','1','4');methods<-expand.grid(method=c('joint','coordinate','warning'),scaling=seq_along(scales),stringsAsFactors=FALSE)
  methods<-rbind(methods,data.frame(method='input_bridge',scaling=2L))
  cache<-lapply(pair,function(p) lapply(scales,function(c) comparison_sets(p[[1]],c)))
  cal<-list();kk<-0L
  for(i in which(jobs$bank=='calibration')) {
    job<-jobs[i,];p<-pair[[as.character(job$case_id)]]
    for(k in seq_len(nrow(methods))) {
      m<-methods[k,];index<-13L
      if(p[[1]]$status=='valid'&&p[[2]]$status=='valid') {
        if(m$method=='input_bridge') {
          grid<-c(0,.0005,.001,.002,.005,.01,.02,.05,.1,.2,.5,1,if(job$stage=='Final')5 else 10)
          index<-which(pair_score(p[[1]],p[[2]])<=grid)[1]-1L
        } else index<-which(vapply(cache[[i]][[m$scaling]],function(s) p[[2]]$action%in%s[[m$method]],TRUE))[1]-1L
      }
      kk<-kk+1L;cal[[kk]]<-data.frame(case_id=job$case_id,layout=job$layout,generator=job$generator,profile=job$scenario,
        method=m$method,scaling=scale_names[m$scaling],index=index)
    }
  }
  cc<-do.call(rbind,cal);key<-function(x) paste(x$layout,x$generator,x$profile,x$method,x$scaling,sep='|')
  thresholds<-tapply(cc$index,key(cc),max);stopifnot(length(thresholds)==300L,all(table(key(cc))==20L))
  results<-list();kk<-0L
  for(i in which(jobs$bank=='assessment')) {
    job<-jobs[i,];p<-pair[[as.character(job$case_id)]];a<-p[[1]];b<-p[[2]]
    full<-action_alphabet(job$stage,1:5);returned<-if(b$status=='valid') b$action else 'failure'
    for(k in seq_len(nrow(methods))) {
      m<-methods[k,];keyrow<-data.frame(layout=job$layout,generator=job$generator,profile=job$scenario,method=m$method,scaling=scale_names[m$scaling])
      index<-unname(thresholds[key(keyrow)])
      labels<-if(a$status!='valid'||index==13L) c(full,'failure') else
        cache[[i]][[m$scaling]][[index+1L]][[if(m$method=='input_bridge') 'joint' else m$method]]
      stable<-a$status=='valid'&&identical(sort(labels),a$action)
      kk<-kk+1L;results[[kk]]<-data.frame(case_id=job$case_id,keyrow,stage=job$stage,calibration_index=index,
        actions=paste(sort(labels),collapse=';'),set_size=length(labels),valid_alphabet=length(full),
        miscoverage=as.integer(!returned%in%labels),normalized_size=if(index==13L||a$status!='valid')1 else (sum(labels%in%full)-1)/(length(full)-1),
        unflagged=as.integer(stable),uninformative=as.integer(all(full%in%labels)),
        changed_failed_unflagged=as.integer(stable&&(b$status!='valid'||a$action!=returned)),anchor_status=a$status,comparison_status=b$status)
    }
  }
  rr<-do.call(rbind,results);expected<-read.csv(file.path(path,'ASSESSMENT_CASES.csv'),stringsAsFactors=FALSE)
  id<-function(x) paste(x$case_id,x$method,x$scaling,sep='|');expected<-expected[match(id(rr),id(expected)),]
  for(n in c('actions','set_size','miscoverage','calibration_index','unflagged','uninformative','changed_failed_unflagged')) stopifnot(identical(as.character(rr[[n]]),as.character(expected[[n]])))
  num<-function(x) vapply(as.character(x),function(s) as.numeric(gmp::as.bigq(s)),0)
  stopifnot(max(abs(rr$normalized_size-num(expected$normalized_size)))<1e-12)
  cal_saved<-read.csv(file.path(path,'CALIBRATION_SCORES.csv'),stringsAsFactors=FALSE)
  stopifnot(identical(cc$index,as.integer(cal_saved$index[match(id(cc),id(cal_saved))])))
  observations<-list();ranges<-list()
  for(k in seq_len(nrow(methods))) {
    m<-methods[k,];r<-rr[rr$method==m$method&rr$scaling==scale_names[m$scaling],]
    for(metric in c('miscoverage','normalized_size','unflagged','uninformative','changed_failed_unflagged')) {
      name<-paste(m$method,scale_names[m$scaling],metric,sep='|');observations[[name]]<-r[[metric]];ranges[[name]]<-c(0,1)
    }
  }
  for(s in scale_names) for(other in c('coordinate','warning')) {
    a<-rr[rr$method==other&rr$scaling==s,];b<-rr[rr$method=='joint'&rr$scaling==s,]
    a<-a[match(b$case_id,a$case_id),];n<-paste0('gain_',other,'_minus_joint|',s)
    observations[[n]]<-a$normalized_size-b$normalized_size;ranges[[n]]<-c(-1,1)
  }
  a<-rr[rr$method=='input_bridge',];b<-rr[rr$method=='joint'&rr$scaling=='1',];a<-a[match(b$case_id,a$case_id),]
  observations[['gain_input_bridge_minus_joint|1']]<-a$normalized_size-b$normalized_size;ranges[['gain_input_bridge_minus_joint|1']]<-c(-1,1)
  stopifnot(length(observations)==57L)
  intervals<-do.call(rbind,lapply(names(observations),function(n) {
    v<-observations[[n]];ran<-ranges[[n]];mu<-mean(v);ci<-kl_interval((mu-ran[1])/diff(ran),length(v),.05/114)
    data.frame(metric=n,n=length(v),mean=mu,lower=ran[1]+diff(ran)*ci[1],upper=ran[1]+diff(ran)*ci[2])
  }))
  expected_ci<-read.csv(file.path(path,'SIMULTANEOUS_INTERVALS.csv'));expected_ci<-expected_ci[match(intervals$metric,expected_ci$metric),]
  stopifnot(max(abs(as.matrix(intervals[c('mean','lower','upper')])-as.matrix(expected_ci[c('mean','lower','upper')])))<1e-10)
  write.csv(rr,file.path(out,'COMPARISON_ASSESSMENT.csv'),row.names=FALSE)
  write.csv(intervals,file.path(out,'COMPARISON_INTERVALS.csv'),row.names=FALSE)
  list(status='PASS',histories=nrow(jobs),assessed=length(unique(rr$case_id)),assessment_rows=nrow(rr),calibration_scores=nrow(cc),intervals=nrow(intervals),
       primary=intervals[intervals$metric=='gain_coordinate_minus_joint|1',])
}

replay_ordinary <- function(path,out) {
  expected<-read.csv(file.path(path,'ASSESSMENT_CASES.csv'),stringsAsFactors=FALSE)
  thresholds<-read.csv(file.path(path,'CALIBRATION_THRESHOLDS.csv'),stringsAsFactors=FALSE)
  numeric_fraction<-function(s) if(s=='infinity') Inf else as.numeric(gmp::as.bigq(s))
  threshold<-vapply(thresholds$threshold,numeric_fraction,0);names(threshold)<-paste(thresholds$layout,thresholds$scenario,sep='|')
  rows<-list()
  for(i in seq_len(nrow(expected))) {
    e<-expected[i,];a<-read_json(file.path(path,'raw',sprintf('case%04d_role1.json',e$case_id)));b<-read_json(file.path(path,'raw',sprintf('case%04d_role2.json',e$case_id)))
    t<-threshold[paste(e$layout,e$scenario,sep='|')];full<-action_alphabet(a$stage,unlist(a$candidates))
    if(!is.finite(t)||a$status!='valid') {joint<-box<-warning<-c(full,'failure')} else {
      z<-return_values(a);budget<-floor(1000*t+1e-10);eta<-if(a$stage=='Final'&&length(z$doses)>1L)min(1,t) else 0
      joint<-frontier(a$stage,z$doses,z$gates,z$q,budget,eta);box<-coordinate_actions(a$stage,z$doses,z$gates,z$q,budget,eta)
      costs<-membership_costs(z$gates);warn<-min(costs$rho)<=budget||(a$stage=='Final'&&length(frontier(a$stage,z$doses,z$gates,z$q,0,eta))>1)
      warning<-if(warn) full else a$action
    }
    returned<-if(b$status=='valid') b$action else 'failure'
    stopifnot(identical(paste(sort(joint),collapse=';'),e$joint_actions),identical(paste(sort(box),collapse=';'),e$box_actions),identical(paste(sort(warning),collapse=';'),e$coarse_actions))
    rows[[i]]<-data.frame(case_id=e$case_id,miscoverage_joint=as.integer(!returned%in%joint),miscoverage_box=as.integer(!returned%in%box),
      reduction_flag=(length(warning)-length(joint))/(length(full)-1),reduction_box=(length(box)-length(joint))/(length(full)-1))
  }
  rr<-do.call(rbind,rows)
  intervals<-do.call(rbind,lapply(setdiff(names(rr),'case_id'),function(n) {
    mu<-mean(rr[[n]]);ci<-kl_interval(mu,nrow(rr),.05/16)
    data.frame(metric=n,mean=mu,lower=ci[1],upper=ci[2])
  }))
  write.csv(rr,file.path(out,'ORDINARY_ASSESSMENT.csv'),row.names=FALSE);write.csv(intervals,file.path(out,'ORDINARY_INTERVALS.csv'),row.names=FALSE)
  list(status='PASS',assessment_histories=nrow(rr),intervals=intervals)
}
