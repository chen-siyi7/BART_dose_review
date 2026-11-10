
ast_save_json <- function(x,path)jsonlite::write_json(x,path,auto_unbox=TRUE,pretty=TRUE,digits=17,na='null',null='null')
ast_sha <- function(path)digest::digest(file=path,algo='sha256')
ast_id <- function(x)gsub('/','',trimws(as.character(x)),fixed=TRUE)
ast_prepare_inputs <- function(raw,out) {
  paths <- file.path(raw,paste0('pone.0232394.s',c('002','003','005'),'.xlsx'))
  stopifnot(all(file.exists(paths)))
  pop <- readxl::read_excel(paths[1],sheet='Populations')
  eff <- readxl::read_excel(paths[2],sheet='Efficacy')
  ae <- readxl::read_excel(paths[3],sheet='AEs')
  pid <- ast_id(pop[['Site/ Subject']]);eid <- ast_id(eff[['Site/Subject']]);aid <- trimws(as.character(ae$Patient))
  stopifnot(!anyDuplicated(pid),!anyDuplicated(eid),all(eid%in%pid),all(aid%in%pid),
    setequal(pid[pop$Safety=='Yes'],eid),all(pop$mITT[match(eid,pid)]==eff$mITT))
  arms <- c('HEC placebo','0.5% SPL7013','1% SPL7013','3% SPL7013')
  conc <- c(0,.5,1,3)
  treated <- pop$Safety=='Yes';ei <- match(pid[treated],eid)
  related <- unique(aid[ae$Relationship%in%c('POSSIBLY RELATED','PROBABLY RELATED','DEFINITELY RELATED')])
  binary <- function(x){stopifnot(all(is.na(x)|x%in%c('Yes','No')));ifelse(is.na(x),NA_integer_,as.integer(x=='Yes'))}
  dat <- data.frame(id=pid[treated],concentration=conc[match(pop[['Treatment Group']][treated],arms)],
    placebo=as.integer(pop[['Treatment Group']][treated]=='HEC placebo'),
    age=as.numeric(sub('/.*','',pop[['Age/ Race']][treated])),
    related_ae=as.integer(pid[treated]%in%related),mITT=pop$mITT[treated]=='Yes',
    cure_toc=binary(eff[['CC TOC']][ei]),cure_eot=binary(eff[['CC EOT']][ei]),stringsAsFactors=FALSE)
  stopifnot(nrow(dat)==125L,!anyNA(dat[c('concentration','placebo','age','related_ae','mITT')]),
    sum(dat$mITT)==111L,sum(dat$related_ae)==30L,
    all(pop[['Treatment Group']][treated]==eff$Treatment[ei]))

  missing_id <- eid[eff$mITT=='Yes' & is.na(eff[['CC TOC']])]
  stopifnot(identical(missing_id,c('10026','10029','14009','14038')))
  missing_index <- match(missing_id,dat$id)
  stopifnot(identical(dat$cure_eot[missing_index],c(1L,1L,1L,0L)))
  bits <- as.matrix(expand.grid(rep(list(0:1),4L)));storage.mode(bits)<-'integer'
  labels <- apply(bits,1L,paste0,collapse='')
  completion <- data.frame(model_id=paste0('efficacy_',labels),bits,LOCF=labels=='1110',stringsAsFactors=FALSE)
  names(completion)[2:5] <- paste0('missing_',missing_id)
  models <- c(list(safety=list(outcome='related_ae',index=seq_len(nrow(dat)),y=dat$related_ae)),
    list(efficacy_primary=list(outcome='cure_toc',index=which(dat$mITT & !is.na(dat$cure_toc)),
      y=dat$cure_toc[dat$mITT & !is.na(dat$cure_toc)])))
  for(i in seq_len(nrow(bits))) {
    y <- dat$cure_toc;y[missing_index] <- bits[i,]
    models[[paste0('efficacy_',labels[i])]] <- list(outcome='cure_toc',index=which(dat$mITT),y=y[dat$mITT])
  }
  stopifnot(length(models)==18L,length(models$efficacy_primary$y)==107L,
    all(vapply(models,function(m)all(m$y%in%0:1),TRUE)))
  ages <- sort(unique(dat$age));weights <- as.numeric(table(factor(dat$age,levels=ages)))/nrow(dat)
  prediction <- do.call(rbind,lapply(c(.5,1,3),function(v)data.frame(concentration=v,placebo=0,age=ages)))
  x <- as.matrix(dat[c('concentration','placebo','age')])
  midpoint <- function(x){v<-sort(unique(x));(v[-1]+v[-length(v)])/2}
  cuts <- lapply(seq_len(ncol(x)),function(j)midpoint(x[,j]));names(cuts)<-colnames(x)
  stopifnot(identical(cuts$concentration,c(.25,.75,2)),identical(cuts$placebo,.5),
    abs(sum(weights)-1)<1e-12)
  arm_summary <- do.call(rbind,lapply(conc,function(v) {
    d<-dat[dat$concentration==v,];e<-d[d$mITT & !is.na(d$cure_toc),]
    data.frame(concentration=v,safety_n=nrow(d),related_AE=sum(d$related_ae),mITT_n=sum(d$mITT),
      observed_primary_n=nrow(e),clinical_cure=sum(e$cure_toc),missing_primary=sum(d$mITT & is.na(d$cure_toc)))
  }))
  stopifnot(identical(arm_summary$safety_n,c(32L,32L,32L,29L)),
    identical(arm_summary$related_AE,c(7L,8L,6L,9L)),
    identical(arm_summary$observed_primary_n,c(26L,30L,26L,25L)),
    identical(arm_summary$clinical_cure,c(3L,7L,12L,7L)))
  input <- list(dat=dat,x=x,models=models,cuts=cuts,prediction=as.matrix(prediction),
    prediction_ages=ages,prediction_weights=weights,candidates=c(.5,1,3),
    missing_ids=missing_id,completions=completion,raw_sha256=setNames(vapply(paths,ast_sha,''),basename(paths)))
  dir.create(out,recursive=TRUE,showWarnings=FALSE)
  saveRDS(input,file.path(out,'INPUTS.rds'));write.csv(dat,file.path(out,'PARTICIPANTS.csv'),row.names=FALSE,na='')
  write.csv(completion,file.path(out,'MISSING_COMPLETIONS.csv'),row.names=FALSE)
  audit <- list(status='PASS',sampling=FALSE,participants=nrow(dat),models=names(models),
    predictor_columns=colnames(x),cuts=cuts,prediction_population='All 125 treated participants, empirical age weights fixed across active concentrations',
    prediction_ages=ages,prediction_weights=weights,missing_id_order=missing_id,LOCF_bits='1110',
    arm_summary=arm_summary,raw_sha256=input$raw_sha256,input_sha256=ast_sha(file.path(out,'INPUTS.rds')))
  ast_save_json(audit,file.path(out,'DATA_AUDIT.json'));invisible(input)
}

ast_fit_one <- function(input,plan,model_id,role) {
  stopifnot(model_id%in%names(input$models),role%in%seq_len(plan$chains),as.character(packageVersion('BART'))=='2.9.9')
  model_index <- match(model_id,names(input$models));seed <- as.integer(plan$seed_base+100L*model_index+role)
  RNGkind('Mersenne-Twister','Inversion','Rejection');set.seed(seed);rng0 <- .Random.seed
  model <- input$models[[model_id]];x <- input$x[model$index,,drop=FALSE]
  xi <- matrix(NA_real_,nrow=ncol(x),ncol=max(lengths(input$cuts)))
  for(j in seq_along(input$cuts))xi[j,seq_along(input$cuts[[j]])]<-input$cuts[[j]]
  warnings <- character();started <- unname(proc.time()[3])
  withCallingHandlers({
    invisible(capture.output(fit <- BART::pbart(x.train=x,y.train=model$y,x.test=input$prediction,
      xinfo=xi,rm.const=FALSE,ntree=as.integer(plan$ntree),k=plan$k,base=plan$base,power=plan$power,
      sparse=FALSE,augment=FALSE,binaryOffset=plan$offset,nskip=as.integer(plan$burnin),
      ndpost=as.integer(plan$retained),nkeeptrain=0L,nkeeptest=as.integer(plan$retained),
      nkeeptreedraws=0L,printevery=1000000L)))
  },warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')})
  stopifnot(identical(dim(fit$prob.test),c(as.integer(plan$retained),nrow(input$prediction))),
    all(is.finite(fit$prob.test)),all(fit$prob.test>=0&fit$prob.test<=1))
  H <- length(input$prediction_ages)
  averages <- vapply(seq_along(input$candidates),function(j) {
    columns <- (j-1L)*H+seq_len(H)
    as.vector(fit$prob.test[,columns,drop=FALSE]%*%input$prediction_weights)
  },numeric(plan$retained));colnames(averages)<-as.character(input$candidates)
  list(status='COMPLETE',model_id=model_id,role=role,seed=seed,initial_rng=rng0,final_rng=.Random.seed,
    training_n=nrow(x),training_events=sum(model$y),draws=averages,
    elapsed_seconds=unname(proc.time()[3])-started,warnings=warnings,
    variable_usage=fit$varcount.mean,BART_version=as.character(packageVersion('BART')))
}

ast_series_diagnostic <- function(x,delta,series,candidate,scenario,risk_limit=NA_real_) {
  safe<-function(fun,y)tryCatch(suppressWarnings(fun(y)),error=function(e)NA_real_)
  mn<-mean(x);rh<-safe(posterior::rhat,x);ess<-safe(posterior::ess_bulk,x);mc<-safe(posterior::mcse_mean,x)
  l<-safe(posterior::mcse_mean,x[,1:2,drop=FALSE]);r<-safe(posterior::mcse_mean,x[,3:4,drop=FALSE]);se<-sqrt(l*l+r*r)
  paired<-if(is.finite(se)&&se>0)abs(mean(x[,1:2])-mean(x[,3:4]))/se else NA_real_
  ok<-min(x)!=max(x)&&all(is.finite(c(rh,ess,mc,paired)))&&rh<=1.05&&ess>=100&&mc>0&&paired<=6
  mult<-qnorm(1-delta/2)
  data.frame(series=series,candidate=candidate,scenario=scenario,risk_limit=risk_limit,mean=mn,
    rhat=rh,ess=ess,mcse=mc,paired_z=paired,qualified=ok,
    lower=if(ok)max(0,mn-mult*mc)else 0,upper=if(ok)min(1,mn+mult*mc)else 1,
    constant_series=min(x)==max(x),stringsAsFactors=FALSE)
}

ast_possible <- function(gates,scores,cutoff=.3) {
  mandatory<-gates[,1]>=cutoff;possible<-gates[,2]>=cutoff
  answer<-if(any(mandatory))character()else 'dose:0'
  for(r in which(possible)) {
    other<-setdiff(which(mandatory),r)
    if(all(scores[r,2]>scores[other,1]|(scores[r,2]==scores[other,1]&r<other)))answer<-c(answer,paste0('dose:',r))
  }
  sort(answer)
}

ast_minimum_report <- function(gates,scores,target,cutoff=.3) {
  J<-nrow(gates);atoms<-unlist(lapply(seq_len(J),function(j)c(paste0('Tgate:',j),paste0('Emean:',j))))
  bands<-setNames(lapply(atoms,function(a){j<-as.integer(sub('.*:','',a));if(startsWith(a,'Tgate:'))gates[j,]else scores[j,]}),atoms)
  possibilities<-ast_possible(gates,scores,cutoff)
  if(length(possibilities)!=1L)return(list(status='UNRESOLVED',possible=possibilities,atoms=atoms,cost=length(atoms),minimum_established=FALSE))
  stopifnot(identical(possibilities,target))
  r<-as.integer(sub('^dose:','',target))
  if(r==0L) {
    stopifnot(all(gates[,2]<cutoff));formula_atoms<-paste0('Tgate:',seq_len(J))
  } else {
    stopifnot(gates[r,1]>=cutoff)
    branches<-list()
    for(retain_target_mean in c(FALSE,TRUE)) {
      lower<-if(retain_target_mean)scores[r,1]else 0
      selected<-c(paste0('Tgate:',r),if(retain_target_mean)paste0('Emean:',r))
      possible_branch<-TRUE
      for(j in setdiff(seq_len(J),r)) {
        options<-list()
        if(gates[j,2]<cutoff)options[[length(options)+1L]]<-paste0('Tgate:',j)
        if(lower==1&&r<j)options[[length(options)+1L]]<-character()
        if(lower>scores[j,2]||(lower==scores[j,2]&&r<j))
          options[[length(options)+1L]]<-paste0('Emean:',j)
        if(!length(options)){possible_branch<-FALSE;break}
        selected<-c(selected,options[[which.min(lengths(options))]])
      }
      if(possible_branch)branches[[length(branches)+1L]]<-selected
    }
    stopifnot(length(branches)>0L)
    formula_atoms<-branches[[which.min(lengths(branches))]]
  }
  formula_atoms<-atoms[atoms%in%formula_atoms]
  subsets<-as.matrix(expand.grid(rep(list(c(FALSE,TRUE)),length(atoms))));best<-Inf;selected<-NULL
  for(i in seq_len(nrow(subsets))) {
    keep<-atoms[subsets[i,]];if(length(keep)>best)next
    gg<-gates;ss<-scores
    for(j in seq_len(J)) {
      if(!paste0('Tgate:',j)%in%keep)gg[j,]<-c(0,1)
      if(!paste0('Emean:',j)%in%keep)ss[j,]<-c(0,1)
    }
    if(identical(ast_possible(gg,ss,cutoff),target)&&length(keep)<best){best<-length(keep);selected<-keep}
  }
  stopifnot(is.finite(best),length(formula_atoms)==best)
  gg<-gates;ss<-scores
  for(j in seq_len(J)) {
    if(!paste0('Tgate:',j)%in%formula_atoms)gg[j,]<-c(0,1)
    if(!paste0('Emean:',j)%in%formula_atoms)ss[j,]<-c(0,1)
  }
  stopifnot(identical(ast_possible(gg,ss,cutoff),target))
  list(status='CERTIFIED',possible=possibilities,atoms=formula_atoms,cost=length(formula_atoms),
    formula_cost=length(formula_atoms),exhaustive_cost=best,exhaustive_atoms=selected,
    minimum_established=TRUE,subsets_evaluated=nrow(subsets))
}

ast_change_costs <- function(counts,means,M,cutoff=.3,eta=0) {
  kappa<-ceiling(M*cutoff);keep<-counts>=kappa
  rho<-ifelse(keep,counts-kappa+1L,kappa-counts)
  lo<-pmax(0,means-eta);hi<-pmin(1,means+eta);J<-length(means)
  costs<-vapply(seq_len(J),function(r) {
    blockers<-setdiff(which(keep&(lo>hi[r]|(lo==hi[r]&seq_len(J)<r))),r)
    as.numeric(if(keep[r])0 else rho[r])+sum(rho[blockers])
  },0)
  data.frame(target=c('dose:0',paste0('dose:',seq_len(J))),eta=eta,
    screening_changes=c(sum(rho[keep]),costs),total_draws=M,screen_count_cutoff=kappa,stringsAsFactors=FALSE)
}

ast_summarize <- function(input,plan,chains,out) {
  delta<-plan$alpha/plan$total_quantities;J<-length(input$candidates)
  stopifnot(length(chains)==length(input$models),J==3L,plan$total_quantities==57L)
  matrix_for<-function(model,j)sapply(chains[[model]],function(z)z$draws[,j])
  diag<-list();k<-0L
  for(j in seq_len(J)) {
    x<-matrix_for('safety',j)
    for(type in c('T_mean','T_gate')) {
      k<-k+1L;diag[[k]]<-ast_series_diagnostic(if(type=='T_mean')x else 1*(x<plan$risk_limit),
        delta,type,j,'shared_safety',plan$risk_limit)
    }
  }
  efficacy_models<-setdiff(names(input$models),'safety')
  for(model in efficacy_models)for(j in seq_len(J)) {
    k<-k+1L;diag[[k]]<-ast_series_diagnostic(matrix_for(model,j),delta,'E_mean',j,model)
  }
  diagnostics<-do.call(rbind,diag);stopifnot(nrow(diagnostics)==57L)
  getrow<-function(type,j,scenario)diagnostics[diagnostics$series==type&diagnostics$candidate==j&diagnostics$scenario==scenario,,drop=FALSE]
  gate<-t(vapply(seq_len(J),function(j) {
    g<-getrow('T_gate',j,'shared_safety');m<-getrow('T_mean',j,'shared_safety')
    bounds<-c(max(0,1-m$upper/plan$risk_limit),min(1,(1-m$lower)/(1-plan$risk_limit)))
    z<-c(max(g$lower,bounds[1]),min(g$upper,bounds[2]));stopifnot(z[1]<=z[2]+1e-12);z
  },numeric(2)))
  gate_means<-vapply(seq_len(J),function(j)getrow('T_gate',j,'shared_safety')$mean,0)
  M<-plan$retained*plan$chains
  counts<-vapply(seq_len(J),function(j)sum(matrix_for('safety',j)<plan$risk_limit),0)
  summaries<-reports<-costs<-list()
  for(i in seq_along(efficacy_models)) {
    model<-efficacy_models[i]
    e<-lapply(seq_len(J),function(j)getrow('E_mean',j,model));means<-vapply(e,function(z)z$mean,0)
    scores<-t(vapply(e,function(z)c(z$lower,z$upper),numeric(2)))
    possible<-ast_possible(gate,scores,plan$screen_cut)
    eligible<-counts>=ceiling(M*plan$screen_cut)
    point<-if(any(eligible))which(eligible)[which.max(means[eligible])]else 0L
    report<-ast_minimum_report(gate,scores,if(length(possible)==1L)possible else NULL,plan$screen_cut)
    report$model_id<-model;report$screening_intervals<-gate;report$efficacy_intervals<-scores
    reports[[model]]<-report
    cc<-ast_change_costs(counts,means,M,plan$screen_cut);cc$model_id<-model;costs[[model]]<-cc
    for(j in seq_len(J))summaries[[length(summaries)+1L]]<-data.frame(model_id=model,candidate=j,
      concentration=input$candidates[j],AE_posterior_mean=getrow('T_mean',j,'shared_safety')$mean,
      AE_lower=getrow('T_mean',j,'shared_safety')$lower,AE_upper=getrow('T_mean',j,'shared_safety')$upper,
      safety_probability=gate_means[j],safety_lower=gate[j,1],safety_upper=gate[j,2],
      cure_posterior_mean=means[j],cure_lower=scores[j,1],cure_upper=scores[j,2],
      fitted_recommendation=paste0('dose:',point),possible_decisions=paste(possible,collapse=';'),
      established=length(possible)==1L,report_length=report$cost,
      LOCF=model=='efficacy_1110',stringsAsFactors=FALSE)
  }
  dir.create(out,recursive=TRUE,showWarnings=FALSE)
  write.csv(diagnostics,file.path(out,'DIAGNOSTICS.csv'),row.names=FALSE,na='')
  write.csv(do.call(rbind,summaries),file.path(out,'POSTERIOR_SUMMARY.csv'),row.names=FALSE)
  write.csv(do.call(rbind,costs),file.path(out,'SCREENING_CHANGE_COSTS.csv'),row.names=FALSE)
  ast_save_json(reports,file.path(out,'INTERVAL_REPORTS.json'))
  result<-list(status='COMPLETE',models=length(input$models),efficacy_scenarios=length(efficacy_models),
    diagnostics_qualified=sum(diagnostics$qualified),diagnostics_total=nrow(diagnostics),
    nominal_family_error=plan$alpha,scalar_error=delta,retained_draws_per_model=M,
    risk_limit=plan$risk_limit,screen_cut=plan$screen_cut,
    scope='Asymptotic MCMC numerical intervals for fixed retrospective model/decision rules; no independent reference or clinical dose-selection validation',
    reports=reports)
  ast_save_json(result,file.path(out,'RESULTS.json'));invisible(result)
}
