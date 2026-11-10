

script <- sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE))
root <- normalizePath(file.path(dirname(script),'..'))
argv <- commandArgs(trailingOnly=TRUE)
output <- if(length(argv))argv[1]else file.path(root,'output/review3_intervals.png')
dir.create(dirname(output),recursive=TRUE,showWarnings=FALSE)
source(file.path(root,'R/reference.R'))
folder <- file.path(root,'data/R_reviewer_comparisons')
plan <- jsonlite::fromJSON(file.path(folder,'PLAN.json'))
stopifnot(plan$budget_scalar_error==.01/60)
rows <- list()
for(stem in c('T4','E1','E4')) {
  d <- readRDS(file.path(folder,'budget',paste0(stem,'.rds')))
  for(b in d$reports) {
    q <- d$q[d$accepted_at<=b$proposals]
    z <- probability_intervals(q,substr(stem,1,1),plan$budget_scalar_error)
    for(n in c('accepted','gate_count','mean','gate_lower','gate_upper','mean_lower','mean_upper'))
      stopifnot(isTRUE(all.equal(as.numeric(z[[n]]),as.numeric(b[[n]]),tolerance=1e-12)))
    rows[[length(rows)+1L]] <- data.frame(endpoint=substr(stem,1,1),dose=as.integer(substr(stem,2,2)),
      proposals=b$proposals,accepted=b$accepted,gate_lower=z$gate_lower,gate_upper=z$gate_upper,
      mean_lower=z$mean_lower,mean_upper=z$mean_upper)
  }
}
all <- do.call(rbind,rows)
a <- all[all$endpoint=='T',];a <- a[order(a$proposals),]
b <- all[all$endpoint=='E' & all$proposals==max(all$proposals),];b <- b[order(b$dose),]
stopifnot(nrow(a)==3L,nrow(b)==2L,all(a$gate_lower<.30 & a$gate_upper>.30),
  max(b$mean_lower)<min(b$mean_upper))
write.csv(all,sub('[.]png$','.csv',output),row.names=FALSE)
library(ggplot2)
source(file.path(root,'R/plot_style.R'))
safety <- ggplot(a,aes(x=proposals/1e6,ymin=gate_lower,ymax=gate_upper)) +
  geom_errorbar(width=.07,linewidth=.55,colour=review_colours['grey']) +
  geom_hline(yintercept=.30,linetype='dashed',linewidth=.4,colour=review_colours['ink']) +
  scale_x_continuous(breaks=a$proposals/1e6) +
  scale_y_continuous(limits=c(0,1),breaks=seq(0,1,.25)) +
  labs(title='(a)',
       x='Proposals per endpoint (millions)',y='Screening probability') + review_theme()
efficacy <- ggplot(b,aes(x=factor(dose,levels=c(4,1),labels=c('Dose 4','Dose 1')),
                       ymin=mean_lower,ymax=mean_upper,colour=factor(dose))) +
  geom_errorbar(width=.12,linewidth=.55) +
  scale_colour_manual(values=c('1'=unname(review_colours['blue']),'4'=unname(review_colours['grey']))) +
  coord_flip() +
  labs(title='(b)',
       x=NULL,y='Posterior mean efficacy') + review_theme() + theme(legend.position='none')
save_review_panels(list(safety,efficacy),sub('[.]png$','',output),6.8,3.35)
stopifnot(file.exists(output),file.info(output)$size>10000)
cat('Figure intervals independently reconstructed from accepted draws\n')
