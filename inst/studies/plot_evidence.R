

script <- sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE))
root <- normalizePath(file.path(dirname(script),'..'))
args <- commandArgs(trailingOnly=TRUE);stopifnot(length(args)==1L)
out <- args[1L];dir.create(out,recursive=TRUE,showWarnings=FALSE)
library(ggplot2)
source(file.path(root,'R/plot_style.R'))


reports <- read.csv(file.path(root,'data/R_history_compression/REPORT_COMPRESSION.csv'))
history <- read.csv(file.path(root,'data/R_history_compression/HISTORY_COMPRESSION.csv'))
reports <- reports[reports$certified,]
groups <- split(reports,reports$history)
compression <- do.call(rbind,lapply(groups,function(d)data.frame(
 history=d$history[1],bank=d$bank[1],reports=nrow(d),
 full_mean=mean(d$full),minimum_mean=mean(d$minimum),
 fraction_saved=mean(d$fraction_saved))))
rownames(compression) <- NULL
check <- history[history$scope=='all',]
idx <- match(compression$history,check$history)
stopifnot(nrow(compression)==32L,all(!is.na(idx)),sum(compression$reports)==223L,
 max(abs(compression$fraction_saved-check$mean_fraction_saved[idx]))<1e-12,
 all(compression$minimum_mean<compression$full_mean))
write.csv(compression,file.path(out,'report_shortening.csv'),row.names=FALSE)
compression$bank <- factor(compression$bank,levels=c('small','deeper','crossed'),
 labels=c('Small forest','Repeated grid','Crossed grids'))
compression$offset <- ave(seq_len(nrow(compression)),compression$bank,
                          FUN=function(x)seq(-.03,.03,length.out=length(x)))
shortening <- ggplot(compression) +
 geom_segment(aes(x=1+offset,xend=2+offset,y=full_mean,yend=minimum_mean),
              linewidth=.35,alpha=.45,colour=review_colours['grey']) +
 geom_point(aes(x=1+offset,y=full_mean),shape=1,size=1.8,stroke=.5,colour=review_colours['grey']) +
 geom_point(aes(x=2+offset,y=minimum_mean),size=1.8,colour=review_colours['blue']) +
 facet_wrap(~bank,nrow=1) +
 scale_x_continuous(breaks=c(1,2),labels=c('Full\nsummary','Shortest\nreport'),limits=c(.7,2.3)) +
 scale_y_continuous(limits=c(0,11),breaks=seq(0,10,2)) +
 labs(x=NULL,y='Mean number of intervals') + review_theme() +
 theme(legend.position='none',panel.grid.major.x=element_blank())
save_review_plot(shortening,file.path(out,'report_shortening'),6.8,3.1)


app <- file.path(root,'data/R_astodrimer_application')
draws <- do.call(rbind,lapply(1:4,function(k)
 readRDS(file.path(app,'chains',paste0('efficacy_primary_chain',k,'.rds')))$draws))
stopifnot(identical(dim(draws),c(128000L,3L)))
summary <- read.csv(file.path(app,'descriptive_context/CURE_CREDIBLE_INTERVALS.csv'))
summary <- summary[summary$model_id=='efficacy_primary',]
diff <- read.csv(file.path(app,'descriptive_context/PAIRED_DIFFERENCES.csv'))
diff <- diff[diff$model_id=='efficacy_primary',]
stopifnot(identical(summary$concentration,c(.5,1,3)),
 identical(diff$comparator_concentration,c(.5,3)))
for(j in 1:3)stopifnot(abs(mean(draws[,j])-summary$posterior_mean[j])<1e-12,
 max(abs(quantile(draws[,j],c(.025,.975),type=7)-
         unlist(summary[j,c('credible_lower','credible_upper')],use.names=FALSE)))<1e-12)
for(j in 1:2) {
 d <- draws[,2]-draws[,c(1,3)[j]]
 stopifnot(abs(mean(d)-diff$posterior_mean_difference[j])<1e-12,
  max(abs(quantile(d,c(.025,.975),type=7)-
          unlist(diff[j,c('credible_lower','credible_upper')],use.names=FALSE)))<1e-12)
}
write.csv(summary,file.path(out,'clinical_cure_summary.csv'),row.names=FALSE)
write.csv(diff,file.path(out,'clinical_cure_differences.csv'),row.names=FALSE)
labels <- c('0.5%','1%','3%')
summary$concentration_label <- factor(labels,levels=labels)
dose_colours <- setNames(unname(review_colours[c('grey','blue','grey')]),labels)
cure <- ggplot(summary,aes(x=concentration_label,y=posterior_mean,colour=concentration_label)) +
 geom_errorbar(aes(ymin=credible_lower,ymax=credible_upper),width=.12,linewidth=.55) +
 geom_point(size=2) + scale_colour_manual(values=dose_colours) +
 scale_y_continuous(limits=c(0,1),breaks=seq(0,1,.25)) +
 labs(title='(a)',
      x='Astodrimer concentration',y='Cure probability') + review_theme() +
 theme(legend.position='none',panel.grid.major.x=element_blank())
diff$comparison <- factor(c('1% vs 0.5%','1% vs 3%'),levels=c('1% vs 3%','1% vs 0.5%'))
comparisons <- ggplot(diff,aes(y=comparison,x=100*posterior_mean_difference)) +
 geom_vline(xintercept=0,linetype='dashed',colour=review_colours['ink'],linewidth=.4) +
 geom_errorbar(aes(xmin=100*credible_lower,xmax=100*credible_upper),
               orientation='y',width=.12,linewidth=.55,colour=review_colours['grey']) +
 geom_point(size=2,colour=review_colours['blue']) +
 scale_x_continuous(breaks=c(-10,0,10,20,30,40),limits=c(-12,43)) +
 labs(title='(b)',
      x='Difference (percentage points)',y=NULL) + review_theme() +
 theme(legend.position='none',panel.grid.major.y=element_blank())
save_review_panels(list(cure,comparisons),file.path(out,'clinical_cure'),6.8,3.45,widths=c(1,1.15))
jsonlite::write_json(list(status='PASS',history_pairs=32L,reports=223L,
 summaries_checked_against_saved_records=TRUE,primary_posterior_draws=nrow(draws),
 posterior_intervals_are_not_Monte_Carlo_intervals=TRUE,plot_engine='ggplot2',
 colour_palette=as.list(review_colours)),file.path(out,'EVIDENCE_FIGURE_AUDIT.json'),
 pretty=TRUE,auto_unbox=TRUE,digits=17)
cat('Report-shortening and clinical figures reconstructed from saved results\n')
