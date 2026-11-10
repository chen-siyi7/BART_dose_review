library(bartDoseReview)

gates <- rep(list(rbind(c(0.7,0.9),c(0.1,0.2))),3)
scores <- rbind(c(0.1,0.2),c(0.7,0.8),c(0.4,0.5))
stopifnot(identical(interval_actions('Final',gates,scores),'dose:2'))
tie <- rbind(c(0.5,0.5),c(0.5,0.5),c(0.1,0.2))
stopifnot(identical(interval_actions('Final',gates,tie),'dose:1'))
failed <- rep(list(rbind(c(0.1,0.2),c(0.7,0.8))),3)
stopifnot(identical(interval_actions('Final',failed,scores),'dose:0'))
uncertain <- rbind(c(0.1,0.9),c(0.1,0.9),c(0.1,0.9))
stopifnot(setequal(interval_actions('Final',gates,uncertain),paste0('dose:',1:3)))

record <- list(stage='Final',target='dose:2',gates=gates,
  scores=lapply(seq_len(nrow(scores)),function(i)scores[i,]))
z <- prepare_certificate(record)
costs <- setNames(seq_along(z$atoms),z$atoms)
best <- minimum_certificate(z,costs)
oracle <- getFromNamespace('exhaustive_certificate','bartDoseReview')(z,costs)
stopifnot(best$status=='CERTIFIED',best$cost==oracle$cost,
  identical(certificate_actions(z,best$atoms),'dose:2'))
for(atom in best$atoms)
  stopifnot(!identical(certificate_actions(z,setdiff(best$atoms,atom)),'dose:2'))
record$scores <- rep(list(c(0,1)),3)
stopifnot(minimum_certificate(prepare_certificate(record))$status=='UNCERTIFIABLE')

counts <- matrix(c(301,301),ncol=1)
cost <- frontier_costs('Final',1:2,counts,c(0.4,0.3))
witness <- frontier_witness(1:2,counts,c(0.4,0.3),2)
stopifnot(cost[['dose:2']]==2,cost[['dose:0']]==4,witness$cost==2,
  sum(abs(witness$counts-counts))==2,
  frontier_costs('Final',1:2,counts,c(0.4,0.3),eta=0.06)[['dose:2']]==0)

X <- matrix(c(0,1),ncol=1)
law <- catalog_law(X,list(0.5),c(5,5),exact=TRUE)
stopifnot(Reduce(`+`,law$probabilities)==gmp::as.bigq(1))
limited <- catalog_law_limited(X,list(0.5),c(5,5),max_states=0L)
stopifnot(limited$status=='RESOURCE_LIMIT',is.null(limited$law))
guarded <- catalog_law_limited(X,list(0.5),c(5,5))
stopifnot(guarded$status=='COMPLETE',
  max(abs(vapply(law$probabilities,as.numeric,0)-unlist(guarded$law$probabilities)))<1e-12)

law <- catalog_law(matrix(0,1,1),list(0.5),12)
a <- prepare_reference(law,1)
logf <- function(x)dnorm(x,sd=1.5,log=TRUE)+7*pnorm(x,log.p=TRUE)+5*pnorm(-x,log.p=TRUE)
shift <- optimize(logf,c(-10,10),maximum=TRUE)$objective
f <- function(x)exp(logf(x)-shift)
den <- integrate(f,-Inf,Inf,rel.tol=1e-10)$value
truth <- integrate(function(x)pnorm(x)*f(x),-Inf,Inf,rel.tol=1e-10)$value/den
for(method in c('direct','aggregated')) {
  set.seed(314L+as.integer(method=='aggregated'))
  draws <- sample_reference(15000,512,2,a,7,12,1,0,method)$q
  stopifnot(length(draws)>100,abs(mean(draws)-truth)<6*sd(draws)/sqrt(length(draws)))
  band <- probability_intervals(draws,'E',0.001)
  stopifnot(truth>=band$mean_lower,truth<=band$mean_upper)
}
