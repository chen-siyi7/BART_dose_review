query_mask_cost <- function(doses,gates,q,target,eta,joint=TRUE) {
  q<-gmp::as.bigq(q);eta<-gmp::as.bigq(eta);J<-length(doses)
  keep<-apply(gates,1L,function(g)all(g>=300))
  rho<-vapply(seq_len(J),function(j)if(keep[j])min(gates[j,]-299)else sum(pmax(0,300-gates[j,])),0)
  masks<-all_masks(J);best<-Inf
  for(i in seq_len(nrow(masks))) {
    b<-masks[i,]
    if(target==0)feasible<-!any(b) else {
      r<-match(target,doses);others<-setdiff(which(b),r)
      upper<-min(gmp::as.bigq(1),q[r]+eta)
      lower<-q[others]-eta;lower[lower<0]<-0
      feasible<-b[r]&&all(upper>lower|(upper==lower&r<others))
    }
    if(feasible)best<-min(best,if(joint)sum(rho[b!=keep])else max(c(0,rho[b!=keep])))
  }
  best
}

audit_query_boundaries <- function() {
  patterns<-list(c(0,0,0),c(1,1,1),c(0,.5,1),c(.5,.5,0),
    c(.8,.2,.5),c(.1,.3,.7))
  checked<-0L;events<-0L
  for(q in patterns)for(bits in seq_len(8L)) {
    gates<-matrix(ifelse(all_masks(3L)[bits,],c(300,301,1000),c(299,0,298)),ncol=1L)
    doses<-c(9,2,7);index<-prepare_score_index(doses,gates,q)
    curves<-frontier_curve(doses,gates,q);events<-events+nrow(curves)
    points<-sort(unique(c(gmp::as.bigq(curves$eta),gmp::as.bigq(c(0,1)))))
    probes<-points
    if(length(points)>1L)probes<-c(probes,(points[-1L]+points[-length(points)])/2)
    for(p in seq_along(probes))for(target in c(0,doses)) {
      eta<-probes[p]
      expected<-query_mask_cost(doses,gates,q,target,eta)
      direct<-frontier_costs('Final',doses,gates,q,eta)[[paste0('dose:',target)]]
      witness<-frontier_witness(doses,gates,q,target,eta)
      stopifnot(indexed_cost(index,target,eta)==expected,curve_cost(curves,target,eta)==expected,
        direct==expected,witness$cost==expected)
      checked<-checked+1L
    }
  }
  list(status='PASS',patterns=length(patterns)*8L,target_queries=checked,stored_events=events)
}

audit_queries <- function(data,out) {
  folder<-file.path(data,'base/query/report')
  raw<-file.path(data,'base/query/source/ordinary_review_benchmark_20261102/raw')
  readtable<-function(n)read.csv(file.path(folder,n),stringsAsFactors=FALSE)
  disposition<-readtable('HISTORY_DISPOSITIONS.csv')
  records<-setNames(lapply(disposition$case_id,function(id) {
    a<-read_json(file.path(raw,sprintf('case%04d_role1.json',id)))
    b<-read_json(file.path(raw,sprintf('case%04d_role2.json',id)))
    list(anchor=a,comparison=b,values=return_values(a))
  }),disposition$case_id)
  for(i in seq_len(nrow(disposition))) {
    d<-disposition[i,];r<-records[[as.character(d$case_id)]]
    stopifnot(r$anchor$status==d$anchor_status,r$comparison$status==d$comparison_status,
      r$anchor$action==d$anchor_action,r$comparison$action==d$returned_action)
  }
  costs<-readtable('COST_QUERIES.csv');events<-readtable('FRONTIER_EVENTS.csv')
  total_events<-0L;event_checks<-0L;cost_checks<-0L;witness_checks<-0L
  for(id in unique(costs$case_id)) {
    v<-records[[as.character(id)]]$values
    index<-prepare_score_index(v$doses,v$gates,v$q)
    curves<-frontier_curve(v$doses,v$gates,v$q)
    total_events<-total_events+nrow(curves)
    label<-function(x)if(x=='none')0 else v$doses[as.integer(x)+1L]
    e<-events[events$case_id==id,,drop=FALSE]
    stopifnot(nrow(curves)==nrow(e))
    for(i in seq_len(nrow(e))) {
      target<-label(e$target[i]);eta<-gmp::as.bigq(e$eta[i])
      row<-curves[curves$target==target&as.logical(gmp::as.bigq(curves$eta)==eta),,drop=FALSE]
      probe<-gmp::as.bigq(e$right_probe[i])
      stopifnot(nrow(row)==1L,row$at_cost==fraction(e$at_cost[i]),
        row$right_cost==fraction(e$right_cost[i]),
        indexed_cost(index,target,eta)==row$at_cost,
        query_mask_cost(v$doses,v$gates,v$q,target,eta)==row$at_cost,
        indexed_cost(index,target,probe)==row$right_cost,
        curve_cost(curves,target,probe)==row$right_cost,
        query_mask_cost(v$doses,v$gates,v$q,target,probe)==row$right_cost)
      event_checks<-event_checks+1L
    }
    rows<-costs[costs$case_id==id,,drop=FALSE]
    for(i in seq_len(nrow(rows))) {
      r<-rows[i,];target<-label(r$target);eta<-gmp::as.bigq(r$eta)
      actual<-c(curve_cost(curves,target,eta),indexed_cost(index,target,eta),
        query_mask_cost(v$doses,v$gates,v$q,target,eta),
        frontier_costs('Final',v$doses,v$gates,v$q,eta)[[paste0('dose:',target)]])
      stopifnot(all(actual==fraction(r$frontier_cost)),
        all(actual==fraction(r$index_cost)),all(actual==fraction(r$scan_cost)),
        query_mask_cost(v$doses,v$gates,v$q,target,eta,FALSE)==fraction(r$coordinate_cost))
      positions<-if(is.na(r$witness_mask)||!nzchar(r$witness_mask))integer()else
        as.integer(strsplit(r$witness_mask,':',fixed=TRUE)[[1L]])+1L
      scores<-gmp::as.bigq(strsplit(r$witness_scores,':',fixed=TRUE)[[1L]])
      winner<-if(length(positions))v$doses[positions[which(scores[positions]==max(scores[positions]))[1L]]] else 0
      z<-membership_costs(v$gates);mask<-seq_along(v$doses)%in%positions
      if(is.na(winner)||winner!=target)stop('Saved witness winner mismatch: case ',id,', target ',target,
        ', eta ',r$eta,', positions ',paste(positions,collapse=':'),', returned ',winner)
      stopifnot(sum(z$rho[mask!=z$keep])==actual[1L],
        all(scores>=0),all(scores<=1),all(abs(scores-gmp::as.bigq(v$q))<=eta))
      cost_checks<-cost_checks+1L;witness_checks<-witness_checks+1L
    }
  }
  replay_sets<-function(rows,name) {
    for(i in seq_len(nrow(rows))) {
      r<-rows[i,];record<-records[[as.character(r$case_id)]];v<-record$values
      budget<-fraction(r$budget);eta<-fraction(r$eta)
      exact<-frontier(record$anchor$stage,v$doses,v$gates,v$q,budget,eta)
      coordinate<-coordinate_actions(record$anchor$stage,v$doses,v$gates,v$q,budget,eta)
      stopifnot(setequal(exact,unlist(jsonlite::fromJSON(r$exact_actions))),
        setequal(coordinate,unlist(jsonlite::fromJSON(r$coordinate_actions))),
        length(exact)==r$exact_size,length(coordinate)==r$coordinate_size,
        length(action_alphabet(record$anchor$stage,v$doses))==r$alphabet_size,
        (record$comparison$action%in%exact)==(r$exact_contains_return=='True'),
        (record$comparison$action%in%coordinate)==(r$coordinate_contains_return=='True'))
    }
    write.csv(rows,file.path(out,name),row.names=FALSE)
    nrow(rows)
  }
  fixed<-readtable('FIXED_ALLOWANCE_RECORDS.csv')
  multiplied<-readtable('CALIBRATED_SENSITIVITY.csv')
  thresholds<-read.csv(file.path(dirname(raw),'CALIBRATION_THRESHOLDS.csv'),stringsAsFactors=FALSE)
  for(i in seq_len(nrow(multiplied))) {
    r<-multiplied[i,];a<-records[[as.character(r$case_id)]]$anchor
    t<-thresholds[thresholds$layout==r$layout&thresholds$scenario==r$scenario,,drop=FALSE]
    stopifnot(nrow(t)==1L)
    allowance<-gmp::as.bigq(t$threshold)*gmp::as.bigq(r$multiplier)
    count<-allowance*1000
    stopifnot(gmp::as.bigq(r$budget)==gmp::numerator(count)%/%gmp::denominator(count),
      gmp::as.bigq(r$eta)==if(a$stage=='Final'&&length(a$candidates)>1L)
        min(gmp::as.bigq(1),allowance)else gmp::as.bigq(0))
  }
  fixed_n<-replay_sets(fixed,'QUERY_FIXED_ALLOWANCE_RECORDS.csv')
  multiplied_n<-replay_sets(multiplied,'QUERY_CALIBRATED_SENSITIVITY.csv')
  summaries<-0L
  summarize<-function(rows,saved,keys) {
    fields<-c('n','exact_singletons','coordinate_singletons','exact_full','coordinate_excess',
      'exact_excludes_return','coordinate_excludes_return','mean_exact_size','mean_coordinate_size')
    for(i in seq_len(nrow(saved))) {
      r<-saved[i,];take<-rep(TRUE,nrow(rows))
      for(k in keys)take<-take&as.character(rows[[k]])==as.character(r[[k]])
      if(r$stage!='all')take<-take&if(r$stage=='Final')rows$stage=='Final'else rows$stage!='Final'
      z<-rows[take,,drop=FALSE]
      actual<-c(nrow(z),sum(z$exact_size==1),sum(z$coordinate_size==1),
        sum(z$exact_size==z$alphabet_size),sum(z$coordinate_size>z$exact_size),
        sum(z$exact_contains_return=='False'),sum(z$coordinate_contains_return=='False'),
        mean(z$exact_size),mean(z$coordinate_size))
      stopifnot(isTRUE(all.equal(actual,as.numeric(r[fields]),tolerance=1e-12)))
    }
    nrow(saved)
  }
  summaries<-summaries+summarize(fixed,readtable('FIXED_FINAL_SUMMARY.csv'),c('level_b','eta'))
  summaries<-summaries+summarize(multiplied,readtable('SENSITIVITY_SUMMARY.csv'),'multiplier')
  write.csv(costs,file.path(out,'QUERY_COSTS.csv'),row.names=FALSE)
  write.csv(events,file.path(out,'QUERY_EVENTS.csv'),row.names=FALSE)
  list(status='PASS',language='R',histories=nrow(disposition),final_histories=length(unique(costs$case_id)),
    target_queries=cost_checks,attaining_witnesses=witness_checks,events=event_checks,
    sparse_events=total_events,fixed_allowance_records=fixed_n,multiplied_allowance_records=multiplied_n,
    reconstructed_summaries=summaries,boundary_tests=audit_query_boundaries(),
    new_fits=FALSE,new_timing=FALSE,
    timing_provenance='Historical query timings are preserved; this replay verifies numerical results only.')
}
