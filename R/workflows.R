.validate_review_bundle <- function(bundle) {
  bundle <- normalizePath(bundle,mustWork=TRUE)
  if(!file.exists(file.path(bundle,'MANIFEST.json')))
    stop('Use the full reproducibility directory containing MANIFEST.json and data/.')
  actual <- setNames(vapply(names(.review_source_hashes),function(n) {
    p <- file.path(bundle,n)
    if(file.exists(p))digest::digest(file=p,algo='sha256')else NA_character_
  },''),names(.review_source_hashes))
  if(!identical(actual,.review_source_hashes) && !identical(actual,.review_original_hashes))
    stop('The reproducibility bundle does not match a recorded source version for this package.')
  verify_manifest(bundle)
  bundle
}

run_analysis <- function(analysis,bundle,out,...) {
  analysis <- match.arg(analysis,c('replay','benchmark','feasibility','reviewer_comparisons',
    'intercept','history_compression','catalog_free','figures','review_intervals','clinical_replay'))
  bundle <- .validate_review_bundle(bundle)
  if(length(out)!=1L || is.na(out) || !nzchar(out) || file.exists(out) || !dir.exists(dirname(out)))
    stop('Choose a new output path whose parent directory exists.')
  parent <- normalizePath(dirname(out))
  installed <- system.file(package='bartDoseReview')
  for(folder in c(bundle,installed))
    if(nzchar(folder) && (parent==folder || startsWith(parent,paste0(folder,'/'))))
      stop('Choose an output path outside the reproducibility bundle and package.')
  out <- file.path(parent,basename(out))
  script <- switch(analysis,replay='replay_all.R',benchmark='scripts/benchmark.R',
    feasibility='scripts/feasibility_study.R',reviewer_comparisons='scripts/reviewer_comparisons.R',
    intercept='scripts/intercept_simulation.R',history_compression='scripts/history_compression.R',
    catalog_free='scripts/catalog_free_extension.R',figures='scripts/plot_evidence.R',
    review_intervals='scripts/plot_review3.R',clinical_replay='scripts/replay_astodrimer.R')
  extra <- unlist(list(...),use.names=FALSE)
  if(!is.null(extra) && (!is.character(extra) || anyNA(extra)))
    stop('Additional script arguments must be character strings.')
  arguments <- if(analysis%in%c('history_compression','catalog_free'))c(bundle,out) else out
  status <- system2(file.path(R.home('bin'),'Rscript'),
    shQuote(c(file.path(bundle,script),arguments,extra)))
  if(status!=0L)stop('The ',analysis,' script failed with status ',status,'.')
  invisible(out)
}

replay_results <- function(bundle,out) {
  path <- run_analysis('replay',bundle,out)
  result <- jsonlite::read_json(file.path(path,'RESULTS.json'))
  stopifnot(identical(result$status,'PASS'),identical(result$new_fits,FALSE),
    identical(result$new_timing,FALSE))
  invisible(result)
}
