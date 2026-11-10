



for (pkg in c('jsonlite', 'digest')) stopifnot(requireNamespace(pkg, quietly = TRUE))
args <- commandArgs(trailingOnly = TRUE); stopifnot(length(args) <= 2L)
script <- sub('^--file=', '', grep('^--file=', commandArgs(), value = TRUE))
app <- normalizePath(if (length(args)) args[1] else file.path(dirname(script), '..'))
out <- if (length(args) > 1L) args[2] else file.path(app, 'descriptive_context')
dir.create(out, showWarnings = FALSE, recursive = TRUE)
sha <- function(p) digest::digest(file = p, algo = 'sha256')
plan_file <- file.path(app, 'PLAN.json'); plan <- jsonlite::fromJSON(plan_file)
input_file <- file.path(app, 'inputs', 'INPUTS.rds'); input <- readRDS(input_file)
stopifnot(identical(plan$input_sha256, sha(input_file)))
plan_sha <- sha(plan_file); candidates <- c(.5, 1, 3)
model_ids <- setdiff(names(input$models), 'safety')
stopifnot(length(model_ids) == 17L, identical(input$candidates, candidates))
cure <- differences <- scenarios <- list(); source_sha <- list()
interval <- function(x) unname(quantile(x, probs = c(.025, .975), type = 7))
for (id in model_ids) {
  draws <- do.call(rbind, lapply(1:4, function(role) {
    f <- file.path(app, 'chains', paste0(id, '_chain', role, '.rds')); z <- readRDS(f)
    stopifnot(z$status == 'COMPLETE', z$model_id == id, z$role == role,
      identical(z$plan_sha256, plan_sha), identical(z$input_sha256, plan$input_sha256),
      identical(dim(z$draws), c(32000L, 3L)), all(is.finite(z$draws)),
      all(z$draws >= 0 & z$draws <= 1))
    source_sha[[basename(f)]] <<- sha(f)
    z$draws
  }))


  largest <- apply(draws, 1L, max)
  tied <- rowSums(draws == largest) > 1L
  winner <- max.col(draws, ties.method = 'first')
  for (j in 1:3) {
    ci <- interval(draws[, j])
    others <- setdiff(1:3, j)
    strict <- draws[, j] > pmax(draws[, others[1]], draws[, others[2]])
    cure[[length(cure) + 1L]] <- data.frame(model_id = id, concentration = candidates[j],
      posterior_mean = mean(draws[, j]), posterior_median = median(draws[, j]),
      credible_lower = ci[1], credible_upper = ci[2],
      probability_best_with_rule_ties = mean(winner == j),
      probability_strictly_best = mean(strict),
      probability_tied_for_best = mean(draws[, j] == largest & tied), stringsAsFactors = FALSE)
  }
  for (j in c(1L, 3L)) {
    delta <- draws[, 2] - draws[, j]; ci <- interval(delta)
    differences[[length(differences) + 1L]] <- data.frame(model_id = id,
      target_concentration = 1, comparator_concentration = candidates[j],
      posterior_mean_difference = mean(delta), posterior_median_difference = median(delta),
      credible_lower = ci[1], credible_upper = ci[2],
      probability_difference_positive = mean(delta > 0), probability_difference_zero = mean(delta == 0),
      stringsAsFactors = FALSE)
  }
  scenarios[[id]] <- list(retained_draws = nrow(draws), posterior_probability_any_best_tie = mean(tied),
    training_n = length(input$models[[id]]$y), LOCF = id == 'efficacy_1110')
}
cure <- do.call(rbind, cure); differences <- do.call(rbind, differences)
primary_cure <- cure[cure$model_id == 'efficacy_primary', ]
primary_differences <- differences[differences$model_id == 'efficacy_primary', ]
completion_cure <- cure[cure$model_id != 'efficacy_primary' & cure$concentration == 1, ]
completion_differences <- differences[differences$model_id != 'efficacy_primary', ]
summary <- list(primary = list(cure = primary_cure, differences = primary_differences),
  completion_sensitivity = list(
    probability_1percent_best_range = range(completion_cure$probability_best_with_rule_ties),
    probability_1percent_strictly_best_range = range(completion_cure$probability_strictly_best),
    cure_1percent_credible_lower_range = range(completion_cure$credible_lower),
    cure_1percent_credible_upper_range = range(completion_cure$credible_upper),
    difference_credible_interval_contains_zero_count = sum(completion_differences$credible_lower <= 0 &
      completion_differences$credible_upper >= 0), total_paired_comparisons = nrow(completion_differences)))
report <- list(status = 'COMPLETE', analysis_type = 'Post hoc descriptive posterior context',
  added_after_frozen_primary_analysis = TRUE,
  scope = 'Posterior uncertainty about model-based age-standardized cure probabilities, conditional on the specified model and data scenario',
  distinct_from_numerical_error = TRUE, in_prespecified_57_quantity_MC_family = FALSE,
  credible_interval = 'Equal-tailed95% marginal posterior intervals, empirical quantile type7; not simultaneous MC error intervals',
  ranking = 'Highest latent age-standardized cure probability among0.5%,1%,3%; screening is not part of this posterior ranking',
  ties = 'probability_best_with_rule_ties assigns ties to ascending concentration; strictly-best and tied-best probabilities are supplied separately',
  standardization_population = 'Common empirical age distribution of125 treated participants; primary cure model uses107 observed mITT outcomes',
  retained_draws_per_model = 128000L, efficacy_scenarios = 17L,
  plan_sha256 = plan_sha, input_sha256 = sha(input_file), source_chain_sha256 = source_sha,
  summary = summary, cure = cure, paired_differences = differences, scenarios = scenarios)
jsonlite::write_json(report, file.path(out, 'POSTERIOR_CONTEXT.json'), auto_unbox = TRUE,
  pretty = TRUE, digits = 17, na = 'null', null = 'null')
write.csv(cure, file.path(out, 'CURE_CREDIBLE_INTERVALS.csv'), row.names = FALSE)
write.csv(differences, file.path(out, 'PAIRED_DIFFERENCES.csv'), row.names = FALSE)
print(primary_cure, row.names = FALSE); print(primary_differences, row.names = FALSE)
cat('1% best-probability range over16 completions:', range(completion_cure$probability_best_with_rule_ties), '\n')
cat('Descriptive posterior context complete; no new fitting and no expanded MC coverage claim.\n')
