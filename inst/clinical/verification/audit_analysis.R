



for (pkg in c('readxl', 'jsonlite', 'digest', 'posterior'))
  stopifnot(requireNamespace(pkg, quietly = TRUE))
args <- commandArgs(trailingOnly = TRUE)
preflight <- '--preflight' %in% args
args <- args[args != '--preflight']
stopifnot(length(args) <= 2L)
script <- sub('^--file=', '', grep('^--file=', commandArgs(), value = TRUE))
app <- normalizePath(if (length(args)) args[1] else file.path(dirname(script), '..'))
out <- if (length(args) >= 2L) args[2] else file.path(app, 'verification')
dir.create(out, recursive = TRUE, showWarnings = FALSE)
outfile <- file.path(out, if (preflight) 'PREFIT_ANALYSIS_AUDIT.json' else 'ANALYSIS_AUDIT.json')
sha <- function(p) digest::digest(file = p, algo = 'sha256')
save_json <- function(x, p) jsonlite::write_json(x, p, auto_unbox = TRUE,
  pretty = TRUE, digits = 17, na = 'null', null = 'null')
checks <- list()
check <- function(ok, label) {
  checks[[label]] <<- isTRUE(ok)
  if (!isTRUE(ok)) stop(label, call. = FALSE)
}
equal <- function(x, y, tolerance = 1e-10) isTRUE(all.equal(x, y,
  tolerance = tolerance, check.attributes = FALSE))
read_json <- function(p) jsonlite::fromJSON(p, simplifyVector = TRUE)
expected_raw <- c('pone.0232394.s002.xlsx' = '11f77074b9008bc5fa67262aedfb128fb92dd2a1ed6adc18f17e97f0c9f26c4a',
  'pone.0232394.s003.xlsx' = '6ae7e90a613ddaabb1859509b70cc9062588d37cb8efd05876a0b035ce5842f6',
  'pone.0232394.s005.xlsx' = 'eda1f419014529668fe3526131a6867b2cc78f540a95b3010ae366b5d72e5370')
result <- tryCatch({
  planfile <- file.path(app, 'PLAN.json'); inputfile <- file.path(app, 'inputs', 'INPUTS.rds')
  plan <- read_json(planfile); input <- readRDS(inputfile)
  plan_sha <- sha(planfile); input_sha <- sha(inputfile)
  check(startsWith(plan$status, 'FROZEN') && identical(plan$sampling_started, FALSE), 'frozen_plan_metadata')
  check(identical(input_sha, plan$input_sha256), 'input_hash')
  source_names <- c('analysis_helpers.R', 'prepare_inputs.R', 'fit_worker.R', 'summarize.R')
  source_sha <- setNames(vapply(file.path(app, source_names), sha, ''), source_names)
  source_root <- normalizePath(file.path(app, '..', '..'))
  source_audit <- new.env(parent = baseenv())
  sys.source(file.path(source_root, 'R/audits.R'), envir = source_audit)
  check(all(vapply(source_names, function(n) source_audit$source_hash_matches(
    file.path(app, n), plan$source_hashes[[n]], source_root), TRUE)), 'frozen_source_hashes')
  check(equal(c(plan$chains, plan$burnin, plan$retained, plan$total_quantities, plan$alpha,
    plan$offset, plan$risk_limit, plan$screen_cut), c(4, 8000, 32000, 57, .05, 0, .45, .3)), 'prespecified_budgets_and_thresholds')
  check(equal(plan$scalar_error, .05 / 57) && identical(plan$BART_version, '2.9.9'), 'error_allocation_and_BART_version')
  raw_options <- c(file.path(app, 'raw'), file.path(app, '..', 'open_data_candidates_20261109', 'astodrimer', 'raw'))
  raw <- raw_options[vapply(raw_options, function(d) all(file.exists(file.path(d, names(expected_raw)))), TRUE)][1]
  check(!is.na(raw), 'original_workbooks_present')
  raw_sha <- setNames(vapply(file.path(raw, names(expected_raw)), sha, ''), names(expected_raw))
  check(identical(unname(raw_sha), unname(expected_raw)) &&
    identical(unname(raw_sha), unname(unlist(plan$raw_sha256[names(expected_raw)]))) &&
    identical(unname(raw_sha), unname(input$raw_sha256)), 'original_workbook_hashes')
  pop <- as.data.frame(readxl::read_excel(file.path(raw, names(expected_raw)[1]), sheet = 'Populations'))
  eff <- as.data.frame(readxl::read_excel(file.path(raw, names(expected_raw)[2]), sheet = 'Efficacy'))
  ae <- as.data.frame(readxl::read_excel(file.path(raw, names(expected_raw)[3]), sheet = 'AEs'))
  key <- function(z) vapply(strsplit(trimws(as.character(z)), '/', fixed = TRUE), paste, collapse = '', FUN.VALUE = '')
  pk <- key(pop[['Site/ Subject']]); ek <- key(eff[['Site/Subject']]); ak <- trimws(as.character(ae$Patient))
  use <- pop$Safety == 'Yes'; ei <- match(pk[use], ek)
  check(!anyDuplicated(pk) && !anyDuplicated(ek) && all(ak %in% pk[use]) &&
    setequal(pk[use], ek) && !anyNA(ei) && all(pop[['Treatment Group']][use] == eff$Treatment[ei]) &&
    all(pop$mITT[use] == eff$mITT[ei]), 'participant_links_and_populations')
  decode <- function(z) { check(all(is.na(z) | z %in% c('Yes', 'No')), 'binary_source_values');
    match(z, c('No', 'Yes')) - 1L }
  arms <- c('HEC placebo', '0.5% SPL7013', '1% SPL7013', '3% SPL7013')
  related <- unique(ak[ae$Relationship %in% c('POSSIBLY RELATED', 'PROBABLY RELATED', 'DEFINITELY RELATED')])
  dat <- data.frame(id = pk[use], concentration = c(0, .5, 1, 3)[match(pop[['Treatment Group']][use], arms)],
    placebo = as.integer(pop[['Treatment Group']][use] == arms[1]),
    age = as.numeric(vapply(strsplit(pop[['Age/ Race']][use], '/', fixed = TRUE), `[`, '', 1L)),
    related_ae = as.integer(pk[use] %in% related), mITT = pop$mITT[use] == 'Yes',
    cure_toc = decode(eff[['CC TOC']][ei]), cure_eot = decode(eff[['CC EOT']][ei]), stringsAsFactors = FALSE)
  check(equal(input$dat, dat), 'independently_rebuilt_participant_inputs')
  check(nrow(dat) == 125L && sum(dat$mITT) == 111L && sum(dat$related_ae) == 30L &&
    sum(dat$mITT & !is.na(dat$cure_toc)) == 107L, 'endpoint_denominators')
  missing_ids <- ek[eff$mITT == 'Yes' & is.na(eff[['CC TOC']])]
  mi <- match(missing_ids, dat$id)
  check(identical(missing_ids, c('10026', '10029', '14009', '14038')) &&
    identical(missing_ids, input$missing_ids) && identical(missing_ids, plan$missing_id_order) &&
    identical(dat$cure_eot[mi], c(1L, 1L, 1L, 0L)), 'missing_order_and_LOCF')
  bits <- as.matrix(expand.grid(rep(list(0:1), 4L), KEEP.OUT.ATTRS = FALSE))
  labels <- apply(bits, 1L, paste0, collapse = '')
  model_ids <- c('safety', 'efficacy_primary', paste0('efficacy_', labels))
  expected_models <- list(safety = list(index = seq_len(nrow(dat)), y = dat$related_ae),
    efficacy_primary = list(index = which(dat$mITT & !is.na(dat$cure_toc)),
      y = dat$cure_toc[dat$mITT & !is.na(dat$cure_toc)]))
  for (i in seq_len(16L)) { y <- dat$cure_toc; y[mi] <- bits[i, ]
    expected_models[[paste0('efficacy_', labels[i])]] <- list(index = which(dat$mITT), y = y[dat$mITT]) }
  check(identical(names(input$models), model_ids) && identical(plan$model_ids, model_ids), 'model_and_completion_order')
  check(all(vapply(model_ids, function(id) equal(input$models[[id]]$index, expected_models[[id]]$index) &&
    equal(input$models[[id]]$y, expected_models[[id]]$y), TRUE)), 'all_18_training_inputs')
  check(equal(as.matrix(input$completions[, 2:5]), bits) &&
    identical(as.character(input$completions$model_id), paste0('efficacy_', labels)) &&
    identical(as.logical(input$completions$LOCF), labels == '1110'), 'all_16_missing_completions')
  ages <- sort(unique(dat$age)); weights <- tabulate(match(dat$age, ages), nbins = length(ages)) / 125
  x <- as.matrix(dat[, c('concentration', 'placebo', 'age')])
  cuts <- lapply(seq_len(3), function(j) { v <- sort(unique(x[, j])); head(v, -1L) + diff(v) / 2 })
  names(cuts) <- colnames(x)
  prediction <- cbind(concentration = rep(c(.5, 1, 3), each = length(ages)), placebo = 0,
    age = rep(ages, 3))
  check(equal(input$x, x) && equal(input$cuts, cuts) && equal(plan$cuts, cuts), 'fixed_covariates_and_split_rules')
  check(equal(input$prediction_ages, ages) && equal(input$prediction_weights, weights) &&
    equal(input$prediction, prediction) && equal(plan$prediction_age_weights, weights) &&
    equal(plan$prediction_ages, ages) && equal(input$candidates, c(.5, 1, 3)), 'common_age_standardization_grid')
  data_audit <- read_json(file.path(app, 'inputs', 'DATA_AUDIT.json'))
  check(identical(data_audit$status, 'PASS') && identical(data_audit$input_sha256, input_sha), 'prepared_data_audit')
  base_report <- list(plan_sha256 = plan_sha, input_sha256 = input_sha,
    source_sha256 = as.list(source_sha), raw_sha256 = as.list(raw_sha),
    original_source_sha256 = plan$source_hashes[source_names],
    source_verification = 'Original hashes or recorded comment-free source transformations',
    primary_article = 'https://doi.org/10.1371/journal.pone.0232394',
    raw_source_urls = as.list(setNames(paste0('https://journals.plos.org/plosone/article/file?id=10.1371/journal.pone.0232394.s',
      c('002', '003', '005'), '&type=supplementary'), names(expected_raw))),
    models = 18L, efficacy_scenarios = 17L, safety_participants = 125L,
    primary_efficacy_observed = 107L, completion_efficacy_participants = 111L,
    scalar_error = .05 / 57, preflight = preflight)
  if (preflight) {
    c(list(status = 'PASS', scope = 'Pre-fit data/code/input checks only; saved-chain/results audit not yet performed',
      checks = checks), base_report)
  } else {


    masks <- as.matrix(expand.grid(rep(list(c(FALSE, TRUE)), 3L), KEEP.OUT.ATTRS = FALSE))
    possible_oracle <- function(gates, scores, cutoff) {
      feasible <- apply(masks, 1L, function(m) all(ifelse(m, gates[, 2] >= cutoff, gates[, 1] < cutoff)))
      decisions <- integer()
      for (i in which(feasible)) {
        eligible <- which(masks[i, ])
        if (!length(eligible)) decisions <- c(decisions, 0L) else for (k in seq_len(nrow(masks))) {
          values <- scores[cbind(seq_len(3), 1L + as.integer(masks[k, ]))]
          decisions <- c(decisions, eligible[which.max(values[eligible])])
        }
      }
      sort(unique(paste0('dose:', decisions)))
    }
    minimum_formula <- function(gates, scores, decision, cutoff) {
      r <- as.integer(sub('dose:', '', decision, fixed = TRUE))
      if (r == 0L) return(if (all(gates[, 2] < cutoff)) 3 else Inf)
      if (gates[r, 1] < cutoff) return(Inf)
      branch_cost <- vapply(c(FALSE, TRUE), function(retain) {
        lower <- if (retain) scores[r, 1] else 0
        other_cost <- vapply(setdiff(seq_len(3), r), function(j) {
          if (lower == 1 && r < j) return(0)
          if (gates[j, 2] < cutoff || lower > scores[j, 2] || (lower == scores[j, 2] && r < j)) return(1)
          Inf
        }, 0)
        1 + as.integer(retain) + sum(other_cost)
      }, 0)
      min(branch_cost)
    }
    check(identical(possible_oracle(matrix(c(.3, 1), 3, 2, byrow = TRUE), matrix(1, 3, 2), .3), 'dose:1') &&
      identical(possible_oracle(matrix(c(0, .29), 3, 2, byrow = TRUE), matrix(c(0, 1), 3, 2, byrow = TRUE), .3), 'dose:0') &&
      length(possible_oracle(matrix(c(0, 1), 3, 2, byrow = TRUE), matrix(c(0, 1), 3, 2, byrow = TRUE), .3)) == 4L,
      'independent_oracle_boundary_cases')
    chain_files <- unlist(lapply(model_ids, function(id) file.path(app, 'chains', paste0(id, '_chain', 1:4, '.rds'))))
    check(all(file.exists(chain_files)) && length(list.files(file.path(app, 'chains'), pattern = '_chain[1-4]\\.rds$')) == 72L,
      'all_72_saved_chains_present')
    chains <- setNames(vector('list', 18L), model_ids); chain_sha <- setNames(character(72L), basename(chain_files))
    seeds <- numeric(72L); elapsed <- numeric(72L); warnings <- list(); k <- 0L
    for (i in seq_along(model_ids)) {
      id <- model_ids[i]; chains[[id]] <- vector('list', 4L)
      for (role in 1:4) {
        k <- k + 1L; file <- chain_files[k]; z <- readRDS(file)
        expected_seed <- as.integer(plan$seed_base + 100 * i + role)
        RNGkind('Mersenne-Twister', 'Inversion', 'Rejection'); set.seed(expected_seed)
        check(identical(z$status, 'COMPLETE') && identical(z$model_id, id) && z$role == role &&
          z$seed == expected_seed && identical(z$initial_rng, .Random.seed), paste0('chain_metadata_', id, '_', role))
        check(is.integer(z$final_rng) && length(z$final_rng) == length(z$initial_rng) &&
          !identical(z$final_rng, z$initial_rng) && identical(z$plan_sha256, plan_sha) &&
          identical(z$input_sha256, input_sha) && identical(z$BART_version, plan$BART_version),
          paste0('chain_provenance_', id, '_', role))
        check(identical(dim(z$draws), c(32000L, 3L)) && all(is.finite(z$draws)) &&
          all(z$draws >= 0 & z$draws <= 1) && identical(colnames(z$draws), c('0.5', '1', '3')) &&
          z$training_n == length(expected_models[[id]]$y) && z$training_events == sum(expected_models[[id]]$y) &&
          is.finite(z$elapsed_seconds) && z$elapsed_seconds >= 0, paste0('chain_dimensions_', id, '_', role))
        chains[[id]][[role]] <- z$draws; seeds[k] <- z$seed; elapsed[k] <- z$elapsed_seconds
        chain_sha[k] <- sha(file); warnings[[basename(file)]] <- z$warnings
      }
    }
    check(length(unique(seeds)) == 72L, 'unique_deterministic_chain_seeds')
    series_matrix <- function(id, j) do.call(cbind, lapply(chains[[id]], function(draws) draws[, j]))
    diagnostic <- function(y, series, candidate, scenario, risk = NA_real_) {
      safe <- function(fun, z) tryCatch(suppressWarnings(fun(z)), error = function(e) NA_real_)
      center <- mean(y); rh <- safe(posterior::rhat, y); ess <- safe(posterior::ess_bulk, y)
      mc <- safe(posterior::mcse_mean, y)
      split_mc <- vapply(list(y[, 1:2, drop = FALSE], y[, 3:4, drop = FALSE]),
        function(z) safe(posterior::mcse_mean, z), 0)
      denom <- sqrt(sum(split_mc^2))
      pair <- if (is.finite(denom) && denom > 0) abs(mean(y[, 1:2]) - mean(y[, 3:4])) / denom else NA_real_
      constant <- length(unique(as.vector(y))) == 1L
      qualified <- !constant && all(is.finite(c(rh, ess, mc, pair))) && rh <= 1.05 && ess >= 100 && mc > 0 && pair <= 6
      width <- qnorm(1 - plan$alpha / (2 * 57)) * mc
      data.frame(series = series, candidate = candidate, scenario = scenario, risk_limit = risk,
        mean = center, rhat = rh, ess = ess, mcse = mc, paired_z = pair, qualified = qualified,
        lower = if (qualified) max(0, center - width) else 0, upper = if (qualified) min(1, center + width) else 1,
        constant_series = constant, stringsAsFactors = FALSE)
    }
    reconstructed <- list()
    for (j in 1:3) { y <- series_matrix('safety', j)
      reconstructed[[length(reconstructed) + 1L]] <- diagnostic(y, 'T_mean', j, 'shared_safety', .45)
      reconstructed[[length(reconstructed) + 1L]] <- diagnostic(as.numeric(y < .45) + 0 * y, 'T_gate', j, 'shared_safety', .45) }
    for (id in model_ids[-1]) for (j in 1:3)
      reconstructed[[length(reconstructed) + 1L]] <- diagnostic(series_matrix(id, j), 'E_mean', j, id)
    diagnostics <- do.call(rbind, reconstructed)
    saved_diagnostics <- read.csv(file.path(app, 'results', 'DIAGNOSTICS.csv'), stringsAsFactors = FALSE)
    check(nrow(diagnostics) == 57L && equal(saved_diagnostics, diagnostics), 'all_57_diagnostics_and_scalar_intervals')
    get_diag <- function(type, j, scenario) diagnostics[diagnostics$series == type & diagnostics$candidate == j & diagnostics$scenario == scenario, ]
    gates <- matrix(NA_real_, 3, 2); counts <- numeric(3)
    for (j in 1:3) {
      m <- get_diag('T_mean', j, 'shared_safety'); g <- get_diag('T_gate', j, 'shared_safety')
      lower_from_mean <- max(0, 1 - m$upper / .45)
      upper_from_mean <- min(1, (1 - m$lower) / .55)
      gates[j, ] <- c(max(g$lower, lower_from_mean), min(g$upper, upper_from_mean))
      counts[j] <- sum(series_matrix('safety', j) < .45)
    }
    check(all(gates[, 1] <= gates[, 2]) && all(gates >= 0 & gates <= 1), 'shared_safety_moment_intersections')
    reports <- jsonlite::fromJSON(file.path(app, 'results', 'INTERVAL_REPORTS.json'), simplifyVector = FALSE)
    summary <- read.csv(file.path(app, 'results', 'POSTERIOR_SUMMARY.csv'), stringsAsFactors = FALSE)
    costs <- read.csv(file.path(app, 'results', 'SCREENING_CHANGE_COSTS.csv'), stringsAsFactors = FALSE)
    totals <- read_json(file.path(app, 'results', 'RESULTS.json'))
    check(nrow(summary) == 51L && nrow(costs) == 68L && identical(names(reports), model_ids[-1]), 'complete_result_dimensions')
    atoms <- as.vector(rbind(paste0('Tgate:', 1:3), paste0('Emean:', 1:3)))
    subsets <- as.matrix(expand.grid(rep(list(c(FALSE, TRUE)), 6L), KEEP.OUT.ATTRS = FALSE))
    subset_audit <- list(); expected_summary <- list(); expected_costs <- list()
    M <- 4 * 32000; kappa <- ceiling(M * .3); fitted_eligible <- counts >= kappa
    flip_cost <- ifelse(fitted_eligible, counts - kappa + 1, kappa - counts)
    for (id in model_ids[-1]) {
      e <- lapply(1:3, function(j) get_diag('E_mean', j, id))
      means <- vapply(e, function(d) d$mean, 0)
      scores <- t(vapply(e, function(d) c(d$lower, d$upper), numeric(2)))
      possible <- possible_oracle(gates, scores, .3)
      point <- if (any(fitted_eligible)) which(fitted_eligible)[which.max(means[fitted_eligible])] else 0L
      check(paste0('dose:', point) %in% possible, paste0('point_contained_', id))
      report <- reports[[id]]
      check(identical(as.character(unlist(report$possible)), possible) && identical(report$model_id, id) &&
        equal(do.call(rbind, lapply(report$screening_intervals, unlist)), gates) &&
        equal(do.call(rbind, lapply(report$efficacy_intervals, unlist)), scores),
        paste0('reported_intervals_and_actions_', id))
      subset_sets <- lapply(seq_len(64L), function(i) {
        gg <- gates; ee <- scores
        for (j in 1:3) { if (!subsets[i, 2 * j - 1L]) gg[j, ] <- c(0, 1)
          if (!subsets[i, 2 * j]) ee[j, ] <- c(0, 1) }
        possible_oracle(gg, ee, .3)
      })
      single <- length(possible) == 1L
      good <- vapply(subset_sets, function(z) single && identical(z, possible), TRUE)
      selected <- as.character(unlist(report$atoms)); report_size <- length(selected)
      check(!anyDuplicated(selected) && all(selected %in% atoms) && report$cost == report_size,
        paste0('reported_atom_names_', id))
      if (single) {
        minimum <- min(rowSums(subsets)[good]); formula <- minimum_formula(gates, scores, possible, .3)
        selected_i <- which(apply(subsets, 1, function(z) identical(atoms[z], selected)))
        check(identical(report$status, 'CERTIFIED') && identical(report$minimum_established, TRUE) &&
          length(selected_i) == 1L && good[selected_i] && report_size == minimum && formula == minimum &&
          report$formula_cost == minimum && report$exhaustive_cost == minimum && report$subsets_evaluated == 64L,
          paste0('formula_and_all_64_subsets_', id))
        ex_atoms <- as.character(unlist(report$exhaustive_atoms))
        ex_i <- which(apply(subsets, 1, function(z) identical(atoms[z], ex_atoms)))
        check(length(ex_i) == 1L && good[ex_i] && length(ex_atoms) == minimum,
          paste0('reported_exhaustive_witness_', id))
      } else {
        minimum <- NA_integer_; formula <- NA_real_
        check(identical(report$status, 'UNRESOLVED') && identical(report$minimum_established, FALSE) &&
          identical(selected, atoms) && report_size == 6L && !any(lengths(subset_sets) == 1L),
          paste0('unresolved_full_report_', id))
      }
      subset_audit[[id]] <- list(possible_decisions = possible, certified = single,
        subsets_checked = 64L, minimum_cost = minimum, formula_cost = formula,
        chosen_atoms = selected, equally_minimal_subsets = if (single) sum(good & rowSums(subsets) == minimum) else 0L)


      state_cost <- apply(masks, 1, function(mask) sum(flip_cost[mask != fitted_eligible]))
      state_action <- apply(masks, 1, function(mask) { eligible <- which(mask)
        if (length(eligible)) eligible[which.max(means[eligible])] else 0L })
      brute_cost <- vapply(0:3, function(r) min(state_cost[state_action == r]), 0)
      expected_costs[[id]] <- data.frame(target = paste0('dose:', 0:3), eta = 0,
        screening_changes = brute_cost, total_draws = M, screen_count_cutoff = kappa, model_id = id,
        stringsAsFactors = FALSE)
      for (j in 1:3) {
        m <- get_diag('T_mean', j, 'shared_safety'); g <- get_diag('T_gate', j, 'shared_safety')
        expected_summary[[length(expected_summary) + 1L]] <- data.frame(model_id = id, candidate = j,
          concentration = c(.5, 1, 3)[j], AE_posterior_mean = m$mean, AE_lower = m$lower, AE_upper = m$upper,
          safety_probability = g$mean, safety_lower = gates[j, 1], safety_upper = gates[j, 2],
          cure_posterior_mean = means[j], cure_lower = scores[j, 1], cure_upper = scores[j, 2],
          fitted_recommendation = paste0('dose:', point), possible_decisions = paste(possible, collapse = ';'),
          established = single, report_length = report_size, LOCF = id == 'efficacy_1110', stringsAsFactors = FALSE)
      }
    }
    check(equal(summary, do.call(rbind, expected_summary)), 'all_51_posterior_and_decision_rows')
    check(equal(costs, do.call(rbind, expected_costs)), 'all_68_screening_change_costs')
    check(identical(totals$status, 'COMPLETE') && totals$models == 18L && totals$efficacy_scenarios == 17L &&
      totals$diagnostics_total == 57L && totals$diagnostics_qualified == sum(diagnostics$qualified) &&
      equal(c(totals$nominal_family_error, totals$scalar_error, totals$retained_draws_per_model,
        totals$risk_limit, totals$screen_cut), c(.05, .05 / 57, M, .45, .3)), 'results_summary_metadata')
    totals_reports <- jsonlite::fromJSON(file.path(app, 'results', 'RESULTS.json'), simplifyVector = FALSE)$reports
    check(identical(totals_reports, reports), 'duplicated_report_payloads_identical')
    c(list(status = 'PASS', scope = 'Independent replay from original workbooks and all saved averaged prediction chains; no native fitting or timing rerun',
      checks = checks, chains = 72L, diagnostics = 57L, diagnostics_qualified = sum(diagnostics$qualified),
      posterior_summary_rows = 51L, screening_cost_rows = 68L, atoms_per_scenario = 6L,
      subsets_per_scenario = 64L, subsets_checked = 17L * 64L,
      feasibility_oracle = 'All8 eligibility masks and all8 lower/upper efficacy corners per subset',
      chain_sha256 = as.list(chain_sha), chain_seeds = as.list(setNames(seeds, basename(chain_files))),
      chain_warning_messages = warnings, scenario_audit = subset_audit,
      limitations = c('The saved chain object contains3 age-standardized concentration means per iteration, not the81 individual prediction probabilities. The prediction grid and aggregation code are hash-verified; the81-to3 aggregation cannot be replayed from those reduced objects.',
        'The57-quantity family is nominal/asymptotic under the stated MCMC assumptions; it is not a finite-sample posterior or clinical validation guarantee.',
        'All16 binary completions address the four missing mITT TOC outcomes; they do not bound population or model uncertainty.',
        'Elapsed seconds are checked as stored nonnegative metadata and are not independently remeasured.')),
      base_report)
  }
}, error = function(e) list(status = 'FAIL', checks = checks, error = conditionMessage(e), preflight = preflight))
save_json(result, outfile)
cat(result$status, ': ', outfile, '\n', sep = '')
if (!identical(result$status, 'PASS')) stop(result$error, call. = FALSE)
