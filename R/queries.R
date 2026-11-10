prepare_score_index <- function(doses, gates, q) {
    gates <- as.matrix(gates)
    q <- gmp::as.bigq(q)
    J <- length(doses)
    stopifnot(J > 0L, length(q) == J, nrow(gates) == J, ncol(gates) > 0L, is.numeric(doses), all(is.finite(doses)), 
        !anyDuplicated(as.character(doses)), all(doses != 0), all(is.finite(gates)), all(gates >= 0), all(gates <= 
            1000), all(gates == floor(gates)), all(q >= 0), all(q <= 1))
    z <- membership_costs(gates)
    eligible <- which(z$keep)
    eligible <- eligible[order(q[eligible], -eligible)]
    positions <- new.env(hash = TRUE, parent = emptyenv())
    for (r in seq_len(J)) assign(as.character(doses[r]), r, positions)
    structure(list(doses = doses, q = q, keep = z$keep, rho = z$rho, eligible = eligible, sorted_scores = q[eligible], 
        prefix = c(0, cumsum(z$rho[eligible])), positions = positions), class = "dose_score_index")
}

indexed_cost <- function(index, target, eta = 0) {
    stopifnot(inherits(index, "dose_score_index"), length(target) == 1L)
    eta <- gmp::as.bigq(eta)
    stopifnot(length(eta) == 1L, eta >= 0, eta <= 1)
    if (target == 0) 
        return(index$prefix[length(index$prefix)])
    key <- as.character(target)
    stopifnot(exists(key, index$positions, inherits = FALSE))
    r <- get(key, index$positions, inherits = FALSE)
    high <- min(gmp::as.bigq(1), index$q[r] + eta)
    n <- length(index$eligible)
    left <- 0L
    right <- n
    while (left < right) {
        mid <- left + (right - left + 1L)%/%2L
        low <- max(gmp::as.bigq(0), index$sorted_scores[mid] - eta)
        blocks <- low > high || (low == high && index$eligible[mid] < r)
        if (blocks) 
            right <- mid - 1L
        else left <- mid
    }
    (if (index$keep[r]) 
        0
    else index$rho[r]) + index$prefix[n + 1L] - index$prefix[left + 1L]
}

curve_cost <- function(curves, target, eta = 0) {
    stopifnot(inherits(curves, "dose_cost_curves"), length(target) == 1L)
    eta <- gmp::as.bigq(eta)
    stopifnot(length(eta) == 1L, eta >= 0, eta <= 1)
    key <- as.character(target)
    lookup <- attr(curves, "lookup")
    stopifnot(exists(key, lookup, inherits = FALSE))
    rows <- get(key, lookup, inherits = FALSE)
    events <- rows$eta
    left <- 1L
    right <- length(events)
    while (left < right) {
        mid <- left + (right - left + 1L)%/%2L
        if (events[mid] <= eta) 
            left <- mid
        else right <- mid - 1L
    }
    if (events[left] == eta) 
        rows$at_cost[left]
    else rows$right_cost[left]
}
