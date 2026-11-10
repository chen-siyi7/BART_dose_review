\name{interval_actions}
\alias{interval_actions}
\title{Dose Recommendations Supported by Numerical Intervals}
\description{Find all actions consistent with the supplied screening and efficacy intervals.}
\usage{
interval_actions(stage, gates, scores, doses = seq_len(length(gates)),
  cutoff = 0.30)
}
\arguments{
\item{stage}{\code{"Final"} for final dose selection; any other stage uses the interim screening rule.}
\item{gates}{A list with one matrix per dose. Each matrix has lower and upper limits in its two columns, with toxicity screening in the first row and efficacy screening in the second.}
\item{scores}{A matrix with one row per dose and lower and upper efficacy-mean limits in its two columns.}
\item{doses}{Dose labels in the order of the supplied matrices.}
\item{cutoff}{Posterior screening probability required for eligibility.}
}
\details{
Final selection uses toxicity screening and ranks eligible doses by their
efficacy means. The earlier dose in the supplied order wins a tie.
Interim actions describe the possible sets passing both screens.
Intervals are treated as a rectangular set of possible numerical values.
}
\value{A sorted character vector. Final actions are labelled \code{"dose:j"}, with \code{"dose:0"} for no eligible dose. Interim sets are labelled \code{"mask:"} followed by eligible dose labels separated by colons.}
\examples{
gates <- rep(list(rbind(c(0.7, 0.9), c(0.1, 0.2))), 3)
scores <- rbind(c(0.1, 0.2), c(0.7, 0.8), c(0.4, 0.5))
interval_actions("Final", gates, scores)
}
