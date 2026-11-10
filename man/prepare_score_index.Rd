\name{prepare_score_index}
\alias{prepare_score_index}
\alias{indexed_cost}
\alias{frontier_curve}
\alias{curve_cost}
\title{Store and Query Costs Across Efficacy Tolerances}
\description{Construct sparse final-selection cost curves or a smaller score index for repeated questions about alternative recommendations.}
\usage{
prepare_score_index(doses, gates, q)
indexed_cost(index, target, eta = 0)
frontier_curve(doses, gates, q)
curve_cost(curves, target, eta = 0)
}
\arguments{
\item{doses}{Distinct nonzero numeric dose labels in preferred tie order.}
\item{gates}{An integer-valued matrix of screening counts between zero and 1,000, with one row per dose. Eligibility requires every count to be at least 300. For the manuscript's final rule, supply the toxicity column.}
\item{q}{Efficacy means between zero and one, in dose order. Numeric values or a \code{gmp::bigq} vector are accepted.}
\item{index}{An object returned by \code{prepare_score_index}.}
\item{curves}{An object returned by \code{frontier_curve}. Retain the complete object, including its lookup attribute, for queries.}
\item{target}{A supplied dose label, or zero for no eligible dose.}
\item{eta}{Allowed absolute efficacy-score change, between zero and one; numeric or a scalar \code{gmp::bigq}.}
}
\details{
Scores and tolerances are represented as exact rationals. Numeric inputs retain
their binary floating-point values; use \code{gmp::as.bigq("1/10")} when an exact
decimal tolerance is intended. Ties favor the earlier candidate.

For \eqn{J} candidates with at most \eqn{L} screening conditions, index preparation
uses \eqn{O(JL + J\log J)} arithmetic operations and linear additional storage.
Hash lookup locates the target and binary search finds its blocking competitors;
an indexed query uses \eqn{O(\log J)} comparisons. Curves use
\eqn{O(JL + J\log J + J^2)} operations and \eqn{O(JL + J^2)} storage.
These counts treat rational arithmetic as one operation; bit costs depend on
the input precision. A curve query uses \eqn{O(\log J)} comparisons.

Each curve records its cost at an event and immediately to its right.
The last right-hand cost extends to tolerance one. The no-dose curve is constant.
Neither construction performs posterior fitting or calibration.
}
\value{
\code{prepare_score_index} returns a \code{dose_score_index} object with sorted
eligible scores, cumulative removal costs and target positions.
\code{frontier_curve} returns a \code{dose_cost_curves} data frame with columns
\code{target}, \code{eta}, \code{at_cost} and \code{right_cost}, and a lookup
attribute for repeated queries. Event tolerances are rational strings.
\code{indexed_cost} and \code{curve_cost} return a scalar minimum total count cost.
}
\examples{
counts <- matrix(c(301, 301, 299), ncol = 1)
index <- prepare_score_index(c(9, 2, 7), counts, c(0.5, 0.5, 0.2))
curves <- frontier_curve(c(9, 2, 7), counts, c(0.5, 0.5, 0.2))
indexed_cost(index, target = 2, eta = 0)
curve_cost(curves, target = 2, eta = gmp::as.bigq("1/100"))
}
