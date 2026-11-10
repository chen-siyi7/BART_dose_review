\name{frontier_costs}
\alias{frontier_costs}
\alias{frontier_witness}
\title{Changes That Permit Alternative Doses}
\description{Compute minimum count changes permitting each action and construct a corresponding final-selection example.}
\usage{
frontier_costs(stage, doses, gates, q, eta = 0)
frontier_witness(doses, gates, q, target, eta = 0)
}
\arguments{
\item{stage}{\code{"Final"} for final selection; another value uses interim screening.}
\item{doses}{Numeric dose labels in selection order.}
\item{gates}{A matrix of screening counts, with one row per dose. Final selection uses the toxicity count column; interim screening uses toxicity and efficacy columns. Eligibility requires counts of at least 300.}
\item{q}{Efficacy means in dose order.}
\item{eta}{Maximum absolute change allowed in each efficacy mean, between zero and one.}
\item{target}{The final dose label to attain, or zero for no eligible dose.}
}
\details{Cost is the sum of absolute screening-count changes. Ranking ties favor the earlier dose in the supplied order.}
\value{\code{frontier_costs} returns a named vector of minimum costs. \code{frontier_witness} returns the target, cost, modified count matrix and modified efficacy scores.}
\examples{
counts <- matrix(c(301, 301), ncol = 1)
frontier_costs("Final", 1:2, counts, c(0.4, 0.3))
frontier_witness(1:2, counts, c(0.4, 0.3), target = 2)
}
