\name{sample_reference}
\alias{sample_reference}
\alias{probability_intervals}
\alias{prepare_tree_prior}
\alias{sample_reference_checkpoints}
\title{Independent Posterior Reference Sampling}
\description{Draw independent posterior samples for a probit BART endpoint and construct finite-sample numerical intervals.}
\usage{
sample_reference(proposals, chunk, trees, a, successes, counts,
  prediction_weights, offset, method)
probability_intervals(q, endpoint, delta)
prepare_tree_prior(X, cuts, counts, minimum = 5L)
sample_reference_checkpoints(checkpoints, chunk, trees, a, successes,
  counts, prediction_weights, offset, method, delta, endpoint)
}
\arguments{
\item{proposals}{Number of independent prior proposals.}
\item{chunk}{Number of proposals processed together.}
\item{trees}{Number of trees in the finite BART ensemble.}
\item{a}{A prepared reference law for \code{"direct"} or \code{"aggregated"}, or a prepared tree prior for \code{"tree_prior"}.}
\item{successes}{Integer endpoint success counts by predictor profile.}
\item{counts}{Integer training counts by predictor profile.}
\item{prediction_weights}{Nonnegative weights by profile, summing to one.}
\item{offset}{Fixed offset on the probit scale.}
\item{method}{\code{"direct"} or \code{"aggregated"}; the checkpoint sampler also accepts \code{"tree_prior"}.}
\item{q}{Accepted posterior probabilities or standardized probability averages, in \eqn{[0,1]}.}
\item{endpoint}{\code{"T"} for toxicity or \code{"E"} for efficacy.}
\item{delta}{Scalar error allocation, between zero and one. Allocate this across all reported quantities when simultaneous coverage is required.}
\item{X}{A predictor matrix with at most 30 listed profiles.}
\item{cuts}{Finite, strictly increasing cutpoint vectors, one per predictor.}
\item{minimum}{Minimum terminal-node training count.}
\item{checkpoints}{Strictly increasing positive integer proposal counts at which to report intervals.}
}
\details{
The ensemble prior has total marginal standard deviation 1.5.
The direct route samples tree partitions individually; the aggregated route
combines leaf contributions conditional on partition counts. The tree-prior
route draws complete admissible trees without enumerating partitions.
Toxicity screening uses \eqn{q<0.45}; efficacy screening uses \eqn{q>0.30}.
}
\value{
\code{sample_reference} returns accepted draws \code{q}, the log likelihood
bound and a comparison acceptance count.
\code{probability_intervals} returns the accepted count, screening count, mean
and numerical bounds. With fewer than two draws the bounds are \eqn{[0,1]}.
\code{prepare_tree_prior} returns split and count arrays.
\code{sample_reference_checkpoints} returns accepted draws, their proposal
indices, checkpoint reports and tree-attempt diagnostics.
}
\examples{
law <- catalog_law(matrix(0, 1, 1), list(0.5), 12)
a <- prepare_reference(law, 1)
set.seed(27)
draws <- sample_reference(2000, 256, 2, a, 7, 12, 1, 0, "aggregated")
probability_intervals(draws$q, "E", 0.01)
}
