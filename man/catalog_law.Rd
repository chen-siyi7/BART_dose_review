\name{catalog_law}
\alias{catalog_law}
\alias{catalog_law_limited}
\alias{prepare_reference}
\title{Prepare the Finite-Point Tree Partition Law}
\description{Enumerate the prior law on partitions induced by complete admissible trees and prepare it for reference sampling.}
\usage{
catalog_law(X, cuts, counts, minimum = 5L, exact = FALSE)
catalog_law_limited(X, cuts, counts, minimum = 5L, max_seconds = 30,
  max_states = 100000L, max_partitions = 10000L,
  max_products = 2000000L)
prepare_reference(law, H)
}
\arguments{
\item{X}{A finite predictor matrix, with one row per listed profile.}
\item{cuts}{A list of finite, strictly increasing cutpoints, one vector per predictor.}
\item{counts}{Nonnegative integer training counts for the listed rows. A zero count identifies a prediction-only profile.}
\item{minimum}{Minimum training count in each terminal node.}
\item{exact}{Use exact rational arithmetic for partition probabilities.}
\item{max_seconds}{Elapsed-time limit for the guarded recursion.}
\item{max_states}{Maximum number of visited node states.}
\item{max_partitions}{Maximum number of partition entries at a node.}
\item{max_products}{Maximum number of child-partition products.}
\item{law}{A completed law returned by \code{catalog_law}, or the \code{law} element of a completed guarded calculation.}
\item{H}{Number of rows in \code{X}.}
}
\details{
The prior split probability at depth \eqn{d} is \eqn{0.95/(d+1)^2}.
Conditioning is on complete admissible trees. Invalid splits are not
renormalized locally. The unguarded implementation accepts at most 12 rows;
the guarded implementation accepts at most 30.
}
\value{
\code{catalog_law} returns normalized partition probabilities, the prior
admissibility probability and the number of stored regions.
\code{catalog_law_limited} returns \code{"COMPLETE"} with a law or
\code{"RESOURCE_LIMIT"} with a reason. A stopped calculation returns no partial law.
\code{prepare_reference} returns the partition and leaf-incidence arrays used by the sampler.
}
\examples{
X <- matrix(c(0, 1), ncol = 1)
law <- catalog_law(X, list(0.5), c(5, 5), exact = TRUE)
reference <- prepare_reference(law, nrow(X))
reference$weights
}
