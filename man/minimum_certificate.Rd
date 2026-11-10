\name{minimum_certificate}
\alias{prepare_certificate}
\alias{certificate_actions}
\alias{minimum_certificate}
\title{Select a Minimum-Cost Supporting Report}
\description{Prepare interval reports, recover the actions they support and select a least-cost subset that determines a target action.}
\usage{
prepare_certificate(record)
certificate_actions(z, selected = z$atoms)
minimum_certificate(z,
  costs = setNames(rep(1, length(z$atoms)), z$atoms))
}
\arguments{
\item{record}{A list containing \code{stage}, \code{gates}, \code{scores} and \code{target}. Gates follow \code{interval_actions}; scores are a list of lower-upper pairs. The target is an action label such as \code{"dose:2"}.}
\item{z}{A report returned by \code{prepare_certificate}.}
\item{selected}{Names of interval assertions to retain. Omitted assertions are replaced by the full interval \code{c(0, 1)}.}
\item{costs}{A named vector of strictly positive, finite assertion costs. Unit costs minimize the number of retained assertions.}
}
\details{
Assertion names identify toxicity screening, efficacy screening or efficacy
means by dose. The selected subset must determine the target for every
numerical value in the remaining rectangular intervals.
}
\value{
\code{prepare_certificate} returns the prepared report with its assertions and bounds.
\code{certificate_actions} returns compatible action labels.
\code{minimum_certificate} returns \code{status}, \code{cost} and \code{atoms}.
Status is \code{"CERTIFIED"} when a supporting subset exists and
\code{"UNCERTIFIABLE"} otherwise.
}
\examples{
record <- list(stage = "Final", target = "dose:2",
  gates = rep(list(rbind(c(0.7, 0.9), c(0.1, 0.2))), 3),
  scores = list(c(0.1, 0.2), c(0.7, 0.8), c(0.4, 0.5)))
z <- prepare_certificate(record)
report <- minimum_certificate(z)
certificate_actions(z, report$atoms)
}
\seealso{\code{\link{interval_actions}}}
