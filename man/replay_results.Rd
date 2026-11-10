\name{replay_results}
\alias{replay_results}
\alias{run_analysis}
\title{Run the Accompanying Study Workflows}
\description{Verify a matching full reproducibility bundle and run saved-result replay, a fresh study or figure generation.}
\usage{
replay_results(bundle, out)
run_analysis(analysis, bundle, out, ...)
}
\arguments{
\item{bundle}{Directory of the matching full reproducibility bundle, containing \code{MANIFEST.json}, \code{replay_all.R} and \code{data/}.}
\item{out}{A new output path whose parent directory exists, outside the bundle and installed package. For \code{"review_intervals"}, supply a PNG filename; other workflows take a directory.}
\item{analysis}{One of \code{"replay"}, \code{"benchmark"}, \code{"feasibility"}, \code{"reviewer_comparisons"}, \code{"intercept"}, \code{"history_compression"}, \code{"catalog_free"}, \code{"figures"}, \code{"review_intervals"} or \code{"clinical_replay"}.}
\item{...}{Additional character arguments passed to the selected script. For \code{"benchmark"}, use \code{"short"}, \code{"workload"} or \code{"all"}.}
}
\details{
The workflows need the separately supplied data and run records.
Source hashes and the complete bundle manifest are verified before execution.
The scripts run in their original bundle layout so their relative input
paths and frozen plans are preserved.
Saved-result replay reconstructs numerical results without repeating model
fitting or timing. Fresh study workflows perform the sampling and fitting
defined by their scripts; fresh native fits require BART 2.9.9.
}
\value{\code{replay_results} invisibly returns the replay result list after checking its status. \code{run_analysis} invisibly returns the output path when its subprocess succeeds.}
\examples{
\dontrun{
bundle <- "/path/to/reproducibility"
replay_results(bundle, "/path/to/new/replay_output")
run_analysis("benchmark", bundle, "/path/to/new/benchmark", "short")
}
}
