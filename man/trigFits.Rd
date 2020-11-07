\name{trigFits}
\alias{trigFits}
\docType{data}
\title{
statistics pertinent to trigonometric regression analysis
of transcription levels of yeast genes
}
\description{
statistics pertinent to trigonometric regression analysis
of transcription levels of yeast genes
}
\usage{data(trigFits)}
\format{
  The format is: \cr
 num [1:6178, 1:3] 0.0922 0.0534 0.0978 0.0735 0.0834 ... \cr
 - attr(*, "dimnames")=List of 2 \cr
  ..$ : chr [1:6178] "YAL001C" "YAL002W" "YAL003W" "YAL004W" ... \cr
  ..$ : chr [1:3] "msep" "amp" "phase"
}
\details{
with expression trajectories from spYCCES in yeastCC package,
a first order trigonometric regression is fit assuming
a period of 64 minutes.  msep is the mean squared error of prediction
at sampling times, amp is estimated amplitude, phase is estimated
phase
}
\source{
yeastCC, also see vignette
}
%\references{
%%  ~~ possibly secondary sources and usages ~~
%}
\examples{
data(trigFits)
library(yeastCC)
data(spYCCES)
alp = spYCCES[,spYCCES$syncmeth=="alpha"]
plot(exprs(alp)["YAL040C",]~alp$time)
abline(h=trigFits["YAL040C", "amp"])
abline(h=-trigFits["YAL040C", "amp"])
}
\keyword{datasets}
