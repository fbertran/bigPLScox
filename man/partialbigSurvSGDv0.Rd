% Generated documentation for partialbigSurvSGDv0 helper
\name{partialbigSurvSGDv0}
\alias{partialbigSurvSGDv0}
\title{Incremental Survival Model Fitting with Pre-Scaled Data}
\description{
Loads a previously scaled design matrix and continues the stochastic gradient
optimisation for a subset of variables.
}
\usage{
partialbigSurvSGDv0(
  name.col,
  datapath,
  ncores,
  resBigscale,
  bigmemory.flag,
  inf.mth
)
}
\arguments{
  \item{name.col}{Character vector containing the column names that should be
  included in the partial fit.}
  \item{datapath}{File system path or connection where the big-memory backing
  file for the scaled design matrix is stored.}
  \item{ncores}{Number of processor cores allocated to the partial fitting
  procedure.}
  \item{resBigscale}{Result object returned by \code{\link{bigscale}}
  containing scaling statistics to be reused.}
  \item{bigmemory.flag}{Logical flag determining whether big-memory backed
  matrices are used when loading and updating the design matrix.}
  \item{inf.mth}{Inference method requested for the partial fit, such as
  none, asymptotic, or bootstrap summaries.}
}
\value{
An updated survival model object that reflects the partial fit using the
selected columns.
}
\seealso{
\code{\link{bigscale}} and \code{\link{bigSurvSGD.na.omit}}.
}
\examples{
\dontrun{
continued <- partialbigSurvSGDv0(
  name.col = c("age", "sex"),
  datapath = tempfile(),
  ncores = 2,
  resBigscale = scaled,
  bigmemory.flag = TRUE,
  inf.mth = "bootstrap"
)
}
}
