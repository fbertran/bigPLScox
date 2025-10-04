% Generated documentation for bigscale utility
\name{bigscale}
\alias{bigscale}
\title{Construct Scaled Design Matrices for Big Survival Models}
\description{
Prepares a large-scale feature matrix for stochastic gradient descent by
applying optional normalisation, stratified sampling, and batching rules.
}
\usage{
bigscale(
  formula,
  data,
  norm.method,
  strata.size,
  batch.size,
  features.mean,
  features.sd,
  parallel.flag,
  num.cores,
  bigmemory.flag,
  num.rows.chunk,
  col.names,
  type
)
}
\arguments{
  \item{formula}{Model formula used to extract the outcome and predictors that
  should be included in the scaled design matrix.}
  \item{data}{Input data source containing the variables referenced in
  \code{formula}.}
  \item{norm.method}{Normalisation strategy (for example centring or
  standardising columns) applied to the feature matrix.}
  \item{strata.size}{Number of observations to retain from each stratum when
  constructing stratified batches.}
  \item{batch.size}{Total size of each mini-batch produced by the scaling
  routine.}
  \item{features.mean}{Optional vector of column means that can be reused to
  normalise multiple data sets in a consistent manner.}
  \item{features.sd}{Optional vector of column standard deviations that pairs
  with \code{features.mean} during scaling.}
  \item{parallel.flag}{Logical flag signalling whether the scaling work should
  be parallelised across cores.}
  \item{num.cores}{Number of processor cores allocated when
  \code{parallel.flag} is \code{TRUE}.}
  \item{bigmemory.flag}{Logical flag specifying whether intermediate results
  should be stored in \pkg{bigmemory}-backed matrices.}
  \item{num.rows.chunk}{Chunk size used when streaming data from on-disk
  objects into memory.}
  \item{col.names}{Optional character vector assigning column names to the
  generated design matrix.}
  \item{type}{Type of model or preprocessing target being prepared, such as
  survival or regression.}
}
\value{
A scaled design matrix along with metadata describing the transformation that
was applied.
}
\seealso{
\code{\link{bigSurvSGD.na.omit}} for fitting models that use the scaled
features.
}
\examples{
\dontrun{
scaled <- bigscale(Surv(time, status) ~ ., data = my_training_data,
  norm.method = "standardize", batch.size = 128)
}
}
