% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/litterDiallel.R
\name{diallelMatrixMakeAndRotate_av}
\alias{diallelMatrixMakeAndRotate_av}
\title{diallelMatrixMakeAndRotate_av}
\usage{
diallelMatrixMakeAndRotate_av(
  data,
  dam.col.name,
  sire.col.name,
  batch.col.name = NULL,
  batch.1.col.name = NULL,
  strains = c("AJ", "B6", "129", "NOD", "NZO", "CAST", "PWK", "WSB"),
  ...
)
}
\arguments{
\item{data}{data frame}

\item{dam.col.name}{dam column name}

\item{sire.col.name}{sire column name}

\item{batch.col.name}{name of batch/random effect column}

\item{batch.1.col.name}{name of additional batch/random effect column}

\item{strains}{character vector with names of inbred parental strains}

\item{...}{additional arguments}
}
\value{
returns diallel incidence matrices, rotated
}
\description{
Make design matrices for diallel, rotate to n-1 space. (a, v) model
}
\examples{
## not run
}
