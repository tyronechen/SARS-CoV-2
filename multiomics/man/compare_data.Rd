% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/multiomics_pipeline.R
\name{compare_data}
\alias{compare_data}
\title{Plot correlation between unimputed and imputed data. This version takes in files as input.}
\usage{
compare_data(
  data_one,
  data_two,
  data_names,
  classes,
  pcomp,
  outfile_path = "compare.pdf"
)
}
\arguments{
\item{data_one}{file path to table 1}

\item{data_two}{file path to table 2}

\item{data_names}{single label for a name}

\item{classes}{file path to classes}

\item{pcomp}{number of principal components to plot}

\item{outfile_path}{path to output file}
}
\description{
Quality control plot to check if data has mutated significantly. This version takes in files as input.
}
\seealso{
\code{\link[=impute_missing]{impute_missing()}}, \code{\link[=replace_missing]{replace_missing()}}, \code{\link[=cor_imputed_unimputed]{cor_imputed_unimputed()}}
}
