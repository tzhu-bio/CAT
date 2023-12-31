% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FootprintAnalysis.R
\name{compareFootprints}
\alias{compareFootprints}
\title{Title   A function to compare the Tn5 signal around the motif in the whole genome.}
\usage{
compareFootprints(
  samples,
  motif,
  smooth_window = 5,
  flanking_length = 200,
  N_cores = 1,
  color = NA,
  save_path = NA,
  file_prefix = NA
)
}
\arguments{
\item{file_prefix}{}
}
\description{
Title   A function to compare the Tn5 signal around the motif in the whole genome.
}
\examples{
  compareFootprints(samples = c("CD8pos_T","Monocytes"),motif = "Arntl_MA0603.1",bindetect_path = "F:/CAT/example/BINDetect/",signal_path = "F:/CAT/example/ATACcorrect/")

}
