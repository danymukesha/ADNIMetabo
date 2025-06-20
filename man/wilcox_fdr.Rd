% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/statistical_tests.R
\name{wilcox_fdr}
\alias{wilcox_fdr}
\title{wilcox_fdr}
\usage{
wilcox_fdr(start, end, group_1, group_2, log, base, comparison)
}
\arguments{
\item{start}{the start of your matrix or data frame}

\item{end}{the end of the matrix or data frame}

\item{group_1}{data frame or matrix for  the first group of your comparison
for the Wilcoxon test}

\item{group_2}{data frame or matrix for  the second group of your comparison
for the Wilcoxon test}

\item{log}{logical TRUE if your data are log transformed}

\item{base}{base for your log transformation}

\item{comparison}{name of your comparison that it will the title
of your output}
}
\value{
the function give a csv file with a mean, the standard deviation,
the p-value of a Wilcoxon test and the foldchange for all variable
in your matrix or data frame
}
\description{
this \code{f()} calculates a Wilcoxon test and the foldchange
}
