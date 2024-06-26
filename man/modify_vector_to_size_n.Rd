% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{modify_vector_to_size_n}
\alias{modify_vector_to_size_n}
\title{Modify the vector to length of n. Remove last elements or add filler elements.}
\usage{
modify_vector_to_size_n(vector, n)
}
\arguments{
\item{vector:}{vector}

\item{n:}{integer, n > 0}
}
\value{
vector with length of n
}
\description{
Modify the vector to length of n. Remove last elements or add filler elements.
}
\examples{
modify_vector_to_size_n(c(1, 2, 3), 4) => c(1, 2, 3, -1)
modify_vector_to_size_n(c(1, 2, 3), 2) => c(1, 2)
}
