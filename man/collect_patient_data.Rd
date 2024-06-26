% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{collect_patient_data}
\alias{collect_patient_data}
\title{Collect data about patients.}
\usage{
collect_patient_data(data, ids, n)
}
\arguments{
\item{data:}{dataframe}

\item{ids:}{numeric vector}

\item{n:}{integer, n > 0 - length of drug trajectory}
}
\value{
list with vectors named "IDs" and "drugs". "IDs" hold information about patients IDs and "drugs" contains n-length drug trajectory per patient
}
\description{
Finds drug trajectory for every patient.
}
