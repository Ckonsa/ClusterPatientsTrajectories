% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{compareDrugTrajectoriesWithDTW}
\alias{compareDrugTrajectoriesWithDTW}
\title{Compare patients drug trajectories with DTW algorithm.}
\usage{
compareDrugTrajectoriesWithDTW(patientsData, punishments)
}
\arguments{
\item{patientsData:}{list with patient IDs and their drug trajectories}

\item{punishments:}{numeric vector}
}
\value{
numeric matrix with distances between patients drug trajectories
}
\description{
Compare patients drug trajectories with DTW algorithm.
}
