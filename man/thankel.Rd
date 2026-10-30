\name{htensor}
\alias{new.trajtens}
\alias{new.unfolds}

\title{Trajectory tensor operations.}

\description{
  A set of routines to create trajectory tensors and its unfoldings stored in
  compact FFT-based form.
}

\usage{
new.trajtens(F, L)
new.unfolds(F, L, rm.repeated = TRUE)
}

\arguments{
  \item{F}{series to construct the trajectory tensor and its unfolds for.}
  \item{L}{the window length, resulting tensor will have `length(L) + 1` modes.}
  \item{rm.repeated}{logical, whether to remove repeating columns in the 
  unfoldings. HOSVD on unfoldings with removed such columns results in HUTD.}
}

\details{
  Fast Fourier Transform provides a very efficient matrix-vector
  multiplication routine for Hankel with Hankel blocks matrices. See the
  paper in 'References' for the details of the algorithm.
}

\references{
  Korobeynikov, A. (2010) \emph{Computation- and space-efficient implementation of
    SSA.} Statistics and Its Interface, Vol. 3, No. 3, Pp. 257-268
}
