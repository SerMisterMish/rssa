\name{tensor}
\alias{n.unfold}
\alias{n.fold}
\alias{n.mode.prod}
\alias{n.mode.prod.list}
\alias{n.mode.fnorm}
\alias{hosvd}
\alias{hooi}

\title{General tensor operations and decompositions}

\description{
  A set of routines to perform general tensor operations (\eqn{n}-mode
  unfolding, folding, \eqn{n}-mode tensor-matrix multiplication, and
  Frobenius norms) as well as higher-order tensor decompositions (HOSVD and
  HOOI).
}

\usage{
n.unfold(tens, n)
n.fold(mat, n, modes)
n.mode.prod(tens, mat, mode)
n.mode.prod.list(tens, mats, modes)
n.mode.fnorm(tens, modes = seq_along(dim(tens)))
hosvd(tens,
      unfolds,
      ranks = NULL,
      svd.method = c("svd", "eigen", "primme"),
      ...,
      report.progress = FALSE,
      calc.est = FALSE)
hooi(tens,
     unfolds,
     ranks,
     svd.method = c("svd", "eigen", "primme"),
     max_iter = 25,
     tol = 1e-05,
     ...,
     report.progress = FALSE)
}

\arguments{
  \item{tens}{a multidimensional array (tensor). For \code{hosvd} and \code{hooi},
    must be a trajectory tensor inheriting from \code{"tensor"} (as produced by
    \code{\link[Rssa:new.trajtens]{new.trajtens}}).}
  \item{mat}{a two-dimensional matrix.}
  \item{mats}{a list of two-dimensional matrices.}
  \item{n, mode}{an integer specifying the target mode (dimension index,
    1-based).}
  \item{modes}{for \code{n.fold}, an integer vector giving the dimensions of the
    resulting tensor. For \code{n.mode.prod.list}, an integer vector of modes
    corresponding to each matrix in \code{mats}. For \code{n.mode.fnorm}, an
    integer vector of modes along which the Frobenius norm is computed.}
  \item{unfolds}{an object of class \code{"unfoldings"} containing precomputed
    FFT-based mode unfoldings (as produced by
    \code{\link[Rssa:new.unfolds]{new.unfolds}}).}
  \item{ranks}{an integer vector (or a single integer) specifying the target
    ranks for each mode. For \code{hosvd}, if \code{NULL}, defaults to full
    ranks (\code{dim(tens)}).}
  \item{svd.method}{SVD implementation to use: standard SVD (\code{"svd"}),
    eigendecomposition of covariance (\code{"eigen"}), or PRIMME eigensolver
    (\code{"primme"}).}
  \item{\dots}{further arguments passed to the underlying SVD solver.}
  \item{report.progress}{logical, whether to display a text progress bar.}
  \item{calc.est}{logical, whether to calculate the reconstructed tensor estimate
    and its residual Frobenius norm in \code{hosvd}.}
  \item{max_iter}{maximum number of alternating least squares iterations in
    \code{hooi}.}
  \item{tol}{stopping tolerance for relative change in residual Frobenius norm
    in \code{hooi}.}
}

\details{
  \bold{General Tensor Operations:}

  \describe{
    \item{\code{n.unfold(tens, n)}}{
      Unfolds (matricizes) an \eqn{N}-mode tensor \code{tens} along mode
      \code{n}. If the tensor has dimensions \eqn{I_1 \times \dots \times I_N}{I_1 x ... x I_N},
      the resulting matrix has dimensions \eqn{I_n \times \prod_{k \neq n} I_k}{I_n x prod_{k != n} I_k}.
    }
    \item{\code{n.fold(mat, n, modes)}}{
      Folds a matrix \code{mat} back into a tensor of dimensions \code{modes}
      along mode \code{n}. This is the exact inverse operation of
      \code{n.unfold}.
    }
    \item{\code{n.mode.prod(tens, mat, mode)}}{
      Computes the \eqn{n}-mode product of tensor \code{tens} with matrix
      \code{mat} (\eqn{\mathcal{A} \times_n \mathbf{M}}{X \times_n M}). The number of
      columns of \code{mat} must equal the size of mode \code{mode} of
      \code{tens}. The resulting tensor has the same dimensions as
      \code{tens}, except that mode \code{mode} has size \code{nrow(mat)}.
    }
    \item{\code{n.mode.prod.list(tens, mats, modes)}}{
      Computes sequential mode products of tensor \code{tens} with a list of
      matrices \code{mats} along modes \code{modes}. If multiple matrices are
      specified for the same mode, they are automatically pre-multiplied
      together for efficiency.
    }
    \item{\code{n.mode.fnorm(tens, modes)}}{
      Calculates the slice-wise Frobenius norm of tensor \code{tens} along each
      mode in \code{modes}.
    }
  }

  \bold{Tensor Decompositions (HOSVD and HOOI):}

  \code{hosvd} implements the Higher-Order Singular Value Decomposition, 
  and \code{hooi} implements the Higher-Order Orthogonal Iteration algorithm 
  (best low-rank Tucker approximation).

  \emph{Important note on inputs:} As of now, the \code{tens} and
  \code{unfolds} arguments in \code{hosvd} and \code{hooi} are expected to have
  the specific structure created from a time series by the
  \code{\link[Rssa:new.trajtens]{new.trajtens}} and
  \code{\link[Rssa:new.unfolds]{new.unfolds}} functions.
  For general-purpose tensor decompositions see, for example,
  \code{\link[rTensor::rTensor]{rTensor}} package.
}

\value{
  For \code{n.unfold}, a matrix.

  For \code{n.fold}, \code{n.mode.prod}, and \code{n.mode.prod.list}, an array
  (tensor).

  For \code{n.mode.fnorm}, a list of numeric vectors containing the Frobenius
  norms along the specified modes.

  For \code{hosvd}, a list with components:
  \item{Z}{core tensor of the decomposition.}
  \item{U}{list of factor matrices (orthonormal mode singular vectors) for each mode.}
  \item{est}{reconstructed tensor estimate (if \code{calc.est = TRUE}).}
  \item{fnorm_resid}{Frobenius norm of the residual (if \code{calc.est = TRUE}).}

  For \code{hooi}, a list with components:
  \item{Z}{core tensor.}
  \item{U}{list of orthogonal factor matrices for each mode.}
  \item{conv}{logical, whether the algorithm converged within \code{max_iter}.}
  \item{iter}{number of iterations performed.}
  \item{est}{reconstructed tensor estimate.}
  \item{norm_percent}{percentage of original tensor Frobenius norm preserved in the estimate.}
  \item{fnorm_resid}{Frobenius norm of the final residual.}
  \item{all_resids}{vector of residual Frobenius norms across iterations.}
}

\references{
  De Lathauwer, L., De Moor, B., and Vandewalle, J. (2000)
  \emph{A multilinear singular value decomposition.}
  SIAM Journal on Matrix Analysis and Applications, 21(4), 1253-1278.

  De Lathauwer, L., De Moor, B., and Vandewalle, J. (2000)
  \emph{On the best rank-1 and rank-(R1, R2, ..., RN) approximation of higher-order tensors.}
  SIAM Journal on Matrix Analysis and Applications, 21(4), 1324-1342.

  Kolda, T. G., and Bader, B. W. (2009)
  \emph{Tensor decompositions and applications.}
  SIAM Review, 51(3), 455-500.
}

\seealso{
  \code{\link[Rssa:new.trajtens]{new.trajtens}},
  \code{\link[Rssa:new.unfolds]{new.unfolds}},
  \code{\link[Rssa:decompose.tssa]{decompose}},
  \code{\link[Rssa:ssa]{ssa}}.
}

\examples{
# --- General tensor operations ---
# Create a 3-mode tensor
A <- array(1:24, dim = c(2, 3, 4))

# Mode-1 unfolding and folding back
A_unf1 <- n.unfold(A, 1)
dim(A_unf1) # 2 x 12
identical(n.fold(A_unf1, 1, c(2, 3, 4)), A) # TRUE

# Mode-2 tensor-matrix product
M <- matrix(1:6, nrow = 2, ncol = 3)
A_prod <- n.mode.prod(A, M, mode = 2)
dim(A_prod) # 2 x 2 x 4

# Sequential mode products
M1 <- matrix(rnorm(4), 2, 2)
M2 <- matrix(rnorm(6), 2, 3)
A_mult <- n.mode.prod.list(A, list(M1, M2), modes = c(1, 2))

# Frobenius norms along each mode
norms <- n.mode.fnorm(A)
length(norms) # 3

# --- Tensor decompositions on trajectory tensors ---
\donttest{
F <- co2
tens <- new.trajtens(F, L = c(234, 119, 117))
unfolds <- new.unfolds(F, L = c(234, 119, 117))

# HOSVD decomposition
h_dec <- Rssa:::hosvd(tens, unfolds, ranks = c(10, 10, 10), calc.est = TRUE)
print(dim(h_dec$Z))

# HOOI decomposition
h_hooi <- Rssa:::hooi(tens, unfolds, ranks = c(10, 10, 10), max_iter = 10)
print(h_hooi$conv)
}
}
