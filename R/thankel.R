#   R package for Singular Spectrum Analysis
#   Copyright (c) 2013 Anton Korobeynikov <anton at korobeynikov dot info>
#
#   This program is free software; you can redistribute it
#   and/or modify it under the terms of the GNU General Public
#   License as published by the Free Software Foundation;
#   either version 2 of the License, or (at your option)
#   any later version.
#
#   This program is distributed in the hope that it will be
#   useful, but WITHOUT ANY WARRANTY; without even the implied
#   warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
#   PURPOSE.  See the GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public
#   License along with this program; if not, write to the
#   Free Software Foundation, Inc., 675 Mass Ave, Cambridge,
#   MA 02139, USA.

.traj.dim.tssa <- function(x) {
  c(x$window, x$length - sum(x$window) + length(x$window))
}

.default.neig.tssa <- function(x, ...) pmin(50, .traj.dim(x))

new.trajtens <- function(F, L) {
  n <- length(L)
  stopifnot(is.numeric(L) && n > 2)
  N <- length(F)
  
  idx <- Reduce(\(l, k) outer(l, 1:L[k], `+`), 3:n, outer(1:L[1], 1:L[2], `+`)) - n + 1
  X <- F[idx]
  dim(X) <- dim(idx)
  class(X) <- "tensor"
  
  X
}

new.unfolds <- function(F, L, rm.repeated = TRUE) {
  N <- length(x)
  n <- length(L)
  stopifnot(is.numeric(L) && n > 2)
  Lsum <- sum(L[-n])
  
  unfolds_r <- list()
  if (!rm.repeated) {
    hn <- hankel(F, Lsum - 1)
    hlast <- hankel(F, sum(L[-(n - 1)]) - 1)
    
    unfolds_r <- lapply(L[-n], \(l) Rssa::new.hbhmat(Re(hn), c(l, 1)))
    unfolds_r[[n]] <- Rssa::new.hbhmat(Re(hlast), c(L[n], 1))
  } else {
    unfolds_r <- lapply(L, FUN = Rssa::new.hmat, F = Re(F))
  }
  
  if (is.complex(F)) {
    unfolds_i <- list()
    
    if (!rm.repeated) {
      unfolds_i <- lapply(L[-n], \(l) Rssa::new.hbhmat(Im(hn), c(l, 1)))
      unfolds_i[[n]] <- Rssa::new.hbhmat(Im(hlast), c(L[n], 1))
    } else {
      unfolds_i <- lapply(L, FUN = Rssa::new.hmat, F = Im(F))
    }
  } else unfolds_i <- NULL
  
  result <- list(Re = unfolds_r, Im = unfolds_i)
  class(result) <- "unfoldings"
  
  result
}

.get.or.create.trajtens <- function(x) {
  .get.or.create(x, "trajtens", new.trajtens(.F(x), L = .traj.dim(x)))
}

.get.or.create.unfolds <- function(x) {
  .get.or.create(x, "unfolds", new.unfolds(.F(x), L = .traj.dim(x),
                                           rm.repeated = x$tensor.decomposition %in% c("hutd", "hutd+hooi")))
}

decompose.tssa <- function(x,
                           neig = NULL,
                           ...,
                           force.continue = FALSE) {
  # Check, whether continuation of decomposition is requested
  if (!force.continue && any(nsigma(x)) > 0)
    stop("Continuation of decomposition is not supported for this method.")
  
  if (is.null(neig))
    neig <- .default.neig(x, ...)
  
  if (length(neig) == 1)
    neig <- pmin(neig, .traj.dim(x))
  
  decomp <- x$tensor.decomposition
  svd.method <- x$svd.method
  
  H <- .get.or.create.trajtens(x)
  unf <- .get.or.create.unfolds(x)
  
  if (identical(decomp, "hosvd") || identical(decomp, "hutd"))
    H.dec <- hosvd(H, unf, neig, svd.method, ...)
  else if (identical(decomp, "hutd+hooi") || identical(decomp, "hooi"))
    H.dec <- hooi(H, unf, neig, svd.method, ...)
  else 
    stop("Unsupported tensor decomposition method")
  
  .set.decomposition(x, Z = H.dec$Z, U = H.dec$U, sigma = n.mode.fnorm(H.dec$Z))
  
  x
}

.init.fragment.tssa <- function(this)
  expression({
    if (any(circular))
      stop("Circular variant of Tensor SSA isn't implemented yet")
    
    if (length(L) < 2)
      stop("Length of 'L' should be at least 2 for Tensor SSA")
    
    x <- as.vector(x)
    N <- length(x)
    
    wmask <- fmask <- weights <- NULL
    
    column.projector <- row.projector <- NULL
  })

.colspan.tssa <- function(x, idx, modes = NULL) {
  U <- .U(x)
  if (is.null(modes)) modes <- seq_along(U)
  if (length(modes) == 1) U[[modes]][, idx, drop = FALSE] else lapply(U[modes], \(u) u[, idx, drop = FALSE])
}

.V.tssa <- function(x) {
  warning("'.V' is the same as '.U' for Tensor SSA. Consider using '.U'")
  .U(x)
}

.rowspan.tssa <- function(x, idx, ...)
  stop("There are no such thing as 'rows' in Tensor SSA. Consider using '.colspan'")

.Z <- function(x, ...)
  UseMethod(".Z")

.Z.default <- function(x)
  x$Z

.Z.tssa <- function(x)
  .decomposition(x, "Z")
