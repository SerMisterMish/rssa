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
  N <- length(F)
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
  if (!force.continue && any(unlist(nsigma(x)) > 0))
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

.expand.group <- function(group, num_modes) {
  n <- length(group)
  if (n == 0L || num_modes == 0L) return(list())
  
  total <- n ^ num_modes
  mat <- matrix(numeric(total * num_modes), nrow = num_modes, ncol = total)
  
  for (i in seq_len(num_modes)) {
    each_rep <- n ^ (i - 1L)
    times_rep <- total / (n * each_rep)
    mat[i, ] <- rep(group, each = each_rep, times = times_rep)
  }
  
  asplit(mat, 2L)
}

.expand.groups <- function(groups, num_modes)
  lapply(groups, .expand.group, num_modes)

reconstruct.tssa <- function(x, groups, ...,
                             drop.attributes = FALSE, cache = TRUE) {
  out <- list();

  dec_modes <- nu(x)
  num_modes <- length(dec_modes)
  if (missing(groups))
    groups = as.list(1:min(dec_modes))
  expanded_groups <- .expand.groups(groups, num_modes)

  # Continue decomposition, if necessary
  .maybe.continue(x, groups = groups, ...)
  
  # Grab indices of pre-cached values
  info <- .get.series.info(x);

  # Do actual reconstruction
  for (i in seq_along(expanded_groups)) {
    group <- expanded_groups[[i]]
    new <- setdiff(group, info)
    cached <- intersect(group, info)
    
    if (length(new)) {
      # Do actual reconstruction
      out[[i]] <- .elseries.tssa(x, new, cache = cache);
      
      # Add pre-cached series
      if (length(cached))
        out[[i]] <- out[[i]] + .get.series(x, cached);
    } else if (length(cached)) {
      out[[i]] <- .get.series(x, cached)
    } else {
      out[[i]] <- 0. * .F(x)
      if (!is.null(x$weights))
        out[[i]][x$weights == 0] <- NA
    }
    
    # Propagate attributes
    out[[i]] <- .apply.attributes(x, out[[i]], fixup = FALSE, drop = drop.attributes)
  }
  
  # Set names
  names(out) <- .group.names(groups)
  
  # Calculate the residuals
  residuals <- .F(x)
  rgroups <- unique(unlist(expanded_groups, recursive = FALSE))
  info <- .get.series.info(x);
  rcached <- intersect(rgroups, info)
  rnew <- setdiff(rgroups, info)
  if (length(rcached))
    residuals <- residuals - .get.series(x, rcached)
  if (length(rnew))
    residuals <- residuals - .elseries.tssa(x, rnew, cache = FALSE)
  
  # Propagate attributes of residuals
  residuals <- .apply.attributes(x, residuals, fixup = FALSE, drop = drop.attributes)
  F <- .apply.attributes(x, .F(x), fixup = FALSE, drop = drop.attributes)
  
  attr(out, "residuals") <- residuals;
  attr(out, "series") <- F;
  
  # Reconstructed series can be pretty huge...
  class(out) <- paste(c(x$kind, "ssa"), "reconstruction", sep = ".")
  invisible(out);
}

.pad0 <- function(x, N) c(x, rep(0, N - length(x)))
.mvpad0 <- function(x, N) rbind(x, matrix(0, nrow = N - nrow(x), ncol = ncol(x)))

.elseries.tssa <- function(x, idx, ..., cache = TRUE) {
  if (any(max(unlist(idx)) > nsigma(x)))
    stop("Too few eigentriples computed for this decomposition")
  
  dec <- .decomposition(x)
  U <- .U(dec)
  Z <- .Z(dec)
  num_modes <- length(U)
  N <- x$length
  
  u_fft <- .get.or.create(x, "tssa.u.fft", 
                          lapply(1:num_modes, \(m) mvfft(.mvpad0(U[[m]], N))))
  
  w <- 1 / .hweights(x) / N
  res <- numeric(prod(x$length));
  for (i in idx) {
    new <- w * do.call(`[`, c(list(Z), as.list(i))) * 
          .hankelize.one.tssa(x, i, u_fft, complex_result = is.complex(.F))
    if (cache) .set.series(x, new, list(i))
    res <- res + new
  }
  
  res;
}

.hankelize.one.tssa <- function(x, idx, u_fft, complex_result = FALSE) {
  res <- fft(Reduce(`*`, lapply(seq_along(u_fft), \(j) u_fft[[j]][, idx[j]])), inverse = TRUE)
  if (complex_result) res else Re(res)
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
