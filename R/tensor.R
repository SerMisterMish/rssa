#   R package for Singular Spectrum Analysis
#   Copyright (c) 2026 Nikita Khromov <hromovn0@gmail.com>
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


n.unfold <- function(tens, n) {
  if (length(n) != 1)
    stop("Only one mode 'n' must be provided")
  
  modes <- dim(tens)
  num_modes <- length(modes)
  if (n > num_modes)
    stop(sprintf("Provided n = %d, but input tensor has only %d modes", n, num_modes))
  
  result <- aperm(tens, c(n, seq(num_modes)[-n]))
  dim(result) <- c(modes[n], prod(modes[-n]))
  
  result
}

n.fold <- function(mat, n, modes) {
  mat_dims <- dim(mat)
  if (prod(modes) != prod(mat_dims))
    stop(sprintf("Product of modes (%d) must be equal to the number of elements in matrix (%d)",
                 prod(modes), prod(mat_dims)))
  
  num_modes <- length(modes)
  if (n < 1 || n > num_modes)
    stop(sprintf("Folding mode 'n' must be between 1 and number of modes (%d). Got n = %d", 
                 num_modes, n))
  
  if (mat_dims[1] != modes[n])
    stop(sprintf("Number of matrix rows (%d) must be equal to modes[n] (%d)", mat_dims[1], modes[n]))

  aperm(array(mat, dim = c(modes[n], modes[-n])), 
              match(1:num_modes, c(n, (1:num_modes)[-n])))
}

.reunfold <- function(mat, tens_modes, old_unfold_mode, new_unfold_mode) {
  if (length(old_unfold_mode) != 1 || length(new_unfold_mode) != 1)
    stop("Both 'old_unfold_mode' and 'new_unfold_mode' must be single integers")
  
  num_modes <- length(tens_modes)
  if (old_unfold_mode < 1 || old_unfold_mode > num_modes || 
      new_unfold_mode < 1 || new_unfold_mode > num_modes) {
    stop(sprintf("Unfolding modes must be between 1 and number of tensor modes (%d)", num_modes))
  }
  
  mat_dims <- dim(mat)
  if (prod(tens_modes) != prod(mat_dims))
    stop(sprintf("Product of 'tens_modes' (%d) must be equal to the number of elements in 'mat' (%d)",
                 prod(tens_modes), prod(mat_dims)))
  
  if (mat_dims[1] != tens_modes[old_unfold_mode])
    stop(sprintf("Number of matrix rows (%d) must be equal to 'tens_modes[old_unfold_mode]' (%d)", 
                 mat_dims[1], tens_modes[old_unfold_mode]))
  
  if (old_unfold_mode == new_unfold_mode)
    return(mat)
  
  old_order <- c(old_unfold_mode, seq_len(num_modes)[-old_unfold_mode])
  new_order <- c(new_unfold_mode, seq_len(num_modes)[-new_unfold_mode])
  
  dim(mat) <- tens_modes[old_order]
  perm <- match(new_order, old_order)
  result <- aperm(mat, perm)
  dim(result) <- c(tens_modes[new_unfold_mode], prod(tens_modes[-new_unfold_mode]))
  
  result
}

.fnorm <- function(x) sqrt(sum(if (is.complex(x)) Re(x)^2 + Im(x)^2 else x^2))

n.mode.fnorm <- function(tens, modes = seq_along(dim(tens))) {
  sq_tens <- if (is.complex(tens)) Re(tens)^2 + Im(tens)^2 else tens^2
  all_modes <- seq_along(dim(tens))
  
  # Permuting the target mode to the front and using 'rowSums' with dims = 1 
  # appears to be faster than just calling 'margin.table', as under the hood it calls apply()
  lapply(modes, function(m) {
    if (m == 1L) sqrt(rowSums(sq_tens, dims = 1L))
    else sqrt(rowSums(aperm(sq_tens, c(m, all_modes[-m])), dims = 1L))
    })
}

n.mode.prod <- function(tens, mat, mode) {
  stopifnot(length(mode) == 1)
  
  mat_dims <- dim(mat)
  if (length(mat_dims) != 2)
    stop("'mat' parameter must be a matrix")
  
  modes <- dim(tens)
  if (mat_dims[2] != modes[mode])
    stop(sprintf("Provided tensor mode (%d) must be equal to the number of columns in 'mat' (%d)", modes[mode], mat_dims[2]))
  num_modes <- length(modes)
  
  result_modes <- modes
  result_modes[mode] <- mat_dims[1]
  
  n.fold(mat %*% n.unfold(tens, mode), mode, result_modes)
}
  
n.mode.prod.list <- function(tens, mats, modes) {
  tens_modes <- dim(tens)
  tens_num_modes <- length(tens_modes)
  
  if (any(modes < 1 | modes > tens_num_modes))
    stop(sprintf("All modes must be between 1 and number of tensor modes (%d)", tens_num_modes))
  if (length(mats) != length(modes))
    stop("'mats' and 'modes' lengths do not match")
  
  if (length(modes) == 1)
    return(n.mode.prod(tens, mats[[1]], mode))
  
  is_duplicated_mode <- duplicated(modes)
  if (any(is_duplicated_mode)) {
    has_duplicate <- is_duplicated_mode | duplicated(modes, fromLast = TRUE)
    duplicated_modes <- modes[has_duplicate]
    duplicated_mode_mats <- mats[has_duplicate]
    unique_dmodes <- unique(duplicated_modes)
    factor_dmodes <- factor(duplicated_modes, levels = unique_dmodes)
    
    grouped_mats <- split(duplicated_mode_mats, factor_dmodes)
    
    pre_multiplied_mats <- unname(lapply(grouped_mats, function(m) 
      Reduce(function(m1, m2) m2 %*% m1, m)))
    
    mats <- c(mats[!has_duplicate], pre_multiplied_mats)
    modes <- c(modes[!has_duplicate], unique_dmodes)
  }
  
  mat <- mats[[1]]
  mode <- modes[1]
  tens_unf <- mat %*% n.unfold(tens, mode)
  old_unfold_mode <- mode
  tens_modes[mode] <- dim(mat)[1]
  for (i in seq_along(modes)[-1]) {
    mat <- mats[[i]]
    mode <- modes[i]
    mat_dims <- dim(mat)
    if (tens_modes[mode] != mat_dims[2])
      stop(sprintf("The %d-mode of the tensor (%d) is not equal to the number of columns in the corresponding matrix (%d)", 
                   mode, tens_modes[mode], mat_dims[2]))
    
    tens_unf <- mat %*% .reunfold(tens_unf, tens_modes, old_unfold_mode, mode)
    old_unfold_mode <- mode
    tens_modes[mode] <- mat_dims[1]
  }
  
  n.fold(tens_unf, old_unfold_mode, tens_modes)
}

hosvd <- function(tens,
                  unfolds,
                  ranks = NULL,
                  svd.method = c("svd", "eigen", "primme"),
                  ...,
                  report.progress = FALSE,
                  calc.est = FALSE)
{
  if (sum(ranks <= 0) != 0)
    stop("ranks must be positive")
  
  svd.method <- match.arg(svd.method)
  if (identical(svd.method, "primme") &&
      !requireNamespace("PRIMME", quietly = TRUE))
    stop("PRIMME package is required for SVD method `primme'")
  
  stopifnot(inherits(tens, "tensor") && inherits(unfolds, "unfoldings"))
  
  if (all(abs(tens) < sqrt(.Machine$double.eps)))
    stop("Zero tensor detected")
  
  modes <- dim(tens)
  num_modes <- length(modes)
  if (is.null(ranks)) {
    ranks <- modes
  }
  else {
    if (!(length(ranks) == 1 || length(ranks) == num_modes))
      stop("Length of the 'ranks' parameter must be either 1 or the same as the number of tensor modes")
    if (length(ranks) == 1)
      ranks <- rep(ranks, num_modes)
    if (any(ranks > modes))
      stop("Each rank must be not greater than the corresponding mode")
  }

  
  if (report.progress)
    pb <- utils::txtProgressBar(min = 0, max = num_modes, style = 3)
  U_list <- vector("list", num_modes)
  for (m in 1:num_modes) {
    R <- unfolds$Re[[m]]
    Ilist <- unfolds$Im
    if (!is.null(Ilist))
      I <- Ilist[[m]]
    else
      I <- NULL
    
    if (identical(svd.method, "svd") || identical(svd.method, "eigen")) {
      if (is.null(I)) {
        unfold <- as.matrix(R)
      } else {
        unfold <- as.matrix(R) + 1i * as.matrix(I)
      }
      
      if (identical(svd.method, "svd"))
        U_list[[m]] <- svd(unfold, nu = min(ranks[m], ncol(unfold)))$u
      else
        U_list[[m]] <- eigen(tcrossprod(unfold, Conj(unfold)), symmetric = TRUE)$vectors[, 1:min(ranks[m], ncol(unfold)), drop = FALSE]
      
    } else if (identical(svd.method, "primme")) {
      matmul <- function(x, y) {
        if (is.matrix(y))
          apply(y, 2, ematmul, emat = x, transposed = FALSE)
        else
          ematmul(x, y, transposed = FALSE)
      }
      
      tmatmul <- function(x, y) {
        if (is.matrix(y))
          apply(y, 2, ematmul, emat = x, transposed = TRUE)
        else
          ematmul(x, y, transposed = TRUE)
      }
      
      if (is.null(I)) {
        A <- function(x, trans) {
          if (identical(trans, "c"))
            tmatmul(R, x)
          else
            matmul(R, x)
        }
        U_list[[m]] <- PRIMME::svds(
          A,
          NSvals = min(ranks[m], ncol(R)),
          m = nrow(R),
          n = ncol(R),
          isreal = TRUE,
          ...
        )$u
      } else {
        A <- function(x, trans) {
          rX <- Re(x)
          iX <- Im(x)
          if (identical(trans, "c")) {
            (tmatmul(R, rX) + tmatmul(I, iX)) +
              1i * (tmatmul(R, iX) - tmatmul(I, rX))
          } else {
            (matmul(R, rX) - matmul(I, iX)) +
              1i * (matmul(R, iX) + matmul(I, rX))
          }
        }
        U_list[[m]] <- PRIMME::svds(
          A,
          NSvals = min(ranks[m], ncol(R)),
          m = nrow(R),
          n = ncol(R),
          isreal = FALSE,
          ...
        )$u
      }
    } else
      stop("Unsupported svd method")
    
    if (report.progress)
      utils::setTxtProgressBar(pb, m)
  }
  if (report.progress)
    close(pb)
  
  Z <- n.mode.prod.list(tens, lapply(U_list, (function(x) Conj(t(x)))), seq_len(num_modes))
  
  if (calc.est) {
    est <- n.mode.prod.list(Z, U_list, seq_len(num_modes))
    resid <- .fnorm(est - tens)
    result <- list(
      Z = Z,
      U = U_list,
      est = est,
      fnorm_resid = resid
    )
  } else {
    result <- list(
      Z = Z,
      U = U_list
    ) 
  }
  
  result
}

hooi <- function(tens,
                 unfolds,
                 ranks,
                 svd.method = c("svd", "eigen", "primme"),
                 max_iter = 25,
                 tol = 1e-05,
                 ...,
                 report.progress = FALSE)
{
  if (sum(ranks <= 0) != 0)
    stop("ranks must be positive")
  
  svd.method <- match.arg(svd.method)
  if (identical(svd.method, "primme") &&
      !requireNamespace("PRIMME", quietly = TRUE))
    stop("PRIMME package is required for SVD method `primme'")
  
  stopifnot(inherits(tens, "tensor") && inherits(unfolds, "unfoldings"))
  
  if (all(abs(tens) < sqrt(.Machine$double.eps)))
    stop("Zero tensor detected")
  
  modes <- dim(tens)
  num_modes <- length(modes)

  if (!(length(ranks) == 1 || length(ranks) == num_modes))
    stop("Length of the 'ranks' parameter must be either 1 or the same as the number of tensor modes")
  if (length(ranks) == 1)
    ranks <- rep(ranks, num_modes)
  if (any(ranks > modes))
    stop("Each rank must be not greater than the corresponding mode")
  
  U_list <- vector("list", num_modes)
  for (m in 1:num_modes) {
    R <- unfolds$Re[[m]]
    Ilist <- unfolds$Im
    if (!is.null(Ilist))
      I <- Ilist[[m]]
    else
      I <- NULL
    
    if (identical(svd.method, "svd") || identical(svd.method, "eigen")) {
      if (is.null(I)) {
        unfold <- as.matrix(R)
      } else {
        unfold <- as.matrix(R) + 1i * as.matrix(I)
      }
      
      if (identical(svd.method, "svd"))
        U_list[[m]] <- svd(unfold, nu = min(ranks[m],  ncol(unfold)))$u
      else 
        U_list[[m]] <- eigen(tcrossprod(unfold, Conj(unfold)), symmetric = TRUE)$vectors[, 1:min(ranks[m], ncol(unfold)), drop = FALSE]
      
    } else if (identical(svd.method, "primme")) {
      matmul <- function(x, y) {
        if (is.matrix(y))
          apply(y, 2, ematmul, emat = x, transposed = FALSE)
        else
          ematmul(x, y, transposed = FALSE)
      }
      
      tmatmul <- function(x, y) {
        if (is.matrix(y))
          apply(y, 2, ematmul, emat = x, transposed = TRUE)
        else
          ematmul(x, y, transposed = TRUE)
      }
      
      if (is.null(I)) {
        A <- function(x, trans) {
          if (identical(trans, "c"))
            tmatmul(R, x)
          else
            matmul(R, x)
        }
        U_list[[m]] <- PRIMME::svds(
          A,
          NSvals = ranks[m],
          m = nrow(R),
          n = ncol(R),
          isreal = TRUE,
          ...
        )$u
      } else {
        A <- function(x, trans) {
          rX <- Re(x)
          iX <- Im(x)
          if (identical(trans, "c")) {
            (tmatmul(R, rX) + tmatmul(I, iX)) +
              1i * (tmatmul(R, iX) - tmatmul(I, rX))
          } else {
            (matmul(R, rX) - matmul(I, iX)) +
              1i * (matmul(R, iX) + matmul(I, rX))
          }
        }
        U_list[[m]] <- PRIMME::svds(
          A,
          NSvals = ranks[m],
          m = nrow(R),
          n = ncol(R),
          isreal = FALSE,
          ...
        )$u
      }
    } else
      stop("Unsupported svd method")
  }
  
  modes_seq <- seq_len(num_modes)
  tens_norm <- .fnorm(tens)
  curr_iter <- 1
  converged <- FALSE
  fnorm_resid <- rep(0, max_iter)
  CHECK_CONV <- function(Z, U_list) {
    est <- n.mode.prod.list(Z, U_list, modes_seq)
    curr_resid <- .fnorm(tens - est)
    fnorm_resid[curr_iter] <<- curr_resid
    if (curr_iter > 1 && abs(curr_resid - fnorm_resid[curr_iter - 1]) / tens_norm < tol)
      return(TRUE)
    else 
      return(FALSE)
  }
  
  if (report.progress)
    pb <- utils::txtProgressBar(min = 0, max = max_iter, style = 3)
  
  while ((curr_iter < max_iter) && (!converged)) {
    for (m in modes_seq) {
      X <- n.mode.prod.list(tens, lapply(U_list[-m], (function(x) Conj(t(x)))), modes_seq[-m])
      
      # TODO: Change this part so it uses primme if svd.method = "primme"
      unfold <- n.unfold(X, m)
      if (identical(svd.method, "svd"))
        U_list[[m]] <- svd(unfold, nu = ranks[m])$u
      else if (identical(svd.method, "eigen") || identical(svd.method, "primme"))
        U_list[[m]] <- eigen(tcrossprod(unfold, Conj(unfold)), symmetric = TRUE)$vectors[, 1:min(ranks[m], ncol(unfold)), drop = FALSE]
    }
    
    Z <- n.mode.prod(X, Conj(t(U_list[[num_modes]])), num_modes)
    
    if (CHECK_CONV(Z, U_list)) {
      converged <- TRUE
      if (report.progress)
        utils::setTxtProgressBar(pb, max_iter)
    }
    else {
      if (report.progress)
        utils::setTxtProgressBar(pb, curr_iter)
      curr_iter <- curr_iter + 1
    }
  }
  if (report.progress)
    close(pb)
  
  fnorm_resid <- fnorm_resid[fnorm_resid != 0]
  norm_percent <- (1 - fnorm_resid[length(fnorm_resid)] / tens_norm) * 100
  est <- n.mode.prod.list(Z, U_list, 1:num_modes)
  invisible(
    list(
      Z = Z,
      U = U_list,
      conv = converged,
      iter = curr_iter,
      est = est,
      norm_percent = norm_percent,
      fnorm_resid = tail(fnorm_resid, 1),
      all_resids = fnorm_resid
    )
  )
}
