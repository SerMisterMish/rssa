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


.plot.ssa.values.tssa <- function(x, ..., numvalues, plot.type = "b") {
  dots <- list(...)
  
  A <- Reduce(c, sapply(numvalues, seq_len))
  B <- Reduce(c, sapply(seq_along(numvalues), \(m) x$sigma[[m]][1:numvalues[m]]))
  mode <- Reduce(c, sapply(seq_along(numvalues), \(m) rep(m, numvalues[m])))
  d <- data.frame(A = A, B = B, mode = mode)
  
  # Provide convenient defaults
  dots <- .defaults(dots,
                    type = plot.type,
                    xlab = "Index",
                    ylab = "norms",
                    main = "Component norms",
                    grid = TRUE,
                    auto.key = list(title = "Mode"),
                    scales = list(y = list(log = TRUE)),
                    par.settings = list(plot.symbol = list(pch = 20)))
  
  do.call("xyplot",
          c(list(x = B ~ A, groups = mode, data = d, ssaobj = x), dots))
}

.contribution.tssa <- function(x, idx, ...) {
  isfcor <- .is.frobenius.orthogonal(x, idx, ...)
  if (!isTRUE(isfcor))
    warning(sprintf("Elementary matrices are not F-orthogonal (max F-cor is %s). Contributions can be irrelevant",
                    format(isfcor, digits = 3)))
  
  lapply(contributions(x, idx), \(cv) round(100 * cv, digits = 2))
}

prepanel.mode.vectors <- function(x, y, ssaobj, mode, symmetric = FALSE) {
  V <- ssaobj$U[[mode]][, y]
  if (is.complex(V)) V <- c(Re(V), Im(V))
  U <- if (identical(x, y)) 1:length(V) else ssaobj$U[[mode]][, x]
  if (is.complex(U)) U <- c(Re(U), Im(U))
  
  res <- prepanel.default.xyplot(U, V)
  if (symmetric) {
    res$ylim <- range(V, -V)
    if (!identical(x, y))
      res$xlim <- range(U, -U)
  }
  
  res
}

panel.mode.vectors <- function(x, y, ssaobj, mode, ..., ref = FALSE) {
  V <- ssaobj$U[[mode]][, y]
  Vcompl <- is.complex(V)
  if (Vcompl) V <- c(Re(V), Im(V))
  U <- if (identical(x, y)) 1:length(V) else ssaobj$U[[mode]][, x]
  Ucompl <- is.complex(U)
  if (Ucompl) U <- c(Re(U), Im(U))
  
  if (ref) {
    panel.abline(h = 0, ..., reference = TRUE)
    if (!identical(x, y))
      panel.abline(v = 0, ..., reference = TRUE)
  }
  
  groups <- if (Vcompl || Ucompl) 2 else 1
  panel.superpose(U, V,
                  panel.groups = panel.xyplot,
                  groups = gl(n = groups, length(V) / groups), ..., subscripts = 1:length(U))
}

.plot.ssa.vectors.tssa <- function(x, ...,
                                   what = 1:num_modes,
                                   plot.contrib, idx, plot.type = "l") {
  dots <- list(...)
  num_modes <- length(.traj.dim(x))
  
  # Fix 'what' passed from plot.ssa
  if (identical(what, c("eigen", "factor")) || 
      identical(what, c("eigen")) || 
      identical(what, c("factor")))
    what <- 1:num_modes
  
  ## FIXME: check for proper lengths
  d <- data.frame(A = idx, B = idx)
  
  # Provide convenient defaults
  dots <- .defaults(dots,
                    type = plot.type,
                    xlab = "",
                    ylab =  "",
                    as.table = TRUE,
                    scales = list(draw = FALSE, relation = "free"),
                    aspect = 1,
                    symmetric = TRUE,
                    ref = TRUE)
  
  if ("main" %in% names(dots))
    main <- if (length(dots$main) == 1) rep(dots$main, length(what)) else dots$main
  else 
    main <- paste(what, "-mode Eigenvectors", sep = "")
    
  contrib_strings <- .contribution(x, idx, ...)
  res_plots <- vector("list", length(what))
  
  for (i in seq_along(what)) {
    m <- what[i]
    plot.formula <-
      A ~ B | factor(A,
                     labels = if (!plot.contrib) A else paste(A, " (", contrib_strings[[m]], "%)", sep = ""))
    
    dots$main <- main[i]
    
    res_plots[[m]] <- do.call("xyplot",
            c(list(x = plot.formula,
                   data = d, ssaobj = x, mode = m,
                   panel = panel.mode.vectors,
                   prepanel = prepanel.mode.vectors),
              dots))
  }
  res_plots
}

.plot.ssa.paired.tssa <- function(x, ...,
                                 what = 1:num_modes,
                                 plot.contrib, idx, idy, plot.type = "l") {
  dots <- list(...)  
  num_modes <- length(.traj.dim(x))
  
  # Fix 'what' passed from plot.ssa
  if (identical(what, c("eigen", "factor")) || 
      identical(what, c("eigen")) || 
      identical(what, c("factor")))
    what <- 1:num_modes
  
  
  # FIXME: check for proper lengths
  d <- data.frame(A = idx, B = idy)
  
  # Provide convenient defaults
  dots <- .defaults(dots,
                    type = plot.type,
                    xlab = "",
                    ylab =  "",
                    as.table = TRUE,
                    scales = list(draw = FALSE, relation = "free"),
                    aspect = 1,
                    symmetric = TRUE,
                    ref = TRUE)
  
  if ("main" %in% names(dots))
    main <- if (length(dots$main) == 1) rep(dots$main, length(what)) else dots$main
  else
    main <- paste("Pairs of ", what, "-mode eigenvectors", sep = "")
  
  contrib_strings_x <- .contribution(x, idx, ...)
  contrib_strings_y <- .contribution(x, idy, ...)
  
  res_plots <- vector("list", length(what))
  
  for (i in seq_along(what)) {
    m <- what[i]
    plot.formula <- A ~ B | factor(A,
                                   labels = if (!plot.contrib) paste(A, "vs", B)
                                   else paste(A, " (", contrib_strings_x[[m]], "%) vs ", B, " (", contrib_strings_y[[m]], "%)", sep = ""))
    
    dots$main <- main[i]
    
    res_plots[[m]] <- do.call("xyplot",
                              c(list(x = plot.formula,
                                     data = d, ssaobj = x, mode = m,
                                     panel = panel.mode.vectors,
                                     prepanel = prepanel.mode.vectors),
                                dots))
  }
  res_plots
}
