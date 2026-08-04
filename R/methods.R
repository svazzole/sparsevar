#' @name print.var
#' @title Print a sparse VAR model
#'
#' @description Compact one-screen summary of a fitted or simulated VAR object.
#' Every figure is read from the object; nothing is re-estimated.
#'
#' @param x an object of class `var`, as returned by [fit_var()] or
#'   [simulate_var()].
#' @param digits number of significant digits used for `lambda` and the MSE.
#' @param ... ignored, present for compatibility with the generic.
#'
#' @return `x`, invisibly.
#'
#' @examples
#' sim <- simulate_var(n = 5, p = 1, nobs = 100)
#' fit <- fit_var(sim$series, p = 1)
#' print(fit)
#'
#' @export
print.var <- function(x, digits = 4, ...) {
  # fit_varx() tags its result "var", not "varx" -- only simulate_varx() sets
  # the varx class -- so a fitted VARX would print as a plain VAR and its
  # exogenous block would go unmentioned. Detect it by structure instead.
  # Fixing the class itself is a behaviour change and is tracked separately.
  if (!is.null(x$B) && !is.null(x$newA)) {
    return(print.varx(x, digits = digits, ...))
  }

  p <- length(x$A)
  nc <- if (p > 0) ncol(x$A[[1]]) else NA_integer_

  cat(sprintf("Sparse VAR(%d)\n", p))
  cat(sprintf(
    "Series: %s    Observations: %s\n",
    format(nc), format(nrow(x$series))
  ))
  cat(print_fit_line(x, digits))
  cat(print_density_line(x$A, "Non-zero coefficients"))

  if (!is.null(x$mse)) {
    cat(sprintf("Mean squared error: %s\n", signif(x$mse, digits)))
  }

  invisible(x)
}

#' @name print.vecm
#' @title Print a sparse VECM model
#'
#' @description Compact one-screen summary of a fitted VECM object. In addition
#' to the VAR fields it reports the rank of the error-correction matrix
#' \eqn{\Pi}, which is what distinguishes a VECM from a VAR in levels.
#'
#' @param x an object of class `vecm`, as returned by [fit_vecm()].
#' @param digits number of significant digits used for `lambda` and the MSE.
#' @param ... ignored, present for compatibility with the generic.
#'
#' @return `x`, invisibly.
#'
#' @export
print.vecm <- function(x, digits = 4, ...) {
  # For a VECM the Gamma matrices number p, whereas the levels-VAR
  # representation in $A has p + 1 of them, so $G is the right thing to count.
  p <- length(x$G)
  nc <- ncol(x$Pi)

  cat(sprintf("Sparse VECM(%d)\n", p))
  cat(sprintf(
    "Series: %s    Observations: %s\n",
    format(nc), format(nrow(x$series))
  ))
  cat(print_fit_line(x, digits))

  pi_rank <- tryCatch(qr(x$Pi)$rank, error = function(e) NA_integer_)
  cat(sprintf(
    "Rank of Pi: %s of %s    non-zero entries: %s of %s\n",
    format(pi_rank), format(nc),
    format(sum(x$Pi != 0)), format(length(x$Pi))
  ))
  cat(print_density_line(x$G, "Non-zero in Gamma"))

  if (!is.null(x$mse)) {
    cat(sprintf("Mean squared error: %s\n", signif(x$mse, digits)))
  }

  invisible(x)
}

#' @name print.varx
#' @title Print a sparse VARX model
#'
#' @description Compact one-screen summary of a fitted VARX object. The
#' endogenous and exogenous coefficient blocks are reported separately, since
#' they are the two things a VARX adds over a VAR.
#'
#' @param x an object of class `varx`, as returned by [fit_varx()].
#' @param digits number of significant digits used for `lambda` and the MSE.
#' @param ... ignored, present for compatibility with the generic.
#'
#' @return `x`, invisibly.
#'
#' @export
print.varx <- function(x, digits = 4, ...) {
  # $A holds the stacked endogenous-plus-exogenous system; $newA and $B are the
  # two blocks a user actually reads, so the dimensions come from those.
  n_endo <- nrow(x$newA[[1]])
  n_exo <- ncol(x$B[[1]])

  cat(sprintf("Sparse VARX(%d)\n", length(x$newA)))
  cat(sprintf(
    "Endogenous series: %s    Exogenous: %s    Observations: %s\n",
    format(n_endo), format(n_exo), format(nrow(x$series))
  ))
  cat(print_fit_line(x, digits))
  cat(print_density_line(x$newA, "Non-zero endogenous"))
  cat(print_density_line(x$B, "Non-zero exogenous"))

  if (!is.null(x$mse)) {
    cat(sprintf("Mean squared error: %s\n", signif(x$mse, digits)))
  }

  invisible(x)
}

#' @name summary.var
#' @title Summarize a sparse VAR model
#'
#' @description The longer view of a fitted VAR: the density of each lag matrix
#' separately rather than pooled, and the residual standard deviation of each
#' equation. Both are what one wants when judging whether the fit has recovered
#' a plausible structure.
#'
#' @param object an object of class `var`.
#' @param ... ignored, present for compatibility with the generic.
#'
#' @return an object of class `summary.var`, a list with elements `p`,
#'   `n_series`, `n_obs`, `penalty`, `method`, `lambda`, `mse`, `density`
#'   (one row per lag) and `residual_sd` (one entry per equation).
#'
#' @details `residual_sd` is computed from `$residuals`, discarding the first
#'   `p` rows, and is *not* the diagonal of `$sigma`. `$sigma` is a shrinkage
#'   estimate whose variances are pulled toward a common target, so its
#'   diagonal is close to constant across equations and does not describe how
#'   the equations differ.
#'
#' @examples
#' sim <- simulate_var(n = 5, p = 2, nobs = 120)
#' fit <- fit_var(sim$series, p = 2)
#' summary(fit)
#'
#' @export
summary.var <- function(object, ...) {
  p <- length(object$A)

  density <- data.frame(
    lag = seq_len(p),
    non_zero = vapply(object$A, function(a) sum(a != 0), integer(1)),
    entries = vapply(object$A, length, integer(1))
  )
  density$density <- density$non_zero / density$entries

  # Deliberately NOT sqrt(diag($sigma)): $sigma is the corpcor::cov.shrink
  # estimate, whose shrinkage pulls the variances toward a common target and
  # here returns them *identical* across equations (0.7106 five times on the
  # example in ?summary.var, against true per-equation values spanning
  # 0.678-0.755). Reporting that as a per-equation figure would be a constant
  # dressed as an estimate.
  #
  # The first p rows of $residuals are not usable either -- they carry the
  # series values rather than residuals -- so they are dropped, which is the
  # same convention the theory threshold uses for its sigma_i (utilsVAR.R).
  residual_sd <- if (!is.null(object$residuals) && nrow(object$residuals) > p) {
    apply(
      object$residuals[(p + 1):nrow(object$residuals), , drop = FALSE],
      2, stats::sd
    )
  } else {
    NULL
  }

  out <- list(
    p = p,
    n_series = if (p > 0) ncol(object$A[[1]]) else NA_integer_,
    n_obs = nrow(object$series),
    penalty = object$penalty,
    method = object$method,
    lambda = object$lambda,
    mse = object$mse,
    density = density,
    residual_sd = residual_sd
  )
  attr(out, "class") <- "summary.var"
  out
}

#' @name print.summary.var
#' @title Print a sparse VAR summary
#'
#' @param x an object of class `summary.var`, as returned by [summary.var()].
#' @param digits number of significant digits used for numeric fields.
#' @param ... ignored, present for compatibility with the generic.
#'
#' @return `x`, invisibly.
#'
#' @export
print.summary.var <- function(x, digits = 4, ...) {
  cat(sprintf("Sparse VAR(%d)\n", x$p))
  cat(sprintf(
    "Series: %s    Observations: %s\n",
    format(x$n_series), format(x$n_obs)
  ))
  cat(print_fit_line(x, digits))

  cat("\nCoefficient density by lag:\n")
  d <- x$density
  for (i in seq_len(nrow(d))) {
    cat(sprintf(
      "  A[%d]: %s of %s (%s%%)\n",
      d$lag[i], format(d$non_zero[i]), format(d$entries[i]),
      formatC(100 * d$density[i], format = "f", digits = 1)
    ))
  }

  if (!is.null(x$residual_sd)) {
    cat("\nResidual standard deviation by equation:\n")
    cat(paste0(
      "  ",
      paste(formatC(x$residual_sd, format = "f", digits = 3), collapse = "  "),
      "\n"
    ))
  }

  if (!is.null(x$mse)) {
    cat(sprintf("\nMean squared error: %s\n", signif(x$mse, digits)))
  }

  invisible(x)
}

# Internal helpers -------------------------------------------------------

# The penalty / selection / lambda line, shared by every print method so the
# three cannot drift apart. Simulated objects carry no penalty, so this returns
# an empty string rather than printing a line of NAs.
print_fit_line <- function(x, digits) {
  if (is.null(x$penalty)) {
    return("")
  }

  lambda_txt <- if (is.null(x$lambda)) {
    ""
  } else {
    sprintf(", lambda = %s", signif(x$lambda, digits))
  }
  method_txt <- if (is.null(x$method)) "" else sprintf(" by %s", x$method)

  sprintf("Penalty: %s (selected%s%s)\n", x$penalty, method_txt, lambda_txt)
}

# The "k of N (p%)" density line for a list of coefficient matrices.
print_density_line <- function(mats, label) {
  if (is.null(mats) || length(mats) == 0) {
    return("")
  }

  nz <- sum(vapply(mats, function(a) sum(a != 0), numeric(1)))
  total <- sum(vapply(mats, length, integer(1)))

  sprintf(
    "%s: %s of %s (%s%%)\n",
    label, format(nz), format(total),
    formatC(100 * nz / total, format = "f", digits = 1)
  )
}
