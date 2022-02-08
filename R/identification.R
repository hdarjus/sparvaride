default_tolerance <- 2 * sqrt(.Machine$double.eps)

#' Functions for working with the GLT framework
#'
#' @param m factor loading matrix with dimensions factors x observations
#' @param tolerance tolerance for numeric errors; if \code{abs(a) < tolerance} then
#' \code{a} is considered to be zero
#' @rdname glt_framework
#' @export
variance_is_identified <- function (m, tolerance = default_tolerance) {
  m_glt <- to_glt(m, tolerance = tolerance)
  delta <- m_glt != 0
  if (NROW(m_glt) < NROW(m) || any(rowSums(delta) == 0)) {  # zero column in m
    FALSE
  } else {
    delta <- delta[, colSums(delta) > 0]
    variance_is_identified_identicator(delta)
  }
}

#' @rdname glt_framework
#' @export
to_glt <- function (m, tolerance) {
  pivot_columns <- compute_rref(m, tolerance)$pivot_columns
  qr <- compute_qr_decomposition(m[, pivot_columns])
  q <- base::qr.Q(qr)
  result <- crossprod(q, m)
  result[abs(result) < tolerance] <- 0
  result
}

compute_qr_decomposition <- base::qr

# apply tolerance after every step and not just at the end
compute_rref <- function (m, tolerance) {
  available_rows <- rep(TRUE, NROW(m))
  column <- 1L
  pivot_columns <- c()
  m[abs(m) < tolerance] <- 0
  while (column <= NCOL(m) && any(available_rows)) {
    if (any(m[available_rows, column] != 0)) {
      row <- which(available_rows)[which.max(abs(m[available_rows, column]))]  # index_of_largest
      m[row, ] <- m[row, ] / m[row, column]
      m[-row, ] <- m[-row, ] - m[-row, column] %o% m[row, ]
      pivot_columns <- c(pivot_columns, column)
      available_rows[row] <- FALSE
    }
    m[abs(m) < tolerance] <- 0
    column <- column + 1L
  }
  list(matrix = m, pivot_columns = pivot_columns)
}

#' @param delta 0-1 identicator of factor loading matrix with dimensions factors x observations
#' @rdname glt_framework
#' @export
is_glt <- function (delta) {
    delta <- delta[rowSums(delta) > 0, colSums(delta) > 0, drop = FALSE]
    if (NROW(delta) == 0) {
        TRUE
    } else {
        cs <- which(colSums(delta) == 1L)
        length(cs) > 0 && is_glt(delta[which(delta[, cs[1]] != 0), -cs[1], drop = FALSE])
    }
}
