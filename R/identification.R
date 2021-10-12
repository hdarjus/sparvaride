default_tolerance <- 2 * sqrt(.Machine$double.eps)

#' @param m factor loading matrix with dimensions factors x observations
#' @param tolerance tolerance for numeric errors; if \code{abs(a) < tolerance} then
#' \code{a} is considered to be zero
#' @export
variance_is_identified <- function (m, tolerance = default_tolerance) {
  m_glt <- to_glt(m, tolerance = tolerance)
  variance_is_identified_identicator(m_glt != 0)
}

compute_qr_decomposition <- base::qr

#' apply tolerance after every step and not just at the end
compute_rref <- function (m, tolerance) {
  column <- 1L
  row <- 1L
  pivot_columns <- c()
  m[abs(m) < tolerance] <- 0
  while (column <= NCOL(m) && row <= NROW(m)) {
    nonzero_indices <- row - 1L + which(m[seq(row, NROW(m)), column] != 0)
    if (length(nonzero_indices) > 0) {
      m[row, ] <- m[row, ] + m[tail(nonzero_indices, 1), ]  # handle the case when m[row, column] == 0
      m[row, ] <- m[row, ] / m[row, column]
      m[-row, ] <- m[-row, ] - m[-row, column] %o% m[row, ]
      pivot_columns <- c(pivot_columns, column)
      row <- row + 1L
    }
    m[abs(m) < tolerance] <- 0
    column <- column + 1L
  }
  list(matrix = m, pivot_columns = pivot_columns)
}

#' @param m factor loading matrix with dimensions factors x observations
#' @param tolerance tolerance for numeric errors; if \code{abs(a) < tolerance} then
#' \code{a} is considered to be zero
to_glt <- function (m, tolerance) {
  pivot_columns <- compute_rref(m, tolerance)$pivot_columns
  qr <- compute_qr_decomposition(m[, pivot_columns])
  q <- base::qr.Q(qr)
  result <- crossprod(q, m)
  result[abs(result) < tolerance] <- 0
  result
}
