test_that("simple varide works", {
  test_cases <- list(
    list(matrix(c(1), 1, 1), FALSE),
    list(matrix(c(1, 1), 2, 1), FALSE),
    list(matrix(c(1, 1, 1), 3, 1), TRUE),
    list(matrix(c(1, 0,
                  1, 0,
                  0, 1,
                  0, 1,
                  1, 1), 5, 2, byrow = TRUE), TRUE),
    list(matrix(c(1, 0,
                  1, 0,
                  0, 1,
                  0, 1,
                  0, 1), 5, 2, byrow = TRUE), FALSE),
    list(matrix(c(1, 0, 0,
                  1, 0, 1,
                  0, 1, 0,
                  0, 1, 1,
                  1, 1, 1), 5, 3, byrow = TRUE), FALSE),
    list(matrix(c(1, 0, 0,
                  1, 0, 1,
                  0, 1, 0,
                  0, 1, 1,
                  0, 1, 1,
                  1, 1, 1,
                  1, 1, 1), 7, 3, byrow = TRUE), TRUE)
  )

  for (test_case in test_cases) {
    expect_equal(variance_is_identified_indicator(t(test_case[[1]])), test_case[[2]])
  }
})

test_that("simulated varide works", {
  # Naive method
  brute_force_varide <- function (delta) {
    r <- NCOL(delta)
    for (q in 1:r) {
      q_subsets <- combn(r, q)  # all possible combinations \binom{r}{q}
      for (q_subset_index in seq_len(NCOL(q_subsets))) {
        q_subset <- q_subsets[, q_subset_index]
        n_nonzero_rows <- sum(rowSums(delta[, q_subset, drop = FALSE]) > 0)
        if (n_nonzero_rows <= 2 * q) {
          return(FALSE)
        }
      }
    }
    return(TRUE)
  }

  simulate_delta <- function (m, r) {
    leading_elements <- sample.int(m, r, replace = FALSE)
    delta <- matrix(0, m, r)
    delta[cbind(leading_elements, 1:r)] <- 1
    for (i in 1:r) {
      if (leading_elements[i] < m) {
        delta[(leading_elements[i] + 1):m, i] <-
          sample(c(0, 1), m - leading_elements[i], replace = TRUE,
                 prob = c(0.2, 1))
      }
    }
    # remove zero rows and columns
    if (all(delta == 0)) {  # entire matrix is zero
      return(matrix(1))
    }
    delta <- delta[rowSums(delta) > 0, , drop = FALSE]
    delta <- delta[, colSums(delta) > 0, drop = FALSE]
    return(delta)
  }

  res <- replicate(100, {
    r <- 5 + sample.int(13, 1)
    m <- 2 * r + sample.int(5, 1)
    delta <- simulate_delta(m, r)
    expect_equal(variance_is_identified_indicator(t(delta)), brute_force_varide(delta))
  })
})
