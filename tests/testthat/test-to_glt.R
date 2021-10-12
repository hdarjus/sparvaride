test_that("GLT is correctly computed", {
  sparse_matrix <- matrix(c(4, 0, 0, 0,
                            0, -2, 0, 0,
                            -1, 1, 0, 0,
                            -3, 0, 0, 0,
                            0, 2, 2, 0,
                            1, -5, 3, 3,
                            1, 1, 1, 1,
                            1, 2, 3, 4,
                            2, 1, 0, -1),
                          ncol = 4, byrow = TRUE)
  sparse_matrix <- t(sparse_matrix)
  rotation <- matrix(c(-0.72883764955166219, 0.64216679939770527, -0.237506809253953011, 0.0028280568103095094,
                       -0.21487437088930778, -0.23757796738246245, 0.028301013720456139, 0.9468816012425926854,
                       -0.44564025072535923, -0.69164998750579221, -0.505624669563513596, -0.2595549175906017858,
                       0.47331749625845398, 0.22977990762676095, -0.828935016578011608, 0.1898380364349491756),
                     nrow = 4, byrow = TRUE)
  rotated_matrix <- rotation %*% sparse_matrix
  expect_equal(to_glt(rotated_matrix, tolerance = default_tolerance), sparse_matrix)
})

test_that("compute_rref works", {
  rref <- matrix(c(2, 3, 1, 0, 0, 0, 0,
                   3, 3, 0, 0, 1, 0, 0,
                   -3, 0, 0, 0, 0, 4, 1),
                 nrow = 3, byrow = TRUE)
  rref <- rref[3:1, 7:1]
  l <- matrix(c(1, 0, 0,
                -1, 3, 0,
                2, 10, -5) / pi,
              nrow = 3, byrow = FALSE)
  l <- t(l)
  m <- l %*% rref
  expect_identical(compute_rref(m, tolerance = default_tolerance)$matrix == 0, rref == 0)
  expect_identical(compute_rref(m, tolerance = default_tolerance)$pivot_columns, c(1L, 3L, 5L))
  expect_equal(compute_rref(m, tolerance = default_tolerance)$matrix, rref)
})
