test_that("multiplication works", {
  expect_false(variance_may_be_identified_indicator(matrix(c(1, 1, 1, 1,
                                                             1, 1, 1, 0),
                                                           2, 4,
                                                           byrow = TRUE)))
  expect_false(variance_may_be_identified_indicator(matrix(c(1, 1),
                                                           1, 2,
                                                           byrow = TRUE)))
  expect_false(variance_may_be_identified_indicator(matrix(c(1, 1, 1, 1, 1, 1,
                                                             1, 1, 1, 0, 0, 0,
                                                             1, 1, 1, 0, 0, 0,
                                                             1, 1, 1, 0, 0, 0),
                                                           4, 6,
                                                           byrow = TRUE)))
  expect_true(variance_may_be_identified_indicator(matrix(c(1, 1, 1),
                                                          1, 3,
                                                          byrow = TRUE)))
  expect_true(variance_may_be_identified_indicator(matrix(c(1, 1, 1, 1, 1,
                                                            1, 1, 1, 0, 0),
                                                          2, 5,
                                                          byrow = TRUE)))
  expect_true(variance_may_be_identified_indicator(matrix(c(1, 1, 1, 1, 1, 1, 1, 1, 1,
                                                            1, 1, 1, 1, 1, 1, 0, 0, 0,
                                                            1, 1, 1, 0, 0, 0, 1, 1, 1,
                                                            1, 1, 1, 0, 0, 1, 1, 1, 0),
                                                          4, 9,
                                                          byrow = TRUE)))
})
