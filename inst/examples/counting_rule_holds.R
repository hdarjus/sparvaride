# Two example matrices
cr_holds <-
  matrix(c(1, 0, 0,
           1, 0, 1,
           0, 1, 0,
           0, 1, 1,
           0, 1, 1,
           1, 1, 1,
           1, 1, 1),
         7, 3, byrow = TRUE)

cr_does_not_hold <-
  matrix(c(1, 0, 0,
           0, 0, 1,
           0, 1, 0,
           0, 1, 0,
           0, 1, 0,
           1, 1, 1,
           1, 1, 1),
         7, 3, byrow = TRUE)

# Check if the counting rule holds
counting_rule_holds(cr_holds)
#> [1] TRUE
counting_rule_holds(cr_does_not_hold)
#> [1] FALSE
