tinytest::exit_if_not(requireNamespace("haven", quietly = TRUE))

library(devtools)
devtools::install_github("chaisemartinPackages/twowayfeweights/R")
library(TwoWayFEWeights)

url = "https://raw.githubusercontent.com/chaisemartinPackages/twowayfeweights/main/wagepan_twfeweights.dta"
wagepan = haven::read_dta(url)

#
## feTR

est_fetr_known = list(
    nr_plus = 820,
    nr_minus = 147,
    nr_weights = 967,
    tot_cells = 1016,
    sum_plus = 1.010529,
    sum_minus = -0.01052899,
    mat = data.frame(
        Coef = -0.13445527172928173,
        SE = 0.07136021078287243,
        `t-stat` = -1.88417705405031,
        Correlation = -0.11825873818614765,
        row.names = "RW_educ",
        check.names = FALSE
    ),
    beta = 0.10662746641838491,
    sensibility = 0.0968691,
    sensibility2 = 3.175859
)

est_fetr = twowayfeweights(
  wagepan,                        # dataset
  "lwage", "nr", "year", "union", # Required Y, G, T, and D args
  type                = "feTR",   # start of optional args...
  summary_measures    = TRUE,
  test_random_weights = "educ"
)

for (i in names(est_fetr_known)) {
    expect_equivalent(
        est_fetr[[i]],
        est_fetr_known[[i]],
        tol = 0.0001,
        info = paste("feTR:", i)
    )
}


#
## fdTR

est_fdtr_known = list(
    nr_plus = 611,
    nr_minus = 405,
    nr_weights = 1016,
    tot_cells = 1016,
    sum_plus = 1.047636,
    sum_minus = -0.04763605,
    mat = data.frame(
        Coef = -0.06649016786275475,
        SE = 0.028938371061793897,
        `t-stat` = -2.297647221427017,
        Correlation = -0.09947980072899107,
        row.names = "RW_educ",
        check.names = FALSE
    ),
    beta = 0.06009597,
    sensibility = 0.03209515,
    sensibility2 = 0.5799135
)

est_fdtr = twowayfeweights(
    wagepan,
    "diff_lwage", "nr", "year", "diff_union", # use differenced versions of Y and D
    type = "fdTR", # changed
    D0 = "union", # added (req'd arg for fdTR type)
    summary_measures = TRUE,
    test_random_weights = "educ"
)

for (i in names(est_fdtr_known)) {
    expect_equivalent(
        est_fdtr[[i]],
        est_fdtr_known[[i]],
        tol = 0.0001,
        info = paste("fdTR:", i)
    )
}
 