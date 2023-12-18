#' twowayfeweights
#'
#' @md
#' @description Estimates the weights attached to the two-way fixed effects
#'   regressions studied in de Chaisemartin & D'Haultfoeuille (2020a), as well
#'   as summary measures of these regressions' robustness to heterogeneous
#'   treatment effects.
#' @param df the data frame for input
#' @param Y the name of Y variable
#' @param G the name of group variable
#' @param T the name of time variable
#' @param D the name of treatment variable
#' @param type the type of estimation. One of "feTR" (the default), "feS",
#'   "fdTR", or "fdS". See section \dQuote{Estimation type} for details.
#' @param D0 the name of the mean of the treatment in group g and at period t.
#'   Note that this is a required argument if type = "fdTR".
#' @param controls the list of names of control variables, empty if not specified
#' @param weights a column of data that replaces the default weight
#' @param other_treatments the list of other treatment variables to include in the regression other than D
#' @param test_random_weights weights when this option is specified, the command
#'   estimates the correlation between each variable in varlist and the weights
#' @section Estimation type: The `type` argument specifies the type of
#'   estimation strategy. It can take one of four different values, defining a
#'   unique combination of regression strategy (either fixed-effects or
#'   first-difference) and inference assumption (either common trends on its own
#'   or common trends plus an additional assumption about treatment stability
#'   over time).
#' - "feTR": fixed-effects regression under the common trends assumption.
#' - "feS": fixed-effects regression under common trends and the assumption that groups' treatment effect does not change over time.
#' - "fdTR": first-difference regression under the common trends assumption.
#' - "fdS": first-difference regression under common trends and the assumption that groups' treatment effect does not change over time.
#' @return twowayfeweights returns data frame that contains the following columns
#'         T: time variable
#'         G: group variable
#'         weight: the result of the weight
#' @export
twowayfeweights <- function(
    df,
    Y,
    G,
    T,
    D,
    type = c("feTR", "feS", "fdTR", "fdS"),
    D0 = NULL,
    summary_measures = NULL,
    controls = c(),
    weights = NULL,
    other_treatments = c(),
    test_random_weights = c(),
    path = NULL
    ) {
  
  # type = match.arg(tolower(type))
  type = match.arg(type)
  if (type == "fdTR" && is.null(D0)) {
    stop("The `D0` argument must also be provided if `type = 'fdTR'`.\n")
  }
  
  if (length(other_treatments) > 0 && type != "feTR") {
    printf("When the other_treatments option is specified, you need to specify the type(feTR) option.")
    return(c())
  }
  
  controls_rename <- get_controls_rename(controls)
  treatments_rename <- get_treatments_rename(other_treatments)
  random_weight_rename <- get_random_weight_rename(test_random_weights)
  
  df_renamed <- twowayfeweights_rename_var(df, Y, G, T, D, D0, controls, other_treatments, test_random_weights)
  df_transformed <- twowayfeweights_transform(df_renamed, controls_rename, weights, treatments_rename)
  df_filtered <- twowayfeweights_filter(df_transformed, Y, G, T, D, D0, type, controls_rename, treatments_rename)

  if (length(other_treatments) == 0) {
    res <- twowayfeweights_calculate(df_filtered, controls_rename, type = type)
    res <- twowayfeweights_result(res$df, res$beta, random_weight_rename)
    df_result <- twowayfeweights_print_results(type, res, D, summary_measures, res$beta, random_weight_rename)
  } else {
    res <- twowayfeweights_calculate_fetr_other_treatment(df_filtered, controls_rename, treatments_rename)
    res <- twowayfeweights_result_other_treatment(res$df, treatments_rename, res$beta, random_weight_rename)
    df_result <- twowayfeweights_print_result_other_treatment(res, treatments_rename, D, res$beta, random_weight_rename)
  }
  
  if (length(path) != 0) {
    write.csv(df_result, path)
  }
  
  return(res)
  
}