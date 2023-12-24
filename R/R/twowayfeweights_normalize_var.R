#' Internal function for normalizing variables
#' 
#' @param df A data frame.
#' @param varname Variable that you want to normalize.
#' @returns A list.
#' @importFrom magrittr %>%
#' @importFrom rlang :=
#' @noRd
twowayfeweights_normalize_var = function(df, varname) {
  
  .data = NULL
  
  var = rlang::sym(varname)
  sdf = df %>%
    dplyr::group_by(.data$G, .data$T) %>%
    dplyr::summarise(tmp_mean_gt = mean(!!var), tmp_sd_gt = stats::sd(!!var))
  
  tmp_sd_gt_sum = sum(sdf$tmp_sd_gt, na.rm=TRUE)
  if (tmp_sd_gt_sum > 0) {
    df = df %>% 
      dplyr::left_join(sdf, by=c("T", "G")) %>%
      dplyr::mutate(!!var := .data$tmp_mean_gt) %>%
      dplyr::select(-.data$tmp_mean_gt) %>%
      dplyr::select(-.data$tmp_sd_gt)
  }
  
  return(list(retcode = (tmp_sd_gt_sum > 0), df = df))
}
