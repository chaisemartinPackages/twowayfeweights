#' @title A print method for twowayfeweights objects
#' @name print.twowayfeweights
#' @description Printed display of twowayfeweights objects.
#' @param x A twowayfeweights object.
#' @param ... Currently ignored.
#' @inherit twowayfeweights examples
#' @export
print.twowayfeweights = function(x, ...) {
  
  treat = if (x$type %in% c("feTR", "fdTR")) {
      "ATT"
    } else if (x$type %in% c("feS", "fdS")) {
      "LATE"
    } else "BLANK"
  
  assumption = if (treat=="ATT") {
    "Under the common trends assumption"
  } else if (treat=="LATE") {
    "Under the common trends, treatment monotonicity, and if groups' treatment effect does not change over time"
  } else "BLANK"
  assumption_string = paste(
    assumption, ", ", 
    sprintf("beta estimates a weighted sum of %d %ss.", x$nr_weights, treat), 
    sep = ""
  )
  # hwidth = min(80, nchar(assumption_string))
  hwidth = nchar(assumption_string)
  
  tot_weights = x$nr_plus + x$nr_minus
  tot_sums = round(x$sum_plus + x$sum_minus, 4)
  
  cat("\n")
  cat(rep("\u2500", hwidth), sep = "")
  cat("\n")
  cat(assumption_string)
  cat("\n")
  cat(sprintf("%d %s receive a positive weight, and %d receive a negative weight.", x$nr_plus, treat, x$nr_minus))
  cat("\n")
  cat(rep("\u2500", hwidth), sep = "")
  cat("\n\n")
  
  # tmat = cbind(
  #   c("Positive weights", "Negative weights", "Total"),
  #   c(round(x$nr_plus, 2), round(x$nr_minus, 2), tot_weights),
  #   c(round(x$sum_plus, 4), round(x$sum_minus, 4), tot_sums)
  # )
  # colnames(tmat) = c(paste0("Treat. var: ", x$params$D), paste0(treat, "s"), paste0("\U03A3", " weights"))
  # rownames(tmat) = rep("", nrow(tmat)) 
  # print(tmat, quote = FALSE, print.gap = 2L, right = TRUE)
  
  tmat = cbind(
    c(round(x$nr_plus, 2), round(x$nr_minus, 2), tot_weights),
    c(round(x$sum_plus, 4), round(x$sum_minus, 4), tot_sums)
  )
  colnames(tmat) = c(paste0(treat, "s"), paste0("\U03A3", " weights"))
  rownames(tmat) = c("Positive weights", "Negative weights", "Total")
  cat(paste0("Treat. var: ", x$params$D), "\n")
  print(tmat, quote = FALSE, print.gap = 2L, right = TRUE)
  
  if (isTRUE(x$summary_measures)) {
    subscr = substr(x$type, 1, 2)
    cat("\nSummary Measures:\n")
    cat(sprintf("  TWFE Coefficient (\U03B2_%s): %.4f", subscr, x$beta), "\n")
    cat(sprintf("  min \U03C3(\U0394) compatible with \U03B2_%s and \U0394_TR = 0: %.4f", subscr, x$sensibility), "\n")
    if (x$sum_minus < 0) {
      cat(sprintf("  min \U03C3(\U0394) compatible with \U03B2_%s and \U0394_TR of a different sign: %.4f", subscr, x$sensibility2), "\n")
    }
  }

  if (!is.null(x$random_weights)) {
    cat("\nTest random weights:\n")
    print(x$mat)
  }
  
  cat("\n")
  cat("The development of this package was funded by the European Union (ERC, REALLYCREDIBLE,GA N. 101043899).", "\n")
  cat("\n")
  
  return(invisible(x))
}


