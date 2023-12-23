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
  weight_string = sprintf("%d %ss receive a positive weight, and %d receive a negative weight.", x$nr_plus, treat, x$nr_minus)
  
  # cat("\n")
  # # cat(rep("\u2500", hwidth), sep = "")
  # cat(cli::rule())
  # cat("\n")
  cat(assumption_string)
  cat("\n")
  cat(weight_string)
  # cat("\n")
  # # cat(rep("\u2500", hwidth), sep = "")
  # cat(cli::rule())
  cat("\n\n")
  
  weight_string = sprintf("%d %s receive a positive weight, and %d receive a negative weight.", x$nr_plus, treat, x$nr_minus)
  
  # header_string = capture.output(
  #   cat(assumption_string, "\n"),
  #   cat(weight_string)
  # )
  # cli::cat_boxx(header_string)
  # cat("\n\n")
  
  tmat = cbind(
    # c(round(x$nr_plus, 2), round(x$nr_minus, 2), tot_weights),
    c(x$nr_plus, x$nr_minus, tot_weights),
    c(round(x$sum_plus, 4), round(x$sum_minus, 4), tot_sums)
  )
  # colnames(tmat) = c(paste0(treat, "s"), paste0("\U03A3", " weights"))
  # rownames(tmat) = c("Positive weights", "Negative weights", "Total")
  # cat(paste0("Treat. var: ", x$params$D), "\n")
  # print(tmat, quote = FALSE, print.gap = 2L, right = TRUE)
  print_treat_matrix(tmat = tmat, tvar = x$params$D, ttype = treat)
  
  if (isTRUE(x$summary_measures)) {
    subscr = substr(x$type, 1, 2)
    # cat("\nSummary Measures:\n")
    cat(cli::style_bold("\nSummary Measures:\n"))
    cat(sprintf("  TWFE Coefficient (\U03B2_%s): %.4f", subscr, x$beta), "\n")
    cat(sprintf("  min \U03C3(\U0394) compatible with \U03B2_%s and \U0394_TR = 0: %.4f", subscr, x$sensibility), "\n")
    if (x$sum_minus < 0) {
      cat(sprintf("  min \U03C3(\U0394) compatible with \U03B2_%s and \U0394_TR of a different sign: %.4f", subscr, x$sensibility2), "\n")
    }
  }

  if (!is.null(x$random_weights)) {
    # cat("\nTest random weights:\n")
    cat(cli::style_bold("\nTest random assignment of weights:\n"))
    print(x$mat)
  }
  
  cat("\n\n")
  # cat("The development of this package was funded by the European Union (ERC, REALLYCREDIBLE,GA N. 101043899).", "\n")
  cat(cli::style_italic("The development of this package was funded by the European Union (ERC, REALLYCREDIBLE,GA N. 101043899).", "\n"))
  cat("\n")
  
  return(invisible(x))
}


## Custom print print method for treatment matrix
print_treat_matrix = function(tmat, tvar, ttype) {

  # Prep matrix (add row and column headers)
  tmat = cbind(
    c("Positive weights", "Negative weights", "Total"),
    tmat
  )
  tmat = rbind(
    c(
      paste0("Treat. var: ", tvar),
      paste0(ttype, "s"),
      paste0("\U03A3", " weights")
    ),
    tmat
  )

  # Calculate the width for each column
  col_widths = apply(tmat, 2, function(x) max(nchar(x)))
  
  # Function to format a single row based on column widths
  format_row = function(row) {
    formatted_elements = mapply(function(x, width, idx) {
      if (idx == 1) {
        sprintf("%-*s", width, x)  # Left-align the first column
      } else {
        sprintf("%*s", width, x)   # Right-align all other columns
      }
    }, row, col_widths, seq_along(row))
    paste(formatted_elements, collapse = "  ")
  }
  
  # Create the header and separator line
  header = format_row(tmat[1, ])
  bold_header = cli::style_bold(header)
  # separator_line = paste(rep("─", nchar(header)), collapse = "")
  # Calculate the separator line width
  separator_line_width = nchar(header)
  
  # Print separator line before header
  # cat(separator_line, "\n", sep = "")
  cat(cli::rule(width = separator_line_width), "\n")
  # Print the bold header
  cat(bold_header, "\n")
  # Print separator line after header
  # cat(separator_line, "\n", sep = "")
  cat(cli::rule(width = separator_line_width), "\n")
  
  # Print the body of the matrix except the last row
  for (i in 2:(nrow(tmat) - 1)) {
    cat(format_row(tmat[i, ]), "\n")
  }
  
  # Print separator line before last row
  # cat(separator_line, "\n", sep = "")
  cat(cli::rule(width = separator_line_width), "\n")
  # Print the last row
  cat(format_row(tmat[nrow(tmat), ]), "\n")
  # Print separator line after last row
  # cat(separator_line, "\n", sep = "")
  cat(cli::rule(width = separator_line_width), "\n")
  
}