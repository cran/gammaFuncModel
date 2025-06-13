## Helper functions
#' Geometric mean function
#' @keywords internal
#' @noRd
geometric_mean <- function(x, na.rm = TRUE) {
  # Remove NA values if na.rm is TRUE
  if (na.rm) {
    x <- x[!is.na(x)]
  }
  # Remove negative values
  x <- x[x >= 0]
  # Compute and return the geometric mean
  if (length(x) > 0) {
    exp(mean(log(x)))
  } else {
    NA  # Return NA if there are no valid values left
  }
}

#' Rounding functions
#' @keywords internal
#' @noRd
round_to_2_sig <- function(x) {
  signif(x, 2)
}

#' @keywords internal
#' @noRd
round_to_3_sig <- function(x) {
  signif(x, 3)
}

#' Anova function that takes care of NULL arguments
#' @importFrom stats anova 
#' @keywords internal
#' @noRd
safe_anova <- function(model1, model2) {
  tryCatch(
    {
      aov_result <- anova(model1, model2)
      return(aov_result$`p-value`[2])  # Return p-value of the second model
    },
    error = function(e) {
      return(NA)  # Return NA if anova fails
    }
  )
}

#' Function that produces a output table
#' @importFrom dplyr %>% rename mutate across where 
#' @importFrom stats p.adjust
#' @keywords internal
#' @noRd
output_table <- function(vecs, coefs, met_vec, REMLstatuses, MLstatuses) {
  # Combine the list of vectors into a data frame
  Result <- do.call(rbind, vecs)
  
  # Create the base column names
  colnames_base <- c(
    paste0(names(coefs), "Effect"), 
    paste0(names(coefs), "Pvalue"), 
    "LRTPvalue1"
  )
  
  # Add "LRTPvalue2" only if there is a second LRTPvalue column
  if (length(vecs[[1]]) > length(colnames_base)) {
    colnames_base <- c(colnames_base, "LRTPvalue2")
  }
  
  
  # Assign the correct column names
  colnames(Result) <- colnames_base #### statuses is a list, convert to vector
  
  
  # Convert Result to a data frame and reorder columns
  Result <- as.data.frame(Result)
  Result$Name <- met_vec
  Result <- Result[, c("Name", setdiff(names(Result), "Name"))]
  
  Result$REMLFitStatus <- REMLstatuses
  
  if(length(MLstatuses) != 0) {
    Result$MLFitStatus <- MLstatuses ### add in the raw p-values
  }
  
  # Adjust p-values and rename columns accordingly
  Result$BHPvalue1 <- p.adjust(Result$LRTPvalue1, method = "BH")
  
  if ("LRTPvalue2" %in% colnames(Result)) {
    Result$BHPvalue2 <- p.adjust(Result$LRTPvalue2, method = "BH")
    Result <- Result %>%
      rename("BHPvalue_total" = "BHPvalue1", "BHPvalue_int" = "BHPvalue2") %>%
      rename("LRTPvalue_total" = "LRTPvalue1","LRTPvalue_int" = "LRTPvalue2") %>%
      mutate(across(where(is.numeric), round_to_2_sig)) 
  } else {
    Result <- Result %>%
      rename("BHPvalue_int" = "BHPvalue1") %>%
      rename("LRTPvalue_int" = "LRTPvalue1") %>%
      mutate(across(where(is.numeric), round_to_2_sig)) 
  }
  
  return(Result)
}

#' Function to precompute f function
#' @importFrom dplyr %>% select distinct
#' @keywords internal
#' @noRd
precompute_f_list <- function(data, model,  grp_name = "Diet", ref = 1) {
  f_list <- list()  
  
  # Loop through each unique combination of Group and ID
  unique_combinations <- data %>% select(ID, !!sym(grp_name)) %>% distinct()
  
  for (i in seq_len(nrow(unique_combinations))) {
    grp_var <- unique_combinations[[grp_name]][i]
    id <- unique_combinations$ID[i]
    data_subset <- data %>% filter(!!sym(grp_name) == grp_var & ID == id)
    f_list[[paste0(id, "-", as.character(grp_var))]] <- generate_f_function(data_subset, model, grp_var,grp_name, id, ref)
  }
  
  return(f_list)
}