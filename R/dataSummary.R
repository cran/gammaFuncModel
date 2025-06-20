### Functions to summarize the data

### Sanity check function
sanityCheck <- function(df, covariates) {
  # Check if the required columns exist
  if (!all(c("ID", "Time") %in% colnames(df))) {
    stop("Error: Data must contain a column 'ID' and a column 'Time'.")
  }
  
  # Check if Time starts at 0 after sorting
  newdat <- df[order(df$Time), ]
  if (newdat$Time[1] == 0) {
    stop("Error: 'Time' should start with 1 instead of 0.")
  }
  
  covariate_cols = c(covariates, "Time")
  # Check for NA values
  if (any(is.na(df[, !(names(df) %in% covariate_cols)]))) {
    stop("Error: Concentration Data should not contain NAs.")
  }
}

### data summary
#' @importFrom dplyr %>% select all_of any_of group_by summarize arrange desc across
#' @importFrom stats median sd
summarizeData <- function(df, covariates) {
  data <- df %>% select(all_of(covariates))
  
  # Separate numeric and categorical columns
  numeric_cols <- names(data)[sapply(data, is.numeric)]
  categorical_cols <- names(data)[sapply(data, function(col) is.factor(col) || is.character(col))]
  
  # Summary for numeric columns 
  numeric_summary <- if(length(numeric_cols) > 0) {
    result <- sapply(data[numeric_cols], function(col) {
      c(
        mean = mean(col, na.rm = TRUE),
        median = median(col, na.rm = TRUE),
        sd = sd(col, na.rm = TRUE),
        num_na = sum(is.na(col)) # Count of NAs in each column
      )
    }) %>% t() %>% as.data.frame()

  } else {
    NULL
  }
  
  
  # Summary for categorical columns (frequency tables)
  categorical_summary <- if (length(categorical_cols) > 0) {
    lapply(data[categorical_cols], function(col) {
      table(col, useNA = "ifany")
    })
  } else {
    NULL
  }
  
  # Distribution of NAs across individuals
  na_distribution <- if ("ID" %in% names(df)) {
    df %>% 
      group_by(ID) %>% 
      summarize(num_na = sum(is.na(across(any_of(numeric_cols))))) %>% 
      arrange(desc(num_na)) %>% 
      as.data.frame()
  } else {
    warning("No 'ID' column found in the data for NA distribution analysis.")
    NULL
  }
  
  # Combine results into a list
  summary_list <- list(
    numeric_summary = numeric_summary,
    categorical_summary = categorical_summary,
    na_distribution = na_distribution
  )
  
  return(summary_list)
}

