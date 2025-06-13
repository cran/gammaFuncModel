#' Implementation of the novel non-linear mixed-effects model based on gamma function form
#' with nested covariance structure where random effects are specified for each Diet level within each subject (ID), capturing within-subject correlation across dietary conditions.
#' to identify metabolites that responds to time differentially across dietary groups
#' 
#' @param data Data frame that contains the 'ID' column along with all covariates as well as concentration column, named 'Concentration', for a single metabolite
#' @param covariates Vector containing the names of the "ID" covariate, grouping covariate and other covariates excluding any "Time" covariates 
#' @param time_terms Vector that contains all additional form of the covariate 'Time" (including the 'Time' covariate), and must contain 'log(Time)', other forms also include I(Time^2) and I(Time^3);
#' @param grp Grouping variable;
#' @param random_formula Random effects formula for the model, nested effects of Diet within ID (could also add random slope for 'Time');
#' @param correlation_formula Correlation formula. Default is autorrgressive but can accomodate other forms such as unstructured covariance or exponential covariance;
#' @param weights specify a variance function that models heteroscedasticity
#' @param time_grp_inter Boolean value that indicates if the model should include interactions terms of 'time_terms' with 'Group';
#' @param include_grp boolean value to indicate whether or not 'grp' should be included in the model construction
#' @param return_ml_model Boolean value that indicates if the model should fit "ML" model as well as "REML" model(default)
#' @return mixed effects models for a single metabolites: one with REML, the other with ML
#' 
#' @importFrom Rdpack reprompt
#' @importFrom nlme lme fixef corSymm varIdent lmeControl
#' @importFrom stats as.formula formula update
#' @references
#' Wickham, H. (2022). dplyr: A Grammar of Data Manipulation. R package version 1.0.10. 
#' Available at: \url{https://CRAN.R-project.org/package=dplyr}
#'
#' Pinheiro, J. C., & Bates, D. M. (2022). nlme: Linear and Nonlinear Mixed Effects Models. R package version 3.1-153.
#' Available at: \url{https://CRAN.R-project.org/package=nlme}
#' 
#' @examples
#' require(gammaFuncModel)
#' require(dplyr)
#' require(nlme)
#' df <- data.frame(
#'   ID = rep(sprintf("%02d", 1:10), each = 9 * 3),  
#'   Time = rep(rep(1:9, each = 3), 10),             
#'   Diet = as.factor(rep(1:3, times = 9 * 10)),    
#'   Age = rep(sample(20:70, 10, replace = TRUE), each = 9 * 3), 
#'   BMI = round(rep(runif(10, 18.5, 35), each = 9 * 3), 1),    
#'   Concentration = round(runif(270, 5, 15), 2)
#' )
#' covariates <- c("ID", "Diet", "Age", "BMI")
#' model <- gammaFunction(
#'   df, 
#'   covariates, 
#'   random_formula = ~ 1  | ID/Diet, 
#'   correlation_formula = corAR1(form = ~ Time | ID/Diet), 
#'   weights = NULL,
#'   include_grp = TRUE)[[1]]
#' summary(model)
#' 
#' @export
gammaFunction <- function(data, covariates, time_terms = c("Time", "log(Time)"), grp = "Diet", 
                          random_formula = ~ 1 + Time | ID/Diet, correlation_formula = corSymm(form = ~ Time | ID/Diet), 
                          weights = varIdent(form = ~ 1 | Time),
                          time_grp_inter = TRUE, return_ml_model = FALSE, include_grp) {
  
    data$Concentration[data$Concentration == 0] <- min(data$Concentration[data$Concentration != 0])/2
    covariates <- setdiff(covariates, "ID")
    
    # Adjust 'covariates' based on include_grp
    if (!include_grp) {
      covariates <- setdiff(covariates, grp)
      time_grp_inter <- FALSE
    }
    
    ### Add the Time terms into the covariates string
    covariates <- c(covariates, time_terms)
    
    # Create a formula string
    formula_string <- paste("log(Concentration) ~ ", paste(covariates, collapse = " + "))
    

    if (time_grp_inter) {
      interaction_terms <- paste(time_terms, "* ", grp, sep = "")
      formula_string <- paste(formula_string, "+ ", paste(interaction_terms, collapse = " + "), sep = "")
    }
    # Convert the formula string to a formula object
    form <- as.formula(formula_string)

    
    REMLfit_status <- "success"  # Default fit status, handle corner cases

    object1 <- tryCatch({
      model <- withCallingHandlers(
        nlme::lme(
          fixed = form,
          random = random_formula,
          correlation = correlation_formula,
          weights = weights,
          data = data,
          control = lmeControl(returnObject = TRUE)
        ),
        warning = function(w) {
          assign("fit_status", "warning", envir = parent.frame())
          invokeRestart("muffleWarning") 
        }
      )
      
      # Fix formula environment *after* fitting
      form <- formula(model)
      environment(form) <- .GlobalEnv
      model <- update(model, fixed = form)
      
      list(model = model, REMLfit_status = if (exists("fit_status")) "warning" else "success")
    }, error = function(e) {
      list(model = NULL, REMLfit_status = "failure")
    })
    


    modelML <- NULL

    MLfit_status <- "Not Fit"  # Default to failure unless successfully updated
    
    if (object1$REMLfit_status != "failure" && return_ml_model) {
      object2 <- tryCatch(
        {
          modelML <- withCallingHandlers(
            update(object1$model, method = "ML"),
            warning = function(w) {
              MLfit_status <<- "warning"  # Set status to warning if a warning occurs
              invokeRestart("muffleWarning")  # Suppress the warning
            }
          )
          # Fix formula environment after fitting
          form <- formula(model)
          environment(form) <- .GlobalEnv
          model <- update(model, fixed = form)
          
          if (MLfit_status != "warning") {
            MLfit_status <- "success"  # Set to success only if no warnings occurred
          }
          list(modelML = modelML, MLfit_status = MLfit_status)
        },
        error = function(e) {
          list(modelML = NULL, MLfit_status = "failure")  # Return failure status on error
        }
      )
      
      # Extract results
      modelML <- object2$modelML
      MLfit_status <- object2$MLfit_status
    }
    

    return(list(
      REML_model = object1$model,
      ML_model = modelML,
      REMLStatus = object1$REMLfit_status,
      MLStatus = MLfit_status
    ))

}



#' Function that produces a summary table for coefficient estimates, their p-values and LRT p-values for every metabolite in the dataframe
#'
#' @param df Data frame containing columns Group(numeric or character); ID(subject ID: character); Time(positive: numeric); other Time terms (numeric); other individual characteristics covariates; as well columns of metabolite concentrations
#' @param met_vec Vector of metabolite names
#' @param covariates Vector containing the names of the "ID" covariate, grouping covariate and other covariates excluding any "Time" covariates
#' @param time_terms Vector that contains all additional form of the covariate 'Time" (including the 'Time' covariate), and must contain 'log(Time)', other forms also include I(Time^2) and I(Time^3);
#' @param grp Grouping variable;
#' @param random_formula Random effects formula for the model, nested effects of Diet within ID (could also add random slope for 'Time');
#' @param correlation_formula Correlation formula. Default is autorgressive but can accommodate other forms such as unstructured covariance or exponential covariance;
#' @param weights specify a variance function that models heteroscedasticity
#' @return Data frame that contains the coefficient estimates, their corresponding p-values; LRT p-values for Time-Group interactions (for every 'Time' term);LRT p-values for Group and Time-Group interactions (for every 'Time' term); as well as the fitted models for each metabolite
#' @importFrom Rdpack reprompt
#' @importFrom dplyr select %>%
#' @importFrom nlme lme fixef corSymm varIdent
#' @references
#' Wickham, H. (2022). dplyr: A Grammar of Data Manipulation. R package version 1.0.10. 
#' Available at: \url{https://CRAN.R-project.org/package=dplyr}
#'
#' Pinheiro, J. C., & Bates, D. M. (2022). nlme: Linear and Nonlinear Mixed Effects Models. R package version 3.1-153.
#' Available at: \url{https://CRAN.R-project.org/package=nlme}
#' @examples
#' \donttest{
#' require(gammaFuncModel)
#' require(dplyr)
#' require(nlme)
#' df <- data.frame(
#'   ID = rep(sprintf("%02d", 1:10), each = 9 * 3), 
#'   Time = rep(rep(1:9, each = 3), 10),             
#'   Diet = as.factor(rep(1:3, times = 9 * 10)),     
#'   Age = rep(sample(20:70, 10, replace = TRUE), each = 9 * 3), 
#'   BMI = round(rep(runif(10, 18.5, 35), each = 9 * 3), 1)     
#' )
#' metvar <- paste0("met", 1:10)
#' concentration_data <- replicate(10, round(runif(270, 5, 15), 2))
#' colnames(concentration_data) <- metvar[1:10]
#' df <- cbind(df, as.data.frame(concentration_data))
#' covariates <- c("ID", "Diet", "Age", "BMI")
#' result <- diffGrpResponse(df, metvar, covariates)[[1]]
#' summary(result)
#' }
#' @export
diffGrpResponse <- function(df, met_vec, covariates, time_terms = c("Time", "log(Time)"), grp = "Diet", 
                         random_formula = ~ 1 + Time | ID/Diet, correlation_formula = corSymm(form = ~ Time | ID/Diet), weights = varIdent(form = ~ 1 | Time)) {
  models <- list()
  vecs <- list()
  REMLstatuses <- numeric()
  MLstatuses <- numeric()
  sanityCheck(df)
  data <- df %>% select(all_of(covariates), "Time")
  covariates <- setdiff(covariates, "ID")
  
  for (i in 1: length(met_vec)) {
    metabolite <- met_vec[i]
    data$Concentration <- df[, metabolite]
    
    models_no_inter_grp <- gammaFunction(data, setdiff(covariates, time_terms), time_terms, grp,  random_formula, correlation_formula, weights, time_grp_inter = FALSE, return_ml_model = TRUE, include_grp = TRUE)
    models_no_inter_no_grp <- gammaFunction(data, setdiff(covariates, time_terms), time_terms, grp,  random_formula, correlation_formula, weights, time_grp_inter = FALSE, return_ml_model = TRUE, include_grp = FALSE)
    models_inter <- gammaFunction(data, setdiff(covariates, time_terms), time_terms, grp, random_formula, correlation_formula, weights, time_grp_inter = TRUE, return_ml_model = TRUE, include_grp = TRUE)


    model1 <- models_inter[[1]]

    ### Stores the model with interaction terms for each metabolites in a list
    models[[metabolite]] <- model1


    model0_grp_ML <- models_no_inter_grp[[2]]
    model0_no_grp_ML <- models_no_inter_no_grp[[2]]
    model1ML <- models_inter[[2]]

    ### LRT for total group effect
    lrt_pvalue1 <- safe_anova(model0_no_grp_ML, model1ML)
    
    ### LRT on for group-time interaction effect
    lrt_pvalue2 <- safe_anova(model0_grp_ML, model1ML)

    if (!is.null(model1)) {
      coefs <- fixef(model1)[-1]
      pvalues <- summary(model1)$tTable[, "p-value"][-1]
    } else {
      coefs <- NA
      pvalues <- NA
    }
    vec <- as.numeric(c(coefs, pvalues, lrt_pvalue1, lrt_pvalue2))
    vecs[[i]] <- vec
    REMLstatuses[i] <- models_inter[[3]]
    MLstatuses[[i]] <- models_inter[[4]]
  }


  Result = output_table(vecs, coefs, met_vec, REMLstatuses, MLstatuses)

  return(list(Result, models))
}


#' Parallelized version of diffGrpResponse()
#'
#' @inheritParams diffGrpResponse
#' @return Data frame that contains the coefficient estimates, their corresponding p-values; LRT p-values for Time-Group interactions (for every 'Time' term);LRT p-values for Group and Time-Group interactions (for every 'Time' term); as well as the fitted models for each metabolite
#' @importFrom Rdpack reprompt
#' @importFrom dplyr select %>%
#' @importFrom nlme lme fixef corSymm varIdent
#' @importFrom future.apply future_lapply
#' @importFrom stats setNames
#' @references
#' Wickham, H. (2022). dplyr: A Grammar of Data Manipulation. R package version 1.0.10. 
#' Available at: \url{https://CRAN.R-project.org/package=dplyr}
#'
#' Pinheiro, J. C., & Bates, D. M. (2022). nlme: Linear and Nonlinear Mixed Effects Models. R package version 3.1-153.
#' Available at: \url{https://CRAN.R-project.org/package=nlme}
#' @note This function uses parallel processing via the `future.apply` package.
#' To enable parallel execution, runs the following before calling this function:
#' 
#'   library(future.apply)
#'   plan(multisession, workers = parallel::detectCores() - 1)
#' 
#' You only need to set the plan once per session.
#' @examples
#' \donttest{
#' require(gammaFuncModel)
#' require(dplyr)
#' require(nlme)
#' df <- data.frame(
#'   ID = rep(sprintf("%02d", 1:10), each = 9 * 3),  
#'   Time = rep(rep(1:9, each = 3), 10),             
#'   Diet = as.factor(rep(1:3, times = 9 * 10)),     
#'   Age = rep(sample(20:70, 10, replace = TRUE), each = 9 * 3), 
#'   BMI = round(rep(runif(10, 18.5, 35), each = 9 * 3), 1)     
#' )
#' metvar <- paste0("met", 1:10)
#' concentration_data <- replicate(10, round(runif(270, 5, 15), 2))
#' colnames(concentration_data) <- metvar[1:10]
#' df <- cbind(df, as.data.frame(concentration_data))
#' covariates <- c("ID", "Diet", "Age", "BMI")
#' library(future.apply)
#' plan(multisession, workers = parallel::detectCores() - 1)
#' result <- diffGrpResponse(df, metvar, covariates)[[1]]
#' summary(result)
#' }
#' @export
diffGrpResponse_parallel <- function(df, met_vec, covariates, 
                              time_terms = c("Time", "log(Time)"), grp = "Diet", 
                              random_formula = ~ 1 + Time | ID/Diet, 
                              correlation_formula = corSymm(form = ~ Time | ID/Diet), 
                              weights = varIdent(form = ~ 1 | Time)) {
  
  sanityCheck(df)  # Still run the sanity check
  
  base_data <- df %>% select(all_of(covariates), "Time")
  covariates <- setdiff(covariates, "ID")
  
  # Define per-metabolite modeling as a function
  run_one_met <- function(metabolite) {
    data <- base_data
    data$Concentration <- df[[metabolite]]
    
    models_no_inter_grp <- gammaFunction(data, setdiff(covariates, time_terms), time_terms, grp,  
                                           random_formula, correlation_formula, weights, 
                                           time_grp_inter = FALSE, return_ml_model = TRUE, include_grp = TRUE)
    
    models_no_inter_no_grp <- gammaFunction(data, setdiff(covariates, time_terms), time_terms, grp,  
                                              random_formula, correlation_formula, weights, 
                                              time_grp_inter = FALSE, return_ml_model = TRUE, include_grp = FALSE)
    
    models_inter <- gammaFunction(data, setdiff(covariates, time_terms), time_terms, grp, 
                                    random_formula, correlation_formula, weights, 
                                    time_grp_inter = TRUE, return_ml_model = TRUE, include_grp = TRUE)
    
    model1 <- models_inter[[1]]
    model0_grp_ML <- models_no_inter_grp[[2]]
    model0_no_grp_ML <- models_no_inter_no_grp[[2]]
    model1ML <- models_inter[[2]]
    
    # LRTs
    lrt_pvalue1 <- safe_anova(model0_no_grp_ML, model1ML)
    lrt_pvalue2 <- safe_anova(model0_grp_ML, model1ML)
    
    if (!is.null(model1)) {
      coefs <- fixef(model1)[-1]
      pvalues <- summary(model1)$tTable[, "p-value"][-1]
    } else {
      coefs <- NA
      pvalues <- NA
    }
    
    vec <- as.numeric(c(coefs, pvalues, lrt_pvalue1, lrt_pvalue2))
    
    return(list(
      vec = vec,
      model = model1,
      REMLstatus = models_inter[[3]],
      MLstatus = models_inter[[4]]
    ))
  }
  
  # Run in parallel
  res_list <- future_lapply(met_vec, run_one_met)
  
  # Extract outputs
  vecs <- lapply(res_list, `[[`, "vec")
  models <- setNames(lapply(res_list, `[[`, "model"), met_vec)
  REMLstatuses <- sapply(res_list, `[[`, "REMLstatus")
  MLstatuses <- sapply(res_list, `[[`, "MLstatus")
  
  # Use the last metabolite's coefs as template (you could also get from any successful model)
  last_coefs <- tryCatch(fixef(res_list[[length(res_list)]]$model)[-1], error = function(e) rep(NA, length(res_list[[1]]$vec)/2 - 2))
  
  # Assemble final output
  Result <- output_table(vecs, last_coefs, met_vec, REMLstatuses, MLstatuses)
  
  return(list(Result, models))
}




#' Function that produces a summary table for coefficient estimates, their p-values and LRT p-values for every metabolite in the dataframe, for a single Group
#' @param df Data frame containing information for a single group, containing columns grp; ID(subject ID: character); Time(positive: numeric); other Time terms (numeric); other individual characteristics covariates; as well columns of metabolite concentrations
#' @param met_vec Vector of metabolite names
#' @param covariates Vector containing the names of the "ID" covariate, grouping covariate and other covariates excluding any "Time" covariates
#' @param time_terms Vector that contains all additional form of the covariate 'Time" (including the 'Time' covariate), and must contain 'log(Time)', other forms also include I(Time^2) and I(Time^3);
#' @param grp Grouping variable (should be a single valued column);
#' @param random_formula Random effects formula for the model, within ID (could also add random slope for 'Time');
#' @param correlation_formula Correlation formula. Default is autorgressive but can accommodate other forms such as unstructured covariance or exponential covariance;
#' @param weights specify a variance function that models heteroscedasticity
#' @return Data frame that contains the coefficient estimates, their corresponding p-values as well as LRT p-values for 'Time' terms
#' @importFrom Rdpack reprompt
#' @importFrom dplyr filter %>%
#' @importFrom nlme lme corAR1 fixef
#' @references
#' Wickham, H. (2022). dplyr: A Grammar of Data Manipulation. R package version 1.0.10. 
#' Available at: \url{https://CRAN.R-project.org/package=dplyr}
#'
#' Pinheiro, J. C., & Bates, D. M. (2022). nlme: Linear and Nonlinear Mixed Effects Models. R package version 3.1-153.
#' Available at: \url{https://CRAN.R-project.org/package=nlme}
#' @examples
#' \donttest{
#' require(gammaFuncModel)
#' require(dplyr)
#' require(nlme)
#' df <- data.frame(
#'   ID = rep(sprintf("%02d", 1:10), each = 9 * 3),  
#'   Time = rep(rep(1:9, each = 3), 10),             
#'   Diet = as.factor(rep(1:3, times = 9 * 10)),     
#'   Age = rep(sample(20:70, 10, replace = TRUE), each = 9 * 3), 
#'   BMI = round(rep(runif(10, 18.5, 35), each = 9 * 3), 1)    
#' )
#' metvar <- paste0("met", 1:10)
#' concentration_data <- replicate(10, round(runif(270, 5, 15), 2))
#' colnames(concentration_data) <- metvar[1:10]
#' df <- cbind(df, as.data.frame(concentration_data))
#' df_single_diet <- subset(df, Diet == 1)
#' covariates <- c("ID","Diet", "Age", "BMI")
#' result_SD <- grpResp2Time(df_single_diet, metvar, covariates)[[1]]
#' summary(result_SD)
#' }
#' @export
grpResp2Time <- function(df, met_vec, covariates, time_terms = c("Time", "log(Time)"), grp = 'Diet',
                         random_formula = ~ 1 | ID, correlation_formula = corAR1(form = ~ Time | ID), weights = NULL) {
  vecs <- list()
  data <- df %>% select(all_of(covariates), "Time")
  covariates <- setdiff(covariates, "ID")
  met_models = list()
  REMLstatuses <- numeric()
  MLstatuses <- numeric()
  
  for (i in 1: length(met_vec) ) {
    metabolite <- met_vec[i]
    data$Concentration <- df[, metabolite]

    models_no_time <- gammaFunction(data, setdiff(covariates, time_terms), time_terms = NULL, grp = grp,  random_formula, correlation_formula, weights, time_grp_inter = FALSE, return_ml_model = TRUE, include_grp = FALSE)
    models <- gammaFunction(data, setdiff(covariates, time_terms), time_terms, grp = grp, random_formula, correlation_formula, weights, time_grp_inter = FALSE, return_ml_model = TRUE, include_grp = FALSE)
    
    model1 <- models[[1]]
    met_models[[metabolite]] <- model1
    
    ### conduct LRT on for Time and log(Time) interaction terms with Group
    model0ML <- models_no_time[[2]]
    model1ML <- models[[2]]
    
    lrt_pvalue <- safe_anova(model0ML, model1ML)
   

    if (!is.null(model1)) {
      coefs <- fixef(model1)[-1]
      pvalues <- summary(model1)$tTable[, "p-value"][-1]
    } else {
      coefs <- NA
      pvalues <- NA
    }
    vec <- as.numeric(c(coefs, pvalues, lrt_pvalue))
    vecs[[i]] <- vec
    REMLstatuses[i] <- models[[3]]
    MLstatuses[i] <- models[[4]]

  }

  Result <- output_table(vecs, coefs, met_vec, REMLstatuses, MLstatuses)
  
  return(list(Result, met_models))
}


#'  Vectorized version of grpRes2Time()
#' 
#' @inheritParams grpResp2Time
#' @return Data frame that contains the coefficient estimates, their corresponding p-values as well as LRT p-values for 'Time' terms
#' @importFrom Rdpack reprompt
#' @importFrom dplyr filter %>%
#' @importFrom nlme lme corAR1 fixef
#' @importFrom stats setNames
#' @references
#' Wickham, H. (2022). dplyr: A Grammar of Data Manipulation. R package version 1.0.10. 
#' Available at: \url{https://CRAN.R-project.org/package=dplyr}
#'
#' Pinheiro, J. C., & Bates, D. M. (2022). nlme: Linear and Nonlinear Mixed Effects Models. R package version 3.1-153.
#' Available at: \url{https://CRAN.R-project.org/package=nlme}
#' @note This function uses parallel processing via the `future.apply` package.
#' To enable parallel execution, runs the following before calling this function:
#' 
#'   library(future.apply)
#'   plan(multisession, workers = parallel::detectCores() - 1)
#' 
#' You only need to set the plan once per session.
#' @examples
#' \donttest{
#' require(gammaFuncModel)
#' require(dplyr)
#' require(nlme)
#' df <- data.frame(
#'   ID = rep(sprintf("%02d", 1:10), each = 9 * 3),  
#'   Time = rep(rep(1:9, each = 3), 10),            
#'   Diet = as.factor(rep(1:3, times = 9 * 10)),    
#'   Age = rep(sample(20:70, 10, replace = TRUE), each = 9 * 3), 
#'   BMI = round(rep(runif(10, 18.5, 35), each = 9 * 3), 1)    
#' )
#' metvar <- paste0("met", 1:10)
#' concentration_data <- replicate(10, round(runif(270, 5, 15), 2))
#' colnames(concentration_data) <- metvar[1:10]
#' df <- cbind(df, as.data.frame(concentration_data))
#' df_single_diet <- subset(df, Diet == 1)
#' covariates <- c("ID","Diet", "Age", "BMI")
#' library(future.apply)
#' plan(multisession, workers = parallel::detectCores() - 1)
#' result_SD <- grpResp2Time_parallel(df_single_diet, metvar, covariates)[[1]]
#' summary(result_SD)
#' }
#' @export
grpResp2Time_parallel <- function(df, met_vec, covariates, time_terms = c("Time", "log(Time)"), grp = 'Diet',
                                  random_formula = ~ 1 | ID, correlation_formula = corAR1(form = ~ Time | ID), weights = NULL) {
  covariates <- setdiff(covariates, "ID")
  data_base <- df %>% select(all_of(covariates), "Time")
  
  results <- future.apply::future_lapply(met_vec, function(metabolite) {
    data <- data_base
    data$Concentration <- df[[metabolite]]
    
    models_no_time <- gammaFunction(data, setdiff(covariates, time_terms), time_terms = NULL, grp = grp,
                                    random_formula, correlation_formula, weights, time_grp_inter = FALSE,
                                    return_ml_model = TRUE, include_grp = FALSE)
    
    models <- gammaFunction(data, setdiff(covariates, time_terms), time_terms, grp = grp,
                            random_formula, correlation_formula, weights, time_grp_inter = FALSE,
                            return_ml_model = TRUE, include_grp = FALSE)
    
    model1 <- models[[1]]
    model0ML <- models_no_time[[2]]
    model1ML <- models[[2]]
    
    lrt_pvalue <- safe_anova(model0ML, model1ML)
    
    if (!is.null(model1)) {
      coefs <- fixef(model1)[-1]
      pvalues <- summary(model1)$tTable[, "p-value"][-1]
    } else {
      coefs <- NA
      pvalues <- NA
    }
    
    vec <- as.numeric(c(coefs, pvalues, lrt_pvalue))
    
    list(vec = vec, model = model1, REML = models[[3]], ML = models[[4]])
  })
  
  vecs <- lapply(results, function(x) x$vec)
  met_models <- setNames(lapply(results, function(x) x$model), met_vec)
  REMLstatuses <- sapply(results, function(x) x$REML)
  MLstatuses <- sapply(results, function(x) x$ML)
  
  coefs <- if (length(vecs) > 0 && !is.null(vecs[[1]])) vecs[[1]][1:(length(vecs[[1]]) - 3)] else NA
  
  Result <- output_table(vecs, coefs, met_vec, REMLstatuses, MLstatuses)
  
  return(list(Result, met_models))
}



#' Function that produces a fitted gamma model for each metabolite
#'
#' @param df Data frame containing columns Group(factor); ID(subject ID: character); Time(positive: numeric); other Time terms (numeric); other individidual characteristics covariates; as well columns of metabolite concentrations
#' @param met_vec the vector of metabolite names
#' @param covariates Vector containing the names of the "ID" covariate, grouping covariate and other covariates excluding any "Time" covariates
#' @param time_terms is the vector that contains all additional form of the covariate 'Time" (including the 'Time' covariate), and must contain 'log(Time)', other forms also include I(Time^2) and I(Time^3);
#' @param grp_name is the grouping variable;
#' @param random_formula is the random effects formula for the model, nested effects of Diet within ID (could also add random slope for 'Time');
#' @param correlation_formula is the correlation formula. Default is autorgressive but can accommodate other forms such as unstructured covariance or exponential covariance;
#' @param weights specify a variance function that models heteroscedasticity;
#' @param graph character string, 'None' by default. If not 'None, in addition to returning models, produces pdf file of graphs based on the specific value of 'graph'. 
#' @param save_path location where the pdf file will be saved; default is NULL, i.e. pdf is saved to a temporary location
#' @return List that contains fitted models for each metabolite and a pdf file for fitted concenration curves.
#' @importFrom Rdpack reprompt
#' @importFrom dplyr select %>% all_of
#' @importFrom nlme lme corSymm varIdent
#' @references
#' Wickham, H. (2022). dplyr: A Grammar of Data Manipulation. R package version 1.0.10. 
#' Available at: \url{https://CRAN.R-project.org/package=dplyr}
#'
#' Pinheiro, J. C., & Bates, D. M. (2022). nlme: Linear and Nonlinear Mixed Effects Models. R package version 3.1-153.
#' Available at: \url{https://CRAN.R-project.org/package=nlme}
#' 
#' @examples
#' \donttest{
#' require(gammaFuncModel)
#' require(dplyr)
#' require(nlme)
#'  df <- data.frame(
#'   ID = rep(sprintf("%02d", 1:10), each = 9 * 3),  
#'   Time = rep(rep(1:9, each = 3), 10),             
#'   Diet = as.factor(rep(1:3, times = 9 * 10)),     
#'   Age = rep(sample(20:70, 10, replace = TRUE), each = 9 * 3), 
#'   BMI = round(rep(runif(10, 18.5, 35), each = 9 * 3), 1)     
#' )
#' metvar <- paste0("met", 1:10)
#' concentration_data <- replicate(10, round(runif(270, 5, 15), 2))
#' colnames(concentration_data) <- metvar[1:10]
#' df <- cbind(df, as.data.frame(concentration_data))
#' covariates <- c("ID", "Diet", "Age", "BMI")
#' mods <- generate_models(
#'   df = df, 
#'  met_vec = metvar, 
#'   covariates = covariates, 
#'   graph = 'None', 
#'   save_path = NULL)
#' }
#' @export
generate_models <- function(df, met_vec, covariates, time_terms = c("Time", "log(Time)"), grp_name = "Diet", 
                         random_formula = ~ 1 + Time| ID/Diet, correlation_formula = corSymm(form = ~ Time | ID/Diet), weights = varIdent(form = ~ 1 | Time),
                         graph = "None", save_path = NULL) {
  met_models <- list()
  vecs <- list()

  data <- df %>% select(all_of(c(covariates, "Time")))
  covariates <- setdiff(covariates, "ID")
  
  
  ### check if this data frame is of multiple groups or a single group
  unique_grouping <- FALSE
  
  if (all(df[[grp_name]] == df[[grp_name]][1])) {
    unique_grouping <- TRUE
  }

  
  for (i in 1: length(met_vec)) {
    metabolite <- met_vec[i]
    data$Concentration <- df[, metabolite]
    

    models_inter <- gammaFunction(data, setdiff(covariates, time_terms), time_terms, grp =  grp_name, random_formula, correlation_formula, weights, time_grp_inter = ifelse(unique_grouping, FALSE, TRUE), return_ml_model = FALSE, include_grp = ifelse(unique_grouping, FALSE, TRUE))

    
    model1 <- models_inter[[1]]
    
    ### Stores the model with interaction terms for each metabolites in a list
    met_models[[metabolite]] <- model1
  }

  
  if (graph != 'None') {  

    generatePlot(graph, df, met_vec, c(covariates, "ID"), grp_name, met_models, save_path)
    
    message("Graph Done")
  }
  
  return(met_models)
}


#' Function that produces Tmax property for a single individual in a particular group, for a specific metabolite
#'
#' @param data Data frame containing columns Group(factor); ID(subject ID: character); Time(positive: numeric); other individiual characteristics covariates(excluding other forms of 'Time')
#' @param model Fitted model for the metabolite in question
#' @param grp_var Value of the grouping variable
#' @param ID Subject ID
#' @param grp_name Name of the grouping variable. Default is 'Diet'
#' @param ref numeric or character. The reference level for the grouping variable, as a factor
#' @return Tmax for this metabolite, in a particular group for a single individual
#' @importFrom Rdpack reprompt
#' @importFrom dplyr select
#' @importFrom nlme lme fixef ranef
#' @importFrom rootSolve uniroot.all
#' @importFrom stats predict
#' @references
#' Wickham, H. (2022). dplyr: A Grammar of Data Manipulation. R package version 1.0.10. 
#' Available at: \url{https://CRAN.R-project.org/package=dplyr}
#'
#' Pinheiro, J. C., & Bates, D. M. (2022). nlme: Linear and Nonlinear Mixed Effects Models. R package version 3.1-153.
#' Available at: \url{https://CRAN.R-project.org/package=nlme}
#' 
#' @examples
#' \donttest{
#' require(gammaFuncModel)
#' require(rootSolve)
#' require(dplyr)
#' require(nlme)
#' df <- data.frame(
#'   ID = rep(sprintf("%02d", 1:10), each = 9 * 3),  
#'   Time = rep(rep(1:9, each = 3), 10),             
#'   Diet = as.factor(rep(1:3, times = 9 * 10)),    
#'   Age = rep(sample(20:70, 10, replace = TRUE), each = 9 * 3), 
#'   BMI = round(rep(runif(10, 18.5, 35), each = 9 * 3), 1),     
#'   Concentration = round(runif(270, 5, 15), 2)
#' )
#' covariates <- c("ID", "Diet", "Age", "BMI")
#' model <- gammaFunction(
#'   df, 
#'   covariates, 
#'   time_grp_inter = FALSE, 
#'   return_ml_model = FALSE, 
#'   include_grp = TRUE
#' )[[1]]
#' test_data = df %>% filter(Diet == 1 & ID == "01") %>% select(-c("Concentration", "ID", "Diet")) 
#' Tmax <- calculate_Tmax(data = test_data, model, grp_var = 1, ID = "01", grp_name = 'Diet', ref = 1)
#' }
#' @export
calculate_Tmax <- function(data, model, grp_var, ID, grp_name = "Diet", ref = 1) {
  fixed_coefs <- fixef(model)
  random_coefs <- ranef(model)
  
  grp <- as.numeric(unique(grp_var))

  d_formula <- function(t) {
    fixed_coefs["log(Time)"] / t + fixed_coefs["Time"]  +
      ifelse(rep(grp, length(t)) != rep(ref, length(t)),
             fixed_coefs[paste0(grp_name, as.character(grp), ":log(Time)")], 0) / t +
      ifelse(rep(grp, length(t)) != rep(ref, length(t)),
             fixed_coefs[paste0(grp_name, as.character(grp), ":Time")], 0) +
      ifelse(rep(!is.null(random_coefs[["ID"]][ID,"Time"]), length(t)),
             random_coefs[["ID"]][ID,"Time"] + random_coefs[[grp_name]][paste0(ID, "/", grp), "Time"], 0)
  }
  
  
  #### Root within boundaries
  root1 <- uniroot.all(d_formula, interval = c(1, 50))
  
  # Add boundary points to the candidates
  Tmax_candidates <- c(root1, min(data$Time), max(data$Time))
  len <- length(Tmax_candidates)
  
  new_data <- data.frame(matrix(ncol = length(colnames(data)), nrow = len))
  colnames(new_data) <- colnames(data)
  
  for(var in setdiff(colnames(data), "Time")) {
    new_data[[var]] <- rep(unique(data[[var]]), len)
  }
  
  new_data[["ID"]] <- rep(unique(ID), len)
  new_data[[grp_name]] <- as.factor(rep(grp, len))
  new_data[["Time"]] <- Tmax_candidates
  
  PredConc_at_Tmax <- predict(model, newdata =  new_data)
  Tmax <- Tmax_candidates[which.max(PredConc_at_Tmax)]
  return(Tmax)
}


#' Function that produces Cmax property for a single individual in a particular group, for a specific metabolite
#'
#' @param data Data frame containing columns Group(factor); ID(subject ID: character); Time(positive: numeric); other individiual characteristics covariates (exlcluding other forms of 'Time')
#' @param model Fitted model for the metabolite in question
#' @param grp_var Value of the grouping variable
#' @param ID Subject ID
#' @param grp_name Name of the grouping variable. Default is 'Diet'
#' @param Tmax for this metabolite, in a particular group for a single individual
#' @return Cmax for this metabolite, in a particular group for a single individual
#' @importFrom Rdpack reprompt
#' @importFrom stats predict setNames
#' @references
#' Wickham, H. (2022). dplyr: A Grammar of Data Manipulation. R package version 1.0.10. 
#' Available at: \url{https://CRAN.R-project.org/package=dplyr}
#'
#' Pinheiro, J. C., & Bates, D. M. (2022). nlme: Linear and Nonlinear Mixed Effects Models. R package version 3.1-153.
#' Available at: \url{https://CRAN.R-project.org/package=nlme}
#' 
#' @examples
#' \donttest{
#' require(gammaFuncModel)
#' require(rootSolve)
#' require(dplyr)
#' require(nlme)
#' df <- data.frame(
#'   ID = rep(sprintf("%02d", 1:10), each = 9 * 3),  
#'   Time = rep(rep(1:9, each = 3), 10),             
#'   Diet = as.factor(rep(1:3, times = 9 * 10)),     
#'   Age = rep(sample(20:70, 10, replace = TRUE), each = 9 * 3),
#'   BMI = round(rep(runif(10, 18.5, 35), each = 9 * 3), 1),     
#'   Concentration = round(runif(270, 5, 15), 2)
#' )
#' covariates <- c("ID", "Diet", "Age", "BMI")
#' model <- gammaFunction(
#'   df, 
#'   covariates, 
#'   time_grp_inter = FALSE, 
#'   return_ml_model = FALSE, 
#'   include_grp = TRUE
#'   )[[1]]
#' test_data = df %>% filter(Diet == 1 & ID == "02") %>% select(-c("Concentration", "ID", "Diet")) 
#' Tmax <- calculate_Tmax(data = test_data, model, grp_var = 1, ID = "02", grp_name = 'Diet', ref = 1)
#' Cmax <- calculate_Cmax(data = test_data, model, grp_var = 1, ID = "02", grp_name = "Diet", Tmax)
#' }
#' @export
calculate_Cmax <- function(data, model, grp_var, ID, grp_name = "Diet",  Tmax) {
  grp <- as.numeric(unique(grp_var))
  
  names <- setdiff(colnames(data), c("ID", "Time", grp_name))
  
  new_data <- data.frame(
    setNames(
      lapply(names, function(var) {
        if (length(unique(data[[var]])) != 1) {
          stop(paste("Covariate", var, "is not unique within the group."))
        }
        unique(data[[var]])
      }),
      names
    )
  )
  new_data[["Time"]] <- Tmax
  new_data[["ID"]] <- unique(ID)
  new_data[[grp_name]] <- as.factor(grp)
  
  Cmax <- exp(predict(model,  new_data))
  
  return(Cmax)
}



#'Function produce predictions from the model 
#'
#' @param data Data frame containing columns Group(factor); ID(subject ID: character); Time(positive: numeric); other individiual characteristics covariates (exlcluding other forms of 'Time')
#' @param model Fitted model for the metabolite in question
#' @param grp_var Value of the grouping variable
#' @param grp_name Name of the grouping variable. Default is 'Diet'
#' @param ID Subject ID
#' @param ref reference group
#' @return f function that produces the prediction from this model for a specific individual in a specific group
#' @importFrom nlme fixef ranef 
#' @examples
#' \donttest{
#' require(gammaFuncModel)
#' require(dplyr)
#' require(nlme)
#' modify.df <- data.frame(
#'   ID = rep(sprintf("%02d", 1:10), each = 9 * 3),
#'   Time = rep(rep(1:9, each = 3), 10),
#'   Diet = as.factor(rep(1:3, times = 9 * 10)),
#'   Age = rep(sample(20:70, 10, replace = TRUE), each = 9 * 3),
#'   BMI = round(rep(runif(10, 18.5, 35), each = 9 * 3), 1),
#'   Concentration = NA
#' )
#' for (i in 1:10) {
#'  for (d in 1:3) {
#'    C0 <- runif(1, 10, 15)    # initial concentration
#'    k <- runif(1, 0.1, 0.3)   # decay rate constant
#'    modify.df$Concentration[modify.df$ID == sprintf("%02d", i) & modify.df$Diet == d] <- 
#'      C0 * exp(-k * modify.df$Time[modify.df$ID == sprintf("%02d", i) & modify.df$Diet == d])
#'  }
#'}
#' covariates <- c("ID", "Diet", "Age", "BMI")
#' model <- gammaFunction(
#'   modify.df, 
#'   covariates, 
#'   time_grp_inter = FALSE, 
#'   return_ml_model = FALSE, 
#'   include_grp = TRUE
#'  )[[1]]
#' test_data = modify.df %>% 
#'   filter(Diet == 1 & ID == "04") %>% 
#'   select(-c("Concentration", "ID", "Diet")) 
#' f_dat = modify.df %>% filter(Diet == 1 & ID == "04") %>% select(-Concentration)
#' f <- generate_f_function(
#'   data = f_dat, 
#'   model = model,  
#'   grp_var = 1, 
#'   grp_name = "Diet", 
#'   ID = "04", 
#'   ref = 1
#'  )
#' }
#' @export
generate_f_function <- function(data, model, grp_var, grp_name = "Diet", ID, ref = 1) {
  # Extract the unique group and ID
  grp <- as.numeric(unique(grp_var))
  ID <- unique(ID)
  
  # Extract fixed and random effects
  fixed_effects <- fixef(model)
  random_effects <- ranef(model)
  
  # Dynamically identify covariates (excluding Time, ID, and grp_name)
  covariate_names <- colnames(data)[
    !grepl("Time", colnames(data)) & !colnames(data) %in% c("ID", grp_name)
  ]
  
  
  # Dynamically extract covariate values
  covariate_values <- sapply(covariate_names, function(var) as.numeric(unique(data[[var]])), simplify = FALSE)
  
  
  f <- function(t) {
    len <- length(t)
    
    # Start with fixed effects intercept and time terms
    predictions <- fixed_effects["(Intercept)"] +
      fixed_effects["log(Time)"] * log(t) +
      fixed_effects["Time"] * t 
    
    
    # Add covariate contributions
    for (var in covariate_names) {
      predictions <- predictions + as.numeric(covariate_values[[var]]) * as.numeric(fixed_effects[var])
    }
    
    # Add group-specific fixed effects
    predictions <- predictions +
      ifelse(rep(grp, len) != rep(ref, len), fixed_effects[paste0(grp_name, as.character(grp))], 0) +
      ifelse(rep(grp, len) != rep(ref, len), fixed_effects[paste0(grp_name, as.character(grp), ":log(Time)")] * log(t), 0) +
      ifelse(paste0(grp_name, as.character(grp), ":Time") %in% names(fixed_effects),
             ifelse(rep(grp, len) != rep(ref, len), fixed_effects[paste0(grp_name, as.character(grp), ":Time")] * t, 0),
             0
      )
    
    
    # Add random effects for ID and group
    predictions <- predictions +
      random_effects$ID[ID, "(Intercept)"] +
      random_effects[[grp_name]][paste0(ID, "/", grp), "(Intercept)"]
    
    
    # Add random time effects safely
    re_time_ID <- if ("Time" %in% colnames(random_effects$ID)) random_effects$ID[ID, "Time"] else 0
    re_time_grp <- if ("Time" %in% colnames(random_effects[[grp_name]])) {
      random_effects[[grp_name]][paste0(ID, "/", grp), "Time"]
    } else 0
    
    
    ## Accounting for Time^2 and Time^3
    if("I(Time^2)" %in% names(fixed_effects)) {
      predictions <- predictions +
        fixed_effects["I(Time^2)"] * t^2 +
        ifelse(paste0("I(Time^2):", grp_name, as.character(grp)) %in% names(fixed_effects),
               ifelse(rep(grp, len) != rep(ref, len), fixed_effects[paste0("I(Time^2):", grp_name, as.character(grp))] * t^2, 0),
               0
        )
    }
    
    if("I(Time^3)" %in% names(fixed_effects)) {
      predictions <- predictions +
        fixed_effects["I(Time^3)"] * t^3 +
        ifelse(paste0("I(Time^3):", grp_name, as.character(grp)) %in% names(fixed_effects),
               ifelse(rep(grp, len) != rep(ref, len), fixed_effects[paste0("I(Time^3):", grp_name, as.character(grp))] * t^3, 0),
               0
        )
    }
    
    return(exp(predictions))
  }
  
  
  return(f)
}




#' Function that produces Half-life property for a single individual in a particular group, for a specific metabolite
#'
#' @param f function that returns the prediction of a metabolite concentration, for a single individual in a particular group
#' @param Tmax Tmax property of a metabolite, for a single individual in a particular group
#' @param Cmax Cmax property of a metabolite, for a single individual in a particular group
#' @return Half-life for this metabolite, in a particular group for a single individual
#' @importFrom Rdpack reprompt
#' @importFrom nlme lme
#' @importFrom stats uniroot
#' @references
#' Wickham, H. (2022). dplyr: A Grammar of Data Manipulation. R package version 1.0.10. 
#' Available at: \url{https://CRAN.R-project.org/package=dplyr}
#'
#' Pinheiro, J. C., & Bates, D. M. (2022). nlme: Linear and Nonlinear Mixed Effects Models. R package version 3.1-153.
#' Available at: \url{https://CRAN.R-project.org/package=nlme}
#' 
#' @examples
#' require(gammaFuncModel)
#' require(rootSolve)
#' require(dplyr)
#' require(nlme)
#' \donttest{
#' modify.df <- data.frame(
#'   ID = rep(sprintf("%02d", 1:10), each = 9 * 3),
#'   Time = rep(rep(1:9, each = 3), 10),
#'   Diet = as.factor(rep(1:3, times = 9 * 10)),
#'   Age = rep(sample(20:70, 10, replace = TRUE), each = 9 * 3),
#'   BMI = round(rep(runif(10, 18.5, 35), each = 9 * 3), 1),
#'   Concentration = NA
#' )
#' for (i in 1:10) {
#'  for (d in 1:3) {
#'    C0 <- runif(1, 10, 15)    # initial concentration
#'    k <- runif(1, 0.1, 0.3)   # decay rate constant
#'    modify.df$Concentration[modify.df$ID == sprintf("%02d", i) & modify.df$Diet == d] <- 
#'      C0 * exp(-k * modify.df$Time[modify.df$ID == sprintf("%02d", i) & modify.df$Diet == d])
#'  }
#'}
#' covariates <- c("ID", "Diet", "Age", "BMI")
#' model <- gammaFunction(
#'   modify.df, 
#'   covariates, 
#'   time_grp_inter = FALSE, 
#'   return_ml_model = FALSE, 
#'   include_grp = TRUE
#'   )[[1]]
#' test_data = modify.df %>% 
#'   filter(Diet == 1 & ID == "03") %>% 
#'   select(-c("Concentration", "ID", "Diet")) 
#' Tmax <- calculate_Tmax(data = test_data, model, grp_var = 1, ID = "03", grp_name = 'Diet', ref = 1)
#' Cmax <- calculate_Cmax(data = test_data, model, grp_var = 1, ID = "03", grp_name = "Diet", Tmax)
#' f_dat = modify.df %>% filter(Diet == 1 & ID == "03") %>% select(-Concentration)
#' f <- generate_f_function(
#'   data = f_dat, 
#'   model = model,  
#'   grp_var = 1, 
#'   grp_name = "Diet", 
#'   ID = "03", 
#'   ref = 1)
#' half_life <- calculate_half_life(f, Tmax, Cmax)
#' }
#' @export
calculate_half_life = function(f, Tmax, Cmax) {
  
  root_finder <- function(t) {
    exp_predictions <- f(t)  
    return(exp_predictions - 0.5 * Cmax)  
  }
  
  result <- tryCatch(
    {
      uniroot(root_finder, interval = c(Tmax, 1000), tol = 1e-8)$root
    },
    error = function(e) {
      if (grepl("f\\(\\) values at end points not of opposite sign", e$message)) {
        return(NA)  
      } else {
        stop(e)  
      }
    }
  )
  
  return(result)
}



#' Function that produces Area Under the Curve(AUC) property for a single individual in a particular group, for a specific metabolite
#'
#' @param f function that returns the prediction of a metabolite concentration, for a single individual in a particular group
#' @param upperbound Numeric value that serves as the upperbound of integration
#' @return AUC for this metabolite, in a particular group for a single individual
#' @importFrom Rdpack reprompt
#' @importFrom nlme lme
#' @importFrom cubature adaptIntegrate
#' @references
#' Wickham, H. (2022). dplyr: A Grammar of Data Manipulation. R package version 1.0.10. 
#' Available at: \url{https://CRAN.R-project.org/package=dplyr}
#'
#' Pinheiro, J. C., & Bates, D. M. (2022). nlme: Linear and Nonlinear Mixed Effects Models. R package version 3.1-153.
#' Available at: \url{https://CRAN.R-project.org/package=nlme}
#' 
#' @examples
#' \donttest{
#' require(gammaFuncModel)
#' require(cubature)
#' require(dplyr)
#' require(nlme)
#' modify.df <- data.frame(
#'   ID = rep(sprintf("%02d", 1:10), each = 9 * 3),
#'   Time = rep(rep(1:9, each = 3), 10),
#'   Diet = as.factor(rep(1:3, times = 9 * 10)),
#'   Age = rep(sample(20:70, 10, replace = TRUE), each = 9 * 3),
#'   BMI = round(rep(runif(10, 18.5, 35), each = 9 * 3), 1),
#'   Concentration = NA
#' )
#' for (i in 1:10) {
#'  for (d in 1:3) {
#'    C0 <- runif(1, 10, 15)    # initial concentration
#'    k <- runif(1, 0.1, 0.3)   # decay rate constant
#'    modify.df$Concentration[modify.df$ID == sprintf("%02d", i) & modify.df$Diet == d] <- 
#'      C0 * exp(-k * modify.df$Time[modify.df$ID == sprintf("%02d", i) & modify.df$Diet == d])
#'  }
#'}
#' covariates <- c("ID", "Diet", "Age", "BMI")
#' model <- gammaFunction(
#'   modify.df, 
#'   covariates, 
#'   time_grp_inter = FALSE, 
#' return_ml_model = FALSE, include_grp = TRUE
#' )[[1]]
#' test_data <- modify.df %>% 
#'   filter(Diet == 1 & ID == "04") %>% 
#'   select(-c("Concentration", "ID", "Diet")) 
#' f_dat = modify.df %>% 
#'   filter(Diet == 1 & ID == "04") %>% 
#'   select(-Concentration)
#' f <- generate_f_function(
#'   data = f_dat, 
#'   model = model,  
#'   grp_var = 1, 
#'   grp_name = "Diet", 
#'   ID = "04", 
#'   ref = 1
#' )
#' AUC <- calculate_AUC(f, 9)
#' AUCInf <- calculate_AUC(f, Inf)
#' }
#' @export
calculate_AUC = function(f, upperbound){
  
  integrate_f <- function(t) {
    f(t)
  }
  
  AUC <- adaptIntegrate(integrate_f, lowerLimit = 1, upperLimit = upperbound, tol = 1e-12)$integral
  return(AUC)
}





#' Function that returns a data frame for Tmax, Cmax, half-life, AUC and AUCInf for metabolites
#'
#' @param df Data frame containing columns Group(factor); ID(subject ID: character); Time(positive: numeric); other individiual characteristics covariates (exlcluding other forms of 'Time')
#' @param met_vec Vector of metabolite names
#' @param models Fitted models for all metabolites of interest
#' @param grp_name Name of the grouping variable
#' @param covariates Vector containing the names of the "ID" covariate, grouping covariate and other covariates excluding any "Time" covariates
#' @param ref reference level for the grouping variable. could be numeric or character
#' @return Data frame with the pharmacokinetic properties of each metabolite
#' @importFrom Rdpack reprompt
#' @importFrom dplyr summarise pick bind_rows mutate %>%
#' @importFrom stats na.omit
#' @references
#' Wickham, H. (2022). dplyr: A Grammar of Data Manipulation. R package version 1.0.10. 
#' Available at: \url{https://CRAN.R-project.org/package=dplyr}
#'
#' Pinheiro, J. C., & Bates, D. M. (2022). nlme: Linear and Nonlinear Mixed Effects Models. R package version 3.1-153.
#' Available at: \url{https://CRAN.R-project.org/package=nlme}
#' 
#' @examples
#' require(gammaFuncModel)
#' require(dplyr)
#' \donttest{
#' df <- data.frame(
#'  ID = rep(sprintf("%02d", 1:10), each = 9 * 3),             
#'  Time = rep(rep(1:9, each = 3), 10),                      
#'  Diet = as.factor(rep(1:3, times = 9 * 10)),             
#'  Age = rep(sample(20:70, 10, replace = TRUE), each = 9 * 3), 
#'  BMI = round(rep(runif(10, 18.5, 35), each = 9 * 3), 1)     
#')
#' metvar <- paste0("met", 1:10)
#' n_rows <- nrow(df)
#' concentration_data <- sapply(1:10, function(m) {
#'  shape <- runif(1, 2, 5)     
#'  scale <- runif(1, 1, 3)     
#'  rgamma(n_rows, shape = shape, scale = scale)
#'})
#' colnames(concentration_data) <- metvar
#' df <- cbind(df, as.data.frame(concentration_data))
#' covariates <- c("ID", "Diet", "Age", "BMI")
#' mods <- generate_models(df = df, met_vec = metvar, covariates = covariates, graph = 'None')
#' result <- pk_calculation(
#'   df = df, 
#'   met_vec = metvar, 
#'   models = mods, 
#'   grp_name = "Diet", 
#'   covariates = covariates
#'  )
#' }
#' @export
pk_calculation <- function(df, met_vec, models, grp_name = "Diet", covariates, ref = 1) {
  names <- setdiff(covariates, c("ID", grp_name))
  names <- c(names, "Time")
  dat = df %>% select(all_of(c(covariates, "Time")))

  
  results <- list()
  for (i in 1:length(met_vec)) {
    metabolite <- met_vec[i]
    model <- models[[metabolite]]
    ref_val <- ref 
    f_funcs <- precompute_f_list(dat, model, grp_name = "Diet", ref = ref_val)
    
    # Calculate Tmax, Cmax, half life, AUC and AUCInf for each individual in each Group
    individual_results <- dat %>%
      group_by(ID, !!sym(grp_name)) %>%
      summarise(
        Tmax = calculate_Tmax(pick(all_of(names)), model,  !!sym(grp_name), ID, ref = ref_val),
        Cmax = calculate_Cmax(pick(all_of(names)), model,  !!sym(grp_name), ID,  grp_name, Tmax),
        half_life = calculate_half_life(f_funcs[[paste0(unique(ID), "-", as.character(unique(!!sym(grp_name))))]], Tmax, Cmax),
        AUC = calculate_AUC(f_funcs[[paste0(unique(ID), "-", as.character(unique(!!sym(grp_name))))]],  upperbound = 8),
        AUCInf = calculate_AUC(f_funcs[[paste0(unique(ID), "-", as.character(unique(!!sym(grp_name))))]], upperbound = Inf),
        .groups = "drop"
      ) %>%
      mutate(Metabolite = metabolite)
    
    results[[metabolite]] <- individual_results
  }
  
 
  final_results <- bind_rows(results) %>%
    group_by(!!sym(grp_name), Metabolite) %>%
    summarize(
      Cmax = round_to_3_sig(geometric_mean(na.omit(.data$Cmax))),
      Tmax = round_to_3_sig(geometric_mean(na.omit(.data$Tmax))) - 1,
      AUC = round_to_3_sig(geometric_mean(na.omit(.data$AUC))),
      AUCInf = round_to_3_sig(geometric_mean(na.omit(.data$AUCInf))),
      Half = round_to_3_sig(geometric_mean(na.omit(.data$half_life))),
      .groups = "drop"
    )
  
  return (final_results)
}