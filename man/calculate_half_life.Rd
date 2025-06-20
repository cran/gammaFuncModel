% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gammaFuncModel.R
\name{calculate_half_life}
\alias{calculate_half_life}
\title{Function that produces Half-life property for a single individual in a particular group, for a specific metabolite}
\usage{
calculate_half_life(f, Tmax, Cmax)
}
\arguments{
\item{f}{function that returns the prediction of a metabolite concentration, for a single individual in a particular group}

\item{Tmax}{Tmax property of a metabolite, for a single individual in a particular group}

\item{Cmax}{Cmax property of a metabolite, for a single individual in a particular group}
}
\value{
Half-life for this metabolite, in a particular group for a single individual
}
\description{
Function that produces Half-life property for a single individual in a particular group, for a specific metabolite
}
\examples{
require(gammaFuncModel)
require(rootSolve)
require(dplyr)
require(nlme)
\donttest{
modify.df <- data.frame(
  ID = rep(sprintf("\%02d", 1:10), each = 9 * 3),
  Time = rep(rep(1:9, each = 3), 10),
  Diet = as.factor(rep(1:3, times = 9 * 10)),
  Age = rep(sample(20:70, 10, replace = TRUE), each = 9 * 3),
  BMI = round(rep(runif(10, 18.5, 35), each = 9 * 3), 1),
  Concentration = NA
)
for (i in 1:10) {
 for (d in 1:3) {
   C0 <- runif(1, 10, 15)    # initial concentration
   k <- runif(1, 0.1, 0.3)   # decay rate constant
   modify.df$Concentration[modify.df$ID == sprintf("\%02d", i) & modify.df$Diet == d] <- 
     C0 * exp(-k * modify.df$Time[modify.df$ID == sprintf("\%02d", i) & modify.df$Diet == d])
 }
}
covariates <- c("ID", "Diet", "Age", "BMI")
model <- gammaFunction(
  modify.df, 
  covariates, 
  time_grp_inter = FALSE, 
  return_ml_model = FALSE, 
  include_grp = TRUE
  )[[1]]
test_data = modify.df \%>\% 
  filter(Diet == 1 & ID == "03") \%>\% 
  select(-c("Concentration", "ID", "Diet")) 
Tmax <- calculate_Tmax(data = test_data, model, grp_var = 1, ID = "03", grp_name = 'Diet', ref = 1)
Cmax <- calculate_Cmax(data = test_data, model, grp_var = 1, ID = "03", grp_name = "Diet", Tmax)
f_dat = modify.df \%>\% filter(Diet == 1 & ID == "03") \%>\% select(-Concentration)
f <- generate_f_function(
  data = f_dat, 
  model = model,  
  grp_var = 1, 
  grp_name = "Diet", 
  ID = "03", 
  ref = 1)
half_life <- calculate_half_life(f, Tmax, Cmax)
}
}
\references{
Wickham, H. (2022). dplyr: A Grammar of Data Manipulation. R package version 1.0.10. 
Available at: \url{https://CRAN.R-project.org/package=dplyr}

Pinheiro, J. C., & Bates, D. M. (2022). nlme: Linear and Nonlinear Mixed Effects Models. R package version 3.1-153.
Available at: \url{https://CRAN.R-project.org/package=nlme}
}
