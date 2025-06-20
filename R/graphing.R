#' Function that generate plots for metabolite models
#' @param graph character string, 'None' by default. If not 'None, in addition to returning models, produces pdf file of graphs based on the specific value of 'graph'. 
#' @param df Data frame containing columns Group(factor); ID(subject ID: character); Time(positive: numeric); other Time terms (numeric); other individidual characteristics covariates; as well columns of metabolite concentrations;
#'    Note: All non-concentration columns must be complete (No missing values); concentration columns can have missing values in the forms of either numeric 0 or 'NA'.
#' @param met_vec the vector of metabolite names
#' @param covariates Vector containing the names of the "ID" covariate, grouping covariate and other covariates excluding any "Time" covariates;
#' @param grp is the grouping variable;
#' @param models a list of fitted non-linear mixed effects metabolite models
#' @param save_path location (file path, not directory) where the pdf file will be saved (must end in '.pdf'); default is NULL, i.e. pdf is saved to a temporary location
#' @return A pdf file for fitted concenration curves that is saved to a user provided file location; otherwise saved to a temporary location
#' @importFrom Rdpack reprompt
#' @importFrom dplyr select all_of group_by summarize %>%
#' @importFrom ggplot2 ggplot aes geom_line labs theme_bw ggtitle element_text guides guide_legend theme margin
#' @importFrom scales hue_pal
#' @importFrom gridExtra grid.arrange
#' @importFrom grid unit
#' @importFrom patchwork wrap_plots plot_annotation plot_layout
#' @importFrom grDevices pdf dev.off
#' @importFrom stats predict 
#' @importFrom rlang sym
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
#' require(patchwork)
#' require(scales)
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
#' mods <- generate_models(df = df, met_vec = metvar, covariates = covariates, graph = 'None')
#' generatePlot(
#'   graph = "individual_separated", 
#'   df = df, 
#'   met_vec = metvar, 
#'   covariates = covariates, 
#'   grp = "Diet", 
#'   models = mods,
#'   save_path = NULL
#'  )
#' }
#' @export
generatePlot <- function(graph, df, met_vec, covariates, grp = "Diet", models, save_path = NULL) {
  data_df <- df %>% dplyr::select(all_of(c(covariates, "Time")))
  
  ### individual plots
  if (graph == "individual_separated") { 
    # Use tempdir() with default name if user doesn't specify path
    if (is.null(save_path)) {
      save_path <- file.path(tempdir(), "Individual_Predictions_separated.pdf")
    }
    
    pdf(save_path, width = 10, height = 10)
    for (i in 1: length(met_vec)){
      metabolite <- met_vec[i]
      mod <- models[[metabolite]]
      data_df$Concentration <- df[, metabolite]
      data_df$pred <- exp(predict(mod, newdata = data_df))
      plots <- list()
      
      for (id in unique(data_df$ID)) {
        dat <- data_df[data_df$ID == id, ]
        
        p <- ggplot(dat, aes(x = Time - 1)) +
          geom_line(aes(y = Concentration, color = !!sym(grp)), linetype = "dashed") +
          geom_line(aes(y = pred, color = !!sym(grp)), linetype = "solid") +
          labs(
            x = "Time (h)",
            y = "Concentration",
            title = paste("Individual", id)
          ) +
          theme_bw() +
          theme(legend.position = "none")  
        
        plots[[id]] <- p
      }
      
      # Use patchwork to combine plots with a single shared legend
      combined_plot <- wrap_plots(plots) +
        plot_layout(guides = "collect") & 
        theme(legend.position = "bottom")
      
      # Add overall title and caption
      combined_plot <- combined_plot +
        plot_annotation(
          title = metabolite,
          caption = "Actual Data (dashed) vs Fitted Data (solid)"
        )
      
      print(combined_plot)
      
    }
    dev.off()
    
    message("PDF saved to: ", save_path, 
            "\nNote: If not specified, file is stored in a temporary directory and may be deleted when the R session ends.")
  }
  
  ### combined plot: all individuals plotted in a single plot instead of one plot
  if (graph == "individual_combined") {
    # Use tempdir() with default name if user doesn't specify path
    if (is.null(save_path)) {
      save_path <- file.path(tempdir(), "Individual_Predictions_combined.pdf")
    }
    
    pdf(save_path, width = 20, height = 20)
    plots_per_page <- 30  
    plot_list <- list()
    
    # Add dynamic color_var column
    data_df$color_var <- if (length(unique(data_df[[grp]])) == 1) {
      as.factor(data_df$ID)
    } else {
      as.factor(data_df[[grp]])
    }
    
    legend_title <- if (length(unique(data_df[[grp]])) == 1) {
      "ID"
    } else {
      grp
    }
    
    for (i in 1:length(met_vec)) {
      metabolite <- met_vec[i]
      mod <- models[[metabolite]]
      data_df$Concentration <- df[, metabolite]
      data_df$pred <- exp(predict(mod, newdata = data_df))
      
      p <- ggplot() +
        geom_line(data = data_df, aes(y = Concentration, x = Time - 1, group = interaction(ID, !!sym(grp)), color = color_var), linetype = "solid") +
        geom_line(data = data_df, aes(y = pred, x = Time - 1, group = interaction(ID, !!sym(grp)), color = color_var), linetype = "dashed", show.legend = FALSE) +
        theme_bw() +
        labs(title = paste0(metabolite), x = "Time(h)", y = "Concentration", color = legend_title)
      
      
      
      plot_list[[length(plot_list) + 1]] <- p
      
      # Print and reset after reaching `plots_per_page`
      if (length(plot_list) == plots_per_page || i == length(met_vec)) {
        do.call(grid.arrange, c(plot_list, ncol = 3))
        plot_list <- list()  # Reset plot list for next page
      }
    }
    dev.off()
    
    message("PDF saved to: ", save_path, 
            "\nNote: If not specified, file is stored in a temporary directory and may be deleted when the R session ends.")
  }
  
  
  
  ### average concentration line across individuals, color by group
  if (graph == "average_predictions") {
    # Use tempdir() with default name if user doesn't specify path
    if (is.null(save_path)) {
      save_path <- file.path(tempdir(), "Average_Predictions.pdf")
    }
    
    pdf(save_path, width = 10, height = 10)
    for (i in  1: length(met_vec)) {
      metabolite <- met_vec[i]
      mod <- models[[metabolite]]
      data_df$Concentration <- df[, metabolite]
      data_df$pred <- exp(predict(mod, newdata = data_df))
      
      avgDat <- data_df %>%
        group_by(!!sym(grp), Time) %>%
        summarize(avgPred = geometric_mean(pred), .groups = "keep")
      
      p <- ggplot() + 
        geom_line(data = data_df, aes(y = Concentration, x = Time - 1, group = interaction(ID, !!sym(grp)), color = !!sym(grp)), linetype = "dashed") + 
        geom_line(data = avgDat, aes(y = avgPred, x = Time -1, color = !!sym(grp)),  linetype = "solid", linewidth = 3, show.legend = FALSE) +
        theme_bw() + 
        
        labs(title = paste0(metabolite),
             x = "Time(h)", 
             y = "Concentration")
      print(p)
    }
    dev.off()
    
    message("PDF saved to: ", save_path, 
            "\nNote: If not specified, file is stored in a temporary directory and may be deleted when the R session ends.")
  }
}
