#' Plotting of histograms of residuals to visualize the normality assumption
#'
#' @param data A data table or data frame with data generated from \code{residual_histogram_data}, where Comparison, residuals, standardized_residuals and more are included. If raw data of type LFDT is the input, it will be converted to residual histo gram data. To allow such transformation set \code{raw} to \code{TRUE}, and use \code{MOR} and \code{groups} argument to control the transformation
#' @param raw Is \code{data} raw, or is it on the appropriate form described in the \code{data} argument
#' @param MOR Is your data table containing mean of replicates? Only need to specify if \code{raw = TRUE}
#' @param groups What are the grouping column names? Only need to specify if \code{raw = TRUE}
#' @param method What method should be used to calculate the optimal bandwidths
#' @param units The units of the analyte considered. These are included in the x-axis and y-axis names of the plots
#' @param studentize Should we studentize the raw residuals. Default is \code{TRUE}. If \code{TRUE}, you may choose if a N(0,1) density curve should be drawn. This may be triggered by setting \code{include_curve = TRUE}
#' @param include_curve Should a density curve of a N(0,1) be displayed over the histogram to visualize the assumed normality. Only relevant if \code{studentize = TRUE}
#' @param save Save plot to disk. Default is \code{FALSE}. If \code{TRUE}, you may specify saving location, file and other with the next four parameters
#' @param save_path Saving location given by a string. Default is the working directory
#' @param save_name Name of the plot used to save. If \code{save_name = 'Auto'}, a generic name is used based on \code{data}
#' @param save_device What file do you want to save the plot as. Valid inputs are 'pdf', 'png', 'jpeg', and more. Using 'png' will produce plots with retina \code{dpi}
#' @param size What should the size of the saved plot be? The number specified is the height of the plot. The width of the plot is calculated based on this and the golden ratio
#'
#' @return A plot (\code{ggplot2} object) based on data
#' @export
#'
#' @examples print(1)

plot_residual_histograms <- function(data, raw = FALSE, MOR = TRUE, groups = "Comparison", method = "SNRR", units = "mmol/L", studentize = TRUE, include_curve = TRUE, save = FALSE, save_path = getwd(), save_name = "Auto", save_device = "pdf", size = 8){
  data <- setDT(data)
  golden_ratio <- (1+sqrt(5)) / 2
  if(raw){
    data <- residual_vs_fitted(data = data, groups = groups, MOR = MOR)
    data <- residual_histogram_data(data = data)
  }
  if(studentize){
    if(include_curve){
      figure <- ggplot(data = data) +
        geom_histogram(mapping = aes(x = studentized_residuals, y = ..density.., fill = "Histogram for residuals"),
                       na.rm = TRUE, alpha = 0.2, color = "black", binwidth = function(x) get_optimal_binwidth(x, method = method)) +
        geom_density(mapping = aes(x = studentized_residuals, color = "Empirical density for residuals"),
                       na.rm = TRUE, alpha = 0.6, size = 1.25) +
        geom_line(mapping = aes(x = x_studentized, y = dnorm(x = x_studentized, mean = 0, sd = 1), color = "Theoretical density for model error terms"), size = 1.25) +
        theme_bw() + facet_wrap(facets = . ~ Comparison, scales = "free") +
        labs(title = "Histograms and density estimates for estimated error terms (residuals)",
             subtitle = "We aim to have histograms with a shape close to the green line, which signify the N(0,1)'s density",
             color = "Which density:",
             fill = "Histogram:") +
        scale_x_continuous(name = "Studentized results", n.breaks = 8) +
        scale_y_continuous(name = "Density", n.breaks = 8) +
        scale_fill_manual(values = c("Histogram for residuals" = "orange")) +
        scale_color_manual(values = c("Empirical density for residuals" = "#55CDEC", "Theoretical density for model error terms" = "#70AD47")) +
        theme(plot.title = element_text(face = "bold", color = "#000000", hjust = 0.5),
              plot.subtitle = element_text(color = "#000000", hjust = 0.5),
              axis.title = element_text(face = "bold", color = "#000000"),
              axis.text.x = element_text(color = "#000000", vjust = 1),
              axis.text.y = element_text(color = "#000000", hjust = 1),
              axis.ticks = element_line(size = 1, lineend = "butt", color = "#55CDEC"),
              legend.background = element_rect(fill = "#000000", color = "#000000"),
              legend.key = element_rect(fill = "#000000", color = "#000000"),
              legend.title = element_text(face = "bold", color = "#FFFFFF"),
              legend.text = element_text(color = "#FFFFFF"),
              strip.background = element_rect(fill = "#000000", color = "#000000"),
              strip.text = element_text(face = "bold", color = "#FFFFFF"))
    }
    else{
      figure <- ggplot(data = data) +
        geom_histogram(mapping = aes(x = studentized_residuals, y = ..density.., fill = "Histogram for residuals"),
                       na.rm = TRUE, alpha = 0.2, color = "black", binwidth = function(x) get_optimal_binwidth(x, method = method)) +
        geom_density(mapping = aes(x = studentized_residuals, color = "Empirical density for residuals"),
                     na.rm = TRUE, alpha = 0.6, size = 1.25) +
        theme_bw() + facet_wrap(facets = . ~ Comparison, scales = "free") +
        labs(title = "Histograms and density estimates for estimated error terms (residuals)",
             subtitle = "We aim to have histograms with a shape close the N(0,1)'s density",
             color = "Which density:",
             fill = "Histogram:") +
        scale_x_continuous(name = "Studentized results", n.breaks = 8) +
        scale_y_continuous(name = "Density", n.breaks = 8) +
        scale_fill_manual(values = c("Histogram for residuals" = "orange")) +
        scale_color_manual(values = c("Empirical density for residuals" = "#55CDEC")) +
        theme(plot.title = element_text(face = "bold", color = "#000000", hjust = 0.5),
              plot.subtitle = element_text(color = "#000000", hjust = 0.5),
              axis.title = element_text(face = "bold", color = "#000000"),
              axis.text.x = element_text(color = "#000000", vjust = 1),
              axis.text.y = element_text(color = "#000000", hjust = 1),
              axis.ticks = element_line(size = 1, lineend = "butt", color = "#55CDEC"),
              legend.background = element_rect(fill = "#000000", color = "#000000"),
              legend.key = element_rect(fill = "#000000", color = "#000000"),
              legend.title = element_text(face = "bold", color = "#FFFFFF"),
              legend.text = element_text(color = "#FFFFFF"),
              strip.background = element_rect(fill = "#000000", color = "#000000"),
              strip.text = element_text(face = "bold", color = "#FFFFFF"))
    }
  }
  else{
    if(include_curve){
      figure <- ggplot(data = data) +
        geom_histogram(mapping = aes(x = residuals, y = ..density.., fill = "Histogram for residuals"),
                       na.rm = TRUE, alpha = 0.2, color = "black", binwidth = function(x) get_optimal_binwidth(x, method = method)) +
        geom_density(mapping = aes(x = residuals, color = "Empirical density for residuals"),
                     na.rm = TRUE, alpha = 0.6, size = 1.25) +
        geom_line(mapping = aes(x = x, y = dnorm(x = x, mean = mean(x), sd = sd(x)), color = "Theoretical density for model error terms"), size = 1.25) +
        theme_bw() + facet_wrap(facets = . ~ Comparison, scales = "free") +
        labs(title = "Histograms and density estimates for estimated error terms (residuals)",
             subtitle = "We aim to have histograms with a shape close to the green line, which signify the N(mu,sigma)'s density",
             color = "Which density:",
             fill = "Histogram:") +
        scale_x_continuous(name = "Raw residuals", n.breaks = 8) +
        scale_y_continuous(name = "Density", n.breaks = 8) +
        scale_fill_manual(values = c("Histogram for residuals" = "orange")) +
        scale_color_manual(values = c("Empirical density for residuals" = "#55CDEC", "Theoretical density for model error terms" = "#70AD47")) +
        theme(plot.title = element_text(face = "bold", color = "#000000", hjust = 0.5),
              plot.subtitle = element_text(color = "#000000", hjust = 0.5),
              axis.title = element_text(face = "bold", color = "#000000"),
              axis.text.x = element_text(color = "#000000", vjust = 1),
              axis.text.y = element_text(color = "#000000", hjust = 1),
              axis.ticks = element_line(size = 1, lineend = "butt", color = "#55CDEC"),
              legend.background = element_rect(fill = "#000000", color = "#000000"),
              legend.key = element_rect(fill = "#000000", color = "#000000"),
              legend.title = element_text(face = "bold", color = "#FFFFFF"),
              legend.text = element_text(color = "#FFFFFF"),
              strip.background = element_rect(fill = "#000000", color = "#000000"),
              strip.text = element_text(face = "bold", color = "#FFFFFF"))
    }
    else{
      figure <- ggplot(data = data) +
        geom_histogram(mapping = aes(x = residuals, y = ..density.., fill = "Histogram for residuals"),
                       na.rm = TRUE, alpha = 0.2, color = "black", binwidth = function(x) get_optimal_binwidth(x, method = method)) +
        geom_density(mapping = aes(x = residuals, color = "Empirical density for residuals"),
                     na.rm = TRUE, alpha = 0.6, size = 1.25) +
        theme_bw() + facet_wrap(facets = . ~ Comparison, scales = "free") +
        labs(title = "Histograms and density estimates for estimated error terms (residuals)",
             subtitle = "We aim to have histograms with a shape close the N(mu,sigma)'s density",
             color = "Which density:",
             fill = "Histogram:") +
        scale_x_continuous(name = "Raw residuals", n.breaks = 8) +
        scale_y_continuous(name = "Density", n.breaks = 8) +
        scale_fill_manual(values = c("Histogram for residuals" = "orange")) +
        scale_color_manual(values = c("Empirical density for residuals" = "#55CDEC")) +
        theme(plot.title = element_text(face = "bold", color = "#000000", hjust = 0.5),
              plot.subtitle = element_text(color = "#000000", hjust = 0.5),
              axis.title = element_text(face = "bold", color = "#000000"),
              axis.text.x = element_text(color = "#000000", vjust = 1),
              axis.text.y = element_text(color = "#000000", hjust = 1),
              axis.ticks = element_line(size = 1, lineend = "butt", color = "#55CDEC"),
              legend.background = element_rect(fill = "#000000", color = "#000000"),
              legend.key = element_rect(fill = "#000000", color = "#000000"),
              legend.title = element_text(face = "bold", color = "#FFFFFF"),
              legend.text = element_text(color = "#FFFFFF"),
              strip.background = element_rect(fill = "#000000", color = "#000000"),
              strip.text = element_text(face = "bold", color = "#FFFFFF"))
    }
  }
  if(save){
    if(tolower(save_name)=="auto"){
      save_name <- paste0("Residual histograms for all MS comparisons",".",save_device)
    }
    else{
      save_name <- paste0(save_name,".",save_device)
    }
    if(save_device == "png"){
      ggsave(filename = save_name, plot = plot, device = save_device, path = save_path, height = size, width = size * golden_ratio, dpi = "retina")
    }
    else{
      ggsave(filename = save_name, plot = plot, device = save_device, path = save_path, height = size, width = size * golden_ratio)
    }
  }
  return(figure)
}
