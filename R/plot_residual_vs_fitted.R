#' Plot residual vs. fitted plots for your data
#'
#' @param data A data table or data frame grouped by Comparison where for each comparison, residuals and fitted values are included. If you use a raw data set of type LFDT you must specify \code{MOR}, \code{raw} and \code{groups} correctly. See \code{residual_vs_fitted()} for more information
#' @param raw Is \code{data} raw, or is it on the appropriate form described in the \code{data} argument
#' @param MOR Is your data table containing mean of replicates? Only need to specify if \code{raw = TRUE}
#' @param groups What are the grouping column names? Only need to specify if \code{raw = TRUE}
#' @param units The units of the analyte considered. These are included in the x-axis and y-axis names of the plots
#' @param studentize Should we studentize the residuals. Default is \code{TRUE}
#' @param include_curve Should a smoothing curve be included? Default is \code{TRUE}
#' @param save Save plot to disk. Default is \code{FALSE}. If \code{TRUE}, you may specify saving location, file and other with the next four parameters
#' @param save_path Saving location given by a string. Default is the working directory
#' @param save_name Name of the plot used to save. If \code{save_name = 'Auto'}, a generic name is used based on \code{data}
#' @param save_device What file do you want to save the plot as. Valid inputs are 'pdf', 'png', 'jpeg', and more. Using 'png' will produce plots with retina \code{dpi}
#' @param size What should the size of the saved plot be? The number specified is the height of the plot. The width of the plot is calculated based on this and the golden ratio
#'
#' @return A plot (\code{ggplot2} object) based on data.
#' @export
#'
#' @examples print(1)

plot_residual_vs_fitted <- function(data, raw = TRUE, MOR = TRUE, groups = "Comparison", units = "mmol/L", studentize = TRUE, include_curve = TRUE, save = FALSE, save_path = getwd(), save_name = "Auto", save_device = "pdf", size = 8){
  data <- setDT(data)
  golden_ratio <- (1+sqrt(5))/2
  if(raw){
    data <- residual_vs_fitted(data = data, groups = groups, MOR = MOR)
  }
  if(studentize){
    data$residuals <- (data$residuals - mean(data$residuals))/sd(data$residuals)
    if(include_curve){
      figure <- ggplot(data = data) +
        geom_smooth(mapping = aes(x = `fitted values`, y = residuals), formula = y ~ x, method = "loess", span = 2, alpha = 0.3) +
        geom_point(mapping = aes(x = `fitted values`, y = residuals), shape = 21, fill = "#6FB2D2", color = "#000000") +
        geom_hline(yintercept = 0, linetype = "dashed") +
        facet_wrap(facets = . ~ Comparison, scales = "free") +
        scale_x_continuous(name = paste0("Fitted values (",units,")"), n.breaks = 8) +
        scale_y_continuous(name = paste0("Residuals (",units,")"), n.breaks = 8) +
        labs(title = "Residuals vs. fitted plots for all comparisons",
             subtitle = paste0("The residuals are ",ifelse(studentize,"studentized","raw")," and the curve displays the smoothed pattern")) +
        theme_bw() +
        theme(plot.title = element_text(hjust = 0.5),
              plot.subtitle = element_text(hjust = 0.5),
              axis.text.x = element_text(vjust = 1),
              axis.text.y = element_text(hjust = 1),
              axis.ticks = element_line(color = "#19282F", size = 1.5, lineend = "round"),
              strip.background = element_rect(fill = "#000000"),
              strip.text = element_text(face = "bold", color = "#FFFFFF"))
    }
    else{
      figure <- ggplot(data = data) +
        geom_point(mapping = aes(x = `fitted values`, y = residuals), shape = 21, fill = "#6FB2D2", color = "#000000") +
        geom_hline(yintercept = 0, linetype = "dashed") +
        facet_wrap(facets = . ~ Comparison, scales = "free") +
        scale_x_continuous(name = paste0("Fitted values (",units,")"), n.breaks = 8) +
        scale_y_continuous(name = paste0("Residuals (",units,")"), n.breaks = 8) +
        labs(title = "Residuals vs. fitted plots for all comparisons",
             subtitle = paste0("The residuals are ",ifelse(studentize,"studentized","raw"))) +
        theme_bw() +
        theme(plot.title = element_text(hjust = 0.5),
              plot.subtitle = element_text(hjust = 0.5),
              axis.text.x = element_text(vjust = 1),
              axis.text.y = element_text(hjust = 1),
              axis.ticks = element_line(color = "#19282F", size = 1.5, lineend = "round"),
              strip.background = element_rect(fill = "#000000"),
              strip.text = element_text(face = "bold", color = "#FFFFFF"))
    }
  }
  else{
    if(include_curve){
      figure <- ggplot(data = data) +
        geom_smooth(mapping = aes(x = `fitted values`, y = residuals), formula = y ~ x, method = "loess", span = 2, alpha = 0.3) +
        geom_point(mapping = aes(x = `fitted values`, y = residuals), shape = 21, fill = "#6FB2D2", color = "#000000") +
        geom_hline(yintercept = 0, linetype = "dashed") +
        facet_wrap(facets = . ~ Comparison, scales = "free") +
        scale_x_continuous(name = paste0("Fitted values (",units,")"), n.breaks = 8) +
        scale_y_continuous(name = paste0("Residuals (",units,")"), n.breaks = 8) +
        labs(title = "Residuals vs. fitted plots for all comparisons",
             subtitle = paste0("The residuals are ",ifelse(studentize,"studentized","raw")," and the curve displays the smoothed pattern")) +
        theme_bw() +
        theme(plot.title = element_text(hjust = 0.5),
              plot.subtitle = element_text(hjust = 0.5),
              axis.text.x = element_text(vjust = 1),
              axis.text.y = element_text(hjust = 1),
              axis.ticks = element_line(color = "#19282F", size = 1.5, lineend = "round"),
              strip.background = element_rect(fill = "#000000"),
              strip.text = element_text(face = "bold", color = "#FFFFFF"))
    }
    else{
      figure <- ggplot(data = data) +
        geom_point(mapping = aes(x = `fitted values`, y = residuals), shape = 21, fill = "#6FB2D2", color = "#000000") +
        geom_hline(yintercept = 0, linetype = "dashed") +
        facet_wrap(facets = . ~ Comparison, scales = "free") +
        scale_x_continuous(name = paste0("Fitted values (",units,")"), n.breaks = 8) +
        scale_y_continuous(name = paste0("Residuals (",units,")"), n.breaks = 8) +
        labs(title = "Residuals vs. fitted plots for all comparisons",
             subtitle = paste0("The residuals are ",ifelse(studentize,"studentized","raw"))) +
        theme_bw() +
        theme(plot.title = element_text(hjust = 0.5),
              plot.subtitle = element_text(hjust = 0.5),
              axis.text.x = element_text(vjust = 1),
              axis.text.y = element_text(hjust = 1),
              axis.ticks = element_line(color = "#19282F", size = 1.5, lineend = "round"),
              strip.background = element_rect(fill = "#000000"),
              strip.text = element_text(face = "bold", color = "#FFFFFF"))
    }
  }
  if(save){
    if(tolower(save_name)=="auto"){
      save_name <- paste0("Residual vs. fitted plots for all MS comparisons",".",save_device)
    }
    else{
      save_name <- paste0(save_name,".",save_device)
    }
    if(save_device == "png"){
      ggsave(filename = save_name, plot = plot, device = save_device, path = save_path, height = size, width = size * golden_ratio, dpi="retina")
    }
    else{
      ggsave(filename = save_name, plot = plot, device = save_device, path = save_path, height = size, width = size * golden_ratio)
    }
  }
  return(figure)

}

