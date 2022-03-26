#' plot precision against concentration
#'
#' @param data A data table with two ID-columns, namely \code{Comparison} and \code{SampleID}
#' @param units The units of the analyte considered. These are included in the x-axis name of the plots
#' @param save Save plot to disk. Default is \code{FALSE}. If \code{TRUE}, you may specify saving location, file and other with the next four parameters
#' @param save_path Saving location given by a string. Default is the working directory
#' @param save_name Name of the plot used to save. If \code{save_name = 'Auto'}, a generic name is used based on \code{data}
#' @param save_device What file do you want to save the plot as. Valid inputs are 'pdf', 'png', 'jpeg', and more. Using 'png' will produce plots with retina \code{dpi}
#' @param size What should the size of the saved plot be? The number specified is the height of the plot. The width of the plot is calculated based on this and the golden ratio
#'
#' @description Plot nice-looking plots of CS-wise precision against concentration. This function is a convenience function that makes plotting much more efficient than doing it manually.
#'
#' @details If the \code{Measure} column in \code{data} is lambda, grouping colors will not be included because no groups are available. If more than one comparison is considered, \code{data} will go through the function \code{facet_wrap()}
#'
#' @return A plot (\code{ggplot2} object) based on data. Use \code{lapply()} to make multiple plots at once. Then combine the plots using the \code{egg} package.
#' @export
#'
#' @examples plot_precision_against_concentration(data = precision_of_replicates(MS_wise(sampled_cs_measurements)))
plot_precision_against_concentration <- function(data, units = "mmol/L", save = FALSE, save_path = getwd(), save_name = "Auto", save_device = "pdf", size = 8){
  if(all(data$Measure=="lambda")){
    if(length(unique(data$Comparison)) > 1){
      plot <- ggplot(data = data) +
        geom_smooth(mapping = aes(x = Concentration, y = Value),
                    method = "lm",
                    formula = y ~ x,
                    alpha = 0.1,
                    se = TRUE,
                    na.rm = TRUE) +
        geom_point(mapping = aes(x = Concentration, y = Value),
                   shape = 21,
                   fill = "blue",
                   size = 2,
                   color = "#000000",
                   alpha = 0.8,
                   na.rm = TRUE) +
        facet_wrap(facets = . ~ Comparison, scales = "free") +
        labs(title = paste0("Estimates for CS-wise ",data$Measure[1]," for ",data$Comparison[1]),
             subtitle = "Are there evident patterns in either MS of given MS comparison") +
        scale_x_continuous(name = paste0("Concentration"," (",units,")"), n.breaks = 10) +
        scale_y_continuous(name = data$Measure[1], n.breaks = 10) + theme_bw() +
        theme(legend.background = element_rect(fill = "#000000", color = "#000000"),
              legend.title = element_text(face = "bold", color = "#FFFFFF"),
              legend.text = element_text(face = "bold", color = "#FFFFFF"),
              legend.position = "top",
              plot.title = element_text(hjust = 0.5),
              plot.subtitle = element_text(hjust = 0.5),
              axis.text.x = element_text(vjust = 1),
              axis.text.y = element_text(hjust = 1),
              axis.ticks = element_line(color = "#19282F", size = 1.5, lineend = "round"),
              strip.background = element_rect(fill = "#000000"),
              strip.text = element_text(face = "bold", color = "#FFFFFF"))
    }
    else{
      plot <- ggplot(data = data) +
        geom_smooth(mapping = aes(x = Concentration, y = Value),
                    method = "lm",
                    formula = y ~ x,
                    alpha = 0.1,
                    se = TRUE,
                    na.rm = TRUE) +
        geom_point(mapping = aes(x = Concentration, y = Value),
                   shape = 21,
                   fill = "blue",
                   size = 2,
                   color = "#000000",
                   alpha = 0.8,
                   na.rm = TRUE) +
        labs(title = paste0("Estimates for CS-wise ",data$Measure[1]," for ",data$Comparison[1]),
             subtitle = "Are there evident patterns in either MS of given MS comparison") +
        scale_x_continuous(name = paste0("Concentration"," (",units,")"), n.breaks = 10) +
        scale_y_continuous(name = data$Measure[1], n.breaks = 10) + theme_bw() +
        theme(legend.background = element_rect(fill = "#000000", color = "#000000"),
              legend.title = element_text(face = "bold", color = "#FFFFFF"),
              legend.text = element_text(face = "bold", color = "#FFFFFF"),
              legend.position = "top",
              plot.title = element_text(hjust = 0.5),
              plot.subtitle = element_text(hjust = 0.5),
              axis.text.x = element_text(vjust = 1),
              axis.text.y = element_text(hjust = 1),
              axis.ticks = element_line(color = "#19282F", size = 1.5, lineend = "round"))
    }

  }
  else{
    plot <- ggplot(data = data) +
      geom_smooth(mapping = aes(x = Concentration, y = Value, color = MS, fill = MS),
                  method = "lm",
                  formula = y ~ x,
                  alpha = 0.1,
                  se = TRUE,
                  na.rm = TRUE) +
      geom_point(mapping = aes(x = Concentration, y = Value, fill = MS),
                 shape = 21,
                 size = 2,
                 color = "#000000",
                 alpha = 0.8,
                 na.rm = TRUE) +
      labs(title = paste0("Estimates for CS-wise ",data$Measure[1]," for ",data$Comparison[1]),
           subtitle = "Are there evident patterns in either MS of given MS comparison",
           fill = "Measurement system",
           color = "Measurement system") +
      scale_x_continuous(name = paste0("Concentration"," (",units,")"), n.breaks = 10) +
      scale_y_continuous(name = data$Measure[1], n.breaks = 10) + theme_bw() +
      scale_fill_manual(values = c("#F76E11","#ED5EDD")) +
      scale_color_manual(values = c("#F76E11","#ED5EDD")) +
      theme(legend.background = element_rect(fill = "#000000", color = "#000000"),
            legend.title = element_text(face = "bold", color = "#FFFFFF"),
            legend.text = element_text(face = "bold", color = "#FFFFFF"),
            legend.position = "top",
            plot.title = element_text(hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.5),
            axis.text.x = element_text(vjust = 1),
            axis.text.y = element_text(hjust = 1),
            axis.ticks = element_line(color = "#19282F", size = 1.5, lineend = "round"))
  }
  if(save){
    if(tolower(save_name)=="auto"){
      save_name <- paste0("Estimates for CS wise ",data$Measure[1]," for ",data$Comparison[1],".",save_device)
    }
    else{
      save_name <- paste0(save_name,".",save_device)
    }
    if(save_device == "png"){
      ggsave(filename = save_name, plot = plot, device = save_device, path = save_path, height = size, width = size * (1+sqrt(5))/2,dpi="retina")
    }
    else{
      ggsave(filename = save_name, plot = plot, device = save_device, path = save_path, height = size, width = size * (1+sqrt(5))/2)
    }
  }
  return(plot)


}
