#' Adding to the generic print function 
#' 
#' @param x The ccobject
#' @param ... passthrough
#' @return a plot
#' @export
#' 

print.ccobject <- function(x, ...) {
  cat(x$Text)
  cat("\n\n")
  print(x$Data, row.names=TRUE)
  cat("\n")
  print(t(x$Estimates))
  cat("\nMean survival time:", round(as.numeric(x$Survival[[5]]),2), "With SE:", round(as.numeric(x$Survival[[6]]),2))
  cat("\n")
}





#' Adding to the generic plot function 
#' 
#' @param x The ccobject
#' @param ... passthrough
#' @return a plot
#' @import ggplot2 tibble forcats
#' @export
#' 

plot.ccobject <- function(x, ...) {
  temp <- x$Estimates %>%  t() %>% as.data.frame() %>% tibble::rownames_to_column(var="Estimator")
  temp$Estimator <- factor(temp$Estimator, labels = )
  
  temp %>% 
    ggplot(aes(y = temp$Estimate, x = fct_relevel(temp$Estimator,c("ZT","BT","CC","AS")), ymax = temp$"0.95UCL", ymin = temp$"0.95LCL")) + 
    geom_point(shape=15, size=5) +  
    geom_errorbar(width = 0.2, size = 1.1) + 
    labs(title="Estimators", x = "", y = "Cost") +
    coord_flip() 
}
