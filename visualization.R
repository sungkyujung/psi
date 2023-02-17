

library(ggplot2)
library(reshape2)

scree.plot <- function(data, title='', intercept = NA){
  
  data.pca = prcomp(data, scale = TRUE)
  Variance = 100 * data.pca$sdev^2 / sum(data.pca$sdev^2) 
  
  cum.Variance = cumsum(Variance)
  
  index <- 1:length(Variance)
  
  ggplot(as.data.frame(x=cbind(index, Variance)),
         aes(x=index, y=Variance)) +
    geom_point() +
    geom_line() +
    labs(y='variance proportion (%)', x='index', title=title) +
    geom_vline(xintercept = intercept + 0.5, alpha = 0.5, 
               linetype = 'dashed', size = 0.5, na.rm = TRUE) +
    annotate("text", na.rm = TRUE, x = intercept + 25, alpha = 0.5, y = 0.5 * max(Variance),
             label = paste0("cut = ", intercept, " (", round(cum.Variance[intercept],1) , "%)" ))

}



heatmap.matrix <- function(data, show_color_bar=TRUE, title='', xlab='', ylab=''){
  
  
  # could use geom_raster or geom_tile -- I think the former is faster
  # could add to theme:
  # panel.border = element_rect(linetype = 1, size=.5, colour = "grey80")
  
  reshaped_data <- as.data.frame(reshape2::melt(data))
  colnames(reshaped_data) <- c('obs', 'var', 'value')
  
  ggplot(data=reshaped_data,
         aes_string(x = 'var', y = 'obs')) +
    geom_raster(aes_string(fill = 'value'), show.legend=show_color_bar) +
    scale_fill_gradient2(low='blue', high='red') +
    theme(panel.background = element_blank(),
          axis.line = element_blank(),
          legend.position="bottom") +
    scale_y_continuous(expand = c(0, 0)) +
    scale_x_continuous(expand = c(0, 0)) +
    labs(title=title,
         x=xlab,
         y=ylab)
  

  
}










validation.plot <- function(lambda, result  , plot.name = '', opt.lambda = NA){
  
  r.test.df = data.frame(lambda * 180 / pi, result$r.test.list , result$r.pen.test.list)
  colnames(r.test.df) <- c("lambda", "Risk", "Risk.Pen")
  
  # ggplot(r.test.df, aes(x=lambda, y=Risk, linetype = "R")) +
  #   geom_line() +
  #   geom_line(data = r.test.df , aes(x=lambda, y=Risk.Pen, color = "R.pen") ) + 
  #   scale_linetype_manual(values = c("R" = "solid", "R.pen" = "solid")) + 
  #   scale_color_manual(values = c("R" = "black","R.pen" = "red")) +
  #   theme_bw()+
  #   labs(title = plot.name) +
  #   theme(legend.position = "none") 
  
  g = ggplot(r.test.df, aes(x=lambda, y=Risk)) +
    geom_line() +
    labs(title = plot.name) +
    labs(x = "λ (°)") + 
    theme(legend.position = "none") 
    
    
  g + geom_vline(xintercept = opt.lambda, alpha = 0.5,  linetype = 'dashed', 
               size = 0.5, na.rm = TRUE) +
    annotate("text", na.rm = TRUE, x = opt.lambda + 45, alpha = 0.7, 
               y = 0.5 * (max(r.test.df$Risk) - min(r.test.df$Risk)) + min(r.test.df$Risk),
               label = paste0("Optimal λ = ", opt.lambda, "°" ))
    
  
}









