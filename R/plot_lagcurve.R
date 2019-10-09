#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_line
#' @importFrom ggplot2 xlab
#' @importFrom ggplot2 ylab
#' @importFrom ggplot2 ggtitle
#' @importFrom tibble tibble
#' @importFrom tidyr gather
#' @importFrom magrittr %>%
#' @export

plot_lagcurve <- function(dlm_object){
  lagdf <- tibble(lag = dlm_object$coefs$lag, 
                  `97.5%` = dlm_object$coefs$upper, 
                  `2.5%` = dlm_object$coefs$lower, 
                  median = dlm_object$coefs$beta)
  plt <- lagdf %>%
    gather(key = 'quantile', value = 'coef', -lag) %>% 
    ggplot(aes(x = lag, y = coef, group = quantile, col = quantile)) + 
    geom_line() + 
    xlab('Lag number') + 
    ylab('Lag coefficient') + 
    ggtitle('Posterior median and credible interval for lag function')
  return(plt)
}