#' confint.glmLimits function
#'
#' This function compute confidence intervals for all parameters estimated with the glmLimits function.
#' @param x an object of class "glmLimits".
#' @param level the confidence level required.
#' @param alpha optional signicance level (ignoring level value).
#' @param dec ignore this... (used to match the quantile of the glmLimits object)
#' @keywords GLM LOD bayesian
#' @export
#' @examples
#' Limits::confint.glmLimits(x, level = .95, alpha = NULL, dec = 4, ...)

confint.glmLimits <- function(x, level = .95, alpha = NULL, dec = 4, ...) {
  if(is.null(alpha) & !is.null(level)) 
    alpha <- 1 - level
  else if(!is.null(alpha) & !is.null(level))
  {
    if (round(alpha, dec) == round(1 - level, dec)){}
    else {
      message("Confint will be calculated using alpha = ", paste0(alpha))
      level <- 1 - alpha
    }
  }
  if(any(!paste0(round(100 * abs(c(0, 1) - alpha / 2),2), "%") %in%  colnames(x$summary$quantiles))) 
    stop("\nFor this level = ",  paste0(level), " value, you need to define the corresponding quantiles = c(",alpha/2, ", ",  1 - alpha/2,") or alpha.corrected = ",alpha, " in glmLimits")
  x$summary$quantiles[, paste0(round(c(alpha / 2, 1 - alpha / 2) * 100, 2),"%")]
}