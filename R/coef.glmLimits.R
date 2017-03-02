#' coef.glmLimits function
#'
#' This function shows estimated coeficients estimated with the glmLimits function.
#' @param x an object of class "glmLimits".
#' @keywords GLM LOD bayesian
#' @export
#' @examples
#' Limits::coef.glmLimits(x, ...)

coef.glmLimits <- function(x, ...) {
  x$summary$statistics[,1:2]
}