#' print.glmLimits function
#'
#' This function shows a short summary of the glmLimits function results.
#' @param x an object of class "glmLimits".
#' @keywords GLM LOD bayesian
#' @export
#' @examples
#' Limits::print.glmLimits()

print.glmLimits <- function(x, ...) {
  print(x$summary)
}