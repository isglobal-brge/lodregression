#' print.summary.rjags function
#'
#' This function shows the summary of the glmLimits function results.
#' @param object an object of class "summary.rjags".
#' @keywords GLM LOD bayesian
#' @export
#' @examples
#' Limits::print.summary.rjags()

print.summary.rjags <- function(object, ...){
  print(object$coefficients)
  cat("\nSignif. codes:  < ", object$alpha, "'*' else ' '  \n")
}