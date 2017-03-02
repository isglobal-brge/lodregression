#' predict.glmLimits function
#'
#' This function compute the predictors estimated with the glmLimits function.
#' @param x an object of class "glmLimits".
#' @param type the type of prediction required. The default ("terms") returns fitted values of each term in the model formula on the linear predictor scale. "binomialresponse" returns binary predictions and "probability" return probabilities if family is "binomial"
#' @param p probability border value in [0,1] to classify binary outcome.
#' @keywords GLM LOD bayesian
#' @export
#' @examples
#' Limits::predict.glmLimits()

predict.glmLimits <- function(x, type = c("terms", "binomialresponse", "probability"), p = 0.5, ...) {
  if(type == "terms") res <- x$pred
  else if(type == "binomialresponse") res <- factor(I(inv.logit(x$pred) > p), labels = c("0", "1"))
  else if(type == "probability"){
    if(nlevels(as.factor(x$response)) > 2) warning("Probability type only for family = 'binomial' models.")
    else res <- inv.logit(x$pred)
  }
  res
}