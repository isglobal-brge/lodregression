#' summary.glmLimits function
#'
#' This function return the summary of the glmLimits function results with class "summary.rjags".
#' @param model the model (of class "jags") extracted from an object of class "glmLimits".
#' @keywords GLM LOD bayesian
#' @export
#' @examples
#' Limits::summary.glmLimits()

summary.glmLimits <- function(model, ...){
  coefs <- coef(model)
  # pred <- model$pred
  betas.sd <- coefs[, "SD"]
  V <- betas.sd^2
  betas <- coefs[, "Mean"] # extract coeficients
  W <- 0.0001
  z <- betas / sqrt(V)
  bayes.factors <- try(as.double(sqrt((V + W) / V) * exp(-(z^2 / 2) * (V / (V + W)))), TRUE)
  conf <- confint(model, alpha = model$alpha.corrected)
  
  betas.sig <- ifelse(conf[, 1] < 0 & conf[, 2] > 0, " ", "*")
  t.values <- betas / betas.sd
  
  if(model$family == "gaussian"){
    
    df <- length(model$response) - nrow(coefs)
    p.values <- 2 * pt(-abs(t.values), df = df)
    bayes <- data.frame(Estimate = round(betas, 3),
                        Sig = betas.sig,
                        Std.Error = round(betas.sd, 3),
                        t.value = round(t.values, 3),
                        p.value = format.pval(p.values),
                        BayesFactor = round(bayes.factors, 4))
    
  } else if(model$family == "binomial"){
    
    p.values <- 2 * pnorm(-abs(t.values))
    bayes <- data.frame(Estimate = round(betas, 3),
                        Sig = betas.sig,
                        Std.Error = round(betas.sd, 3),
                        z.value = round(t.values, 3),
                        p.value = format.pval(p.values),
                        BayesFactor = round(bayes.factors, 4))
    
  }
  
  res <- list(coefficients = bayes, alpha = model$alpha.corrected)
  class(res) <- "summary.rjags"
  return(res)
} 