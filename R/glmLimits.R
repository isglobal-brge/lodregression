#' glmLimits function
#'
#' This function allows you to adjust multiple Bayesian GLM models with predictors subjected to a limit of detection (LOD) specified with a BUGS-language description of the prior distribution, and a set of data. 
#' @param formula an object of class "formula":  a symbolic description of the model to be fitted. 
#' @param data a data.frame containing the variables in the model.
#' @param LimitVariables character string of names of the variables subjected to a limit.
#' @param lod numeric vector with the same length of LimitVariables containing the lower limit values for the variables with LOD (with the same order). 
#' @param loq optional numeric vector with the same length of LimitVariables containing the upper limit values for the variables with limits (with the same order)
#' @param family character with a description of the error distribution to be used in the model. Possible values are "gaussian" (default) for continuous outcomes or "binomial" for binary outcomes.
#' @param logLimit a logical indicating if limit variables would be assumed distributed as log-normal.
#' @param n.iter number of iterations of the Markov chain to run to be used in rjags::update(), rjags::coda.samples() and rjags::dic.samples().
#' @param n.chains the number of parallel chains for the model to be used in rjags::jags.model().
#' @param alpha.corrected significance level to be used for the confident intervals.
#' @param quantiles numeric vector containing the percentiles in [0,1] to obtain of the estimated parameters.
#' @param progress.bar type of progress bar. Possible values are "text" (default), "gui", and "none" to be used in rjags::update().
#' @param dic.type type of penalty to use in rjags::dic.samples().
#' @param n.adapt the number of iterations for adaptation to be used in rjags::jags.model(). If n.adapt = 0 then no adaptation takes place.
#' @return It returns an object of class "glmLimits" containing a model of class "jags".
#' @keywords GLM LOD bayesian
#' @export
#' @examples
#' glmLimits(formula, data, LimitVariables, lod = NULL, loq = NULL, family = "gaussian", logLimit = FALSE, n.iter = 1000, n.chains = 1, alpha.corrected = 0.05, quantiles = unique(c(alpha.corrected/2, 0.025, 0.975, 1-(alpha.corrected/2))), quiet = FALSE, progress.bar = "text", dic.type = "popt", n.adapt = 1000, ...)

glmLimits <- function(formula, data, LimitVariables, lod = NULL, loq = NULL, 
                        family = "gaussian", logLimit = FALSE, 
                        n.iter = 1000, 
                        n.chains = 1,
                        alpha.corrected = 0.05,
                        quantiles = unique(c(alpha.corrected/2, 0.025, 0.975, 1-(alpha.corrected/2))),
                        quiet = FALSE,
                        # zlevels = 5,
                        progress.bar = "text",
                        dic.type = "popt",
                        n.adapt = 1000,
                        ...) {
  
  
  fam <- match(family, c("gaussian", "binomial"), nomatch = NA)
  if (is.na(fam))
    stop (" 'family' argument should be gaussian or binomial")
  rm(fam)
  ## keep only the arguments which should go into the model frame
  mf <- match.call(expand.dots = FALSE)
  m0 <- match(c("formula", "data", "subset", "na.action"), names(mf), nomatch = 0L)
  mf <- mf[c(1L, m0)]
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  ## 1) allow model.frame to update the terms object before saving it.
  mt <- attr(mf, "terms")
  n <- nrow(data)
  
  y <- mf[,attr(mt,"response")]
  xnames <- attr(mt,"term.labels")
  X <- mf[,xnames, drop = FALSE]
  if(!exists("LimitVariables"))
    stop("'LimitVariables' argument is needed.\n It's necesary to specify the name of the variables that are subjected to a limit.")
  x <- X[,match(x = LimitVariables, table = colnames(X))]
  
  # Number of variables subjected to a limit
  nlimitvars <- length(LimitVariables)
  # 1! variable subjected to a limit
  
  if(nlimitvars == 1)
    {
    
    if(!is.null(LimitVariables) & is.null(lod) & is.null(loq))
      stop("Minimum one of both arguments ('lod' and 'loq') should not be NULL")  
    
    if(!is.null(LimitVariables) & !is.null(lod))
      {
      ind1 <- x <= lod
      if(is.null(loq)) ind <- ind1  
    } 
    
    if(!is.null(LimitVariables)  & !is.null(loq)){
      ind2 <- x >= loq 
      if(is.null(lod)) ind <- ind2
    } 
    
    if(!is.null(LimitVariables) & !is.null(lod) & !is.null(loq))
      ind <- ind1 | ind2 # x <= lod | x >= loq
    
    
    x[ind] <- NA
    if(!is.null(lod) & !is.null(loq))
      {
        censored <- rep(1, n)
        censored[ind1] <- 0
        censored[ind2] <- 2
      } else if(is.null(loq))
      censored <- as.integer(ifelse(is.na(x), 0, 1))
    else if(is.null(lod))
      censored <- as.integer(ifelse(is.na(x), 2, 1))
    } 
  
 
  
  
  # > 1 variables subjected to a LOD
  
  if(nlimitvars > 1)
    {
    if(!is.null(LimitVariables) & !is.null(lod) & is.null(loq))
    {
      if(length(lod) != nlimitvars)
        stop("length(lod) must be equal to nlimitvars and in respective order")
      else
        ind <- as.data.frame(mapply(function(x, lod) x <= lod, x = as.list(x), lod = as.list(lod)))
    } else if(!is.null(LimitVariables) & is.null(lod) & !is.null(loq))
      stop("Multiple variables subjected to a limit only for LOD")
    else if(!is.null(LimitVariables) & !is.null(lod) & !is.null(loq))
      stop("Multiple variables subjected to a limit only for LOD")
    else
      stop("'lod' argument should not be NULL")  
    x <- mapply(function(x, ind) {x[ind] <- NA; x}, x = as.list(x), ind = as.list(ind))
    censored <- apply(x, 2, function(x) as.integer(ifelse(is.na(x), 0, 1)))
    colnames(censored) <- paste0("censored", 1:nlimitvars)
    # x <- as.list(as.data.frame(x))
    # names(x) <- paste0("x", 1:nlimitvars)
  }
  

  covariates <- ncol(X) > nlimitvars
  
  if(covariates)
  {
    if(ncol(X) == nlimitvars + 1)
    {
      Z <- list(z1 = X[,!xnames %in% LimitVariables])
      names(Z) <- xnames[!xnames %in% LimitVariables]
    } else
      Z <- as.list(X[,!xnames %in% LimitVariables])
    znames <- names(Z)
    names(Z) <- paste0("z", 1:length(Z))
  } else {
    Z <- NULL
    znames <- NULL
    Zcat <- Zcont <-NULL
    catcov <- contcov <- NULL
  }
  
  # Si hi ha covariables, difreneciem continues de categoriques i creem dummys
  if(covariates)
  {
    # Search for possible factors from Z data
    nfactors <- function(z)
    {
      nlev <- nlevels(as.factor(z))
      if( is.factor(z))
        res <- nlev
      else 
        res <- 0
      res
    }
    nfactors <- unlist(lapply(Z, nfactors))
    is.z.factor <- unlist(lapply(Z, is.factor))
   
    catcov <- any(is.z.factor)
    contcov <- any(!is.z.factor)
    
    # Transform Z due to factors if factors
    if(catcov){
      faux <- function(z, Zname){
        namelevels <- levels(as.factor(z))
        Zaux <- as.data.frame(matrix(rep(NA, length(z) * length(namelevels[-1])), 
                                     ncol = length(namelevels[-1]), 
                                     dimnames = list(x = NULL, y = paste0(Zname, "", namelevels[-1]))))
        
        for(i in 1:(nlevels(as.factor(z)) - 1)){
          Zaux[,i] <- as.integer(I(z == namelevels[i + 1]))
        }
        Zaux
      }
      if(sum(is.z.factor) == 1)
        Zaux <- faux(z = Z[is.z.factor][[1]], Zname = znames[is.z.factor])
      else
        Zaux <- mapply(faux, z = Z[is.z.factor], Zname = znames[is.z.factor], 
                       USE.NAMES = FALSE)
      
      if(is.list(Zaux)){
        if(is.data.frame(Zaux))
        {Zcat <- as.list(Zaux)}
        else if(length(Zaux) > 1)
        {Zcat <- as.list(do.call(cbind, Zaux))}
        else
        {Zcat <- Zaux}
      } else {
        if(is.data.frame(Zaux))
        {Zcat <- as.list(Zaux)}
        else if(length(Zaux) == length(y)){
          Zcat <- list(Zaux)
          names(Zcat) <- paste0(znames[is.z.factor], levels(Z[is.z.factor][[1]])[-1])
        }
        else{
          Zcat <- NULL
          zCatnames <- NULL
          catcov <- FALSE
        }
      }
      rm(Zaux)
      if(!is.null(Zcat)){
        zCatnames <- names(Zcat) 
        names(Zcat) <- paste0("zCat", 1:length(Zcat))
      }
    } else {
      Zcat<- NULL
      zCatnames <- NULL
      catcov <- FALSE
    }
    
    
    # Transform Z due to continuous variables
    if(contcov){
      Zcont <- Z[nfactors == 0] # Si >1 continua, sera llista
      zContnames <- znames[!is.z.factor]
      if(!is.list(Zcont)){
        if(length(Zcont) == length(y)){
          Zcont <- list(Zcont)
          names(Zcont) <- paste0("zCont", 1:length(Zcont))
        }
        else{
          Zcont <- NULL
          zContnames <- NULL
          contcov <- FALSE
        }
      }
      if(!is.null(Zcont)){
        zContnames <- znames[!is.z.factor] # names(Zcont) 
        names(Zcont) <- paste0("zCont", 1:length(Zcont))
      }
    }
    else {
      Zcont<- NULL
      zContnames <- NULL
      contcov <- FALSE
    }
  } else {
    catcov <- contcov <- FALSE
    Zcat <- Zcont <- NULL
    zCatnames <- zContnames <- NULL
  }
  rm(Z)
  
  # Limits matrix
  if(is.null(lod) & !is.null(loq)) {
    LOD <- 0; LOQ <- loq
    } else if(!is.null(lod) & is.null(loq) & nlimitvars == 1) {
    LOQ <- max(X[, LimitVariables], na.rm = TRUE) + 
      abs(diff(range(X[,LimitVariables], na.rm = TRUE)))/2 * 100 
    LOD <- lod
  } else if(!is.null(lod) & !is.null(loq) & nlimitvars == 1)
  {
    LOD <- lod
    LOQ <- loq
  }
  
  if(nlimitvars == 1) 
    Lim <- as.matrix(data.frame(lower = rep(LOD, n), upper = rep(LOQ, n)))
  else {
    LOD <- lod
    xmax <- apply(x, 2, function(x) max(x, na.rm = TRUE))
    xdiffrange <- apply(x, 2, function(x) abs(diff(range(x, na.rm = TRUE))))
    LOQ <- xmax + xdiffrange / 2 * 100
    names(lod) <- names(LOD) <- names(LOQ)
    Lim <- mapply(function(lod, loq, n) 
      as.matrix(data.frame(lower = rep(lod, n), upper = rep(loq, n))), 
                  lod = as.list(LOD), 
                  loq = as.list(LOQ), 
                  n = n, SIMPLIFY = FALSE)
    names(Lim) <- paste0("Lim", 1:length(lod))
  }
  
  # bugsData for the model with 1 limited variable:
  if(nlimitvars == 1){
    if (logLimit) {  
      bugsData <- c(list(N = n,
                         censored = censored,
                         y = y, 
                         x = x, 
                         Lim = Lim, tau.a = 1, tau.b = 2), 
                    Zcat, Zcont)
    }   
    else {
      bugsData <- c(list(N = n,
                         censored = censored,
                         y = y, 
                         x = x, 
                         Lim = Lim), 
                    Zcat, Zcont)
    }
    
    # Initialization values needed
    ## set the missing values to random variates from U(0, lod) and U(loq, Inf)
    bugsInits <- list(list(x = as.numeric(rep(NA, length(x)))))
    
    
    if(!is.null(LimitVariables) & !is.null(loq)){
      nMissing <- sum(ind2, na.rm = TRUE)
      upper <- as.numeric(loq + abs(diff(range(x, na.rm = TRUE)))/2 * 10)
      bugsInits[[1]]$x[ind2] <- as.numeric(runif(nMissing, loq, upper))
      # x >= loq 
      # censored <- !censored
    }
         
    
    if(!is.null(LimitVariables) & !is.null(lod)){
      nMissing <- sum(ind1, na.rm = TRUE)
      bugsInits[[1]]$x[ind1] <- as.numeric(runif(nMissing, 0, lod))
      # x <= lod
    }    
  }
  
  
  # bugsData for the model with > 1 limited variable:
  if(nlimitvars > 1){
    x <- as.list(as.data.frame(x))
    names(x) <- paste0("x", 1:nlimitvars)
    
    if (logLimit) {  
      bugsData <- c(list(N = n,
                         y = y, 
                         tau.a = 1, tau.b = 2), 
                    x, Lim, as.list(as.data.frame(censored)),
                    Zcat, Zcont)
    }   
    else {
      bugsData <- c(list(N = n,
                         y = y), 
                    x, Lim, as.list(as.data.frame(censored)),
                    Zcat, Zcont)
    }
    
    
    # buugs <<- bugsData # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    # Initialization values needed
    ## set the missing values to random variates from U(0, lod) and U(loq, Inf)
    aux <- as.list(as.data.frame(replicate(nlimitvars, as.numeric(rep(NA, n)))))
    
    names(aux) <- names(x)
    bugsInits <- list(aux)
    rm(aux)
    
    nMissing <- apply(censored, 2, function(x) sum(!x))

    if(!is.null(LimitVariables) & !is.null(loq)){
      warning("For multiple arguments in 'lod', 'loq' argument will not be used!")

    if(!is.null(LimitVariables) & !is.null(lod)){
      for(i in 1:nlimitvars)
        bugsInits[[1]][paste0("x", i)][X[,LimitVariables[i]] <= lod[i]] <- as.numeric(runif(nMissing[i], 0, lod[i]))
    }
  }

}
  if(n.chains == 2){
    bugsInits[[n.chains]] <- bugsInits[[1]]
  } else if(n.chains > 1){
    bugsInits[2:n.chains] <- replicate(n.chains - 1, bugsInits[[1]], simplify = FALSE)
  }
  require(rjags)


  if (family == "gaussian" & !logLimit) {
    
    modelString <- paste("
                         model{
                         for( i in 1:N ) {
                         y[i] ~ dnorm( y.hat[i] , tau )
                         y.hat[i] <- intercept + ",  
                         if(nlimitvars == 1) {paste0("beta", LimitVariables, " * x[i] ")} else 
                           paste0("beta", LimitVariables, " * x", 1:nlimitvars, "[i]", collapse = " + "), 
                         if(contcov) paste("+", paste0(zContnames, " * zCont", 1:length(Zcont), "[i]", collapse = " + ")), 
                         if(catcov) paste("+", paste0(zCatnames, " * zCat", 1:length(Zcat), "[i]", collapse = " + ")), 
                         if(nlimitvars == 1) {
                           paste(paste0("
                          censored[i] ~ dinterval(x[i], Lim[i,])
                                  "),
                           paste0("
                          x[i] ~ dnorm(mux, taux)"), sep = "")
                           } else if(nlimitvars > 1) {
                           paste(paste0("
                                  censored", 1:nlimitvars, 
                                  "[i] ~ dinterval(x", 1:nlimitvars, "[i], Lim", 1:nlimitvars, "[i,])
                                  ", collapse = "\n"), 
                           paste0("
                                  x", 1:nlimitvars, 
                                  "[i] ~ dnorm(mux", 1:nlimitvars, 
                                  ", taux", 1:nlimitvars, ")", collapse = "\n"), sep = "")
                         }, "
                         }
                         
                         tau <- pow(sigma, -2)
                         ## prior
                         sigma ~ dunif(0, 100)
                         intercept ~ dnorm(0 , .0001)
                         ", paste0("beta", LimitVariables, " ~ dnorm(0 , .0001)", collapse = " \n"), "
                         ", if(contcov) paste0(zContnames," ~ dnorm(0, .0001)", collapse = "\n"), "
                         ", if(catcov) paste0(zCatnames," ~ dnorm(0, .0001)", collapse = "\n"), 
                         if(nlimitvars == 1){
                           paste("
                         mux ~ dnorm(0, .0001)
                         taux <- pow(sigma2, -2)
                         sigma2 ~ dunif(0, 100)")
                         } else if(nlimitvars > 1){
                           paste(
                             paste0("
                                    mux", 1:nlimitvars, " ~ dnorm(0, .0001)", collapse = "\n"),
                             paste0("
                                    taux", 1:nlimitvars, " <- pow(sigmax", 1:nlimitvars, ", -2)", collapse = "\n"),
                             paste0("
                                    sigmax", 1:nlimitvars, " ~ dunif(0, 100)", collapse = "\n"),
                             sep = "")
                         }, "
                         }", sep = "")
  
  } else if (family == "gaussian" & logLimit) {
    
    modelString <- paste("
                         model{
                         for( i in 1:N ) {
                         y[i] ~ dnorm( y.hat[i] , tau )
                         y.hat[i] <- intercept + ",  
                         if(nlimitvars == 1) {paste0("beta", LimitVariables, " * x[i] ")} else 
                           paste0("beta", LimitVariables, " * x", 1:nlimitvars, "[i]", collapse = " + "), 
                         if(contcov) paste("+", paste0(zContnames, " * zCont", 1:length(Zcont), "[i]", collapse = " + ")), 
                         if(catcov) paste("+", paste0(zCatnames, " * zCat", 1:length(Zcat), "[i]", collapse = " + ")), 
                         if(nlimitvars == 1) {
                           paste(paste0("
                                        censored[i] ~ dinterval(x[i], Lim[i,])
                                        "),
                                 paste0("
                                        x[i] ~ dlnorm(mux, taux)"), sep = "")
                         } else if(nlimitvars > 1) {
                           paste(paste0("
                                        censored", 1:nlimitvars, 
                                        "[i] ~ dinterval(x", 1:nlimitvars, "[i], Lim", 1:nlimitvars, "[i,])
                                        ", collapse = "\n"), 
                                 paste0("
                                        x", 1:nlimitvars, 
                                        "[i] ~ dlnorm(mux", 1:nlimitvars, 
                                        ", taux", 1:nlimitvars, ")", collapse = "\n"), sep = "")
                         }, "
                         }
                         
                         tau <- pow(sigma, -2)
                         sigma ~ dunif(0, 100)
                         intercept ~ dnorm(0 , .0001)
                         ", paste0("beta", LimitVariables, " ~ dnorm(0 , .0001)", collapse = " \n"), " 
                         ", if(contcov) paste0(zContnames," ~ dnorm(0, .0001)", collapse = "\n"), "
                         ", if(catcov) paste0(zCatnames," ~ dnorm(0, .0001)", collapse = "\n"), if(nlimitvars == 1){
                           paste("
                         mux ~ dnorm(0, .0001)
                         taux ~ dgamma(tau.a, tau.b)
                         sigmax ~ dunif(0, 100)")
                         } else if(nlimitvars > 1){
                           paste(
                             paste0("
                                    mux", 1:nlimitvars, " ~ dnorm(0, .0001)", collapse = "\n"),
                             paste0("
                                    taux", 1:nlimitvars, " ~ dgamma(tau.a, tau.b)", collapse = "\n"),
                             paste0("
                                    sigmax", 1:nlimitvars, " <- 1 / sqrt(taux", 1:nlimitvars, ")", collapse = "\n"),
                             sep = "")
                         }, "
                         }", sep = "")
  
  } else if (family == "binomial") {
    # http://www.johnmyleswhite.com/notebook/2010/08/20/using-jags-in-r-with-the-rjags-package/
    # 1/(1+exp(-y.hat[i]))
    load.module("glm")
    modelString <- paste("
                         model{
                         for( i in 1:N ) {
                         y[i] ~ dbern(  ilogit(y.hat[i]) ) 
                         y.hat[i] <- intercept + ",  
                         if(nlimitvars == 1) {paste0("beta", LimitVariables, " * x[i] ")} else 
                           paste0("beta", LimitVariables, " * x", 1:nlimitvars, "[i]", collapse = " + "), 
                         if(contcov) paste("+", paste0(zContnames, " * zCont", 1:length(Zcont), "[i]", collapse = " + ")), 
                         if(catcov) paste("+", paste0(zCatnames, " * zCat", 1:length(Zcat), "[i]", collapse = " + ")), 
                         if(nlimitvars == 1) {
                           paste(paste0("
                                        censored[i] ~ dinterval(x[i], Lim[i,])
                                        "),
                                 paste0("
                                        x[i] ~ dnorm(mux, taux)"), sep = "")
                         } else if(nlimitvars > 1) {
                           paste(paste0("
                                        censored", 1:nlimitvars, 
                                        "[i] ~ dinterval(x", 1:nlimitvars, "[i], Lim", 1:nlimitvars, "[i,])
                                        ", collapse = "\n"), 
                                 paste0("
                                        x", 1:nlimitvars, 
                                        "[i] ~ dnorm(mux", 1:nlimitvars, 
                                        ", taux", 1:nlimitvars, ")", collapse = "\n"), sep = "")
                         }, "
                         }

                         ## prior
                         intercept ~ dnorm(0 , .0001)
                         ", paste0("beta", LimitVariables, " ~ dnorm(0 , .0001)", collapse = " \n"), "
                         ", if(contcov) paste0(zContnames," ~ dnorm(0, .0001)", collapse = "\n"), "
                         ", if(catcov) paste0(zCatnames," ~ dnorm(0, .0001)", collapse = "\n"), 
                         if(nlimitvars == 1){
                           paste("
                                 mux ~ dnorm(0, .0001)
                                 taux <- pow(sigma2, -2)
                                 sigma2 ~ dunif(0, 100)")
                         } else if(nlimitvars > 1){
                           paste(
                             paste0("
                                    mux", 1:nlimitvars, " ~ dnorm(0, .0001)", collapse = "\n"),
                             paste0("
                                    taux", 1:nlimitvars, " <- pow(sigmax", 1:nlimitvars, ", -2)", collapse = "\n"),
                             paste0("
                                    sigmax", 1:nlimitvars, " ~ dunif(0, 100)", collapse = "\n"),
                             sep = "")
                         }, "
                         }", sep = "")
  }
  
  writeLines(modelString, con = "model.bug")
  gc(reset = TRUE, verbose = FALSE)
  
  # buugsinits <<- bugsInits # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  # buugsdata <<- bugsData # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  
  jagsModel <- jags.model(file = "model.bug", data = bugsData, 
                          inits = bugsInits, quiet = quiet, n.chains = n.chains,
                          n.adapt = n.adapt)
  gc(reset = TRUE, verbose = FALSE)
  
  update(jagsModel, n.iter = n.iter, progress.bar = progress.bar)
  gc(reset = TRUE, verbose = FALSE)
  
  if(n.chains > 1)
  {
    load.module("dic") # necessary for pD and deviance monitor
    dicSamples <- try(dic.samples(jagsModel, n.iter = n.iter, type = dic.type), TRUE)
    gc(reset = TRUE, verbose = FALSE)
  } else dicSamples <- NULL
  
  
  parameters <- c("intercept", paste0("beta", LimitVariables), zContnames, zCatnames)
  simSamples <- coda.samples(jagsModel, variable.names = parameters, n.iter = n.iter)
  gc(reset = TRUE, verbose = FALSE)
  
  
  yhat <-  summary(coda.samples(jagsModel, variable.names = c("y.hat"), n.iter = n.iter), 
                   quantiles = c(0.5))$quantiles[paste0("y.hat[", 1:n, "]")]
  gc(reset = TRUE, verbose = FALSE)
  
  stats <- summary(simSamples, quantiles = quantiles)
  ans <- list(simSamples = simSamples, summary = stats, model = jagsModel, dicSamples = dicSamples, 
              pred = as.vector(yhat), response = y, family = family, alpha.corrected = alpha.corrected)
  class(ans) <- "glmLimits"
  class(ans$simSamples) <- "mcmc.list"
  class(ans$model) <- "jags"
  gc(reset = TRUE)
  ans
  }

