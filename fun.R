#functions

#########----------For the probabilities------------###########
ANODev <- function(model, m){
  anova.model <- anova(model)
  dev <- anova.model$`Resid. Dev`  
  pvalues <- rep(0, m)
  rownames <- rep(0, m)
  likelihoodstats.simple <- rep(0, m)
  
  for (i in 1:m){
    
    pvalues[i] <- 1 - pchisq(dev[i*2 - 1] - dev[i*2 + 1], 2)
    rownames[i] <- paste("Harmonic:", i)
    likelihoodstats.simple[i] <- dev[i*2 - 1] - dev[i*2 + 1]
    
  }
  cbind(rownames, likelihoodstats.simple, pvalues)
}

normalresid <- function(model){
  pres <- residuals(model, type = "pearson")
  dres <- residuals(model, type = "deviance")
  normalresid <- rep(0, length(pres))
  
  normalresid <- dres + dres^(-1)*log(pres/dres)
  
  normalresid
}

#########----------For the amounts------------###########

FDev <- function(model, m){
  anova.model <- anova(model)
  dev <- anova.model$`Resid. Dev`  
  pvalues <- rep(0, m)
  rownames <- rep(0, m)
  F.models <- rep(0, m)
  dispersion <- summary(model)$dispersion
  degf <- summary(model)$df[2]
  
  for (i in 1:m){
    
    F.models[i] <- dev[i*2 - 1] - dev[i*2 + 1]
    pvalues[i] <- 1 - pf(F.models[i]/dispersion, 2, degf)
    rownames[i] <- paste("Harmonic:", i)
    
    
  }
  cbind(rownames, F.models, pvalues)
}