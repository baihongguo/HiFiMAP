glmmkin2randomvec <- function(obj, Z = NULL, N.randomvec = 1000, group.idx=NULL, cluster.idx=NULL, robust = FALSE) {
  
  if(class(obj) != "glmmkin") stop("Error: \"obj\" must be a class glmmkin object.")
  N <- length(obj$id_include)
  random.vectors <- matrix(rnorm(N*N.randomvec), nrow=N, ncol=N.randomvec)
  
  if(!is.null(obj$P) && !robust) {
    eig <- eigen(obj$P, symmetric = TRUE)
    random.vectors <- tcrossprod(eig$vectors, t(random.vectors * sqrt(pmax(eig$values, 0))))
    rm(eig)
  } else {
    if(obj$n.groups != 1 && (is.null(group.idx) || !all.equal(seq_len(obj$n.groups), sort(unique(group.idx))))) stop("Error: heteroscedastic linear mixed models should include a valid group.idx argument.")
    if(is.null(group.idx)) group.idx <- rep(1, N)
    
    # -----------------------------------------------------------------------
    # UNIVERSAL GLMM WORKING VARIANCE CALCULATION
    # -----------------------------------------------------------------------
    # Extract the family from obj directly, or reconstruct it from obj$call
    fit_family <- obj$family
    if (is.null(fit_family) && !is.null(obj$call$family)) {
      fit_family <- eval(obj$call$family)
    }
    
    # Fallback to Gaussian if the family is completely missing from the object
    if (is.null(fit_family)) {
      fit_family <- gaussian()
    }
    
    # Universal formula: works automatically for Gaussian, Binomial, Poisson, etc.
    eta <- obj$linear.predictors
    mu <- obj$fitted.values
    mu.eta <- fit_family$mu.eta(eta)
    sqrtW <- mu.eta / sqrt(fit_family$variance(mu))
    # -----------------------------------------------------------------------
    
    if(!robust) {
      # Apply the working variance weights (sqrtW) to scale the random vectors
      random.vectors <- (sqrt(obj$theta[group.idx]) / sqrtW) * random.vectors
    } else {
      res <- as.numeric(obj$Y - tcrossprod(obj$X, t(obj$coefficient)))
      if(is.null(cluster.idx)) random.vectors <- random.vectors * res
      else random.vectors <- random.vectors[match(cluster.idx,unique(cluster.idx)),] * res
    }
    
    if(!is.null(Z)) {
      if(class(Z) != "list") stop("Error: \"Z\" must be a list of matrices.")
      if(length(Z) != length(obj$theta) - obj$n.groups) stop("Error: number of matrices in \"Z\" does not match the number of variance components in \"obj\".")
      for(i in 1:length(Z)) {
        if(nrow(Z[[i]]) != N) stop("Error: \"Z\" matrix ", i, " is not compatible in sample size with \"obj\".")
        p <- ncol(Z[[i]])
        if(obj$theta[i+obj$n.groups] < 0) stop("Error: negative variance component estimates are not allowed.")
        if(obj$theta[i+obj$n.groups] == 0) next
        random.vectors2 <- matrix(rnorm(p*N.randomvec), nrow=N.randomvec, ncol=p)
        random.vectors <- random.vectors + sqrt(obj$theta[i+obj$n.groups]) * tcrossprod(Z[[i]], random.vectors2)
      }
    }
    
    if(!is.null(obj$P)) random.vectors <- crossprod(obj$P, random.vectors)
    else random.vectors <- crossprod(obj$Sigma_i, random.vectors) - tcrossprod(obj$Sigma_iX, tcrossprod(crossprod(random.vectors, obj$Sigma_iX), obj$cov))
  }
  
  out <- list(
    theta = obj$theta, 
    scaled.residuals = obj$scaled.residuals, 
    random.vectors = as.matrix(random.vectors), 
    ncovar = ncol(obj$X), 
    RSigma_R = as.matrix(forceSymmetric(crossprod(random.vectors, solve(obj$Sigma_i, random.vectors)))), 
    id_include = obj$id_include
  )
  class(out) <- "glmmkin.randomvec"
  return(out)
}

