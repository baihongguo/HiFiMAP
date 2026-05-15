#' @keywords internal
.Q_pval <- function(Q, lambda, method = "davies") {
  if(method == "davies") {
    tmp <- suppressWarnings(CompQuadForm::davies(q = Q, lambda = lambda, acc = 1e-6))
    pval <- tmp$Qq
    if((tmp$ifault > 0) | (pval <= 1e-5) | (pval >= 1)) method <- "kuonen"
  }
  if(method == "kuonen") {
    pval <- .pKuonen(x = Q, lambda = lambda)
    if(is.na(pval)) method <- "liu"
  }
  if(method == "liu") pval <- CompQuadForm::liu(q = Q, lambda = lambda)
  return(pval)
}

#when davies method loses accuarcy 
#' @keywords internal
.pKuonen <- function (x, lambda, delta = rep(0, length(lambda)), df = rep(1, length(lambda)))
{
  delta <- delta[lambda != 0]
  df <- df[lambda != 0]
  lambda <- lambda[lambda != 0]
  if(length(lambda) != length(delta)) stop("Error: inconsistent length in lambda and delta!")
  if(length(lambda) != length(df)) stop("Error: inconsistent length in lambda and df!")
  if (length(lambda) == 1) {
    pchisq(x/lambda, df = df, ncp = delta, lower.tail = FALSE)
  }
  d <- max(lambda)
  lambda <- lambda/d
  x <- x/d
  k0 <- function(zeta) {
    -sum(df * log(1 - 2 * zeta * lambda))/2 + sum((delta * lambda *
                                                     zeta)/(1 - 2 * zeta * lambda))
  }
  kprime0 <- function(zeta) {
    sapply(zeta, function(zz) {
      sum(((delta + df) * lambda)/(1 - 2 * zz * lambda) + 2 * (delta *
                                                                 zz * lambda^2)/(1 - 2 * zz * lambda)^2)
    })
  }
  kpprime0 <- function(zeta) {
    sum((2 * (2 * delta + df) * lambda^2)/(1 - 2 * zeta * lambda)^2 + 8 *
          delta * zeta * lambda^3/(1 - 2 * zeta * lambda)^3)
  }
  if (any(lambda < 0)) {
    lmin <- max(1/(2 * lambda[lambda < 0])) * 0.99999
  }
  else if (x > sum((df+delta)*lambda)) {
    lmin <- -0.01
  }
  else {
    lmin <- -length(lambda)*max(df+delta)/(2 * x)
  }
  lmax <- min(1/(2 * lambda[lambda > 0])) * 0.99999
  hatzeta <- uniroot(function(zeta) kprime0(zeta) - x, lower = lmin,
                     upper = lmax, tol = 1e-08)$root
  w <- sign(hatzeta) * sqrt(2 * (hatzeta * x - k0(hatzeta)))
  v <- hatzeta * sqrt(kpprime0(hatzeta))
  if (abs(hatzeta) < 1e-04)
    NA
  else pnorm(w + log(v/w)/w, lower.tail = FALSE)
}


