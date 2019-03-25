############################################################################
## Relevant functions required for the Probability of Agreement analysis  ##
############################################################################

## Author: Nathaniel T. Stevens
##         Assistant Professor, University of Waterloo
##   Date: March 2019 

#############################
## Log-Likelihood Function ##
#############################
log.likelihood <- function(param, n, r1vec, r2vec, data){
  # Calculate the log-likelihood contribution for each subject and then sum them all
  mu <- param[1]
  alpha <- param[2]
  beta <- param[3]
  sp <- param[4]
  s1 <- param[5]
  s2 <- param[6]
  l <- rep(0, n)
  for(i in 1:n){
    r1 <- r1vec[i]
    r2 <- r2vec[i]
    Yi1 <- subset(data, subset = (data$Subject==i & data$MS==1))$Measurements
    Yi2 <- subset(data, subset = (data$Subject==i & data$MS==2))$Measurements
    A <- sum((Yi1-rep(mu, r1))^2)
    B <- sum(Yi1-rep(mu, r1))^2
    C <- sum((Yi2-rep(alpha + beta*mu, r2))^2)
    D <- sum(Yi2-rep(alpha + beta*mu, r2))^2
    E <- sum(Yi1-rep(mu, r1)) * sum(Yi2-rep(alpha + beta*mu, r2))
    # Log-likelihood contribution for subject i
    l[i] <- -(r1 / 2 + r2 / 2) * log(2 * pi) - r1 * log(s1) - r2 * log(s2) - log1p(sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) / 2 - A / s1 ^ 2 / 2 + sp ^ 2 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) * B / s1 ^ 4 / 2 - C / s2 ^ 2 / 2 + sp ^ 2 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) * beta ^ 2 * D / s2 ^ 4 / 2 + sp ^ 2 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) * beta * E / s1 ^ 2 / s2 ^ 2
  }

  # The log-likelihood function:
  l <- sum(l)
  return(-l) # Return negative because optim() is a minimizer
}

##############################
## Log-Likelihood Gradients ##
##############################
gradient <- function(param, n, r1vec, r2vec, data){
  # Calculate the gradients for each parameter / subject combination and sum over all subjects
  # for the gradients for each parameter.
  mu <- param[1]
  alpha <- param[2]
  beta <- param[3]
  sp <- param[4]
  s1 <- param[5]
  s2 <- param[6]
  g <- matrix(0, nrow = 6, ncol = n)
  for(i in 1:n){
    r1 <- r1vec[i]
    r2 <- r2vec[i]
    Yi1 <- subset(data, subset = (data$Subject==i & data$MS==1))$Measurements
    Yi2 <- subset(data, subset = (data$Subject==i & data$MS==2))$Measurements
    A <- sum((Yi1-rep(mu, r1))^2)
    C <- sum((Yi2-rep(alpha + beta*mu, r2))^2)
    F <- sum(Yi1-rep(mu, r1))
    G <- sum(Yi2-rep(alpha + beta*mu, r2))
    B <- F^2
    D <- G^2
    E <- F*G
    
    # Derivative with respect to mu:
    g[1,i] <- F / s1 ^ 2 - r1 * sp ^ 2 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) * F / s1 ^ 4 + beta * G / s2 ^ 2 - r2 * sp ^ 2 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) * beta ^ 3 * G / s2 ^ 4 - r1 * sp ^ 2 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) * beta * G / s1 ^ 2 / s2 ^ 2 - r2 * sp ^ 2 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) * beta ^ 2 * F / s1 ^ 2 / s2 ^ 2
    # Derivative with respect to alpha:
    g[2,i] <- G / s2 ^ 2 - r2 * sp ^ 2 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) * beta ^ 2 * G / s2 ^ 4 - r2 * sp ^ 2 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) * beta * F / s1 ^ 2 / s2 ^ 2
    # Derivative with respect to beta:
    g[3,i] <- -beta * r2 * sp ^ 2 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) / s2 ^ 2 - beta * r2 * sp ^ 4 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) ^ 2 / s2 ^ 2 * B / s1 ^ 4 + mu * G / s2 ^ 2 + (-2 * G * sp ^ 2 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) * beta ^ 2 * mu * r2 - 2 * D * beta ^ 3 * r2 * sp ^ 4 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) ^ 2 / s2 ^ 2 + 2 * D * sp ^ 2 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) * beta) / s2 ^ 4 / 2 + (-F * sp ^ 2 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) * beta * mu * r2 - 2 * E * beta ^ 2 * r2 * sp ^ 4 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) ^ 2 / s2 ^ 2 + E * sp ^ 2 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2))) / s1 ^ 2 / s2 ^ 2
    # Derivative with respect to sp:
    g[4,i] <- -sp * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2) / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) + sp / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) * B / s1 ^ 4 - sp ^ 3 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) ^ 2 * B / s1 ^ 4 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2) + sp / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) * beta ^ 2 * D / s2 ^ 4 - sp ^ 3 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) ^ 2 * beta ^ 2 * D / s2 ^ 4 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2) + 2 * sp / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) * beta * E / s1 ^ 2 / s2 ^ 2 - 2 * sp ^ 3 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) ^ 2 * beta * E / s1 ^ 2 / s2 ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)
    # Derivative with respect to s1:
    g[5,i] <- -r1 / s1 + sp ^ 2 * r1 / s1 ^ 3 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) + A / s1 ^ 3 + sp ^ 4 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) ^ 2 * B / s1 ^ 7 * r1 - 2 * sp ^ 2 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) * B / s1 ^ 5 + sp ^ 4 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) ^ 2 * beta ^ 2 * D / s2 ^ 4 * r1 / s1 ^ 3 + 2 * sp ^ 4 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) ^ 2 * beta * E / s1 ^ 5 / s2 ^ 2 * r1 - 2 * sp ^ 2 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) * beta * E / s1 ^ 3 / s2 ^ 2
    # Derivative with respect to s2:
    g[6,i] <- -r2 / s2 + sp ^ 2 * r2 * beta ^ 2 / s2 ^ 3 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) + sp ^ 4 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) ^ 2 * B / s1 ^ 4 * r2 * beta ^ 2 / s2 ^ 3 + C / s2 ^ 3 + sp ^ 4 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) ^ 2 * beta ^ 4 * D / s2 ^ 7 * r2 - 2 * sp ^ 2 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) * beta ^ 2 * D / s2 ^ 5 + 2 * sp ^ 4 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) ^ 2 * beta ^ 3 * E / s1 ^ 2 / s2 ^ 5 * r2 - 2 * sp ^ 2 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) * beta * E / s1 ^ 2 / s2 ^ 3
  }
  
  # The gradients:
  g <- apply(X = g, MARGIN = 1, FUN = sum)
  return(-g) # Return negative because optim() is a minimizer
}

###############################
## Fisher Information Matrix ##
###############################
fisher.info <- function(param, n, r1vec, r2vec){
  # Calculate the Fisher Information Matrix associated with the 6 parameters
  # Expected values of the various sums of squares
  mu <- param[1]
  alpha <- param[2]
  beta <- param[3]
  sp <- param[4]
  s1 <- param[5]
  s2 <- param[6]
  dldmu2 <- rep(0, n); dldmualpha <- rep(0, n); dldmubeta <- rep(0, n); dldmusp <- rep(0, n); dldmus1 <- rep(0, n); dldmus2 <- rep(0, n); dldalpha2 <- rep(0, n); dldalphabeta <- rep(0, n); dldalphasp <- rep(0, n); dldalphas1 <- rep(0, n); dldalphas2 <- rep(0, n); dldbeta2 <- rep(0, n); dldbetasp <- rep(0, n); dldbetas1 <- rep(0, n); dldbetas2 <- rep(0, n); dldsp2 <- rep(0, n); dldsps1 <- rep(0, n); dldsps2 <- rep(0, n); dlds12 <- rep(0, n); dlds1s2 <- rep(0, n); dlds22 <- rep(0, n);
  for(i in 1:n){
    r1 <- r1vec[i]
    r2 <- r2vec[i]
    A <- r1*(sp^2+s1^2)
    B <- r1*(r1*sp^2+s1^2)
    C <- r2*((beta^2)*(sp^2)+s2^2)
    D <- r2*(r2*(beta^2)*(sp^2)+s2^2)
    E <- (r1*r2)*beta*sp^2
    F <- 0
    G <- 0
    # Second partial derivatives for each subject
    dldmu2[i] <- -r1 / s1 ^ 2 + r1 ^ 2 * sp ^ 2 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) / s1 ^ 4 - r2 * beta ^ 2 / s2 ^ 2 + r2 ^ 2 * sp ^ 2 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) * beta ^ 4 / s2 ^ 4 + 2 * r1 * r2 * sp ^ 2 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) * beta ^ 2 / s1 ^ 2 / s2 ^ 2
    dldmualpha[i] <- -r2 * beta / s2 ^ 2 + r2 ^ 2 * sp ^ 2 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) * beta ^ 3 / s2 ^ 4 + r1 * r2 * sp ^ 2 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) * beta / s1 ^ 2 / s2 ^ 2
    dldmubeta[i] <- 2 * r1 * beta * r2 * sp ^ 4 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) ^ 2 / s2 ^ 2 * F / s1 ^ 4 + (-beta * mu * r2 + G) / s2 ^ 2 - (-sp ^ 2 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) * beta ^ 3 * mu * r2 ^ 2 - 2 * G * beta ^ 4 * r2 ^ 2 * sp ^ 4 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) ^ 2 / s2 ^ 2 + 3 * G * sp ^ 2 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) * beta ^ 2 * r2) / s2 ^ 4 - (-sp ^ 2 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) * beta * mu * r1 * r2 - 2 * G * beta ^ 2 * r2 * sp ^ 4 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) ^ 2 / s2 ^ 2 * r1 + G * sp ^ 2 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) * r1) / s1 ^ 2 / s2 ^ 2 - (-2 * F * beta ^ 3 * r2 ^ 2 * sp ^ 4 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) ^ 2 / s2 ^ 2 + 2 * F * sp ^ 2 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) * beta * r2) / s1 ^ 2 / s2 ^ 2
    dldmusp[i] <- -2 * r1 * sp / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) * F / s1 ^ 4 + 2 * r1 * sp ^ 3 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) ^ 2 * F / s1 ^ 4 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2) - 2 * r2 * sp / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) * beta ^ 3 * G / s2 ^ 4 + 2 * r2 * sp ^ 3 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) ^ 2 * beta ^ 3 * G / s2 ^ 4 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2) - 2 * r1 * sp / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) * beta * G / s1 ^ 2 / s2 ^ 2 + 2 * r1 * sp ^ 3 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) ^ 2 * beta * G / s1 ^ 2 / s2 ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2) - 2 * r2 * sp / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) * beta ^ 2 * F / s1 ^ 2 / s2 ^ 2 + 2 * r2 * sp ^ 3 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) ^ 2 * beta ^ 2 * F / s1 ^ 2 / s2 ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)
    dldmus1[i] <- -2 * F / s1 ^ 3 - 2 * r1 ^ 2 * sp ^ 4 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) ^ 2 * F / s1 ^ 7 + 4 * r1 * sp ^ 2 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) * F / s1 ^ 5 - 2 * r2 * sp ^ 4 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) ^ 2 * beta ^ 3 * G / s2 ^ 4 * r1 / s1 ^ 3 - 2 * r1 ^ 2 * sp ^ 4 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) ^ 2 * beta * G / s1 ^ 5 / s2 ^ 2 + 2 * r1 * sp ^ 2 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) * beta * G / s1 ^ 3 / s2 ^ 2 - 2 * r2 * sp ^ 4 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) ^ 2 * beta ^ 2 * F / s1 ^ 5 / s2 ^ 2 * r1 + 2 * r2 * sp ^ 2 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) * beta ^ 2 * F / s1 ^ 3 / s2 ^ 2
    dldmus2[i] <- -2 * r1 * sp ^ 4 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) ^ 2 * F / s1 ^ 4 * r2 * beta ^ 2 / s2 ^ 3 - 2 * beta * G / s2 ^ 3 - 2 * r2 ^ 2 * sp ^ 4 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) ^ 2 * beta ^ 5 * G / s2 ^ 7 + 4 * r2 * sp ^ 2 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) * beta ^ 3 * G / s2 ^ 5 - 2 * r1 * sp ^ 4 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) ^ 2 * beta ^ 3 * G / s1 ^ 2 / s2 ^ 5 * r2 + 2 * r1 * sp ^ 2 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) * beta * G / s1 ^ 2 / s2 ^ 3 - 2 * r2 ^ 2 * sp ^ 4 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) ^ 2 * beta ^ 4 * F / s1 ^ 2 / s2 ^ 5 + 2 * r2 * sp ^ 2 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) * beta ^ 2 * F / s1 ^ 2 / s2 ^ 3
    dldalpha2[i] <- -r2 / s2 ^ 2 + r2 ^ 2 * sp ^ 2 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) * beta ^ 2 / s2 ^ 4
    dldalphabeta[i] <- -r2 * mu / s2 ^ 2 - (-sp ^ 2 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) * beta ^ 2 * mu * r2 ^ 2 - 2 * G * beta ^ 3 * r2 ^ 2 * sp ^ 4 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) ^ 2 / s2 ^ 2 + 2 * G * sp ^ 2 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) * beta * r2) / s2 ^ 4 - (-2 * F * beta ^ 2 * r2 ^ 2 * sp ^ 4 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) ^ 2 / s2 ^ 2 + F * sp ^ 2 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) * r2) / s1 ^ 2 / s2 ^ 2
    dldalphasp[i] <- -2 * r2 * sp / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) * beta ^ 2 * G / s2 ^ 4 + 2 * r2 * sp ^ 3 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) ^ 2 * beta ^ 2 * G / s2 ^ 4 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2) - 2 * r2 * sp / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) * beta * F / s1 ^ 2 / s2 ^ 2 + 2 * r2 * sp ^ 3 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) ^ 2 * beta * F / s1 ^ 2 / s2 ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)
    dldalphas1[i] <- -2 * r2 * sp ^ 4 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) ^ 2 * beta ^ 2 * G / s2 ^ 4 * r1 / s1 ^ 3 - 2 * r2 * sp ^ 4 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) ^ 2 * beta * F / s1 ^ 5 / s2 ^ 2 * r1 + 2 * r2 * sp ^ 2 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) * beta * F / s1 ^ 3 / s2 ^ 2
    dldalphas2[i] <- -2 * G / s2 ^ 3 - 2 * r2 ^ 2 * sp ^ 4 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) ^ 2 * beta ^ 4 * G / s2 ^ 7 + 4 * r2 * sp ^ 2 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) * beta ^ 2 * G / s2 ^ 5 - 2 * r2 ^ 2 * sp ^ 4 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) ^ 2 * beta ^ 3 * F / s1 ^ 2 / s2 ^ 5 + 2 * r2 * sp ^ 2 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) * beta * F / s1 ^ 2 / s2 ^ 3
    dldbeta2[i] <- -(-2 * beta ^ 2 * r2 ^ 2 * sp ^ 4 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) ^ 2 / s2 ^ 2 + sp ^ 2 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) * r2) / s2 ^ 2 - r2 * sp ^ 2 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) * (-4 * beta ^ 2 * r2 * sp ^ 4 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) ^ 2 / s2 ^ 2 + sp ^ 2 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2))) / s2 ^ 2 * B / s1 ^ 4 - r2 * mu ^ 2 / s2 ^ 2 + (8 * G * beta ^ 3 * r2 ^ 2 * sp ^ 4 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) ^ 2 / s2 ^ 2 * mu - 2 * D * beta ^ 2 * r2 * sp ^ 2 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) * (-4 * beta ^ 2 * r2 * sp ^ 4 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) ^ 2 / s2 ^ 2 + sp ^ 2 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2))) / s2 ^ 2 - 8 * D * beta ^ 2 * r2 * sp ^ 4 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) ^ 2 / s2 ^ 2 - 8 * G * sp ^ 2 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) * beta * mu * r2 + 2 * D * sp ^ 2 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) + 2 * sp ^ 2 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) * beta ^ 2 * mu ^ 2 * r2 ^ 2) / s2 ^ 4 / 2 + (4 * F * beta ^ 2 * r2 ^ 2 * sp ^ 4 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) ^ 2 / s2 ^ 2 * mu - 2 * E * beta * r2 * sp ^ 2 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) * (-4 * beta ^ 2 * r2 * sp ^ 4 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) ^ 2 / s2 ^ 2 + sp ^ 2 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2))) / s2 ^ 2 - 4 * E * beta * r2 * sp ^ 4 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) ^ 2 / s2 ^ 2 - 2 * r2 * sp ^ 2 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) * mu * F) / s1 ^ 2 / s2 ^ 2
    dldbetasp[i] <- -2 * beta * r2 * sp / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) / s2 ^ 2 + 2 * beta * r2 * sp ^ 3 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) ^ 2 / s2 ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2) - 4 * beta * r2 * sp ^ 3 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) ^ 2 / s2 ^ 2 * B / s1 ^ 4 + 4 * beta * r2 * sp ^ 5 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) ^ 3 / s2 ^ 2 * B / s1 ^ 4 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2) + (-4 * G * sp / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) * beta ^ 2 * mu * r2 + 4 * G * sp ^ 3 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) ^ 2 * beta ^ 2 * mu * r2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2) - 8 * D * beta ^ 3 * r2 * sp ^ 3 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) ^ 2 / s2 ^ 2 + 8 * D * beta ^ 3 * r2 * sp ^ 5 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) ^ 3 / s2 ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2) + 4 * D * sp / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) * beta - 4 * D * sp ^ 3 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) ^ 2 * beta * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) / s2 ^ 4 / 2 + (-2 * F * sp / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) * beta * mu * r2 + 2 * F * sp ^ 3 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) ^ 2 * beta * mu * r2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2) - 8 * E * beta ^ 2 * r2 * sp ^ 3 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) ^ 2 / s2 ^ 2 + 8 * E * beta ^ 2 * r2 * sp ^ 5 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) ^ 3 / s2 ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2) + 2 * E * sp / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) - 2 * E * sp ^ 3 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) / s1 ^ 2 / s2 ^ 2
    dldbetas1[i] <- -2 * beta * r2 * sp ^ 4 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) ^ 2 / s2 ^ 2 * r1 / s1 ^ 3 - 4 * beta * r2 * sp ^ 6 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) ^ 3 / s2 ^ 2 * B / s1 ^ 7 * r1 + 4 * beta * r2 * sp ^ 4 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) ^ 2 / s2 ^ 2 * B / s1 ^ 5 + (-4 * G * sp ^ 4 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) ^ 2 * beta ^ 2 * mu * r2 * r1 / s1 ^ 3 - 8 * D * beta ^ 3 * r2 * sp ^ 6 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) ^ 3 / s2 ^ 2 * r1 / s1 ^ 3 + 4 * D * sp ^ 4 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) ^ 2 * beta * r1 / s1 ^ 3) / s2 ^ 4 / 2 + (-2 * F * sp ^ 4 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) ^ 2 * beta * mu * r2 * r1 / s1 ^ 3 - 8 * E * beta ^ 2 * r2 * sp ^ 6 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) ^ 3 / s2 ^ 2 * r1 / s1 ^ 3 + 2 * E * sp ^ 4 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) ^ 2 * r1 / s1 ^ 3) / s1 ^ 2 / s2 ^ 2 - 2 * (-F * sp ^ 2 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) * beta * mu * r2 - 2 * E * beta ^ 2 * r2 * sp ^ 4 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) ^ 2 / s2 ^ 2 + E * sp ^ 2 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2))) / s1 ^ 3 / s2 ^ 2
    dldbetas2[i] <- -2 * beta ^ 3 * r2 ^ 2 * sp ^ 4 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) ^ 2 / s2 ^ 5 + 2 * beta * r2 * sp ^ 2 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) / s2 ^ 3 - 4 * beta ^ 3 * r2 ^ 2 * sp ^ 6 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) ^ 3 / s2 ^ 5 * B / s1 ^ 4 + 2 * beta * r2 * sp ^ 4 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) ^ 2 / s2 ^ 3 * B / s1 ^ 4 - 2 * mu * G / s2 ^ 3 + (-4 * G * sp ^ 4 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) ^ 2 * beta ^ 4 * mu * r2 ^ 2 / s2 ^ 3 - 8 * D * beta ^ 5 * r2 ^ 2 * sp ^ 6 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) ^ 3 / s2 ^ 5 + 8 * D * beta ^ 3 * r2 * sp ^ 4 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) ^ 2 / s2 ^ 3) / s2 ^ 4 / 2 - 2 * (-2 * G * sp ^ 2 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) * beta ^ 2 * mu * r2 - 2 * D * beta ^ 3 * r2 * sp ^ 4 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) ^ 2 / s2 ^ 2 + 2 * D * sp ^ 2 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) * beta) / s2 ^ 5 + (-2 * F * sp ^ 4 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) ^ 2 * beta ^ 3 * mu * r2 ^ 2 / s2 ^ 3 - 8 * E * beta ^ 4 * r2 ^ 2 * sp ^ 6 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) ^ 3 / s2 ^ 5 + 6 * E * beta ^ 2 * r2 * sp ^ 4 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) ^ 2 / s2 ^ 3) / s1 ^ 2 / s2 ^ 2 - 2 * (-F * sp ^ 2 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) * beta * mu * r2 - 2 * E * beta ^ 2 * r2 * sp ^ 4 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) ^ 2 / s2 ^ 2 + E * sp ^ 2 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2))) / s1 ^ 2 / s2 ^ 3
    dldsp2[i] <- -(r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2) / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) + 2 * sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2) ^ 2 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) ^ 2 + 1 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) * B / s1 ^ 4 - 5 * sp ^ 2 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) ^ 2 * B / s1 ^ 4 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2) + 4 * sp ^ 4 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) ^ 3 * B / s1 ^ 4 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2) ^ 2 + 1 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) * beta ^ 2 * D / s2 ^ 4 - 5 * sp ^ 2 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) ^ 2 * beta ^ 2 * D / s2 ^ 4 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2) + 4 * sp ^ 4 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) ^ 3 * beta ^ 2 * D / s2 ^ 4 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2) ^ 2 + 2 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) * beta * E / s1 ^ 2 / s2 ^ 2 - 10 * sp ^ 2 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) ^ 2 * beta * E / s1 ^ 2 / s2 ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2) + 8 * sp ^ 4 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) ^ 3 * beta * E / s1 ^ 2 / s2 ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2) ^ 2
    dldsps1[i] <- 2 * sp * r1 / s1 ^ 3 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) - 2 * sp ^ 3 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2) / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) ^ 2 * r1 / s1 ^ 3 + 4 * sp ^ 3 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) ^ 2 * B / s1 ^ 7 * r1 - 4 * sp / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) * B / s1 ^ 5 - 4 * sp ^ 5 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) ^ 3 * B / s1 ^ 7 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2) * r1 + 4 * sp ^ 3 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) ^ 2 * B / s1 ^ 5 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2) + 4 * sp ^ 3 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) ^ 2 * beta ^ 2 * D / s2 ^ 4 * r1 / s1 ^ 3 - 4 * sp ^ 5 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) ^ 3 * beta ^ 2 * D / s2 ^ 4 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2) * r1 / s1 ^ 3 + 8 * sp ^ 3 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) ^ 2 * beta * E / s1 ^ 5 / s2 ^ 2 * r1 - 4 * sp / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) * beta * E / s1 ^ 3 / s2 ^ 2 - 8 * sp ^ 5 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) ^ 3 * beta * E / s1 ^ 5 / s2 ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2) * r1 + 4 * sp ^ 3 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) ^ 2 * beta * E / s1 ^ 3 / s2 ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)
    dldsps2[i] <- 2 * sp * r2 * beta ^ 2 / s2 ^ 3 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) - 2 * sp ^ 3 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2) / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) ^ 2 * r2 * beta ^ 2 / s2 ^ 3 + 4 * sp ^ 3 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) ^ 2 * B / s1 ^ 4 * r2 * beta ^ 2 / s2 ^ 3 - 4 * sp ^ 5 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) ^ 3 * B / s1 ^ 4 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2) * r2 * beta ^ 2 / s2 ^ 3 + 4 * sp ^ 3 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) ^ 2 * beta ^ 4 * D / s2 ^ 7 * r2 - 4 * sp / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) * beta ^ 2 * D / s2 ^ 5 - 4 * sp ^ 5 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) ^ 3 * beta ^ 4 * D / s2 ^ 7 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2) * r2 + 4 * sp ^ 3 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) ^ 2 * beta ^ 2 * D / s2 ^ 5 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2) + 8 * sp ^ 3 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) ^ 2 * beta ^ 3 * E / s1 ^ 2 / s2 ^ 5 * r2 - 4 * sp / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) * beta * E / s1 ^ 2 / s2 ^ 3 - 8 * sp ^ 5 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) ^ 3 * beta ^ 3 * E / s1 ^ 2 / s2 ^ 5 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2) * r2 + 4 * sp ^ 3 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) ^ 2 * beta * E / s1 ^ 2 / s2 ^ 3 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)
    dlds12[i] <- r1 / s1 ^ 2 - 3 * sp ^ 2 * r1 / s1 ^ 4 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) + 2 * sp ^ 4 * r1 ^ 2 / s1 ^ 6 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) ^ 2 - 3 * A / s1 ^ 4 + 4 * sp ^ 6 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) ^ 3 * B / s1 ^ 10 * r1 ^ 2 - 11 * sp ^ 4 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) ^ 2 * B / s1 ^ 8 * r1 + 10 * sp ^ 2 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) * B / s1 ^ 6 + 4 * sp ^ 6 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) ^ 3 * beta ^ 2 * D / s2 ^ 4 * r1 ^ 2 / s1 ^ 6 - 3 * sp ^ 4 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) ^ 2 * beta ^ 2 * D / s2 ^ 4 * r1 / s1 ^ 4 + 8 * sp ^ 6 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) ^ 3 * beta * E / s1 ^ 8 / s2 ^ 2 * r1 ^ 2 - 14 * sp ^ 4 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) ^ 2 * beta * E / s1 ^ 6 / s2 ^ 2 * r1 + 6 * sp ^ 2 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) * beta * E / s1 ^ 4 / s2 ^ 2
    dlds1s2[i] <- 2 * sp ^ 4 * r1 / s1 ^ 3 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) ^ 2 * r2 * beta ^ 2 / s2 ^ 3 + 4 * sp ^ 6 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) ^ 3 * B / s1 ^ 7 * r1 * r2 * beta ^ 2 / s2 ^ 3 - 4 * sp ^ 4 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) ^ 2 * B / s1 ^ 5 * r2 * beta ^ 2 / s2 ^ 3 + 4 * sp ^ 6 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) ^ 3 * beta ^ 4 * D / s2 ^ 7 * r1 / s1 ^ 3 * r2 - 4 * sp ^ 4 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) ^ 2 * beta ^ 2 * D / s2 ^ 5 * r1 / s1 ^ 3 + 8 * sp ^ 6 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) ^ 3 * beta ^ 3 * E / s1 ^ 5 / s2 ^ 5 * r1 * r2 - 4 * sp ^ 4 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) ^ 2 * beta * E / s1 ^ 5 / s2 ^ 3 * r1 - 4 * sp ^ 4 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) ^ 2 * beta ^ 3 * E / s1 ^ 3 / s2 ^ 5 * r2 + 4 * sp ^ 2 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) * beta * E / s1 ^ 3 / s2 ^ 3
    dlds22[i] <- r2 / s2 ^ 2 - 3 * sp ^ 2 * r2 * beta ^ 2 / s2 ^ 4 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) + 2 * sp ^ 4 * r2 ^ 2 * beta ^ 4 / s2 ^ 6 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) ^ 2 + 4 * sp ^ 6 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) ^ 3 * B / s1 ^ 4 * r2 ^ 2 * beta ^ 4 / s2 ^ 6 - 3 * sp ^ 4 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) ^ 2 * B / s1 ^ 4 * r2 * beta ^ 2 / s2 ^ 4 - 3 * C / s2 ^ 4 + 4 * sp ^ 6 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) ^ 3 * beta ^ 6 * D / s2 ^ 10 * r2 ^ 2 - 11 * sp ^ 4 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) ^ 2 * beta ^ 4 * D / s2 ^ 8 * r2 + 10 * sp ^ 2 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) * beta ^ 2 * D / s2 ^ 6 + 8 * sp ^ 6 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) ^ 3 * beta ^ 5 * E / s1 ^ 2 / s2 ^ 8 * r2 ^ 2 - 14 * sp ^ 4 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) ^ 2 * beta ^ 3 * E / s1 ^ 2 / s2 ^ 6 * r2 + 6 * sp ^ 2 / (1 + sp ^ 2 * (r1 / s1 ^ 2 + r2 * beta ^ 2 / s2 ^ 2)) * beta * E / s1 ^ 2 / s2 ^ 4
  }

  # Second partial derivatives
  dldmu2 <- sum(dldmu2)
  dldmualpha <- sum(dldmualpha)
  dldmubeta <- sum(dldmubeta)
  dldmusp <- sum(dldmusp) 
  dldmus1 <- sum(dldmus1)
  dldmus2 <- sum(dldmus2)
  dldalpha2 <- sum(dldalpha2)
  dldalphabeta <- sum(dldalphabeta)
  dldalphasp <- sum(dldalphasp)
  dldalphas1 <- sum(dldalphas1)
  dldalphas2 <- sum(dldalphas2)
  dldbeta2 <- sum(dldbeta2)
  dldbetasp <- sum(dldbetasp)
  dldbetas1 <- sum(dldbetas1)
  dldbetas2 <- sum(dldbetas2)
  dldsp2 <- sum(dldsp2)
  dldsps1 <- sum(dldsps1)
  dldsps2 <- sum(dldsps2)
  dlds12 <- sum(dlds12)
  dlds1s2 <- sum(dlds1s2)
  dlds22 <-sum(dlds22)
  
  # Create the matrix
  J <- matrix(0, nrow = 6, ncol = 6)
  J[1,] <- c(dldmu2, dldmualpha, dldmubeta, dldmusp, dldmus1, dldmus2)
  J[2,] <- c(dldmualpha, dldalpha2, dldalphabeta, dldalphasp, dldalphas1, dldalphas2)
  J[3,] <- c(dldmubeta, dldalphabeta, dldbeta2, dldbetasp, dldbetas1, dldbetas2)
  J[4,] <- c(dldmusp, dldalphasp, dldbetasp, dldsp2, dldsps1, dldsps2)
  J[5,] <- c(dldmus1, dldalphas1, dldbetas1, dldsps1, dlds12, dlds1s2)
  J[6,] <- c(dldmus2, dldalphas2, dldbetas2, dldsps2, dlds1s2, dlds22)
  
  return(-J) #change the sign since we want the negative of the expected value instead of the expected value
}

############################
## P(agreement) Gradients ##
############################
delta.method <- function(delta, p, alpha, beta, s1, s2){
  # Function that calculates the gradient vector for theta(p)
  D <- matrix(0, nrow = 6)
  # Derivative with respect to mu
  D[1,] <- 0
  # Derivative with respect to alpha
  D[2,] <- -pi ^ (-0.1e1 / 0.2e1) * exp(-(p * beta + alpha - delta - p) ^ 2 / (s1 ^ 2 + s2 ^ 2) / 2) * sqrt(2) * (s1 ^ 2 + s2 ^ 2) ^ (-0.1e1 / 0.2e1) / 2 + pi ^ (-0.1e1 / 0.2e1) * exp(-(p * beta + alpha + delta - p) ^ 2 / (s1 ^ 2 + s2 ^ 2) / 2) * sqrt(2) * (s1 ^ 2 + s2 ^ 2) ^ (-0.1e1 / 0.2e1) / 2
  # Derivative with respect to beta
  D[3,] <- -pi ^ (-0.1e1 / 0.2e1) * exp(-(p * beta + alpha - delta - p) ^ 2 / (s1 ^ 2 + s2 ^ 2) / 2) * sqrt(2) * p * (s1 ^ 2 + s2 ^ 2) ^ (-0.1e1 / 0.2e1) / 2 + pi ^ (-0.1e1 / 0.2e1) * exp(-(p * beta + alpha + delta - p) ^ 2 / (s1 ^ 2 + s2 ^ 2) / 2) * sqrt(2) * p * (s1 ^ 2 + s2 ^ 2) ^ (-0.1e1 / 0.2e1) / 2
  # Derivative with respect to sp
  D[4,] <- 0
  # Derivative with respect to s1
  D[5,] <- pi ^ (-0.1e1 / 0.2e1) * exp(-(p * beta + alpha - delta - p) ^ 2 / (s1 ^ 2 + s2 ^ 2) / 2) * sqrt(2) * (p * beta + alpha - delta - p) * (s1 ^ 2 + s2 ^ 2) ^ (-0.3e1 / 0.2e1) * s1 / 2 - pi ^ (-0.1e1 / 0.2e1) * exp(-(p * beta + alpha + delta - p) ^ 2 / (s1 ^ 2 + s2 ^ 2) / 2) * sqrt(2) * (p * beta + alpha + delta - p) * (s1 ^ 2 + s2 ^ 2) ^ (-0.3e1 / 0.2e1) * s1 / 2
  # Derivative with respect to s2
  D[6,] <- pi ^ (-0.1e1 / 0.2e1) * exp(-(p * beta + alpha - delta - p) ^ 2 / (s1 ^ 2 + s2 ^ 2) / 2) * sqrt(2) * (p * beta + alpha - delta - p) * (s1 ^ 2 + s2 ^ 2) ^ (-0.3e1 / 0.2e1) * s2 / 2 - pi ^ (-0.1e1 / 0.2e1) * exp(-(p * beta + alpha + delta - p) ^ 2 / (s1 ^ 2 + s2 ^ 2) / 2) * sqrt(2) * (p * beta + alpha + delta - p) * (s1 ^ 2 + s2 ^ 2) ^ (-0.3e1 / 0.2e1) * s2 / 2
  return(D)
}

######################
## Modified QQ-Plot ##
######################
mod.qqplot <- function(data, n){ 
  ## This function creates a modified QQ-plot as described in Stevens et al. (2015) but for unbalanced repeated measurements
  ## Inputs:
  #  data = a dataframe with four columns (Subject, MS, Repeat, Measurements) and N = sum{r_ij} rows. The entries of Measurements
  #         correspond to the repeated measurements by the given system on each of the n subjects. Note that the observations 
  #         are ordered by subject with repeated measurements appearing consecutively with MS1 measurements appearing first and MS2 second
  #     n = the number of subjects used in the study
  
  ## Sample call: mod.qqplot(data = rr, n = 21)
  
  mm1 <- rep(0, n) # vector of subject-means (MS1)
  mm2 <- rep(0, n) # vector of subject-means (MS2)
  for(i in 1:n){
    mm1[i] <- mean(subset(data, (subset = data$Subject==i & data$MS==1))$Measurements)
    mm2[i] <- mean(subset(data, (subset = data$Subject==i & data$MS==2))$Measurements) 
  }
  mu1 <- mean(mm1)
  mu2 <- mean(mm2)
  sd1 <- sd(mm1)
  sd2 <- sd(mm2)
  rand.norm1 <- matrix(0, nrow = n, ncol = 50)
  rand.norm2 <- matrix(0, nrow = n, ncol = 50)
  for (i in 1:50){
    rand.norm1[,i] <- rnorm(n, mu1, sd1)
    rand.norm2[,i] <- rnorm(n, mu2, sd2)
  }
  par(mfrow=c(1,2))
  # MS1
  blackpts1 <- qqnorm(mm1, ylim = range(mm1, rand.norm1), main = "QQ-Plot of MS1 Part Averages vs Standard Normal", pch = 19, cex = 0.75)
  for (i in 1:50){
    graypts1 <- qqnorm(rand.norm1[,i], plot.it = FALSE)
    points(graypts1, col = "grey87", cex = 0.75)
  }
  points(blackpts1, col = "black", pch = 19, cex = 0.75) 
  qqline(mm1, col = "red")
  # MS2
  blackpts2 <- qqnorm(mm2, ylim = range(mm2, rand.norm2), main = "QQ-Plot of MS2 Part Averages vs Standard Normal", pch = 19, cex = 0.75)
  for (i in 1:50){
    graypts2 <- qqnorm(rand.norm2[,i], plot.it = FALSE)
    points(graypts2, col = "grey87", cex = 0.75)
  }
  points(blackpts2, col = "black", pch = 19, cex = 0.75) 
  qqline(mm2, col = "red")
}

########################
## Repeatabiltiy Plot ##
########################
repeat.plot <- function(data, n){
  ## This function creates a repeatability plot as described in Stevens et al. (2015) but for unbalanced repeated measurements
  ## Inputs:
  #  data = a dataframe with four columns (Subject, MS, Repeat, Measurements) and N = sum{r_ij} rows. The entries of Measurements
  #         correspond to the repeated measurements by the given system on each of the n subjects. Note that the observations 
  #         are ordered by subject with repeated measurements appearing consecutively with MS1 measurements appearing first and MS2 second
  #     n = the number of subjects used in the study
  
  ## Sample call: repeat.plot(data = rr, n = 21)
  
  mm1 <- rep(0, n) # vector of subject-means (MS1)
  mm2 <- rep(0, n) # vector of subject-means (MS2)
  r1vec <- rep(0, n) # vector of repeated measurement numbers (MS1)
  r2vec <- rep(0, n) # vector of repeated measurement numbers (MS2)
  for(i in 1:n){
    mm1[i] <- mean(subset(data, subset = (data$Subject==i & data$MS==1))$Measurements)
    mm2[i] <- mean(subset(data, subset = (data$Subject==i & data$MS==2))$Measurements)
    r1vec[i] <- max(subset(data, subset = (data$Subject==i & data$MS==1))$Repeat)
    r2vec[i] <- max(subset(data, subset = (data$Subject==i & data$MS==2))$Repeat)
  }
  mm1 <- rep(mm1, r1vec) # vector of MS1 subject-means repeated r_i1 times
  res1 <- data[data$MS==1,]$Measurements - mm1
  mm2 <- rep(mm2, r2vec) # vector of MS1 subject-means repeated r_i2 times
  res2 <- data[data$MS==2,]$Measurements - mm2
  par(mfrow=c(1,2))
  plot(mm1, res1, pch = 19, cex = 0.5, main = "MS1 Residual Measurements vs Subject-Average", ylab = "Observed Residuals", xlab = "Subject Averages", xlim = c(min(mm1,mm2),max(mm1,mm2)), ylim = c(min(res1,res2),max(res1,res2)))
  abline(h = 0, col = "red")
  plot(mm2, res2, pch = 19, cex = 0.5, main = "MS2 Residual Measurements vs Subject-Average", ylab = "Observed Residuals", xlab = "Subject Averages", xlim = c(min(mm1,mm2),max(mm1,mm2)), ylim = c(min(res1,res2),max(res1,res2)))
  abline(h = 0, col = "red")
}

####################################
## Repeated Measures Scatter Plot ##
####################################
rep.scatter.plot <- function(data, n){
  Y1 <- subset(data, subset = data$MS==1)$Measurements
  Y2 <- subset(data, subset = data$MS==2)$Measurements
  sub_means1 <- rep(0, n)
  sub_means2 <- rep(0, n)
  for(i in 1:n){
    sub_i1 <- subset(data, subset = (data$Subject==i & data$MS==1))
    sub_means1[i] <- mean(sub_i1$Measurements)
    sub_i2 <- subset(data, subset = (data$Subject==i & data$MS==2))
    sub_means2[i] <- mean(sub_i2$Measurements)
  }
  par(mfrow = c(1,1))
  plot(sub_means1, sub_means2, ylim = c(min(sub_means1, sub_means2), max(sub_means1, sub_means2)), xlim = c(min(sub_means1, sub_means2), max(sub_means1, sub_means2)), xlab = "MS1 Subject Averages", ylab = "MS2 Subject Averages", main = "Repeated Measures Scatterplot", pch = 16)
  abline(a = 0, b = 1, col = "red", lty = 2)
}

###########################
## Ordinary Scatter Plot ##     (this only works when number of measurements per subject is the same across devices)
###########################
scatter.plot <- function(data, n){
  Y1 <- subset(data, subset = data$MS==1)$Measurements
  Y2 <- subset(data, subset = data$MS==2)$Measurements
  par(mfrow = c(1,1))
  plot(Y1, Y2, ylim = c(min(Y1, Y2), max(Y1, Y2)), xlim = c(min(Y1, Y2), max(Y1, Y2)), xlab = "MS1 Readings", ylab = "MS2 Readings", main = "Scatterplot", pch = 16)
  abline(a = 0, b = 1, col = "red", lty = 2)
}


###########################
## P(agreement) Function ##
###########################
prob.agree <- function(n, delta, data){
  ## This function conducts the probability of agreement analysis for the 
  ## comparison of two homoscedastic measurement systems with possibly
  ## unbalanced measurements.
  
  ## Sample call: prob.agree(n = 21, delta = 5, data = rr)
  
  ## Inputs:
  #     n = the number of parts/subjects being measured by each system
  # delta = the the maximum allowable difference between measurement by two systems
  #  data = a dataframe with four columns (Subject, MS, Repeat, Measurements) and N = sum{r_ij} rows. The entries of Measurements
  #         correspond to the repeated measurements by the given system on each of the n subjects. Note that the observations 
  #         are ordered by subject with repeated measurements appearing consecutively with MS1 measurements appearing first and MS2 second
  
  
  # Make sure the correct inputs are inputted
  if(is.null(n) == TRUE){
    print("You must enter the number of subjects, n.")
  }else if(is.null(delta) == TRUE){
    print("You must define the clinically acceptable difference by inputing delta.")
  }else if(is.null(data) == TRUE){
    print("You must enter data for both measurement systems.")
  }else{ # Function inputs are correct, therefore we can proceed
    
    # Get initial guesses for parameters (needed for likelihood maximization)
    Y1 <- subset(data, subset = data$MS==1)$Measurements
    Y2 <- subset(data, subset = data$MS==2)$Measurements
    sub_means1 <- rep(0, n)
    sub_means2 <- rep(0, n)
    sub_vars1 <- rep(0, n)
    sub_vars2 <- rep(0, n)
    r1vec <- rep(0, n)
    r2vec <- rep(0, n)
    for(i in 1:n){
      sub_i1 <- subset(data, subset = (data$Subject==i & data$MS==1))
      sub_means1[i] <- mean(sub_i1$Measurements)
      sub_vars1[i] <- var(sub_i1$Measurements)
      r1vec[i] <- max(sub_i1$Repeat)
      sub_i2 <- subset(data, subset = (data$Subject==i & data$MS==2))
      sub_means2[i] <- mean(sub_i2$Measurements)
      sub_vars2[i] <- var(sub_i2$Measurements)
      r2vec[i] <- max(sub_i2$Repeat)
    }
    muinit <- mean(c(Y1,Y2)) # overall mean for both MS observations
    alphainit <- mean(Y2) - mean(Y1)
    betainit <- 1
    s1init <- sqrt(mean(sub_vars1)) # repeatability for MS1 (pooled across parts)
    s2init <- sqrt(mean(sub_vars2)) # repeatability for MS2 (pooled across parts)
    spinit <- mean(c(sqrt(var(Y1) - s1init^2), sqrt(var(Y2) - s2init^2)))
    
    # Saving initial parameters in matrix
    init.param <- t(as.matrix(c(muinit,alphainit,betainit,spinit,s1init,s2init)))
    
    # Maximize the Likelihood
    options(warn = -1)
    maximize <- optim(init.param, fn = log.likelihood, gr = gradient, n, r1vec, r2vec, data, method = "BFGS", hessian = T)
    
    # Get MLEs and standard errors 
    param.hat <- as.matrix(t(maximize$par)) # save parameters estimates in matrix
    muhat <- param.hat[1]
    alphahat <- param.hat[2]
    betahat <- param.hat[3]
    sphat <- param.hat[4]
    s1hat <- param.hat[5]
    s2hat <- param.hat[6]
    
    #Get Inverse of Fisher Information to get Assymptotic Variances
    I <- fisher.info(param.hat, n, r1vec, r2vec)
    inverseI <- solve(I) 
    
    # Square root of diagonal of Inverse Fisher information matrix: find Standard Errors
    standard.errors <- as.matrix(sqrt(diag(inverseI)))   
    # standard.errors <- as.matrix(sqrt(diag(solve(maximize$hessian))))   
    
    # Construct estimate of P(agreement) and get its standard error for CIs
    p <- seq(muhat-3*sphat, muhat+3*sphat, 0.1)
    thetahat <- matrix(0, nrow = length(p))
    SEtheta <- matrix(0, nrow = length(p))
    for(i in 1:length(p)){
      k1hat <- (-delta-alphahat-p[i]*(betahat-1))/(sqrt((s1hat^2)+(s2hat^2)))
      k2hat <- (delta-alphahat-p[i]*(betahat-1))/(sqrt((s1hat^2)+(s2hat^2)))
      thetahat[i] <- pnorm(k2hat)-pnorm(k1hat) # estimated theta (evaluated at MLEs)
      D <- delta.method(delta = delta,p = p[i], alpha = alphahat, beta = betahat, s1 = s1hat,s2 = s2hat)
      SEtheta[i] <- sqrt(t(D) %*% inverseI %*% D) # Standard error of Probability of Agreement
    }
    LCL <- thetahat-1.96*SEtheta
    LCL[LCL < 0] = 0
    UCL <- thetahat+1.96*SEtheta
    LCL[LCL > 1] = 1
    # Probability of Agreement plot
    par(mfrow=c(1,1))
    plot(p, thetahat, main = bquote(plain("Probability of Agreement Plot with ") ~ delta == .(delta)), ylim = c(0,1), col = "red", type = "l", ylab = "θ(s)", xlab = "s")
    points(p, thetahat-1.96*SEtheta, type = "l", col = "blue", lty = 2)
    points(p, thetahat+1.96*SEtheta, type = "l", col = "blue", lty = 2)
    labels <- c("θ(s)","95% CI")
    legend("bottomright", inset = 0.02, labels, col = c("red", "blue"), lty = c(7,2))
  }
  
  # Order of return: (mu, alpha, beta, sp, s1, s2)
  return(list(Estimates = param.hat, St.Errors = standard.errors, PAgree = data.frame(True.Val = p, PA = thetahat)))
}


