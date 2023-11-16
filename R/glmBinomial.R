glmBinomial <- function(myObject){
  Y = myObject$response
  baseLine = myObject$baseline_variable
  Y_star = Y < baseLine
  X = myObject$design_matrix
  q = myObject$q

  tol = myObject$tol
  epsilon = myObject$epsilon
  ite_max = myObject$ite_max
  ite = myObject$ite

  beta = rep(0,q)

  while (epsilon > tol & ite <= ite_max){
    eta = X %*% beta
    mu = exp(eta)/(1+exp(eta))
    nu = mu*(1-mu)
    V = diag(x = as.vector(nu))
    Z = eta + solve(V) %*% (Y_star-mu)
    beta_new = solve(t(X) %*% V %*% X) %*% t(X) %*% V %*% Z
    epsilon = sqrt(t(beta_new-beta)%*%(beta_new-beta))
    beta = beta_new
    ite = ite + 1
    beta_t = t(beta)
  }
  Var_beta = solve(t(X) %*% V %*% X)
  rownames(Var_beta)[1] = c("Intercept")
  colnames(Var_beta)[1] = c("Intercept")
  Var_dis<-format(Var_beta, digits = 4, scientific = 4)

  Wald_CI<-data.frame(Estimate=beta, Std_Deviation=sqrt(diag(Var_beta)))
  rownames(Wald_CI)[1] = c("Intercept")
  Wald_CI$CI_Lower = Wald_CI$Estimate-1.96*Wald_CI$Std_Deviation
  Wald_CI$CI_Upper = Wald_CI$Estimate+1.96*Wald_CI$Std_Deviation
  Wald_CI = format(Wald_CI, digits= 7)

  r_p = (Y-mu)/sqrt(mu)

  H=sqrt(V) %*% X %*% solve(t(X) %*% V %*% X) %*% t(X) %*% sqrt(V)
  Leverage=diag(H)

  r_PS = r_p/sqrt(1-Leverage)

  CookD = Leverage/(q*(1-Leverage))*(r_PS)^2

  newList = list("variance covariance matrix" = Var_dis,
                 "estimate, standard deviation, and confidence intervals" = Wald_CI,
                 "Pearson Residual" = r_p,
                 "Leverage" = Leverage,
                 "Standarized Pearson Residual" = r_PS,
                 "Cook's Distance" = CookD)

  return(newList)
}
