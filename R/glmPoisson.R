glmPoisson <- function(myObject){
  Y = myObject$response
  log_Y = myObject$transformed_response
  X = myObject$design_matrix
  n = myObject$n
  q = myObject$q

  tol = myObject$tol
  epsilon = myObject$epsilon
  ite_max = myObject$ite_max
  ite = myObject$ite

  beta = solve(t(X) %*% X) %*% t(X) %*% log_Y
  tol = 0.0001

  while (epsilon > tol & ite <= ite_max){
    eta = X %*% beta
    mu = exp(eta)
    nu = exp(eta)
    V = diag(x = as.vector(nu))
    Z = eta + solve(V) %*% (Y-mu)
    beta_new = solve(t(X) %*% V %*% X) %*% t(X) %*% V %*% Z
    epsilon = sqrt(t(beta_new-beta)%*%(beta_new-beta))
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

  D_i = 2*(Y*(log(Y/mu)-(Y-mu)))
  r_D = sign(Y-mu)*sqrt(abs(D_i))

  H=sqrt(V) %*% X %*% solve(t(X) %*% V %*% X) %*% t(X) %*% sqrt(V)
  Leverage=diag(H)

  r_PS = r_p/sqrt(1-Leverage)
  r_DS = r_D/sqrt(1-Leverage)

  CookD = Leverage/(q*(1-Leverage))*(r_PS)^2

  D = sum(D_i)
  chiP = sum(r_p^2)
  res = data.frame(D = D, chiP = chiP, chi = qchisq(0.95,(n-q)))
  res = format(res, digits= 3)
  colnames(res) = c("Deviance", "chi^2_P", "chi^2_{n-p,0.95}")

  newList = list("variance covariance matrix" = Var_dis,
                 "estimate, standard deviation, and confidence intervals" = Wald_CI,
                 "Pearson Residual" = r_p,
                 "Deviance Residual" = r_D,
                 "Leverage" = Leverage,
                 "Standarized Pearson Residual" = r_PS,
                 "Standarized Deviance Residual" = r_DS,
                 "Cook's Distance" = CookD,
                 "GoF of the model" = res)

  return(newList)
}
