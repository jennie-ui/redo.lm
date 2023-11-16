# Iteratively Reweighted Least Squares (IRWLS)
# main effects model with no interaction terms
glmSetup <- function(res_var, base_var, input_data, tolerance){
  df = data.frame(input_data)
  n = nrow(df)
  baseline_position = which(colnames(input_data) == base_var)
  y_position = which(colnames(input_data) == res_var)
  baseLine = input_data[,baseline_position]
  Y = input_data[,y_position]
  log_Y = log(Y + 0.1)
  X = as.matrix(cbind(1,df[-c(y_position)]))
  q = ncol(X)

  tol = tolerance
  epsilon = 99
  ite_max = 25
  ite = 0

  distribution_plot = hist(Y,
                           main = paste0("Histogram of ",colnames(df)[1]),
                           xlab = colnames(df)[1])

  myObject = list("histogram" = distribution_plot,
                  "response" = Y,
                  "transformed_response" = log_Y,
                  "baseline_variable" = baseLine,
                  "design_matrix" = X,
                  "n" = n,
                  "q" = q,
                  "tol" = tol,
                  "epsilon" = epsilon,
                  "ite_max" = ite_max,
                  "ite" = ite)

  return(myObject)
}
