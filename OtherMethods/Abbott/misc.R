


# Choosing the GP length-scale prior
# Abbott & Funk use an inverse-gamma prior distribution parameterised such that
# ... P(X<l) = 0.01 and P(X<u) = 0.99 for some user-determined l and u. They use
# ... l = 14 and u = 90 in their analysis. This function allows the user to find
# ... the appropriate 

library(pracma)

findLengthscalePriorParams = function(l, u, plot=FALSE) {
  
  # Define the target function to optimise. y is our guess, theta is the target 0.01 and 0.99 quantiles
  tail_delta = function(y, target_quantiles) {
    
    # Force alpha to be positive and smaller than 1e16
    alpha = ifelse(y[1] <= 0, 1e-6, ifelse(y[1] >= 1e6, 1e6, y[1]))
    beta = ifelse(y[2] <= 0, 1e-6, ifelse(y[2] >= 1e6, 1e6, y[2]))
    
    error = c(
      0.01 - pgamma(1/target_quantiles[1], shape = alpha, rate = beta, lower.tail=FALSE),
      0.99 - pgamma(1/target_quantiles[2], shape = alpha, rate = beta, lower.tail=FALSE)
    )
    return(error)
  }
  
  # Come up with initial guesses (may need to edit this)
  y_guess = c(5, 100)
  y = fsolve( function(y) tail_delta(y, c(l, u)), y_guess )$x
  
  if (plot) {
    inv_gamma_pdf <- function(x, alpha, beta) (beta^alpha / gamma(alpha)) * x^(-alpha - 1) * exp(-beta / x)
    curve(inv_gamma_pdf(x, y[1], y[2]), from = 0.01, to = 100, col = "blue", lwd = 2, ylab = "Density", main = "Inverse Gamma PDF")
    abline(v = c(l, u), col = "red", lty = 2)
  }
  
  return(y)
  
}
