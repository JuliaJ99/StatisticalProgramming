# ------- Group member names --------------------------------------------------
# Julia Kaczmarczyk - s1977294 (Julia K)
# Akira Ishiyama - s2445245
# Julia Jose - s2421229 (Julia J)

# ------- Github repo --------------------------------------------------------
# https://github.com/JuliaJ99/StatisticalProgramming

# ------- Individual Contributions --------------------------------------------
# 

# Rough proportions: Julia K - , Akira - , Julia J - 

# ------- Outline of project --------------------------------------------------
# Code to simulate Newton's method based numerical optimisation.



## ------ Code ----------------------------------------------------------------

newt <- function(theta, func, grad, hess = NULL, ..., tol = 1e-8, fscale = 1, 
                 maxit = 100, max.half = 20, eps = 1e-6){
#   newt is a function that runs numerical optimisation, with the underlying 
# theory based on Newton's method. It searches for an optimal set of thetas 
# that minimises an objective function given as an input argument. The 
# optimality condition is defined by the following criteria; function gradient 
# of zero and positive definite Hessian at theta. newt takes some initial 
# thetas as one of its arguments as a start point to carry on a set of
# iterations to search for thetas that meets the aforementioned criteria.
#   An inner function 'define_hessian' is provided to approximate the Hessian 
# matrix if one is not provided with the function argument.
#   newt also provides warnings to the function user if the specific condition
# is deemed impossible to compute the optimal thetas.
  
  # Finite difference approximation for the Hessian matrix if the 'hess'
  # argument is missing
  define_hessian <- function(theta){
    Hfd <- matrix(0, length(theta), length(theta))
    # generate a matrix of 0s to later replace the values with the 
    # second derivatives
    for (i in 1:length(theta)){
      new_theta <- theta
      new_theta[i] <- new_theta[i] + eps
      Hfd[i,] <- (grad(new_theta, ...) - grad(theta, ...)) / eps
      # this step essentially computes the second derivatives. 
      # second derivatives can be obtained by applying the first principle to
      # 'grad', the first derivatives.
      }
    (t(Hfd) + Hfd) / 2
    # make the hessian matrix symmetric
  }
  
  if(missing(hess)){
    hess <- define_hessian
  }
    
  # Break the function if any one of the objectives or derivatives are not 
  # finite at the initial theta
  if ((abs(func(theta)) == Inf) || 
      any(abs(grad(theta)) == Inf) || 
      any(c(abs(hess(theta))) == Inf)) {
    stop('Objective function or derivatives not finite at the initial theta')
  }

  new_theta <- theta
  iter <- 0 # the number of iteration done in the while loop
  
  # Loop to search for the theta that minimises the objective function
  while (iter <= maxit){
    hessian <- hess(new_theta)
    # second derivatives
    while(inherits(try(chol(hessian), silent = TRUE), "try-error")){
      hessian <- hessian + diag(dim(hessian)[1])
    }
    # perturb Hessian until cholesky decomposition is possible
    gradient <- grad(new_theta)
    # first derivatives
    D <- func(new_theta)
    # objective function values
    chol_hessian <- chol(hessian)
    # Cholesky decomposition is a preparatory step towards using 'chol2inv' in
    # the next line
    delta <- -1 * (Matrix::chol2inv(chol_hessian)) %*% gradient
    # use 'Matrix::chol2inv(chol_hessian)' to find the inverse of 'hessian' to 
    # compute 'delta' that minimises D_delta (objective function value) below, 
    # based on Taylor's theorem
    D_delta <- func(new_theta + delta)
    # objective function values updated with delta
    
    step <- 0 # the number of times delta is halved
    while (D_delta > D){
      step <- step + 1
      
      # Break the function if the number of steps halved goes above the 
      # number set by 'max.half'
      if (step > max.half){
        stop('Step fails to reduce')
      }
      
      delta <- delta / 2
      D_delta <- func(new_theta + delta)
      # update objective function values to test whether it is now smaller 
      # than D, the condition set in the while loop, after delta being halved
    }
    
    new_theta <- new_theta + delta
    # optimal theta obtained by the loop above that minimises the objective 
    # function value

    # Examine if the gradient of the objective function at new_theta is 
    # approximately 0, after taking into account of the tolerance set by the 
    # function user. if so, the convergence is reached
    if ((max(abs(grad(new_theta))) < (tol * (abs(D_delta) + fscale)))){
      
      # Break the function if the hessian is not positive definite at 
      # convergence
      if (inherits(try(chol(hess(new_theta)), silent = TRUE), "try-error")){
        stop('Hessian not positive definite at convergence')
      }
      
      # Convergence achieved, return the results
      cat('Convergence reached')
      return(list('f' = func(new_theta), 'theta' = new_theta, 
                  'iter' = iter, 'g' = grad(new_theta), 
                  'Hi' = chol2inv(chol(hess(new_theta)))))
      
    }
    iter <- iter + 1
  }
  
  # If the iteration count reaches to be greater than or equal to the number 
  # set by 'maxit', notify the user that convergence was not reached
  if (iter >= maxit){
    cat('Maximum number of iterations reached without convergence')
  }
  
}


## Testing---------------------------------------------------------------------
rb <- function(th,k=2) {
  k*(th[2]-th[1]^2)^2 + (1-th[1])^2
}
gb <- function(th,k=2) {
  c(-2*(1-th[1])-k*4*th[1]*(th[2]-th[1]^2),k*2*(th[2]-th[1]^2))
}
hb <- function(th,k=2) {
  h <- matrix(0,2,2)
  h[1,1] <- 2-k*2*(2*(th[2]-th[1]^2) - 4*th[1]^2)
  h[2,2] <- 2*k
  h[1,2] <- h[2,1] <- -4*k*th[1]
  h
}

newt(theta=c(1,3),rb,gb,hb)
newt(theta=c(1,3),rb,gb)
