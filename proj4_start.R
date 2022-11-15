newt <- function(theta, func, grad, hess = NULL, ..., tol = 1e-8, fscale = 1, 
                 maxit = 100, max.half = 20, eps = 1e-6){
#!! Explain what newt does
  
  # Finite difference approximation for the Hessian matrix if the 'hess'
  # argument is missing
  if (missing(hess)){
    hess <- function(grad, theta, eps, ...){
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
      #!! how should i explain this step
    }
  }
  
  # Break the function if any one of the objectives or derivatives are not 
  # finite at the initial theta
  if ((abs(func(theta)) == Inf) || 
      any(abs(grad(theta)) == Inf) || 
      any(c(abs(hess(theta))) == Inf)) {
    stop('Objective function or derivatives not finite at the initial theta')
  }
  
  #!! How should i explain this step
  pd_fix<- function(hess_th){
    while(inherits(try(chol(hess_th), silent = TRUE), "try-error")){
      hess_th <- hess_th + diag(dim(hess_th)[1])
    }
    hess_th
  }
  
  new_theta <- theta
  iter <- 0 # the number of iteration done in the while loop
  
  # Loop to search for the theta that minimises the objective function
  while (iter <= maxit){
    hess_th <- hess(new_theta)
    # second derivatives
    hess_th <- pd_fix(hess_th)
    grad_th <- grad(new_theta)
    # first derivatives
    D <- func(new_theta)
    # objective function values
    chol_hess_th <- chol(hess_th)
    # Cholesky decomposition is a preparatory step towards using chol2inv in
    # the next line. it computes the inverse much more efficiently compared to
    # 'solve' function
    delta <- -1 * (Matrix::chol2inv(chol_hess_th)) %*% grad_th
    # find the inverse of 'hess_th', the Hessian given the theta, and multiply
    # it with -1 and 'grad_th', the first derivative, to compute 'delta' that 
    # minimises the D_delta below, based on Taylor's theorem
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
    convergence <- tol * (abs(D_delta) + fscale)
    #!! need explaination?
    
    
    # Examine if the gradient of the objective function at new_theta is 
    # approximately 0, after taking into account of the tolerance set by the 
    # function user.
    if (all(abs(grad(new_theta)) < rep(convergence,
                                       times = length(grad(new_theta))))){
      
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

th <-c(1,1)

newt(theta=c(1,3),rb,gb,hb)

