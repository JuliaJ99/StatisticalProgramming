# ------- Group member names --------------------------------------------------
# Julia Kaczmarczyk - s1977294 (Julia K)
# Akira Ishiyama - s2445245
# Julia Jose - s2421229 (Julia J)

# ------- Github repo --------------------------------------------------------
# https://github.com/JuliaJ99/StatisticalProgramming

# ------- Individual Contributions --------------------------------------------
# Julia K - Coding the whole main framework
# Akira Ishiyama - code editing, commenting
# Julia J - code editing, code debugging, code testing

# Rough proportions: Julia K - 36%, Akira - 32%, Julia J - 32%

# ------- Outline of project --------------------------------------------------

#         Simulate Newton's method based numerical optimisation

#   The aim of the function to be coded is to find a minimiser of a given 
# objective function. Assuming that an objective function is smooth with three
# bounded derivatives, one can derive its gradient and Hessian matrix. For a
# parameter to a function minimiser, all gradients at the parameter value has to 
# be 0, and the Hessian matrix has to be positive definite.
#   The Newton's optimisation method starts with an initial guess, and 
# finds its objective function value, gradient, and Hessian. It then minimises
# the function, and finds an improved guess of the parameter. This iteration
# is carried on until the gradient of the function at a parameter value reaches 
# approximately 0. To ensure the gradient converges to 0, each iteration
# has to check whether a newly proposed parameter actually reduces the function
# value. 
#   The method is added with two modification in this project to guarantee 
# convergence. Firstly, if the Hessian is not positive definite, the 
# function newt perturbs it to be so. Secondly, to prevent divergence, 
# newt halves the 'delta' value (the minimiser of an objective function based on 
# Taylor's theorem) until the function value is reduced for an iteration. 


#   Algorithm: A brief walk-through of the newt function is as follows;
# 1. Approximate the Hessian if not provided
# 2. Check if the objectives and derivatives are finite at the initial theta
# 3. Iteration 1, while the iteration count is less than the limit set by
#    'maxit', compute all the necessary information to run optimisation such as
#    the gradient and the Hessian at the given theta for each iteration. 
#    improve the theta parameter in each iteration by computing 'delta',
#    function minimiser based on Taylor's theorem, using existing information.
# 4. Iteration 2, in case the 'delta' parameter diverges the gradient, while 
#    the function value after taking into account of 'delta' is greater than
#    the original function value, half the 'delta' parameter. Stop the newt 
#    function if unable to improve the minimisation, after halving 'delta' to
#    the limit 'max.half' set by the function user.
# 5. Check the improved gradient at each iteration. If the gradient has
#    converged to a value less than the benchmark set by the function user,
#    report that the convergence is reached. 
# 6. Keep going through the two iteration until a convergence is reached, or 
#    stop the function newt if the iteration count for iteration 1 goes above
#    the limit 'maxit' set by the function user.



## ------ Code ----------------------------------------------------------------

newt <- function(theta, func, grad, hess = NULL, ..., tol = 1e-8, fscale = 1, 
                 maxit = 100, max.half = 20, eps = 1e-6){
#   newt is a function that runs numerical optimisation, with the underlying 
# theory based on Newton's method. It searches for an optimal set of theta 
# values that minimises an objective function given as an input argument. The 
# optimality condition is defined by the following criteria: function gradient 
# is zero and Hessian is positive definite at theta. Gradient of zero is tested 
# via its convergence to approximately zero. newt takes some initial thetas as a 
# starting point to carry on a set of iterations to search for thetas that meet 
# the aforementioned criteria.
#
#   newt also provides warnings to the function user if the specific condition
# is deemed impossible to compute the optimal thetas.
#
# Input: 
#   theta: initial values to start the optimisation, 
#   func: objective function to be minimised, 
#   grad: gradient of the function, 
#   hess: Hessian matrix of the function, 
#         newt will compute the hessian using finite differencing if hess = NULL
#   ...: additional arguments of func, grad, and hess,
#   tol: the convergence tolerance to test if the gradient is close enough to 0, 
#   fscale: a rough estimate of the magnitude of func near the optimum, 
#   maxit: the maximum number of Newton iterations to try, 
#   max.half: the maximum number of times 'delta', the minimiser of an
#             objective function based on Taylor's theorem, should be halved, 
#   eps: the finite difference intervals to use for Hessian approximation
# Output: 
#   f: the value of the objective function at the minimum,
#   theta: the value of the parameters at the minimum,
#   iter: the number of iterations taken to reach the minimum,
#   g: the gradient at the minimum,
#   Hi: the inverse of the Hessian matrix at the minimum.
  

  # Finite difference approximation for the Hessian matrix if the 'hess'
  # argument is missing
  define_hessian <- function(theta,...){
    Hfd <- matrix(0, length(theta), length(theta))
    # generate a matrix of 0s to later replace the values with the 
    # second derivatives
    grad_old <- grad(theta, ...)
    for (i in 1:length(theta)){
      new_theta <- theta
      new_theta[i] <- new_theta[i] + eps
      Hfd[i,] <- (grad(new_theta, ...) - grad_old) / eps
      # this step computes the second derivatives, Hessian. 
      # second derivatives can be obtained by applying the first principle to
      # 'grad', the gradient function.
      }
    (t(Hfd) + Hfd) / 2
    # make the hessian matrix symmetric
  }
  
  if(missing(hess) | is.null(hess)){
    # if argument 'hess' is missing from newt, define a Hessian by calling the
    # 'define_hessian' function to approximate.
    hess <- define_hessian
  }
    
  # Break the function if any one of the objectives or derivatives are not 
  # finite at the initial theta
  if ((abs(func(theta,...)) == Inf) || 
      any(abs(grad(theta,...)) == Inf) || 
      any(c(abs(hess(theta,...))) == Inf)) {
    stop('Objective function or derivatives not finite at the initial theta')
  }

  new_theta <- theta
  iter <- 0 # the number of iteration done in the while loop
  
  # Loop to search for the theta that minimises the objective function
  while (iter <= maxit){
    hessian <- hess(new_theta,...)
    # second derivatives
    
    counter_posdef<-0
    while(inherits(try(chol(hessian), silent = TRUE), "try-error")){
      # perturb Hessian until cholesky decomposition is possible by adding a 
      # multiple of the identity matrix
      
      a<-dim(hessian)[1]
      hessian <- hessian + diag(a)*1e-6*norm(hessian)*(10^counter_posdef)
      counter_posdef<-counter_posdef+1
      if(counter_posdef>1000000){
        stop('Hessian cannot be pertubed to be positive definite')
      }
    }
    
    gradient <- grad(new_theta,...)
    # first derivatives
    D <- func(new_theta,...)
    # objective function values
    chol_hessian <- chol(hessian)
    # Cholesky decomposition is a preparatory step towards using 'chol2inv' in
    # the next line
    delta <- -1 * (Matrix::chol2inv(chol_hessian)) %*% gradient
    # use 'Matrix::chol2inv(chol_hessian)' to find the inverse of 'hessian' to 
    # compute 'delta' that minimises D_delta (objective function value) below, 
    # based on Taylor's theorem
    D_delta <- func(new_theta + delta,...)
    # objective function values updated with delta
    
    step <- 0 # the number of times delta is halved
    while (D_delta > D){
      # while the updated objective function value 'D_delta' is greater
      # than the original function value 'D', keeps halving the 'delta'
      # parameter
      step <- step + 1
      
      # Break the function if the number of steps halved goes above the 
      # number set by 'max.half'
      if (step > max.half){
        stop('Fails to reduce the function value by halving delta')
      }
      
      delta <- delta / 2
      D_delta <- func(new_theta + delta,...)
      # update objective function values to test whether it is now smaller 
      # than D, the condition set in the while loop, after delta being halved
    }
    
    new_theta <- new_theta + delta
    # optimal theta obtained by the loop above that minimises the objective 
    # function value

    # Examine if the gradient of the objective function at new_theta is 
    # approximately 0, after taking into account of the tolerance set by the 
    # function user. if so, the convergence is reached
    if ((max(abs(grad(new_theta,...))) < (tol * (abs(D_delta) + fscale)))){
      
      # Break the function if the hessian is not positive definite at 
      # convergence
      if (inherits(try(chol(hess(new_theta,...)), silent = TRUE), "try-error")){
        warning('Hessian not positive definite at convergence')
        Hi<-NULL
      }else{
        Hi<-chol2inv(chol(hess(new_theta,...)))
      }
      
      
      
      # Convergence achieved, return the results
      return(list('f' = func(new_theta,...), 'theta' = new_theta, 
                  'iter' = iter, 'g' = grad(new_theta,...), 
                  'Hi' = Hi))
      
    }
    iter <- iter + 1
  }
  
  # If the iteration count reaches to be greater than or equal to the number 
  # set by 'maxit', notify the user that convergence was not reached
  if (iter >= maxit){
    stop('Maximum number of iterations reached without convergence')
  }
  
}