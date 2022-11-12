newt <- function(theta,func,grad,hess=NULL,...,tol=1e-8,fscale=1,maxit=100,max.half=20,eps=1e-6){
  
  # checking if objective or derivatives are finite at the initial theta 
  if (abs(D(theta))==Inf|abs(grad(theta))==Inf|abs(hess(theta))==Inf) {
    stop('Objective function or derivatives not finite at the initial theta!') #stop or warning?
  }
  
  # if hess not supplied use the finite difference approximation
  if (hess=NULL){
    hess <- function(grad,theta,eps,...){
      grad_th <- grad(theta)
      Hfd <- matrix(0,length(theta),length(theta))
      for (i in 1:length(theta)){
        new_theta <- theta
        new_theta[i] <- new_theta[i] + eps
        grad_new <- grad(new_theta)
        Hfd[i,] <- (grad_new - grad_th)/eps
      }
      (t(Hfd)+Hfd)/2
    }
  }
  
  
  new_theta <- theta
  iter <- 0
  
  while (iter<=maxit){
    hess_th <- hess(new_theta)
    grad_th <- grad(new_theta)
    old_D <- D(new_theta)
    new_D <- old_D + 1
    delta <- -chlol2inv(hess_th)%*%grad_th # inversion could be done nicer
    step <- 0 # number of times we halve delta
    while (new_D < old_D){
      if (step == max.half){
        stop('Step fails to reduce!')
      }
      new_D <- old_D+t(delta)%*%grad_th+1/2*t(delta)%*%hess_th%*%delta #or new_D = D(theta+delta)
      delta <- delta/2
      step <- step + 1
    }
    new_theta <- new_theta + delta
    convergence <- tol*(abs(new_D)+fscale)#brackets or no brackets???
    if (abs(grad) < convergence){
      if (class(pd_check(hess(new_theta)))=='try-error'){
        stop('Hessian not positive definite at convergence!') # stop or warning?
      }
      #converged
    }
    iter <- iter + 1
  }
  
  # max number of iterations
  if (iter == maxit){
    stop('Maximum number of iterations reached withot convergence!')
  }
  
  #check convergance
  
  
  
  
} #end of the newt function


# checking if hessian is positive definite
pd_check<- function(hess_th){
  if (class(try(chol(hess_th),silent=TRUE))=='try-error'){
    stop('Hessian not positive definite!')
  }
}

# fixing a non-pd hessian
# this fixes the hessian, but idk how we should adjust the values of theta/grad/D
# to the new hessian
pd_fix<- function(hess_th){
  while(class(try(chol(hess_th),silent=TRUE))[1]=='try-error'){
    hess_th <- hess_th + diag(dim(hess_th)[1]) # diagonal should be multipled 
                              #by some value, i chose to do a while loop for the
                              #multiple but idk if it works for small values for example
  }
  hess_th
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

hb(th)
A <-hess(gb,th,eps=1e-6,th=th,k=2)




