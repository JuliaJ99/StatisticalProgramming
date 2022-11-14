
# ** Julia J, 14/11 morning updates:
# Made a few changes to pd check and pd fix
# Made sure the 'max iter without convergence' message appears
# Still not sure of the convergence formula. By tweaking it, either we reach convergence too early around 1.9-1.3, 
# or we don't converge at all

# Please don't remove any of the print statements yet. I know it's annoying but it's useful to observe the process


newt <- function(theta,func,grad,hess=NULL,...,tol=1e-8,fscale=1,maxit=100,max.half=20,eps=1e-6){
  # if hess not supplied use the finite difference approximation
  
  if (missing(hess)){ #*** Julia J: replaced ==Null with missing
    hess <- function(grad,theta,eps,...){
      grad_th <- grad(theta)
      Hfd <- matrix(0,length(theta),length(theta))
      
      cat('entered hessian define func')
      for (i in 1:length(theta)){
        new_theta <- theta 
        new_theta[i] <- new_theta[i] + eps
        grad_new <- grad(new_theta)
        Hfd[i,] <- (grad_new - grad_th)/eps
      }
      
      (t(Hfd)+Hfd)/2
    }
  }
  
  # checking if objective or derivatives are finite at the initial theta 
  if ((abs(func(theta))==Inf) || 
      any(abs(grad(theta))==Inf) || 
      any(c(abs(hess(theta)))==Inf)) { # should hess be included?
    stop('Objective function or derivatives not finite at the initial theta!') #stop or warning?
  }
  
  new_theta <- theta
  iter <- 0
  
  while (iter <= maxit){
    cat('\n\nNew Iteration',iter)
    cat('\n Theta start of loop',new_theta)
    hess_th <- hess(new_theta)
    hess_th<-pd_fix(hess_th)
    grad_th <- grad(new_theta)
    old_D <- func(new_theta)
    chol_hess_th <- chol(hess_th)
    cat('\n func value at old_theta',old_D)
    delta <- -1*(chol2inv(chol_hess_th))%*%grad_th
    # forwardsolve finds the the component to be backsolved with chol_hess_th. 
    
    new_D <- old_D+ (t(delta)%*%grad_th)+ (0.5*t(delta)%*%hess_th%*%delta) #or new_D = D(theta+delta)
    
    # new_D<-func(new_theta+delta)
    
    step <- 0 # number of times we halve delta
    while (new_D > old_D){
      step <- step + 1
      if (step > max.half){
        stop('Step fails to reduce!')
      }
      delta <- delta/2
      new_D <- old_D+ (t(delta)%*%grad_th )+ (0.5*t(delta)%*%hess_th%*%delta) #or new_D = D(theta+delta)
      # new_D<-func(new_theta+delta)
      cat('\nEntered step halving, new_D with new halved delta is ',new_D)
    }
    
    cat('\n func value at new_theta',new_D)
    
    new_theta <- new_theta + delta
    cat('\n theta + delta',new_theta)
    convergence <- tol*(abs(new_D)+fscale)#brackets or no brackets???
    cat('\n Convergence value',convergence)
    cat('\n Grad value',grad(new_theta))
    if (all(abs(grad(new_theta)) < rep(convergence,times=length(grad(new_theta))))){
      
      # ** Julia J: The documentation for 'try' suggests using inherits for this kind of situation
      
      if (inherits(try(chol(hess(new_theta)),silent=TRUE),"try-error")){
        stop('Hessian not positive definite at convergence!') # stop or warning?
      }
      
      #converged - output the optimized theta
      cat('\n Convergence reached!')
      return(list('f'=func(new_theta),'theta'=new_theta))
    }
    
    iter <- iter + 1
  }
  
  # max number of iterations
    if (iter >= maxit){
      cat('Maximum number of iterations reached without convergence!')
    }
  
  #check convergance
  
  
} #end of the newt function

# 
# # checking if hessian is positive definite
# pd_check<- function(hess_th){
#  
#   if (class(try(chol(hess_th),silent=TRUE))=='try-error'){
#     stop('Hessian not positive definite!')
#   }
# }

# fixing a non-pd hessian
# this fixes the hessian, but idk how we should adjust the values of theta/grad/D
# to the new hessian
pd_fix<- function(hess_th){
  
  # ** Julia J: The documentation for 'try' suggests using inherits for this kind of situation
  
  while(inherits(try(chol(hess_th),silent=TRUE),"try-error")){
    hess_th <- hess_th + diag(dim(hess_th)[1]) # diagonal should be multipled 
                              # by some value, i chose to do a while loop for the
                              # multiple but idk if it works for small values for example
  }
  return(hess_th)
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

# hb(th)
# A <-hess(gb,th,eps=1e-6,th=th,k=2)

newt(theta=c(1,3),rb,gb,hb)

# fOR TOMO:
# library(debug)
# mtrace(newt)
# newt(theta=c(1,3),rb,gb,hb)
# mtrace(newt)
# Q
# all(c(1,1)<c(2,2))

A<- matrix(1:4,nrow=2,ncol=2)
length(class(try(chol(A),silent=TRUE))=='try-error')
inherits(try(chol(A),silent=TRUE),"try-error")

