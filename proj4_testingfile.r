#setwd("C:\Users\kaklu\StatisticalProgramming\proj4.r")
setwd("C:\\Users\\julia\\Downloads")
source('proj4.r') 

#*************** Simon's example from practical file***********

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



newt(theta=c(1,2),rb,gb,hb)
newt(theta=c(1,2),rb,gb,fscale=2)

#********** Functions custome written ****************************

# Functions
# optim/nlm maximizes when control$fnscale is negative. Will try testing that.

poly1<-function(xy=c(1,1)){
  # const=c(...)
  x<-xy[1];y<-xy[2]
  2*x^3 + (6*x*y^2)-(3*y^3)-150*x
}

poly1_grad<-function(xy=c(1,1)){
  # const=c(...)
  x<-xy[1];y<-xy[2]
  wrt_x<-6*x^2 + 6*y^2 -150
  wrt_y<-12*x*y-9*y^2
  c(wrt_x,wrt_y)
}
# baseline
optim_start<-Sys.time()
optim_result<-optim(c(1,8),fn=poly1,gr=poly1_grad)
optim_end<-Sys.time()
optim_result$par
cat("\nOptim time ",optim_end-optim_start)

# nlm_result<-nlm(p=c(2,3),f=poly1)
# nlm_result$estimate

ours_start<-Sys.time()
our_result<-newt(c(1,8),func=poly1,grad=poly1_grad)
ours_end<-Sys.time()
our_result$thetas
cat('\nOur newt time ',ours_end-ours_start)

newt()

trig<-function(thetas=c(1,1)){
  x<-thetas[1]
  y<-thetas[2]
  sin(x)*cos(y) + sin(x)
}

trig_grad<-function(thetas=c(1,1)){
  x<-thetas[1]
  y<-thetas[2]
  wrt_x<-cos(x)*cos(y) + cos(x)
  wrt_y<-sin(x)*-1*sin(y)
  c(wrt_x,wrt_y)
}

optim_result<-optim(c(1,1),fn=trig,gr=trig_grad,control=list('fnscale' = -1))
optim_result$par
our_result<-newt(theta=c(1,1),func=trig,grad=trig_grad,fscale=-1,maxit = 200)
our_result$theta

#************ Poisson example from notes ****************

t80 <- 1:13 ## years since 1980
y <- c(12,14,33,50,67,74,123,141,165,204,253,246,240) ## AIDS cases
gll <- function(theta,t,y) {
  ## grad of -ve log lik of Poisson AIDS early epidemic model
  alpha <- theta[1];beta <- theta[2] ## enhances readability
  ebt <- exp(beta*t) ## avoid computing twice
  -c(sum(y)/alpha - sum(ebt),sum(y*t) - alpha*sum(t*ebt)) ## -dl/dbeta
} ## gll
nll <- function(theta,t,y) {
  ## -ve log likelihood for AIDS model y_i ~ Poi(alpha*exp(beta*t_i))
  ## theta = (alpha,beta)
  mu <- theta[1] * exp(theta[2] * t) ## mu = E(y)
  -sum(dpois(y,mu,log=TRUE)) ## the negative log likelihood
} ## nll

result<-newt(theta=c(10,.1),func=nll,grad=gll,hess=NULL,y=y,t=t80)
result

optim_result<-optim(c(10,0.1),fn=nll,gr=gll,y=y,t=t80)
optim_result







# define_hessian <- function(theta,...){
#   grad<-gll
#   Hfd <- matrix(0, length(theta), length(theta))
#   # generate a matrix of 0s to later replace the values with the 
#   # second derivatives
#   grad_old <- grad(theta, ...)
#   for (i in 1:length(theta)){
#     new_theta <- theta
#     new_theta[i] <- new_theta[i] + eps
#     Hfd[i,] <- (grad(new_theta, ...) - grad_old) / eps
#     # this step computes the second derivatives, Hessian. 
#     # second derivatives can be obtained by applying the first principle to
#     # 'grad', the gradient function.
#   }
#   (t(Hfd) + Hfd) / 2
#   # make the hessian matrix symmetric
# }
# define_hessian(c(10,.1),t=t80,y=y,eps=)
