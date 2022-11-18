setwd("C:/Users/julia/Desktop/Julia_Git_repo/StatisticalProgramming")
source('proj4.r') 

# Simon's example:
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
newt(theta=c(1,2),rb,gb)

# Functions
# optim/nlm maximizes when control$fnscale is negative. Will try testing that.

poly1<-function(xy=c(1,1),...){
  const=c(...)
  x<-xy[1];y<-xy[2]
  2*x^3 + (6*x*y^2)-(3*y^3)+const*x
}

poly1_grad<-function(xy=c(1,1),...){
  const=c(...)
  x<-xy[1];y<-xy[2]
  wrt_x<-6*x^2 + 6*y^2 +const
  wrt_y<-12*x*y-9*y^2
  c(wrt_x,wrt_y)
}
# baseline
optim_start<-Sys.time()
optim_result<-optim(c(20,8),fn=poly1,gr=poly1_grad,const=-150)
optim_end<-Sys.time()
optim_result$par
cat("\nOptim time ",optim_end-optim_start)

# nlm_result<-nlm(p=c(2,3),f=poly1)
# nlm_result$estimate

ours_start<-Sys.time()
our_result<-newt(theta=c(20,8),func=poly1,grad=poly1_grad,hess=NULL,const=-150)
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

