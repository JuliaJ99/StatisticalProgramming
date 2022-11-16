

source("proj4.r") 

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

newt(theta=c(1,3),rb,gb,hb)
newt(theta=c(1,3),rb,gb)

# Functions
# optim/nlm maximizes when control$fnscale is negative. Will try testing that.

poly1<-function(xy=c(1,1)){
  x<-xy[1];y<-xy[2]
  2*x^3 + (6*x*y^2)-(3*y^3)-150*x
}

poly1_grad<-function(xy=c(1,1)){
  x<-xy[1];y<-xy[2]
  wrt_x<-6*x^2 + 6*y^2 -150
  wrt_y<-12*x*y-9*y^2
  c(wrt_x,wrt_y)
}
# baseline
optim_result<-optim(c(2,3),fn=poly1,gr=poly1_grad)
optim_result$par
nlm_result<-nlm(p=c(2,3),f=poly1)
nlm_result$estimate
our_result<-newt(theta=c(2,3),func=poly1,grad=poly1_grad)
our_result$theta

