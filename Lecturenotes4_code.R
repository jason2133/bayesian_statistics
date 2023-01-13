### Slide 6 ###

a<-3.06 ; b<-2.56; n<-20; x<-12 # n<-200; x<-120
an <- a+x; bn <- b+n-x; S<-10000; theta_post <- rbeta(S,an,bn)
median(theta_post); qbeta(0.5,an,bn)
posterior_exploss_grid = function(delta){
  grid <- seq(0,1,.001)
  risk <- sum(abs(delta-grid)*dbeta(grid,an,bn))*0.001}
delta = seq(0, 1, by = 0.005);
post_exploss_abs<-apply(as.matrix(delta),1,posterior_exploss_grid)
plot(delta, post_exploss_abs, type = 'l', col='red',
     lwd = 3, ylab ='posterior expected loss')
mtext(expression(paste("Absolute loss"," ",
                       italic(abs(theta-delta)),sep="")), side=3,line=.2)
delta[which.min(post_exploss_abs)]

### Slide 8 ###

loss_function = function(theta,delta){
  if (delta < theta){
    return(9*abs(theta - delta))} else{
      return(abs(theta - delta))}
}
posterior_exploss_mc = function(delta, S = 10000){
  theta = rbeta(S, 3.06+12, 20-12+2.56)
  loss <- apply(as.matrix(theta),1,loss_function,delta)
  risk = mean(loss)}
delta = seq(0, 1, by = 0.005); set.seed(11)
post_exploss <- apply(as.matrix(delta),1,posterior_exploss_mc)
plot(delta, post_exploss, type = 'l', col='blue',
     lwd = 3, ylab ='posterior expected loss')
delta[which.min(post_exploss)]; qbeta(0.9,3.06+12, 20-12+2.56)

### Slide 11 ###

a<-3.06 ; b<-2.56; n<-20; x<-12 # n<-200; x<-120
qbeta(c(.025,.975),an,bn); qbeta(c(.025,.975),a,b)
quantile(theta_post,c(.025,.975))
theta.grid<-seq(0,1,length=100)
plot(theta.grid,dbeta(theta.grid,an,bn),type="l",col="red",
     xlab=expression(theta),
     ylab=expression(paste(italic("p("),theta,"|x)")))
abline(v=qbeta(c(.025,.975), an,bn),lty=3)
segments(qbeta(0.025,an,bn),-0.14,qbeta(0.975,an,bn),-0.14,lwd=4)
text(0.58, 0.8, "0.95", cex = 1.4)
text(0.58, 0.5, "Posterior Probability", cex = 1.0)
mtext("Quantile-based Credible Interval", side=3,line=.2)
legend("topleft",c("Posterior","Credible Interval"),
       lwd=c(2,4),lty=c(1,1),col=c("red",1), cex=0.8)

### Slide 14 (WinBUGS) ###

model{
  x ~ dbin(theta,10)
  theta ~ dbeta(1,1)
}

### Slide 15 (WinBUGS) ###

model{
  x ~ dbin(theta,n) #likelihood
  theta ~ dbeta(alpha,beta) #prior
  prob <- step(theta - 0.396) * step(0.766 - theta) }
list(x=12,n=20,alpha=3.06,beta=2.56) #data
list(theta=0.1) #starting/initial value