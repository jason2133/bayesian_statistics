
### Slide 7 ###

library(mnormt); rho=0.99; Sig=matrix(c(1,rho,rho,1),2,2)
theta1=seq(-3,3,0.1); theta2=seq(-3,3,0.1);
n=length(theta1); z=matrix(0,n,n)
for (i in 1:n){z[i,]=dmnorm(cbind(theta1[i],theta2),c(0,0),Sig)}
par(cex=1.2); par(tcl=-0.25); par(mgp=c(4,0.6,0))
contour(theta1,theta2,z,levels=seq(0.005,0.005),
        labels=" ",lwd=2,axes=F)
axis(1,at=c(-4,-3,-1.5,0,1.5,3),labels=rep("",6))
axis(2,at=c(-4,-3,-1.5,0,1.5,3),labels=rep("",6))
mtext(expression(theta[2]),side=2,cex=1.5,line=0.5)
mtext(expression(theta[1]),side=1,cex=1.5,line=0.5)

iter=101; theta1=numeric(iter); theta2=numeric(iter)
for (i in 2:iter){
  theta1[i]=rnorm(1,rho*theta2[i-1],sqrt(1-rho^2))
  theta2[i]=rnorm(1,rho*theta1[i],sqrt(1-rho^2))}
points(theta1,theta2,pch=19,lwd=1.5)
for (i in 1:iter){
  segments(theta1[i],theta2[i],theta1[i+1],theta2[i],col="red")
  segments(theta1[i+1],theta2[i],theta1[i+1],theta2[i+1],col="red")}

hist(theta1,nclass=30,freq=F); curve(dnorm(x,0,1),-3,3,add=T)

### Detailed contour plot ###

theta1=seq(-3,3,0.1); theta2=seq(-3,3,0.1);
n=length(theta1); z=matrix(0,n,n)
for (i in 1:n){z[i,]=dmnorm(cbind(theta1[i],theta2),c(0,0),Sig)}
par(cex=1.2); par(tcl=-0.25); par(mgp=c(4,0.6,0))
contour(theta1,theta2,z,levels=c(0.005,seq(0.2,1,0.2)),
        labels=" ",lwd=2,axes=F)
axis(1,at=c(-4,-3,-1.5,0,1.5,3),labels=rep("",6))
axis(2,at=c(-4,-3,-1.5,0,1.5,3),labels=rep("",6))
mtext(expression(theta[2]),side=2,cex=1.5,line=0.5)
mtext(expression(theta[1]),side=1,cex=1.5,line=0.5)
iter=101; theta1=numeric(iter); theta2=numeric(iter)
for (i in 2:iter){
  theta1[i]=rnorm(1,rho*theta2[i-1],sqrt(1-rho^2))
  theta2[i]=rnorm(1,rho*theta1[i],sqrt(1-rho^2))}
points(theta1,theta2,pch=19,lwd=1.5)
for (i in 1:iter){
  segments(theta1[i],theta2[i],theta1[i+1],theta2[i],col="red")
  segments(theta1[i+1],theta2[i],theta1[i+1],theta2[i+1],col="red")}

### Slide 12 ###

x1 <- 1; x2 <- -1; rho <- 0.9; x <- c(x1,x2)
Sigma = matrix(c(1,rho,rho,1),2,2); nDraws <- 10000

library(mnormt); MC.samples <- rmnorm(nDraws, x, Sigma)
gibbs.samples <- matrix(0,nDraws,2); theta2 <- 0
for (i in 1:nDraws){
  theta1 <- rnorm(1, x1 + rho*(theta2-x2), sqrt(1-rho^2))
  gibbs.samples[i,1] <- theta1
  theta2 <- rnorm(1, x2 + rho*(theta1-x1), sqrt(1-rho^2))
  gibbs.samples[i,2] <- theta2
}
apply(MC.samples[9000:10000,],2,mean);
apply(gibbs.samples[9000:10000,],2,mean); x
cov(MC.samples[9000:10000,]);cov(gibbs.samples[9000:10000,]);Sigma

index <- seq(3000,10000,10)
apply(gibbs.samples[index,],2,mean); x
cov(gibbs.samples[index,]); Sigma

### Slide 13 ###

par(mfrow=c(2,4))
plot(MC.samples[,1], type="l", main='MC draws',
     xlab="Iteration number", ylab=expression(theta[1]))
hist(MC.samples[,1], freq = FALSE, main='MC draws',
     ylim = c(0,0.5), xlab=expression(theta[1]))
curve(dnorm(x,x1,1),-2,4,add=T,col="red",lwd=3)
plot(cumsum(MC.samples[,1])/seq(1,nDraws),type="l",main='MC draws',
     xlab='Iteration number', ylab='cumulative mean of theta')
lines(seq(1,nDraws),1*matrix(1,1,nDraws),col="red",lwd=3)
acf(MC.samples[,1], main='MC draws', lag.max = 20)
plot(gibbs.samples[,1], type="l", main='Gibbs draws',
     xlab='Iteration number', ylab=expression(theta[1]))
hist(gibbs.samples[,1], freq = FALSE, main='Gibbs draws',
     ylim = c(0,0.5), xlab=expression(theta[1]))
curve(dnorm(x,x1,1),-2,4,add=T,col="red",lwd=3)
plot(cumsum(gibbs.samples[,1])/seq(1,nDraws),type="l",
     main='Gibbs draws',xlab='Iteration number',
     ylab='cumulative mean of theta')
lines(seq(1,nDraws),1*matrix(1,1,nDraws),col="red",lwd=3)
acf(gibbs.samples[,1], main='Gibbs draws', lag.max = 20)

### Slide 16 ###

x = c(4,5,4,1,0,4,3,4,0,6,3,3,4,0,2,6,3,3,5,4,5,3,1,4,4,1,5,5,
      3,4,2,5,2,2,3,4,2,1,3,2,2,1,1,1,1,3,0,0,1,0,1,1,0,0,3,1,0,3,2,
      2,0,1,1,1,0,1,0,1,0,0,0,2,1,0,0,0,1,1,0,2,3,3,1,1,2,1,1,1,1,2,
      4,2,0,0,0,1,4,0,0,0,1,0,0,0,0,0,1,0,0,1,0,1)
year=1851:1962; n = length(x)
plot(year[1:41],x[1:41],col=3,pch=1,xlim=c(1850,1962),
     ylim=c(-1,10),xlab="years",ylab="counts")
par(new=T)
plot(year[42:112],x[42:112],col=4,pch=1,xlim=c(1850,1962),
     ylim=c(-1,10),xlab="years",ylab="counts"); abline(v=1891,col=2)
text(1922,8.5, "change point : around 1891")
text(1910,4.7,paste("Mean before 1891 = ",
                    round(mean(x[1:41]),3),sep="")); abline(h=3.908,lty=3)
text(1910,1.4,paste("Mean after 1891 = ",
                    round(mean(x[42:n]),3),sep="")); abline(h=0.901,lty=3)

### Slide 18 ###

alpha = 0.001; beta = 0.001; delta = 0.001; gamma = 0.001; m = 10
set.seed(123456); M=2000; postsamp = NULL; par(mfrow=c(3,2))
fullcm1 = function(m,lambda,phi,x,n,alpha,beta,gamma,delta){
  lambda^(alpha-1+sum(x[1:m]))*exp(-(beta+m)*lambda)*
    phi^(gamma-1+sum(x)-sum(x[1:m]))*exp(-(delta+n-m)*phi)}
for (i in 1:M){
  lambda = rgamma(1,sum(x[1:m])+alpha,m+beta)
  phi = rgamma(1,sum(x)-sum(x[1:m])+gamma,n-m+delta)
  fullcm = NULL
  for (j in 1:n) {
    fullcm=c(fullcm,fullcm1(j,lambda,phi,x,n,alpha,beta,gamma,delta))
  }
  fullcm=fullcm/sum(fullcm); m = sample(1:n,size=1,prob=fullcm)
  postsamp = rbind(postsamp,c(lambda,phi,m))
}
post.samp <- data.frame(postsamp);names <- c("lambda","phi","m")
summary(post.samp)
for (i in 1:2){
  ts.plot(post.samp[,i],xlab="iteration",ylab=names[i])
  hist(post.samp[,i],xlab="",main=names[i],prob=T)}
ts.plot(post.samp[,3],xlab="iteration",ylab=names[3])
plot(table(post.samp[,3]),type="h",xlab="",main=names[3])

### Slide 19 ###

library(rjags)
jags_code = "model{
for(year in 1 : n){ D[year] ~ dpois(lambda[year])
log(lambda[year]) <- b[1]+step(year - changeyear)*b[2]}
for (j in 1:2) {b[j] ~ dnorm(0.0,1.0E-6)}
changeyear ~ dunif(1,n)
lambda1 <- exp(b[1])
phi1 <- exp(b[1]+b[2])}"
xdata=list(n=112, D=c(4,5,4,1,0,4,3,4,0,6,3,3,4,0,2,6,3,3,5,
                      4,5,3,1,4,4,1,5,5,3,4,2,5,2,2,3,4,2,1,3,2,1,1,1,1,1,3,0,0,1,0,1,1,
                      0,0,3,1,0,3,2,2,0,1,1,1,0,1,0,1,0,0,0,2,1,0,0,0,1,1,0,2,2,3,1,1,2,
                      1,1,1,1,2,4,2,0,0,0,1,4,0,0,0,1,0,0,0,0,0,1,0,0,1,0,0))
init_value = function() {list(b = c(0, 0) , changeyear = 50)}
model<-jags.model(textConnection(jags_code),n.chains=1,data=xdata)
update(model, 10000, progress.bar="none")
samp <- coda.samples(model, variable.names=c("b","changeyear"),
                     n.iter=20000, progress.bar="none")
summary(samp); plot(samp)