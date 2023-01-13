### Slide 5 ###

mu<-c(-3,0,3); s2<-c(.33,.33,.33); w<-c(.45,.1,.45)
ths<-seq(-5,5,length=100); set.seed(1); S<-2000
d<-sample(1:3,S, prob=w,replace=TRUE)
th<-rnorm(S,mu[d],sqrt(s2[d]))
THD.MC<-cbind(th,d)
par(mfrow=c(1,3))
plot(THD.MC[,1], type="l", main='MC draws', xlab='Iteration number')
hist(THD.MC[,1], freq = FALSE, main=' ',
     ylim = c(0,0.5), xlab=expression(theta[1]))
ths<-seq(-6,6,length=1000)
lines( ths, w[1]*dnorm(ths,mu[1],sqrt(s2[1])) +
         w[2]*dnorm(ths,mu[2],sqrt(s2[2])) +
         w[3]*dnorm(ths,mu[3],sqrt(s2[3])),lwd=2 )
acf(THD.MC[,1],main=' ',lag.max = 20)

### Slide 6 ###

th<-0; THD.MCMC<-NULL; S<-10000; set.seed(1)
for(s in 1:S) {
  d<-sample(1:3 ,1,prob= w*dnorm(th,mu,sqrt(s2)))
  th<-rnorm(1,mu[d],sqrt(s2[d]) )
  THD.MCMC<-rbind(THD.MCMC,c(th,d) ) }
par(mfrow=c(1,3))
plot(THD.MCMC[,1], type="l", main='Gibbs draws',
     xlab='Iteration number')
hist(THD.MCMC[,1], freq = FALSE, main=' ',ylim = c(0,0.5))
ths<-seq(-6,6,length=1000)
lines( ths, w[1]*dnorm(ths,mu[1],sqrt(s2[1])) +
         w[2]*dnorm(ths,mu[2],sqrt(s2[2])) +
         w[3]*dnorm(ths,mu[3],sqrt(s2[3])),lwd=2 )
acf(THD.MCMC[,1], main=' ',lag.max = 20)

### Slide 11 ###

library(rjags)
x=c(5,1,5,14,3,19,1,1,4,22); k=length(x)
t=c(94.320,15.72,62.88,125.76,5.24,31.44,1.048,1.048,2.096,10.48)
model_string <- textConnection("model{
for(i in 1:k){
x[i] ~ dpois(theta[i]*t[i]) # Likelihood
}
for(i in 1:k){
theta[i] ~ dgamma(1,beta) # Priors
}
beta ~ dgamma(0.1, 1)
}")
inits<-list(theta=rgamma(k,1,1),beta=1); data<-list(x=x,k=k,t=t)
model <- jags.model(model_string,data=data,inits=inits)
update(model, 10000); params <- c("theta","beta")
samples <- coda.samples(model, variable.names=params,
                        n.iter=20000, progress.bar="none")
summary(samples); plot(samples)

### Slide 12 ###

x=c(5,1,5,14,3,19,1,1,4,22); k=length(x)
t=c(94.320,15.72,62.88,125.76,5.24,31.44,1.048,1.048,2.096,10.48)
S <- 25000; post_samp <- matrix(NA,S,k+1); c <- 0.1; d <- 1
theta <- x/t; beta <- 1 #initial values
for(s in 1:S){ # Gibbs sampling
  for(i in 1:k){
    theta[i] <- rgamma(1,x[i]+1,t[i]+beta)
  }
  beta <- rgamma(1,k+c,d+sum(theta))
  post_samp[s,] <- c(beta,theta)
}
boxplot(post_samp[,2:11],outline=FALSE,ylab=expression(theta))
par(mfrow=c(1,2))
plot(post_samp[,1],type="l",xlab="Iteration",ylab=expression(beta))
hist(post_samp[,1], breaks=30, freq=F, xlab=expression(beta))
lines(density(post_samp[,1]),lwd=2,col="red")

### Slide 14 ###

ar <- function(n,theta) {
  y <- numeric(n); y[1] <- rnorm(1)
  for (i in 2:n) { y[i] <- rnorm(1,theta*y[i-1],1)}
  return(y)}
library(coda); set.seed(1); n<-1000; theta <- c(0.0,0.5,0.75,0.99)
par(mfrow=c(2,2))
for (i in 1:4) { y <- ar(n,theta[i])
plot(1:n,y,type='l',xlab="t",
     main=bquote(theta==.(theta[i])~
                   (n[eff]==.(round(effectiveSize(y))))))}

### Slide 15 ###

library(coda)
mu0<-1.9 ; tau02<-0.95^2; alpha <- 0.5; beta <- 0.005 # priors
x<-c(1.64,1.70,1.72,1.74,1.82,1.82,1.82,1.90,2.08) # data
n<-length(x) ; mean.x<-mean(x) ; var.x<-var(x)
set.seed(1)
S<-1000
PHI<-matrix(nrow=S,ncol=2)
PHI[1,]<-phi<-c( mean.x, 1/var.x)
### Gibbs sampling
for(s in 2:S) {
  # generate a new theta value from its full conditional
  mun<- ( mu0/tau02 + n*mean.x*phi[2] ) / ( 1/tau02 + n*phi[2] )
  sig2n<- 1/( 1/tau02 + n*phi[2] )
  phi[1]<-rnorm(1, mun, sqrt(sig2n) )
  # generate a new sigma^2 value from its full conditional
  alphan<- alpha+n/2
  betan<- beta + ((n-1)*var.x + n*(mean.x-phi[1])^2)/2
  phi[2]<- rgamma(1, alphan, betan)
  PHI[s,]<-phi}

acftheta <- acf(PHI[,1])
acfsigma2 <- acf(1/PHI[,2])
effectiveSize( PHI[,1])
effectiveSize(1/PHI[,2])

### Slide 16 ###

library(R2jags); set.seed(3)
y<-c(1.64,1.70,1.72,1.74,1.82,1.82,1.82,1.82,1.90,2.08)
mu0 <- 0; sig20<-100; alpha<-0.01; beta<-0.01; n<-length(y)
model.script = cat("model{
for(i in 1:n){
y[i] ~ dnorm(theta,prec.var)}
theta ~ dnorm(mu0,prec.var0)
prec.var0 <- 1/sig20
prec.var ~ dgamma(a, b)
sig2 <- 1/prec.var
}",file="semiconjugate.txt")
semiconjugate_res <- jags(data = list("y"=y,"n"=n,"mu0"=mu0,
                                      "sig20"=sig20,"a"=alpha,"b"=beta),
                          parameters.to.save=c("theta","sig2"),
                          model.file="semiconjugate.txt",n.chains = 3, n.burnin=500,
                          n.thin=1, n.iter=1000)
traceplot(semiconjugate_res)
semiconjugate_res

### Slide 17 ###

library(coda); par(mfrow=c(1,3)); effectiveSize(THD.MCMC[,1])
S<-100000; set.seed(1); index=seq(10000,100000,100)
for(s in 1:S) { d<-sample(1:3 ,1,prob= w*dnorm(th,mu,sqrt(s2)))
th<-rnorm(1,mu[d],sqrt(s2[d]));THD.MCMC<-rbind(THD.MCMC,c(th,d))}
effectiveSize(THD.MCMC[index,1])
plot(THD.MCMC[index,1], type="l", main='Gibbs draws')
hist(THD.MCMC[index,1], freq = FALSE, main=' ',ylim = c(0,0.5))
ths<-seq(-6,6,length=1000)
lines(ths,w[1]*dnorm(ths,mu[1],sqrt(s2[1]))+w[2]*dnorm(ths,mu[2],
                                                       sqrt(s2[2]))+w[3]*dnorm(ths,mu[3],sqrt(s2[3])),lwd=2)
acf(THD.MCMC[index,1], main=' ',lag.max = 20)