### Slide 5 ###

library(rjags)
jags_code ="model{
for( i in 1:n) { logage[i] <- log(x[i])
y[i] ~ dnorm(mu[i] , tau); mu[i] <- beta0+ beta1*logage[i] }
beta0 ~ dnorm(0, 0.001); beta1 ~ dnorm(0, 0.001)
tau ~ dgamma(0.01, 0.01); sigma <- 1/sqrt(tau)}"
dugong = list(x=c(1.0, 1.5, 1.5, 1.5, 2.5, 4.0, 5.0, 5.0, 7.0,
                  8.0, 8.5, 9.0, 9.5, 9.5, 10.0, 12.0, 12.0, 13.0,
                  13.0, 14.5, 15.5, 15.5, 16.5, 17.0, 22.5, 29.0, 31.5),
              y = c(1.80, 1.85, 1.87, 1.77, 2.02, 2.27, 2.15, 2.26, 2.47,
                    2.19, 2.26, 2.40, 2.39, 2.41, 2.50, 2.32, 2.32, 2.43,
                    2.47, 2.56, 2.65, 2.47, 2.64, 2.56, 2.70, 2.72, 2.57), n=27)
jags_reg = jags.model(textConnection(jags_code), data=dugong)
update(jags_reg, 10000) #progress.bar="none")
samp <- coda.samples(jags_reg,variable.names=
                       c("beta0","beta1","sigma"), n.iter=10000)
summary(samp); plot(samp)
lm(dugong$y ~ log(dugong$x)) #classical OLS

### Slide 12 ###

y=c(8.4,9.5,11.8,10.4,13.3,14.8,13.2,14.7,16.4,16.5,18.9,18.5)
x=c(20,22,24,26,28,30,32,34,36,38,40,42)
plot(x,y); OLS <- lm(y ~ x); summary(OLS)
library(MASS); p<-2; n<-length(y); S<-10^4; X<-cbind(rep(1,n),x)
beta.hat <- solve(t(X)%*%X)%*%t(X)%*%y
s2 <- t(y)%*%(diag(n)-X%*%solve(t(X)%*%X)%*%t(X))%*%y/(n-p)
vec.b <- matrix(0,nrow=p,ncol=S); vec.sigma2 <- rep(0,S)
for (s in 1:S){
  tau <- rgamma(1,(n-2)/2,(n-2)*s2/2)
  vec.sigma2[s] <- 1/tau
  b <- mvrnorm(1,beta.hat,solve(t(X)%*%X)/tau)
  vec.b[,s] <- b}
apply(vec.b,1,mean);beta.hat;mean(1/vec.sigma2);1/s2

vec.new <- rep(0,S); xnew <- 35
for (s in 1:S){
  tau <- rgamma(1,(n-2)/2,(n-2)*s2/2)
  b <- mvrnorm(1,beta.hat,solve(t(X)%*%X)/tau)
  vec.new[s] <- rnorm(1, mean=b[1]+b[2]*xnew, sd=sqrt(1/tau))}
mean(vec.new); predict(OLS,newdata=data.frame(x=xnew))
quantile(vec.new,c(0.025,0.975))

### Slide 15 ###

y=c(8.4,9.5,11.8,10.4,13.3,14.8,13.2,14.7,16.4,16.5,18.9,18.5)
x=c(20,22,24,26,28,30,32,34,36,38,40,42); OLS<-lm(y~x)
n=length(y); n.iters <- 30000; postsamp <- matrix(0,n.iters,3)
colnames(postsamp)<-c("beta0","beta1","sigma2")
beta0<-OLS$coef[1]; beta1<-OLS$coef[2]; s2<-var(OLS$residuals)
postsamp[1,] <- c(beta0,beta1,s2)
mu0 <- 0; s20 <- 1000; a0 <- 0.01; b0 <- 0.01
for(i in 2:n.iters){
  V <- n/s2+mu0/s20; M <- sum(y-x*beta1)/s2+1/s20
  beta0 <- rnorm(1,M/V,1/sqrt(V)) # posterior sample of beta0
  V <- sum(x^2)/s2+mu0/s20; M <- sum(x*(y-beta0))/s2+1/s20
  beta1 <- rnorm(1,M/V,1/sqrt(V)) # posterior sample of beta1
  an <- n/2 + a0; bn <- sum((y-beta0-x*beta1)^2)/2 + b0
  s2 <- 1/rgamma(1,an,bn) #posterior sample of sigma2
  postsamp[i,] <- c(beta0,beta1,s2)
}
pairs(postsamp); res<-matrix(0,3,4)
rownames(res)<-c("Intercept","Slope","sigma2")
colnames(res) <- c("Mean","SD","Q025","Q975")
res[,1]<-apply(postsamp,2,mean); res[,2]<-apply(postsamp,2,sd)
res[,3]<-apply(postsamp,2,quantile,0.025)
res[,4]<-apply(postsamp,2,quantile,0.975)
print(res,digits=3)

### Slide 20 ###

library(LearnBayes); data(donner)
attach(donner); X=cbind(1,age,male)
fit=glm(survival~X-1,family=binomial(link=probit))
summary(fit)

y<-survival; p<-3; n<-length(y)
beta0<-rep(0,p); Sigma0<-diag(100,p)
beta<-beta0; nsim<-3000; nwarm<-1000;
beta.post<-matrix(0, nsim, p)
Sigma.beta<-solve(t(X)%*%X+solve(Sigma0)) #posterior variance
z<-rep(0, n)

library(MASS)
install.packages("msm")
library(msm) # a package with the truncated normal distribution
mvrnorm; rtnorm

### Slide 21 ###

for (isim in 1:nsim+nwarm){
  for (i in 1:n){ # generate latent variables z
    if(y[i]==1) z[i]<-rtnorm(1,t(X[i,])%*%beta, 1, 0, Inf)
    else
      z[i]<-rtnorm(1,t(X[i,])%*%beta, 1, -Inf, 0)
  }
  Mu.beta<-Sigma.beta%*%(t(X)%*%z+solve(Sigma0)%*%beta0)
  beta<-mvrnorm(1, Mu.beta, Sigma.beta)
  if(isim>nwarm) beta.post[isim-nwarm,]<-beta}

par(mfrow=c(1,3))
plot(beta.post[,1], type="l"); plot(beta.post[,2], type="l")
plot(beta.post[,3], type="l")
acf(beta.post[,1]); acf(beta.post[,2]); acf(beta.post[,3])

beta.est<-c(0,0)
for(i in 1:p) beta.est[i]<-mean(beta.post[1:nsim, i])
beta.est
plot(density(beta.post[1:nsim, 1]), xlab="beta_1", main="")
abline(v=beta.est[1], lty=2)
plot(density(beta.post[1:nsim, 2]), xlab="beta_2", main="")
abline(v=beta.est[2], lty=2)
plot(density(beta.post[1:nsim, 3]), xlab="beta_3", main="")
abline(v=beta.est[3], lty=2)

### Slide 22 ###

library(LearnBayes); data(donner)
attach(donner); X=cbind(1,age,male)
fit=glm(survival~X-1,family=binomial(link=probit))
summary(fit)

probit_model <- "model{
for(i in 1:n){
Y[i] ~ dbin(p[i],1)
p[i] <- phi(beta[1] + beta[2]*x1[i]+beta[3]*x2[i])}
for(j in 1:3){ beta[j] ~ dnorm(0,0.01)}
}"
library(rjags)
dat <- list(Y=survival,n=length(survival),x1=age,x2=male)
model1 <- jags.model(textConnection(probit_model),
                     data = dat,n.chains=3, quiet=TRUE)
update(model1, 10000, progress.bar="none")
samp <- coda.samples(model1, variable.names=c("beta"),
                     n.iter=20000, progress.bar="none")
summary(samp); plot(samp); effectiveSize(samp)