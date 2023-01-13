### Slide 5 ###

n=m=100;x=10; nsim<-10^5; set.seed(13); a=1; b=1
theta=rbeta(nsim,a+x,b+n-x) # method of composition
tildex=rbinom(nsim,m,theta)

freq=table(tildex);x.pred=as.integer(names(freq))
predp=freq/sum(freq)
plot(x.pred,predp,type="h",ylab="Posterior Predictive PMF")
sum((tildex>=20)*1)/nsim

pbetabinom=function(a,b,x,n,m,predx){
  predp=0*predx; ax=a+x; bx=b+n-x
  lcon=lgamma(m+1)-lgamma(predx+1)-lgamma(m-predx+1)
  predp=exp(lcon+lbeta(predx+ax,m-predx+bx)-lbeta(ax,bx))
  return(predp)
}
predx <- 20:100; predp_exact<-pbetabinom(1,1,10,100,100,predx)
sum(predp_exact)

### Slide 6 (OpenBUGS) ###

model{ #betabinom.txt
  x ~ dbin(theta, n) # Model the data
  xtilde ~ dbin(theta, m) # Prediction of future binomial
  theta ~ dbeta(a, b) # The prior
  prob <- step(xtilde - 20) # Pred prob that xtilde >= 20
}
list(n = 100, m = 100, x = 10, a = 1, b = 1) #data
list(theta = 0.5, xtilde =10) #starting/initial values

list(theta = 0.1, xtilde =9) #starting/initial values for chain 2
list(theta = 0.9, xtilde =11) #starting/initial values for chain 3

### Slide 8 ###

library(R2OpenBUGS)
data <- list(n = 100, m = 100, x = 10, a = 1, b = 1)
inits <- function(){list(theta = 0.5, xtilde =10)}

betabinom.sim <- bugs(data, inits, model.file="betabinom.txt",
                      parameters=c("theta","xtilde"), n.chains=1, n.iter=20000,debug=T)
print(betabinom.sim)
plot(betabinom.sim)

betabinom.sim <- bugs(data, inits, model.file="betabinom.txt",
                      parameters=c("theta","xtilde"), n.chains=1, n.iter=20000,codaPkg=T)

codaobject <- read.bugs(betabinom.sim)
plot(codaobject)

### Slide 12 ###

library(rjags)
n <- 20; x <- 12; alpha <- 1; beta <- 1
model_string <- "model{
x ~ dbinom(theta,n) # Likelihood (sum of observed Bernoulli)
theta ~ dbeta(alpha, beta) # Prior
xtilde ~ dbern(theta) # Prediction of future Bernoulli
prob <- equals(xtilde,1) #Predictive probability
}"
betabin_model <- jags.model(textConnection(model_string),
                            data = list(x=x,n=n,alpha=alpha,beta=beta))
update(betabin_model, 10000); # Burnin for 10000 samples
postsamp <- coda.samples(betabin_model,
                         variable.names=c("theta","xtilde","prob"),n.iter=10000)
summary(postsamp); plot(postsamp)

### Slide 15 ###

a<-0.5; b<-0.5; n<-20; x<-sum(rpois(n,0.5)); grid<-seq(0.01,2,.01)
like<-dpois(x,n*grid); like<-like/sum(like)
prior<-dgamma(grid,a,b); prior<-prior/sum(prior)
post<-like*prior; post<-post/sum(post) #post<-dgamma(grid,x+a,n+b)
plot(grid,like,type="l",lty=2, col=1, xlab="theta",ylab="Density")
lines(grid,prior,col="blue"); lines(grid,post,lwd=2,col="red")
legend("topright",c("Likelihood","Prior","Posterior"),
       lwd=c(1,1,2),lty=c(2,1,1),col=c(1,"blue","red"))
(x+a)/(n+b); sqrt(x+a)/(n+b) #posterior mean & s.d.
qgamma(c(0.025,0.975), x+a, n+b) #posterior quantiles

### Slide 16 ###

n.samples <- 10000; set.seed(10)
postsamp <- rgamma(n.samples,x+a,n+b)

hist(postsamp,breaks=25,xlim=0:1,main="Posterior density",freq=F)
lines(density(postsamp),lty=2,col="red")
lines(grid,dgamma(grid,x+a,n+b), lty=1, col="green")
legend("topright",c("Density estimate","Exact density"),
       lwd=c(1,1),lty=c(2,1),col=c("red","green"))

mean(postsamp); sd(postsamp)
quantile(postsamp,c(0.025,0.975))

### Slide 17 ###

gamma.sim <- function(alpha) {
  fx <- function(x) x^(alpha-1)*exp(-x)/gamma(alpha)
  g <- function(x) (1/alpha)*exp(-x/alpha)
  #M <- 3^(3/2)/sqrt(2*pi*exp(1))
  M <-optimize(f=function(x){fx(x)/g(x)},
               maximum=T,interval=c(0,100))$objective
  while (TRUE) {
    Y <- -log(runif(1))*alpha; V <- runif(1, 0, M*g(Y))
    if (V < fx(Y)) return(Y)} }
nsim<-10^4; alpha=3/2; gamma.rng<-replicate(nsim,gamma.sim(alpha))
hist(gamma.rng,nclass=20,freq=F,ylim=c(0,0.5),col="gray")
curve(dgamma(x,alpha,1),0,10,add=TRUE, col="red4",lwd=2)
#postsamp<-replicate(nsim,gamma.sim(x+a))/(n+b)