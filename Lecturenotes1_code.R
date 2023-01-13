### Slide 5 ###

pmf <- c(0.2, 0.3, 0.1, 0.15, 0.05, 0.2)
rng.discrete1 <- function(pmf) {
  cdf <- cumsum(pmf); u <- runif(1); x <- sum(u > cdf)
  return(x) }
nsim <- 10^4; set.seed(10)
rng1.X <- replicate(nsim, rng.discrete1(pmf))
table(rng1.X)/nsim #sample(0:5,100,prob=pmf,rep=T)

rng.discrete2 <- function(pmf) {
  repeat{ y<-runif(1,max=length(pmf)); v <- runif(1)
  if (v < pmf[floor(y)+1]/max(pmf))
    break }
  return(floor(y))}
nsim <- 10^4; set.seed(10)
rng2.X <- replicate(nsim, rng.discrete2(pmf))
table(rng2.X)/nsim #sample(0:5,100,prob=pmf,rep=T)

### Slide 6 ###

n <- 10^3; u <- runif(n); x <- u^(1/3)
hist(x,prob=T,main=" "); curve(3*x^2,0,1,col="red",lwd=2,add=T)

nsim=2500;M=3
mv=runif(nsim,max=M); y=runif(nsim)
x=y[mv<3*y^2] #g(y)=1
plot(y,mv,xlab="y",ylab="MV g(y)",xlim=c(0,1),ylim=c(0,3))
u1=mv[mv<3*y^2]; points(x,u1,pch=19,col="red4")
z=seq(0,1,0.01);lines(z,dbeta(z,3,1),lwd=2,col="red4");abline(h=M)
rng.beta1 <-function(alpha) { # alpha>1
  while (TRUE) {
    mv = runif(1,max=alpha); y=runif(1)
    if (mv< alpha*y^(alpha-1)) return(y)}}
nsim <- 10^3; set.seed(1); x <- replicate(nsim,rng.beta1(3))
hist(x,nclass=10,prob=T); curve(3*x^2,0,1,col="red",lwd=2,add=T)

### Slide 12 ###

alpha <- c(0,0,1,0,0,0)
Pmat <- matrix(c(.5,.5,0,0, 0,0,.25,.5,.25,0,
                 0,0,0,.25,.5,.25,0,0,0,0,.25,.5,.25,0,0,0,0,
                 .25,.5,.25,0,0,0,0,.5,.5),nrow=6, ncol=6, byrow=TRUE)
Pmat; print(alpha %*% Pmat, digits = 5)
Markov.sim <- function(n, Pmat, x0) {
  sim <- numeric(n); m <- ncol(Pmat); temp <-x0
  for (i in 1:n) {
    nextstate <- sample(1:m,1,prob=Pmat[temp,])
    sim[i] <- nextstate; temp <-sim[i]}
  return(sim)}
Markov.sim(1,Pmat,3)
nsim <-10^4; Mat.sim <- replicate(nsim,Markov.sim(1,Pmat,3))
table(Mat.sim)/nsim

### Slide 13 ###

powermat <- function(mat,k) {
  if (k == 0) return (diag(dim(mat)[1]))
  if (k == 1) return(mat)
  if (k > 1) return( mat %*% powermat(mat, k-1))
}
print(alpha %*% powermat(Pmat,5),digits=5)
table(replicate(nsim,Markov.sim(5,Pmat,3)[5]))/nsim
powermat(Pmat,10)
powermat(Pmat,100)
powermat(Pmat,200)

### Slide 15 ###

pmf <- c(0.2, 0.3, 0.1, 0.15, 0.05, 0.2); set.seed(1)
nburnin = 2000; nsample = 10^4; nsim = nburnin + nsample
x = numeric(nsim); x[1] = sample(0:5, 1) # initialize the chain
for (i in 2:nsim) {
  y = sample(0:5, 1); rho = min(pmf[y+1]/pmf[x[i-1]+1], 1)
  x[i] = x[i-1] + (y-x[i-1]) * (runif(1)<rho)
}
x = x[(nburnin): nsim]; table(x)/nsample
library(markovchain); markovchainFit(data=x)
Rt <- markovchainFit(data=x)$estimate
Pmat <- round(Rt@transitionMatrix, 3)
powermat(Pmat,10); powermat(Pmat,100)

### Slide 16 ###

a=3; b=1; nburnin=2000; nsample = 1000; nsim=nburnin+nsample
x=rep(runif(1),nsim) # initialize the chain
for (i in 2:nsim){
  y=runif(1); rho=dbeta(y,a,b)/dbeta(x[i-1],a,b)
  x[i]=x[i-1] + (y-x[i-1])*(runif(1)<rho)
}
par(mfrow=c(1,2), mar=c(2,2,1,1))
hist(x[(nburnin+1):nsim], breaks=20,
     col="blue", cex.main=0.5, freq=FALSE)
curve(dbeta(x,a,b), col="sienna", lwd=2, add=TRUE)
xind = replicate(nsample,rng.beta1(3)) #rbeta(nsample,a,b)
hist(xind, breaks=20,col="grey", cex.main=0.5, freq=FALSE)
curve(dbeta(x,a,b), col="sienna", lwd=2, add=TRUE)
plot(x[(nburnin+1):nsim],type="l"); plot(xind, type="l")
acf(x[(nburnin+1):nsim]); acf(xind)