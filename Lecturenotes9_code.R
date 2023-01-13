### Slide 3 ###


xobs<-c(46,58,40,47,47,53,43,48,50,
         55,49,50,52,56,49,54,51,50,52,50)
xbar<-mean(xobs); SSE<-sum((xobs-xbar)^2); n<-length(xobs)
S <- 10^4; sigma2 <- SSE/rchisq(S, n-1)
theta <- rnorm(S, mean=xbar, sd=sqrt(sigma2)/sqrt(n))
posts = data.frame(theta, sigma2)

require(MASS)
empi<-kde2d(posts[,1],posts[,2],lims=c(48,52,0,40))

th.grid=seq(48,52,len=1000); sig2.grid=seq(0.01,40,len=1000)
joint.post = function(th,sig2){
 a=(SSE/2)^((n-1)/2)/gamma((n-1)/2)*sqrt(n/(2*pi))
 a*sig2^(-1-n/2)*exp(-SSE/(2*sig2))*exp(-n*(th-xbar)^2/(2*sig2))}
post.grid = outer(th.grid,sig2.grid,joint.post)

plot(posts,ylim=c(2,42),xlim=c(46,54))

contour(th.grid,sig2.grid,post.grid,xlab=expression(theta),
   ylab=expression(sigma^2),lwd=3,add=TRUE,col="blue")

contour(empi,lwd=3,col="red",add=TRUE)


xtilde=rnorm(S,theta,sqrt(sigma2)); mean(xtilde) #xbar
tval <- (xtilde-xbar)/(sqrt((1+1/n)*var(xobs)))
hist(tval,prob=T,nclass=30,ylim=c(0,0.5))
curve(dt(x,n-1),-5,5,add=T,col="red")


### Slide 4 ###

set.seed(21); w<-57; S<-10^4; sigma2 <- SSE/rchisq(S, n-1)
theta <- rnorm(S, mean=xbar, sd=sqrt(sigma2)/sqrt(n))
xtilde=rnorm(S,theta,sqrt(sigma2))
tval <- (xtilde-xbar)/(sqrt((1+1/n)*var(xobs)))
a=(w-xbar)/sqrt((1+1/n)*var(xobs)); 1-pt(a,n-1); mean(xtilde > w)

### Slide 7 ###

mu0<-1.9 ; tau02<-0.95^2; alpha <- 0.5; beta <- 0.005 # priors
x<-c(1.64,1.70,1.72,1.74,1.82,1.82,1.82,1.90,2.08) # data
n<-length(x) ; mean.x<-mean(x) ; var.x<-var(x)

G<-100 ; H<-100 # discrete approximation
mean.grid<-seq(1.505,2.00,len=G); prec.grid<-seq(1.75,175,len=H)
post.grid<-matrix(nrow=G,ncol=H)
for(g in 1:G){ for(h in 1:H) {
  post.grid[g,h]<- dnorm(mean.grid[g], mu0, sqrt(tau02)) *
    dgamma(prec.grid[h],alpha,beta) *
    prod( dnorm(x,mean.grid[g],1/sqrt(prec.grid[h])) )} }
post.grid<-post.grid/sum(post.grid)
mean.post<-apply(post.grid,1,sum)
prec.post<-apply(post.grid,2,sum)

par(mfrow=c(1,3),mar=c(2.75,2.75,.5,.5),mgp=c(1.70,.70,0))
image( mean.grid,prec.grid,post.grid,col=gray( (10:0)/10 ),
       xlab=expression(theta), ylab=expression(tilde(sigma)^2) )
plot(mean.grid,mean.post,type="l",xlab=expression(theta),
     ylab=expression( paste(italic("p("),
                            theta,"|",italic(y[1]),"...",italic(y[n]),")",sep="")))
plot(prec.grid,prec.post,type="l",xlab=expression(tilde(sigma)^2),
     ylab=expression( paste(italic("p("),
                            tilde(sigma)^2,"|",italic(y[1]),"...",italic(y[n]),")",sep="")))

### Slide 9 ###

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

### Slide 10 ###

par(mfrow=c(1,3),mar=c(2.75,2.75,.5,.5),mgp=c(1.70,.70,0))
m1<-5
plot(PHI[1:m1,],type="l",xlim=range(PHI[1:100,1]),
     ylim=range(PHI[1:100,2]),
     lty=1,col="gray",xlab=expression(theta),
     ylab=expression(tilde(sigma)^2))
text( PHI[1:m1,1], PHI[1:m1,2], c(1:m1) )
m1<-15
plot(PHI[1:m1,],type="l",xlim=range(PHI[1:100,1]),
     ylim=range(PHI[1:100,2]),
     lty=1,col="gray",xlab=expression(theta),
     ylab=expression(tilde(sigma)^2))
text( PHI[1:m1,1], PHI[1:m1,2], c(1:m1) )
m1<-100
plot(PHI[1:m1,],type="l",xlim=range(PHI[1:100,1]),
     ylim=range(PHI[1:100,2]),
     lty=1,col="gray",xlab=expression(theta),
     ylab=expression(tilde(sigma)^2))
text(PHI[1:m1,1], PHI[1:m1,2], c(1:m1))

### Slide 11 ###

par(mfrow=c(1,3),mar=c(2.75,2.75,.5,.5),mgp=c(1.70,.70,0))
sseq<-1:1000
image( mean.grid,prec.grid,post.grid,col=gray( (10:0)/10 ),
       xlab=expression(theta), ylab=expression(tilde(sigma)^2) ,
       xlim=range(PHI[,1]),ylim=range(PHI[,2]) )
points(PHI[sseq,1],PHI[sseq,2],pch=".",cex=3 )
plot(density(PHI[,1],adj=2), xlab=expression(theta),main="",
     xlim=c(1.55,2.05),
     ylab=expression( paste(italic("p("),
                            theta,"|",italic(y[1]),"...",italic(y[n]),")",sep="")))
plot(density(PHI[,2],adj=2), xlab=expression(tilde(sigma)^2),
     main="",ylab=expression( paste(italic("p("),
                                    tilde(sigma)^2,"|",italic(y[1]),"...",italic(y[n]),")",sep="")))
quantile(PHI[,1],c(.025,.5,.975))
quantile(PHI[,2],c(.025,.5,.975))

### Slide 12 ###

library(rjags)
y<-c(1.64, 1.70, 1.72, 1.74, 1.82, 1.82, 1.82, 1.82, 1.90, 2.08)
mu0<-0; tau02<-0.95^2; alpha<-0.01; beta<-0.01; n=length(y)
model_string <- "model{
for(i in 1:n){
y[i] ~ dnorm(theta,prec.var)
}
theta ~ dnorm(mu0,prec.var0)
prec.var0 <- 1/tau02
prec.var ~ dgamma(a, b)
sig2 <- 1/prec.var
}"
model <- jags.model(textConnection(model_string),n.chains = 3,
                    data = list(y=y,n=n,mu0=mu0,tau02=tau02,a=alpha,b=beta))
update(model,10000,progress.bar="none") #Burnin for 10000 samples
samp <- coda.samples(model, variable.names=c("theta","sig2"),
                     n.iter=20000, progress.bar="none")
summary(samp); plot(samp)

### Slide 16 ###


m <- 500; k <- 15; B <- 5
invcdfx<-function(u,y,B){ log(1-u*(1-exp(-B*y)))/(-y)}

rng_gibbs <- function(k,B){
  x.tmp <- y.tmp <- matrix(nrow=k+1)
  y.tmp[1] <- runif(1,0,5); u1<-runif(1,0,1)
  x.tmp[1] <- invcdfx(u1,y.tmp[1],B)
  for (j in 2:(k+1)) {
    u2<-runif(1,0,1); y.tmp[j] <- invcdfx(u2,x.tmp[j-1],B)
    u3<-runif(1,0,1); x.tmp[j] <- invcdfx(u3,y.tmp[j],B)
  }
  return(c(x.tmp[k+1],y.tmp[k+1]))
}


xy <- replicate(m, rng_gibbs(k,B));xy <- t(xy)
x <- xy[,1];y <- xy[,2]
hist(x,breaks=30,prob=TRUE)

require(pracma)
mar.x = function(x){
  fun <- function(x,y) exp(-x*y)
  const=1/integral2(fun,0,B,0,B)$Q
  const*(1/x*(1-exp(-B*x)))}
  
curve(mar.x(x),0,5,add=T,col="red")  
