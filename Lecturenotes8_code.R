### Slide 3 ###

n <- 10; x <- c(3.5, 4.4, 3.4, 3.8, 2.1, 3.4, 3.9, 4.4, 2.7, 2.6)
theta <- 3.5; a <- 0.1; b <- 0.1; like <-1; grid <- seq(0.01,5,.01)
for(i in 1:n){like <- like*dnorm(x[i],theta,1/sqrt(grid)) }
like <- like/sum(like)
prior <- dgamma(grid,a,b); prior <- prior/sum(prior)
post <- like*prior; post <- post/sum(post)
plot(grid,post,type="l",lwd=2, xlab="tau",ylab="Density", col="red")
lines(grid,prior, col="blue"); lines(grid,like,lty=2)
legend("topright",c("Likelihood","Prior","Posterior"),
       lwd=c(1,1,2),lty=c(2,1,1),col=c(1,"blue","red"))

### Slide 4 ###

SSE <- sum((x-theta)^2)
postm <- (n/2+a)/(SSE/2+b); postvar <- (n/2+a)/(SSE/2+b)^2
qgamma(c(0.025,0.975),n/2+a,SSE/2+b)

library(rjags)
model_string <- "model{
for(i in 1:n){
y[i] ~ dnorm(theta,tau)
}
tau ~ dgamma(a,b)
sigma <- pow(tau,-0.5)
}"
model <- jags.model(textConnection(model_string),
                    data = list(y=x,n=n,theta=theta,a=a,b=b))
update(model, 10000, progress.bar="none")
postsamp <- coda.samples(model, variable.names=c("tau","sigma"),
                         n.iter=20000, progress.bar="none")
summary(postsamp)
plot(postsamp)

### Slide 9 ###

mu0<-1.9 ; k0<-1; alpha<-0.5; beta<-.005 # prior
x<-c(1.64,1.70,1.72,1.74,1.82,1.82,1.82,1.90,2.08) # data
n<-length(x) ; xbar<-mean(x) ; s2<-var(x)

### Slide 10 ###

alphan<-alpha+n/2; mun<-(k0*mu0+n*xbar)/(k0+n)
betan<- beta+1/2*(n-1)*s2+ n*k0*(xbar-mu0)^2/(2*(n+k0))
dinvgamma<-function(x,a,b){
  ld<-a*log(b)-lgamma(a)-(a+1)*log(x)-b/x
  exp(ld)}
gs<-100; theta<-seq(1.6,2.0,length=gs)
is2<-seq(15,160 ,length=gs); s2g<-seq(.001,.045,length=gs)
ld.th.is2<-ld.th.s2<-matrix(0,gs,gs)
for(i in 1:gs) { for(j in 1:gs) {
  ld.th.is2[i,j]<-dnorm(theta[i],mun,1/sqrt(is2[j]*10),log=TRUE)+
    dgamma(is2[j],alphan,betan, log=TRUE )
  ld.th.s2[i,j]<-dnorm(theta[i],mun,sqrt(s2g[j]/10),log=TRUE)+
    log(dinvgamma(s2g[j],alphan,betan))
}}
par(mfrow=c(1,2),mgp=c(1.75,.75,0));gr<-gray((10:0)/10)
image(theta,is2,exp(ld.th.is2),col=gr,xlab=expression(theta),
      ylab=expression(tau))
image(theta,s2g,exp(ld.th.s2), col=gr,xlab=expression(theta),
      ylab=expression(sigma^2) )

### Slide 11 ###

set.seed(1); N<-10000
sig2.post<-1/rgamma(N, alphan, betan)
theta.post<-rnorm(N, mun, sqrt(sig2.post/(k0+n)))
quantile(theta.post,c(.025,.975)); quantile(sig2.post,c(.025,.975))
layout(matrix(c(1,1,2,3),2,2,byrow=T))
par(mar=c(3,3,1,1),mgp=c(1.75,.75,0))
image(theta,s2g,exp(ld.th.s2), col=gr,xlab=expression(theta),
      ylab=expression(sigma^2),xlim=c(1.60,2.0),ylim=c(.001,.07))
points(theta.post[1:5000], sig2.post[1:5000],pch=".",
       xlab=expression(theta),ylab=expression(sigma^2),
       xlim=c(1.65,1.95),ylim=c(.005,.07))
plot(density(sig2.post,adjust=3),main="",xlab=expression(sigma^2),
     xlim=c(0,.075),ylab=expression( paste(italic("p("),
                                           sigma^2,"|",italic(x[1]),"...",italic(x[n]),")",sep="")))
plot(density(theta.post,adjust=3),main="",xlab=expression(theta),
     xlim=c(1.60,2.0),ylab=expression( paste(italic("p("),
                                             theta,"|",italic(x[1]),"...",italic(x[n]),")",sep="")))