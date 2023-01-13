### Slide 8 ###

sig = 40; n=1; xobs = 850 # observed data with n=1
m0A=900; k0A=16; t0A = sqrt(sig^2/k0A)
wnA=n/(n+k0A); mnA=wnA*xobs+(1-wnA)*m0A; tnA=sqrt(sig^2/(n+k0A))
m0B=900; k0B=1/16; t0B = sqrt(sig^2/k0B)
wnB=n/(n+k0B); mnB=wnB*xobs+(1-wnB)*m0B; tnB=sqrt(sig^2/(n+k0B))
curve(dnorm(x,m0A,t0A),500,1100,col="blue",xlab="theta",ylab=" ")
curve(dnorm(x,mean=xobs,sd=sig),500,1100,add=T)
points(xobs,0,pch=16,cex=2)
curve(dnorm(x,mnA,tnA),500,1100,col="blue",add=TRUE,lty=2,lwd=2)
curve(dnorm(x,m0B,t0B),500,1100,col="red",add=TRUE)
curve(dnorm(x,mnB,tnB),500,1100,col="red",add=TRUE,lty=2,lwd=2)

### Slide 10 ###

sig = 40; n=1; xobs = 850 # observed data with n=1
mu0=900; tau0 <- 100 #k0A=16; t0A = sqrt(sig^2/k0A)
grid <- seq(500,1100,.5)
like <- dnorm(xobs,grid,sig); like <- like/sum(like)
prior <- dnorm(grid,mu0,tau0); prior <- prior/sum(prior)
post <- like*prior; post <- post/sum(post)
plot(grid,like,type="l",lty=2,ylim=c(0,0.01),
     col=1, xlab="theta",ylab=""); points(xobs,0,pch=16,cex=2)
lines(grid,prior,col="blue"); lines(grid,post,lwd=2,col="red")
legend("topright",c("Likelihood","Prior","Posterior"),
       lwd=c(1,1,2),lty=c(2,1,1),col=c(1,"blue","red"))

### Slide 11 ###

sig = 40; n=1; xobs = 850 # observed data with n=1
mu0=900; tau0 <- 100 #k0A=16; t0A = sqrt(sig^2/k0A)
taun2 <- 1/(sig^{-2}+tau0^{-2}); taun <- sqrt(taun2)
mun <- taun2*(xobs/sig^2+mu0/tau0^2)

c(mun+qnorm(0.025)*taun, mun+qnorm(0.975)*taun)

library(rjags)
model_string <- "model{
x ~ dnorm(theta,prec.samp)
prec.samp <- pow(sigma,-2)
theta ~ dnorm(mu0,prior.prec)
prior.prec <- pow(tau0,-2)
}"
model <- jags.model(textConnection(model_string),
                    data = list(x=xobs,sigma=sig,mu0=mu0,tau0=tau0))
update(model, 10000, progress.bar="none")
samp <- coda.samples(model, variable.names=c("theta"),
                     n.iter=20000, progress.bar="none")
summary(samp); plot(samp)

### Slide 14 ###

x=c(1.64,1.70,1.72,1.74,1.82,1.82,1.82,1.90,2.08); n<-length(x)
mu0<-1.9;tau0<-(.5*1.9);tau02<-tau0^2;xbar<-mean(x);sigma2<-var(x)
mun<-(mu0/(tau02)+n*xbar/sigma2)/(1/tau02+n/sigma2)
taun2<-1/(1/tau02+n/sigma2)
qnorm(c(.025,.975),mun,sqrt(taun2)); 1-pnorm(1.9,mun,sqrt(taun2))
set.seed(2); S<-10^5; postsamp<-rnorm(S,mun,sqrt(taun2))
quantile(postsamp,c(.025,.975)); mean(postsamp>1.9)
xs<-seq(0,mu0*2,length=500)
plot(xs,dnorm(xs,mun,sqrt(taun2)),type="l",xlab=expression(theta),
     ylab=expression(paste(italic("p("),theta,"|",italic(x[1]),
                           "...",italic(x[n]),",",sigma^2==0.017,")",sep="")),lwd=2)
lines(xs,dnorm(xs,mu0,sqrt(tau02)),type="l",col="gray",lwd=2)

pred_theta_sim<-rnorm(1, mun, sqrt(taun2))
pred_x_sim1<-rnorm(1, pred_theta_sim, sqrt(sigma2)); pred_x_sim1
S <- 100 # 100 predictions
pred_theta_sim <- rnorm(S, mun, sqrt(taun2))
pred_x_sim <- rnorm(S, pred_theta_sim, sqrt(sigma2))


### Slide 15 ###

library(rjags)
x=c(1.64,1.70,1.72,1.74,1.82,1.82,1.82,1.90,2.08); n<-length(x)
mu0<-1.9;tau0<-(.5*1.9);tau02<-tau0^2;xbar<-mean(x);sigma2<-var(x)

model_string <- "model{
 for (i in 1:n) {
  x[i] ~ dnorm(theta,1/sigma2) #likelihood
	}
  #xbar ~ dnorm(theta,n/sigma2) 
  theta ~ dnorm(mu0,1/tau02) #prior
  xtilde ~ dnorm(theta,1/sigma2) #prediction
}"


normnorm_model <- jags.model(textConnection(model_string),
 data = list(x=x,n=n,sigma2=sigma2,mu0=mu0,tau02=tau02))

update(normnorm_model, 5000); # Burnin for 5000 samples

postsamp <- coda.samples(normnorm_model,
   variable.names=c("theta","xtilde"),n.iter=10000)

summary(postsamp)
plot(postsamp)
