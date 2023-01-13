### Slide 11 ###

library(tidyverse)
betapars<-matrix(c(0.5,0.5,5,1,1,3,2,2,2,5,1,1),6,2, byrow=TRUE)
theta <- seq(.001, .999, length.out = 100); BETA <- NULL

for (j in 1:6){
df <- data.frame(theta = theta, 
pdf=dbeta(theta, betapars[j,1],betapars[j,2]))
df$Type<-paste("Beta(",betapars[j,1],",",betapars[j,2],")",sep ="") 
BETA <- rbind(BETA, df)
}

ggplot(BETA,aes(theta,pdf))+geom_line(aes(color=Type,linetype=Type)
,size = 1)+facet_wrap(~Type,scale = "free")+
  theme(axis.text.y=element_blank())+
  scale_x_continuous(breaks=seq(0, 1, 0.5))+
  scale_color_manual(values=c("black","olivedrab","red",
                              "hotpink4","lightpink","azure4"))+
  scale_linetype_manual(values=c("solid","longdash","dotdash",
                                 "solid","dashed","longdash"))


### Slide 14 ###

n <- 30; x <- 6; a <- 3.1; b <- 1.5
grid <- seq(0,1,.01)
like <- dbinom(x,n,grid); like <- like/sum(like) #scaled
prior <- dbeta(grid,a,b); prior <- prior/sum(prior)
post <- like*prior; post <- post/sum(post)
#post <- dbeta(grid,x+a,n-x+b); post <- post/sum(post)

plot(grid,like,type="l",lty=2,col=1,xlab="theta",ylab="Density")
lines(grid,prior,col="blue")
lines(grid,post,lwd=2,col="red")
legend("topright",c("Likelihood","Prior","Posterior"),
       lwd=c(1,1,2),lty=c(2,1,1),col=c(1,"blue","red"))

### Slide 15 ###

n.samples <- 10000; postsamp <- rbeta(n.samples,x+a,n-x+b)
# postsamp <- replicate(n.samples, beta.rej1(x+a,n-x+b))

hist(postsamp,breaks=25,xlim=0:1,main="Posterior density",freq=F)
lines(density(postsamp),lty=2,col="red")
lines(grid,dbeta(grid,x+a,n-x+b), lty=1, col="green")
legend("topright",c("Density estimate","Exact density"),
       lwd=c(1,1),lty=c(2,1),col=c("red","green"))
mean(postsamp); quantile(postsamp,c(0.05,0.95))
(x+a)/(n+a+b); qbeta(c(0.05,0.95),x+a,n-x+b)

### Slide 17 ###

theta <- seq(0,1,length=1000)
a<-3.06 ; b<-2.56; n<-20; x<-12 # n<-200; x<-120
plot(theta,dbeta(theta,a+x,b+n-x),type="l", col="red",ylab=
       expression(paste(italic("p("),theta,"|x)",sep="")),
     xlab=expression(theta), ylim=c(0,12), lwd=2)
mtext(expression(paste("Beta(3.06,2.56) prior,",italic("n"),"=20,",
                       italic(sum(x[i])),"=12",sep="")),side=3,line=.1)
abline(v=c((a+x-1)/(a+b+n-2),x/n,a/(a+b)),
       col=c("green","red","blue"),lwd=c(2,2))
lines(theta,dbeta(theta,a,b),type="l",col="blue",lwd=2)
legend("topleft",c("Prior","Posterior","MLE"),
       lwd=c(2,2,2),lty=c(1,1,1),col=c("blue","red","green"))
theta[which.max(dbeta(theta,a,b))] #numerical prior mode
theta[which.max(dbeta(theta,a+x,b+n-x))] #numerical posterior mode

set.seed(13); postsamp <- replicate(10^4, beta.rej1(x+a,n-x+b))
mode_estimate <- function(x) {
  d.val <- density(x); d.val$x[which.max(d.val$y)]
}
mode_estimate(postsamp)


