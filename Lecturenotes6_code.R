### Slide 2 ###

x <- seq(0,5,by=0.00001)
plot(x,dgamma(x,1,rate=1),xlim=c(0,5),ylim=c(0,1.5),lwd=1.5,
     type="l",col="green3",xlab=expression(theta),ylab="density")
lines(x,dgamma(x,2,rate=2),col="red2",lwd=1.5,type="l")
lines(x,dgamma(x,5,rate=2),col="blue2",lwd=1.5,type="l")
lines(x,dgamma(x,0.001,rate=0.001),col="orange3",lwd=1.5,type="l")
legend("topright", legend=c("gamma(1,1)","gamma(2,2)",
                            "gamma(5,2)","gamma(0.001,0.001)"),
       col=c("green3","red2","blue2","orange3"), lwd=2, cex=1)

### Slide 7 ###

a<-2 ; b<-1; # prior parameters
n1<-111 ; sx1<-217 # group 1
n2<-44 ; sx2<-66 # group 2
(a+sx1)/(b+n1); (a+sx2)/(b+n2) # posterior mean
(a+sx1-1)/(b+n1); (a+sx2-1)/(b+n2) # posterior mode
qgamma(c(.025, .975), a+sx1, b+n1) # posterior 95% CI
qgamma(c(.025, .975), a+sx2, b+n2) # posterior 95% CI
x <- seq(0,6,by=0.00001)
plot(x,dgamma(x,shape=2,rate=1),xlim=c(0,6),ylim=c(0,3),lwd=1.5,
     type="l",col="green3",xlab=expression(theta),ylab="density")
lines(x,dgamma(x,shape=219,rate=112),col="red2",lwd=1.5,type="l")
lines(x,dgamma(x,shape=68,rate=45),col="blue2",lwd=1.5,type="l")
legend("topright",legend=c("Prior","Post.(No Bach)","Post.(Bach)"),
       col=c("green3","red2","blue2"), lwd=2, cex=1)


### Slide 8 ###

th1_mc<-rgamma(10^5,a+sx1,b+n1); th2_mc<-rgamma(10^5,a+sx2,b+n2)
mean(th1_mc>th2_mc)

nbach=dnbinom(0:10,size=(a+sx1),mu=(a+sx1)/(b+n1))
bach=dnbinom(0:10,size=(a+sx2),mu=(a+sx2)/(b+n2))
barplot(rbind(nbach,bach),xlab="Children",
        names.arg=c("0","1","2","3","4","5","6","7","8","9","10"),
        col=c("orange3","green3"),beside=TRUE,
        main="Posterior predictive distributions")
legend("topright",c("No Bachelor","Bachelor"),
       bty="n",fill=c("orange3","green3"))
x1tilde_mc<-rpois(10^5,th1_mc); x2tilde_mc<-rpois(10^5,th2_mc)
mean(x1tilde_mc>x2tilde_mc)

### Slide 10 ###

require(R2jags)
Poigam_jags = "model{
sx1 ~ dpois(theta1); theta1<-n1*th1 #Group 1
sx2 ~ dpois(theta2); theta2<-n2*th2 #Group 2
xtilde1 ~ dpois(th1) # Prediction of future Poisson from Group 1
xtilde2 ~ dpois(th2) # Prediction of future Poisson from Group 2
th1 ~ dgamma(a, b); th2 ~ dgamma(a, b) #gamma prior
postprob <- step(th1-th2) # Post prob that th1 > th2
postpredprob <- step(xtilde1-xtilde2)-equals(xtilde1,xtilde2)
}"
datalist= list(sx1=217, sx2=66, n1=111, n2=44, a = 2, b = 1)
initlist= function() {list(th1 = 1, th2=1, xtilde1=10, xtilde2=5)}
params = c("th1", "th2", "postprob", "postpredprob")
jags_out = jags(model=textConnection(Poigam_jags), data =datalist,
                inits=initlist, n.chains = 3, n.iter = 15000, parameters = params)
jags_out; traceplot(jags_out)



