### Slide 3 ###

jointpmf<-matrix(c(3,3,2,1,3,3,3,2,2,3,3,3,1,2,3,3)/40,4,4,byrow=T)
dimnames(jointpmf)[[1]]<-1:4;dimnames(jointpmf)[[2]]<-1:4;jointpmf
rng.comp <- function(jointp, nsim){
  jointsamp<-matrix(0,nsim,2);nY<-dim(jointp)[1];nX<-dim(jointp)[2]
  for (k in 1:nsim) {
    ypmf <- apply(jointp,1,sum)
    y <- sample(1:nY,1,prob=ypmf,replace=TRUE)
    x <- sample(1:nX,1,prob=jointp[y,]/ypmf[y],replace=TRUE)
    jointsamp[k,] <- c(x,y)}
  return(jointsamp)}
set.seed(1);nsim<-10^4;xy.rng<-data.frame(rng.comp(jointpmf,nsim))
names(xy.rng) <- c("Y", "X"); table(xy.rng)/nsim

### Slide 4 ###

rng.comp1 <- function(ntrial){
  y <- runif(1,0,1)
  x <- rbinom(1,ntrial,y)
  return(x)
}
set.seed(3); nsim<-10^4
x.rng<- replicate(nsim,rng.comp1(10))
table(x.rng)/nsim

### Slide 7 ###

trans.mat <- function(A) {
  n <- length(A[,1]); temp <- A
  for (i in 1:n){ temp[i,] <- A[i,]/sum(A[i,])}
  return(temp)}
gibbs_discrete <- function(pmat1, pmat2, i = 1, iter){
  jointsamp<-matrix(0,iter,2);nX<-dim(pmat1)[1];nY<-dim(pmat2)[1]
  for(k in 1:iter){
    j <- sample(1:nX, 1, prob = pmat1[i, ])
    i <- sample(1:nY, 1, prob = pmat2[j, ])
    jointsamp[k, ] <- c(i, j)}
  return(jointsamp)}
powermat <- function(mat,k) {
  if (k == 0) return (diag(dim(mat)[1]))
  if (k == 1) return(mat)
  if (k > 1) return( mat %*% powermat(mat, k-1))
}
Pmat.YX <- trans.mat(jointpmf); Pmat.XY <- trans.mat(t(jointpmf))
nburnin=2000; nSample = 10^4; niter=nburnin+nSample; set.seed(11)
gibbs_samp <- data.frame(gibbs_discrete(Pmat.YX,Pmat.XY,1,niter))
names(gibbs_samp)<-c("Y","X")
table(gibbs_samp[(nburnin+1):(niter),])/nSample
Pmat.X <- Pmat.YX %*% Pmat.XY
powermat(Pmat.X,10); powermat(Pmat.X,100); apply(jointpmf,1,sum)

### Slide 8 ###

gibbs_unifbin <- function(ntrial, y = 0.5, iter){
  jointsamp <- matrix(0, iter, 2)
  for(k in 1:iter){
    x <- rbinom(1, size = ntrial, prob = y)
    y <- rbeta(1, x+1, ntrial - x+1)
    jointsamp[k, ] <- c(x, y)
  }
  return(jointsamp)
}
nburnin=2000; nSample = 10^4; niter=nburnin+nSample; set.seed(11)
gibbs_samp <- data.frame(gibbs_unifbin(10,0.5,niter))
names(gibbs_samp)<-c("Y","X")
table(gibbs_samp[(nburnin+1):(niter),1])/nSample
hist(gibbs_samp[(nburnin+1):(niter),2],freq=F,main=" ")

### Slide 9 ###

y<-seq(0.001,0.999,0.001); fy <-dbeta(y,3,5); max(fy)
(gamma(8)/(gamma(3)*gamma(5)))*(1/3)^2*(2/3)^4

### Slide 10 ###

beta.rej1<-function(a,b) {
  M1<-(a-1)^(a-1)*(b-1)^(b-1)*(a+b-2)^(2-a-b); success<-FALSE
  while (!success) {
    W<-runif(1); V<-runif(1); accept_ratio<-W^(a-1)*(1-W)^(b-1)/M1
    success<-(V<accept_ratio)}
  Y<-W; return(W)}
a<-3; b<-5; nreps<-10^4
beta.rand <- replicate(nreps, beta.rej1(a,b))