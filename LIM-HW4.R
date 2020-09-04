### 1. (a)
load("arsenic.rdata")
X=arsenic[,-1]
y=log(arsenic[,1])
X$gender=as.numeric(X$gender)-1 # Female:1, Male:0
X$age=as.numeric(X$age)
X=cbind(1, X)
colnames(X)=c("intercept","arsenic.water","age","gender")

beta.hat=qr.coef(qr(x=X), y=y)
beta.hat

n=nrow(X)
p=ncol(X)
#beta.hat=as.matrix(beta.hat, nrow=3, ncol=1)
X=as.matrix(X)
yhat=X%*%beta.hat
sigma.sq=sum((y-yhat)^2)/(n-p)
sigma.sq

#The error variance is 0.2459478.

### 1. (c)
xnew=c(1,0,30,0)
alpha=0.01
est.mean=xnew%*%beta.hat
tperc=qt(1-alpha/2, n-p)
sigma=sqrt(sigma.sq)
sqrtquadform=sqrt(1+t(xnew)%*%solve(t(X)%*%X)%*%xnew)
moe=tperc*sigma*sqrtquadform
exp(c(est.mean-moe, est.mean+moe))

#99% prediction interval for yet-to-be observed value of arsenic.toenail for 30-year-old male subject with 
#0 arsenic.water is (0.04239963, 0.95417323). \ 

### 1. (d)
set.seed(5701)
n=nrow(X)
p=ncol(X)
reps=1e4
xnew=c(1,0,30,0)
qrX.mat=qr(x=X)
sqrtquadform=sqrt(1+t(xnew)%*%qr.solve(crossprod(X))%*%xnew)[1]
t.perc=qt(0.975, n-p)
captured.list=numeric(reps)
for(r in 1:reps){
  Z=sigma*sqrt(12)*(runif(n)-0.5)
  y.vec=X%*%beta.hat+Z
  beta.hat.real=qr.coef(qr=qrX.mat, y=y.vec)
  sE=sqrt(sum((y.vec-X%*%beta.hat.real)^2)/(n-p))
  est.mean=sum(beta.hat.real*xnew)
  moe=t.perc*sE*sqrtquadform
  left.pt=exp(est.mean-moe)
  right.pt=exp(est.mean+moe)
  ynew=xnew %*% beta.hat+sE*sqrt(12)*(runif(1)-0.5)
  ynew=exp(ynew)
  captured.list[r]=1*(left.pt<=ynew)&(ynew<=right.pt)
}
prop.test(x=sum(captured.list), n=reps, conf.level = 0.99, correct=FALSE)$conf.int[1:2]

# Since the confidence level is not laid between simulation-based confidence interval, 
#we can say that the coverage probability of this random prediction interval is not equal to 95%. 

### 2. (a)
n=10
mu=5
x.list=rnorm(n,mean=mu, sd=sqrt(mu))
mu.list=seq(1,10,0.001)
neg.loglk=function(mu.list, x.list){
  n=length(x.list)
  n.mu=length(mu.list)
  neg.list=numeric(n.mu)
  for(i in 1:n.mu){
    neg.list[i]=(n/2)*log(2*pi)+(n/2)*log(mu.list[i])+sum((x.list-mu.list[i])^2)/(2*mu.list[i])
  }
  return(neg.list)
}
plot(mu.list, neg.loglk(mu.list, x.list), type="l", xlab="mu.list",ylab="-loglikelihood")


### 2. (b)
## Dichotomous search for finding MLE of mu
## Minimize a univariate strictly quasiconvex function
## over the interval [a0,b0]
##
##  Arguments 
##    g, the function to minimize, where g(u, ...) is
##       the function evaluated at u.
##    a0, left endpoint of the initial interval of uncertainty.
##    b0, right endpoint of the initial interval of uncertainty.
##    L, the maximum length of the final interval of uncertainty.
##    eps, search parameter, must be less than L/2.
##    quiet, should the function stay powuiet?
##    ..., additional argument specifications for g
##  
##  Returns the midpoint of the final interval of uncertainty
dsearch=function(g, a0, b0,  L=1e-7, eps=(L/4), quiet=FALSE, ...)
{
  mm=mean(c(a0,b0))
  while(b0-a0 > L){
    tL = mm-eps
    tR = mm+eps
    g.at.tL=g(tL, ...)
    g.at.tR= g(tR, ...)
    if(g.at.tL < g.at.tR){
      b0=tR
    } 
    else if(g.at.tL > g.at.tR){
      a0=tL
    } 
    else{
      b0=tR
      a0=tL
    }
    if(!quiet) cat("new interval is", a0, b0, "\n")
    mm=mean(c(a0,b0))
  }
  return(mm)
}


set.seed(5701)
n=10
mu=5
x.list=rnorm(n,mean=mu, sd=sqrt(mu))
mu.list=seq(1,10,0.001)
## Objective function
neg.loglk=function(mu.list, x.list){
  n=length(x.list)
  n.mu=length(mu.list)
  neg.list=numeric(n.mu)
  for(i in 1:n.mu){
    neg.list[i]=(n/2)*log(2*pi)+(n/2)*log(mu.list[i])+sum((x.list-mu.list[i])^2)/(2*mu.list[i])
  }
  return(neg.list)
}

## Minimize g with Dichotomous search
mu.MLE=dsearch(g=neg.loglk, a0=min(x.list), b0=max(x.list), 
               quiet=TRUE, x.list=x.list)
xbar=mean(x.list)
ss=var(x.list)
mu.MLE;xbar;ss

#Maximum likelihood estimate of $\mu_*$ is 4.816727. \  
#The realization of $\bar X$ is 4.913975. \ 
#The realization of $S^2$ is 4.300483. \  



### 2. (c)
set.seed(5701)
n=10
mu.list=c(5,10,50)
reps=1e4
for(i in 1:length(mu.list)){
  mu=mu.list[i]
  xbar.diff=numeric(reps)
  ss.diff=numeric(reps)
  mu.MLE.diff=numeric(reps)
  xbar.mle.diff=numeric(reps)
  for(r in 1:reps){
    x.list=rnorm(n=n, mean=mu, sd=sqrt(mu))
    xbar=mean(x.list)
    ss=var(x.list)
    mu.MLE=dsearch(g=neg.loglk, a0=min(x.list), b0=max(x.list), quiet=TRUE, x.list=x.list)
    xbar.diff[r]=abs(xbar-mu)
    ss.diff[r]=abs(ss-mu)
    mu.MLE.diff[r]=abs(mu.MLE-mu)
    xbar.mle.diff[r]=xbar.diff[r]-mu.MLE.diff[r]
  }
  cat("mu=",mu,", 99% approx. simulation-based CI for","\n")
  CI.xbar.diff=t.test(xbar.diff, conf.level = 0.99)$conf.int[1:2] 
  cat("E(|X_bar - mu|) is", CI.xbar.diff,"\n")
  CI.ss.diff=t.test(ss.diff, conf.level = 0.99)$conf.int[1:2] 
  cat("E(|Sample Variance - mu|) is", CI.ss.diff,"\n")
  CI.mu.MLE.diff=t.test(mu.MLE.diff, conf.level = 0.99)$conf.int[1:2] 
  cat("E(|MLE - mu|) is", CI.mu.MLE.diff,"\n")
  CI.xbar.mle.diff=t.test(xbar.mle.diff, conf.level = 0.99)$conf.int[1:2] 
  cat("E(|X_bar - mu|-|MLE - mu|) is", CI.xbar.mle.diff,"\n")
}
# According to the result, MLE is the best since it has lowest mean absolute error and narrowest confidence interval between three estimates.
