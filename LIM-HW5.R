### 1. (a)
beta.hat.lam=function(X, y, lam){
  n=nrow(X)
  p=ncol(X)
  xbar=apply(X[,-1],2,mean)
  ybar=mean(y)
  Xtilde=scale(X[,-1], center=xbar, scale=FALSE)
  ytilde=y-ybar
  q=min(c(n-1, p-1))
  out=svd(Xtilde, nu=q, nv=q)
  H=diag(out$d[1:q]/(out$d[1:q]^2+lam))
  bhatm1=out$v%*%H%*%t(out$u)%*%ytilde
  bhat1=ybar-sum(xbar*bhatm1)
  beta.hat=c(bhat1, bhatm1)
  return(beta.hat)
}

### 1. (b)
set.seed(5701)
n=50
p=5
sigma.star=1
beta.star=rnorm(p, mean=0, sd=sqrt(0.05))
### design matrix X
X=matrix(NA, nrow=n, ncol=p)
X[,1]=1
X[,-1]=rnorm(n*(p-1))
### y
y=X%*%beta.star+rnorm(n)*sigma.star
lambda=2
### beta.hat when lambda=2
beta.hat=beta.hat.lam(X,y,lambda)
### Getting gradient of objective function of beta.hat
firstterm=-2*crossprod(X, y-X%*%beta.hat)
secondterm=2*lambda*beta.hat
secondterm[1]=0
gradient=firstterm+secondterm
gradient
sum(gradient)

### 1. (c)
### a function that performs K-fold CV that finds best lambda minimizing prediction error
find.best.lam=function(X, y, lam.vec, K){
  n=nrow(X)
  p=ncol(X)
  cv.err=rep(0, length(lam.vec))
  idx=sample(n)
  for(k in 1:K){
    fold=which(idx==k)
    X.train=X[-fold,]
    X.test=X[fold,]
    y.train=y[-fold]
    y.test=y[fold]
    
    Xbar.tr=apply(X.train[,-1],2,mean)
    ybar.tr=mean(y.train)
    Xtilde.tr=scale(X.train[,-1], center=Xbar.tr, scale=F)
    ytilde.tr=y.train-ybar.tr
    q=min(c(nrow(X.train)-1, p-1))
    out=svd(Xtilde.tr, nu=q, nv=q)
    for(j in 1:length(lam.vec)){
      H=diag(out$d[1:q]/(out$d[1:q]^2+lam.vec[j]))
      bhatm1=out$v%*%H%*%t(out$u)%*%ytilde.tr
      bhat1=ybar.tr-sum(Xbar.tr*bhatm1)
      beta.hat.train=c(bhat1, bhatm1)
      err=sum((y.test-X.test%*%beta.hat.train))
      cv.err[j]=cv.err[j]+err
    }
  }
  result=list()
  best.lambda=lam.vec[which.min(cv.err)]
  result$b=beta.hat.lam(X, y, best.lambda)
  result$best.lam=best.lambda
  result$cv.error=cv.err
  return(result)
}



### 1. (d). 
library(MASS)
set.seed(5701)
### parameter setting
n=100; p=50; theta=0.95; lam.vec=10^seq(-8,8,by=0.5); sigma.star=1; K.vec=c(0,5,10); reps=200
### Set Design matrix and beta.star
Sigma=matrix(0, nrow=p-1, ncol=p-1)
for(i in 1:(p-1)){
  for(j in 1:(p-1)){
    Sigma[i,j]=theta^abs(i-j)
  }
}
X=cbind(1, mvrnorm(n, mu=rep(0, p-1), Sigma=Sigma))
mat=qr.solve(crossprod(X))%*%t(X)
beta.star=rnorm(p, mean=0, sd=sqrt(p^(-1)))

for(r in 1:reps){
  ### Generate y
  y=X%*%beta.star+rnorm(n)*sigma.star
  loss.mat.1=loss.mat.2= matrix(0, nrow=reps, ncol=length(K.vec))
  ### Compute the realizations of beta.hat(0),beta.hat(5),beta.hat(10)
  
  for(i in 1:length(K.vec)){
    k=K.vec[i]
    if(k==0){
      ### save beta.hat for each replication
      beta.hat=mat%*%y
    }
    else {
      result=find.best.lam(X, y, lam.vec, k)
      beta.hat=result$b
    }
    ### save Loss for three estimates' loss
    loss.mat.1[r,i]=sum((beta.hat-beta.star)^2)
    loss.mat.2[r,i]=sum((X%*%beta.hat-X%*%beta.star)^2)
  }
} 
cbind(apply(loss.mat.1,2,mean),apply(loss.mat.1,2,sd)/sqrt(reps))
cbind(apply(loss.mat.2,2,mean),apply(loss.mat.2,2,sd)/sqrt(reps))



### 1. (d). 6 Estimates
### Designate matrix that saves the best lambdas for each K
### size: 200 x 3
loss.mat=matrix(0, nrow=reps, ncol=6)
mat=qr.solve(crossprod(X))%*%t(X)
for(r in 1:reps){
  y=X%*%beta.star+rnorm(n)*sigma.star
  result=numeric(length(K.vec))
  for(i in 1:length(K.vec)){
    k=K.vec[i]
    if(k==0){
      beta.hat=mat%*%y
    }
    else {
      result=find.best.lam(X, y, lam.vec, k)
      beta.hat=result$b
    }
    loss.mat[r,i]=sum((beta.hat-beta.star)^2)
    loss.mat[r,i+3]=sum((X%*%beta.hat-X%*%beta.star)^2)
  }
}
loss.mat=apply(loss.mat, 2, mean)
estimates=sqrt(loss.mat)/reps
### The realizations of beta.hat(lam=0), beta.hat(lam_with 5-CV), beta.hat(lam_with 10-CV)
estimates


### 2. (b)
ridge.penal.lgst=function(X, y, n.list, lam, t.list, m=NULL, tol, maxit){
  p=ncol(X)
  if(is.null(m))
    m=c(0, rep(1, p-1))
  b=rep(0,p)
  lam.m=lam*m
  lam.M=diag(lam.m)
  X.t.n.list.y=crossprod(X, n.list*y)
  k=0
  add=tol+1
  while((k<=maxit)&(sum(abs(add))>tol)){
    k=k+1
    pi.t=ilogit(as.numeric(X%*%b))
    W=diag(n.list*pi.t*(1-pi.t))
    minusGrad=X.t.n.list.y-crossprod(X, n.list*pi.t)-lam.m*b+lam.m*t.list
    Hess=crossprod(X, W%*%X)+ lam.M
    add=qr.solve(Hess, minusGrad)
    b=b+add
    cat("k=",k,"b=",b,"\n")
  }
  b=as.numeric(b)
  return(list(b=b, total.iterations=k))
}


### 2. (c)
ilogit=function(u){
  return(exp(u)/(1+exp(u)))
}

set.seed(5701)
n=10; p=5
X=cbind(1, matrix(rnorm(n*(p-1)), nrow=n, ncol=(p-1)))
beta.star=rnorm(p, mean=0, sd=1/sqrt(p))
n.list=rep(30,n)
pi.list=ilogit(as.numeric(X%*%beta.star))
y=rbinom(n=n, size=n.list, prob=pi.list)/n.list
### Fit binomial logistic regression model
lam=10^4
t.list=c(0, beta.star[-1])

fit=ridge.penal.lgst(X, y, n.list, lam, t.list, tol=1e-7, maxit=100)
fit

### The gradient with final iterate
m=c(0, rep(1, p-1))
pi.hat=ilogit(as.numeric(X%*%fit$b))
Grad=crossprod(X, n.list*(pi.hat-y))+lam*fit$b*m-lam*t.list
Grad


### 3. (b)
ridge.penal.lgst=function(X, y, n.list, lam, m=NULL, tol=1e-7, maxit=100, b.start=NULL){
  p=ncol(X)
  if(is.null(m))
    m=c(0, rep(1, p-1))
  lam.m=lam*m
  lam.M=diag(lam.m)
  X.t.n.list.y=crossprod(X, n.list*y)
  if(is.null(b.start))
    b=rep(0,p)
  else
    b=b.start
  k=0
  add=tol+1
  while((k<=maxit)&(sum(abs(add))>tol)){
    k=k+1
    pi.t=ilogit(as.numeric(X%*%b))
    W=diag(n.list*pi.t*(1-pi.t))
    minusGrad=X.t.n.list.y-crossprod(X, n.list*pi.t)-lam.m*b
    Hess=crossprod(X, W%*%X)+ lam.M
    add=qr.solve(Hess, minusGrad)
    b=b+add
    #cat("k=",k,"b=",b,"\n")
  }
  b=as.numeric(b)
  return(list(b=b, total.iterations=k))
}

bridge.binl.lsg=function(X, y, n.list, lam, alpha, tol=1e-7, maxit=100){
  p=ncol(X)
  b=ridge.penal.lgst(X=X, y=y, n.list=n.list, lam=lam)$b
  k=0
  dummy=tol+1
  while((k<=maxit)&(dummy>tol)){
    k=k+1
    m=abs(b)^(alpha-2)
    m[1]=0
    newb=ridge.penal.lgst(X=X, y=y, n.list=n.list, m=m, lam=lam, b.start=b)$b
    dummy=sum(abs(newb-b))
    b=newb
    cat("k=",k,"b=",b,"\n")
  }
  b=as.numeric(b)
  return(list(b=b, total.iterations=k))
}


### 3. (c)
set.seed(5701)
n=10; p=5
X=cbind(1, matrix(rnorm(n*(p-1)), nrow=n, ncol=(p-1)))
beta.star=rnorm(p, mean=0, sd=1/sqrt(p))
n.list=rep(30,n)
pi.list=ilogit(as.numeric(X%*%beta.star))
y=rbinom(n=n, size=n.list, prob=pi.list)/n.list
lam=5
alpha=1.1

fit=bridge.binl.lsg(X=X, y=y, n.list=n.list, lam=lam, alpha=alpha, tol=1e-7, maxit=100)
m=c(0, rep(1,p-1))
pi.hat=ilogit(as.numeric(X%*%fit$b))
Grad=crossprod(X, n.list*(pi.hat-y))+lam*abs(fit$b)^(alpha-1)*sign(fit$b)*m
Grad
