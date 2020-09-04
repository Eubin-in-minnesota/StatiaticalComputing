### 1.
library(MASS)
# function that makes covariance matrix 
covar.mat=function(n, sigma, rho){
  Sigma=matrix(NA, nrow=n, ncol=n)
  for(i in 1:n){
    for(j in 1:n){
      if(i==j){
        Sigma[i,j]=sigma^2
      }
      else
        Sigma[i,j]=sigma^2*rho
    }
  }
  return(Sigma)
}

set.seed(5701)
reps=1e4
sigma.star=1
n.vec=c(10,50,200,500,1000)
alpha.vec=c(0.01,0.05,0.1)
beta.star=c(1,1,2,1,-2)
p=5
xnew=c(1,0,1,0,1)
true.mean=sum(beta.star*xnew)
captured.list=numeric(reps)
# result matrix
res.mat=matrix(0,length(n.vec)*length(alpha.vec),5)
k=1
for(n in n.vec){
  for(alpha in alpha.vec){
    for(r in 1:reps){
      res.mat[k,1]=n
      res.mat[k,2]=alpha
      Sigma=covar.mat(n=p-1, sigma=sigma.star, rho=0.6)
      # create the design matrix X
      Xtilde=mvrnorm(n=n, mu=rep(0, p-1), Sigma = Sigma)
      X.mat=cbind(1, Xtilde)
      # create epsilon vector
      epsilon.new=rexp(1,1)-1
      epsilon.vec=rexp(n,1)-1
      
      y.vec=X.mat%*%beta.star+epsilon.vec
      
      # compute the relization of betahat
      qrX.mat=qr(x=X.mat)
      beta.hat=qr.coef(qr=qrX.mat, y=y.vec)
      sE=sqrt(sum((y.vec-X.mat%*%beta.hat)^2)/(n-p))
      
      # compute useful quantities
      sqrtquadform=sqrt(1+t(xnew)%*%qr.solve(crossprod(X.mat))%*%xnew)[1]
      t.perc=qt(1-alpha/2, n-p)
      
      # compute the 1-alpha % prediction interval for Y_new
      est.mean=sum(beta.hat*xnew)
      moe=t.perc*sE*sqrtquadform
      left.pt=est.mean-moe
      right.pt=est.mean+moe
      # generate the realization of Y_new
      ynew=sum(beta.star*xnew)+epsilon.new
      # check if ynew was captured
      captured.list[r]=1*(left.pt<=ynew)&(ynew<=right.pt)
    }
    res.mat[k,3:4]=prop.test(x=sum(captured.list), n=reps, conf.level=0.99, correct=FALSE)$conf.int[1:2]
    # indication to check whether simulation-based 99% score approximate CI covered 1-alpha
    if(res.mat[k,3]<1-alpha&&res.mat[4]>1-alpha)
      res.mat[k,5]=1 # captured=1
    else
      res.mat[k,5]=0 # failed to capture=0
    k=k+1
  }
}
colnames(res.mat)=c("n","alpha","lower bound","upper bound","captured 1-alpha?")
res.mat

### 2.
# arguments
# X: nxp design matrix
# y: n size response vector
# K: number of folds >= 2

olscv=function(X,y,K){
  n=length(y)
  # make random indice from 1 to n
  idx=sample(1:n)
  # shuffle the data(design matrix)
  X=X[idx,]
  
  # Sum of Square Error
  SE=0
  
  folds=rep_len(1:K, nrow(X))
  for(k in 1:K){
    fold=which(folds==k)
    train.X=X[-fold,]
    test.X=X[fold,]
    train.y=y[-fold]
    test.y=y[fold]
    
    # train model
    m.lm=lm.fit(x=train.X, y=train.y)
    beta=m.lm$coefficients
    
    # y prediction
    y.predict=test.X%*%beta
    
    # update 
    SE=SE+sum((test.y-y.predict)^2)
  }
  # return MSE
  return(SE/n)
}

### 3. (a). i)
library(MASS)
# H0: beta_4 = 0, H1: otherwise
power.calculate=function(n, beta1=1, beta2, beta3, beta4, rho, alpha, reps) {
  beta.star=c(beta1=1, beta2, beta3, beta4)
  sigma.star=1
  p=4
  d=1
  
  # create Sigma 
  Sigma=covar.mat(n=p-1, sigma=sigma.star, rho=rho)
  # create the design matrix X
  Xtilde=mvrnorm(n=n, mu=rep(0, p-1), Sigma = Sigma)
  X.mat=cbind(1, Xtilde)
  
  # under H_0
  qrXnull=qr(X.mat[, 1:3])
  # under H_1
  qrX=qr(X.mat)
  
  # to save p-values
  pval.list=rep(0, reps)
  for(r in 1:reps){
    # generate y
    y=X.mat%*%beta.star+sigma.star*rnorm(n)
    
    # rss0(under H_0)
    beta.null=qr.coef(qrXnull, y)
    rss0=sum((y-X.mat[, 1:3]%*%beta.null)^2)
    
    # rss1(under H_1)
    beta=qr.coef(qrX, y)
    rss1=sum((y-X.mat%*%beta)^2)
    
    # calc F statistic and store current p-value
    f = ((rss0 - rss1) / d) / (rss1 / (n - p))
    pval.list[r] = 1 - pf(f, d, n - p)
  }
  return(mean(pval.list < alpha))
}


# parameter settings
alpha=0.01
reps=1e3
beta1=beta2=beta3=1
beta4.vec=c(0.5, 1)
rho.vec=c(0.3, 0.9)
n.vec=seq(from=10, by=5, length.out=100)

res.mat = matrix(0, nrow=length(beta4.vec)*length(rho.vec), ncol=4)
colnames(res.mat) = c("beta4", "rho", "n", "power")

set.seed(5701)
k=1
for(beta4 in beta4.vec){
  for(rho in rho.vec){
    dum=2
    final.n=0
    final.power=0
    for(n in n.vec){
      power=power.calculate(n, beta1, beta2, beta3, beta4, rho, alpha, reps)
      if(abs(power-0.8)<dum){
        final.n=n
        final.power=power
        dum=abs(power-0.8)
      }
    }
    res.mat[k,] = c(beta4, rho, final.n, final.power)
    k=k+1
  }
}
res.mat

### 3. (a). ii) (1)
set.seed(5701)
get.ic=function(X, y){
  n=dim(X)[1]
  p=dim(X)[2]
  beta.hat=lm.fit(x=X, y=y)$coefficients
  rss=sum((y-X%*%beta.hat)^2)
  common=n*log(2*pi)+n*log(rss/n) + n
  aic=common + 2*(p+1)
  bic=common + (p+1)*log(n)
  return(c(aic, bic))
}

# arguments
# which: 1-aic 2-bic
prob.fullmodel.IC.is.smaller=function(n, beta4, rho, reps, which){
  beta.star=c(1, 1, 1, beta4)
  sigma.star=1
  p=4
  
  # create Sigma 
  Sigma=covar.mat(n=p-1, sigma=sigma.star, rho=rho)
  # create the design matrix X
  Xtilde=mvrnorm(n=n, mu=rep(0, p-1), Sigma = Sigma)
  X.mat=cbind(1, Xtilde)
  
  # indicator=1 if Full model's IC < null model's IC
  count=numeric(reps)
  for(r in 1:reps){
    # generate y
    y=X.mat%*%beta.star+sigma.star*rnorm(n)
    
    # IC under H_0
    ic0=get.ic(X=X.mat[,1:3],y)[which]
    
    # IC under H_1
    ic1=get.ic(X=X.mat,y)[which]
    
    if(ic0>ic1){
      count[r]=1
    }
  }
  # return the probability that Full model's IC < null model's IC
  return(sum(count)/reps)
}

beta4.vec=c(0.5, 1)
rho.vec=c(0.3, 0.9)
# parameter set from i)
param.mat1 = matrix(0, nrow=4, ncol=4)
colnames(param.mat1) = c("beta.4", "rho", "chosen n", "Prob. that full model's AIC is smaller")
param.mat1[1,] = c(0.5, 0.3, 60, 0)
param.mat1[2,] = c(0.5, 0.9, 350, 0)
param.mat1[3,] = c(1, 0.3, 20, 0)
param.mat1[4,] = c(1, 0.9, 110, 0)

for(i in 1:4){
  beta4=param.mat1[i,1]
  rho=param.mat1[i,2]
  n=param.mat1[i,3]
  beta.star=c(1,1,1,beta4)
  sigma.star=1
  reps=1e4
  param.mat1[i,4]=prob.fullmodel.IC.is.smaller(n=n, beta4=beta4, rho=rho, reps=reps, which=1)
}
param.mat1

### 3. (a). ii) (2)
param.mat2 = matrix(0, nrow=4, ncol=4)
colnames(param.mat2) = c("beta.4", "rho", "chosen n", "Prob. that full model's BIC is smaller")
param.mat2[1,] = c(0.5, 0.3, 60, 0)
param.mat2[2,] = c(0.5, 0.9, 350, 0)
param.mat2[3,] = c(1, 0.3, 20, 0)
param.mat2[4,] = c(1, 0.9, 110, 0)

for(i in 1:4){
  beta4=param.mat2[i,1]
  rho=param.mat2[i,2]
  n=param.mat2[i,3]
  beta.star=c(1,1,1,beta4)
  sigma.star=1
  reps=1e4
  param.mat2[i,4]=prob.fullmodel.IC.is.smaller(n=n, beta4=beta4, rho=rho, reps=reps, which=2)
}
param.mat2

### 3. (a). ii) (3)
set.seed(5701)
param.mat3=matrix(0,nrow=4,ncol=4)
colnames(param.mat3) = c("beta.4", "rho", "chosen n", "Prob. that CV-based MSE of Full model is smaller")
param.mat3[1,] = c(0.5, 0.3, 60, 0)
param.mat3[2,] = c(0.5, 0.9, 350, 0)
param.mat3[3,] = c(1, 0.3, 20, 0)
param.mat3[4,] = c(1, 0.9, 110, 0)

for(i in 1:4){
  beta4=param.mat3[i,1]
  rho=param.mat3[i,2]
  n=param.mat3[i,3]
  beta.star=c(1,1,1,beta4)
  sigma.star=1
  reps=1e4
  K=5
  p=4
  # create Sigma 
  Sigma=covar.mat(n=p-1, sigma=sigma.star, rho=rho)
  # create the design matrix X
  Xtilde=mvrnorm(n=n, mu=rep(0, p-1), Sigma = Sigma)
  X.mat=cbind(1, Xtilde)
  # create y
  y=X.mat%*%beta.star+sigma.star*rnorm(n)
  count=0
  for(r in 1:reps){
    if(olscv(X.mat,y,K)<olscv(X.mat[,1:3],y,K)){
      count=count+1
    }
    param.mat3[i,4]=count/reps
  }
  
}
param.mat3

### 3. (a) ii. Rank the methods performance and compare them to the power of LRT/F test at 1% significance level.
final=cbind(res.mat, param.mat1[,4], param.mat2[,4], param.mat3[,4])
colnames(final)[5:7]=c("P(AIC makes correct choice)","P(BIC makes correct choice)","P(5-fold CV makes correct choice)")
final

### 3. (a) iii. (1)
set.seed(5701)
n=100
beta4=0
rho=0.3
library(MASS)
# H0: beta_4 = 0, H1: otherwise
typeIerror.calculate=function(n, beta4, rho, alpha, reps){
  beta.star=c(1, 1, 1, beta4)
  sigma.star=1
  p=4
  d=1
  
  # create Sigma 
  Sigma=covar.mat(n=p-1, sigma=sigma.star, rho=rho)
  # create the design matrix X
  Xtilde=mvrnorm(n=n, mu=rep(0, p-1), Sigma = Sigma)
  X.mat=cbind(1, Xtilde)
  
  # under H_0
  qrXnull=qr(X.mat[, 1:3])
  # under H_1
  qrX=qr(X.mat)
  
  # to save p-values
  pval.list=numeric(reps)
  for(r in 1:reps){
    # generate y
    y=X.mat%*%beta.star+sigma.star*rnorm(n)
    
    # rss0(under H_0)
    beta.null=qr.coef(qrXnull, y)
    rss0=sum((y-X.mat[, 1:3]%*%beta.null)^2)
    
    # rss1(under H_1)
    beta=qr.coef(qrX, y)
    rss1=sum((y-X.mat%*%beta)^2)
    
    # calc F statistic and store current p-value
    f = ((rss0 - rss1) / d) / (rss1 / (n - p))
    pval.list[r] = 1 - pf(f, d, n - p)
  }
  return(mean(pval.list < alpha))
} 
p1=typeIerror.calculate(n=n, beta4=0, rho=0.3, alpha=0.01, reps=1e4)
p1

### 3. (a) iii. (2)
set.seed(5701)
reps=1e4
n=100
beta4=0
rho=0.3
p2=prob.fullmodel.IC.is.smaller(n=n, beta4=beta4, rho=rho, reps=reps, which=1)
p2

### 3. (a) iii. (3)
set.seed(5701)
p3=prob.fullmodel.IC.is.smaller(n=n, beta4=beta4, rho=rho, reps=reps, which=2)
p3

### 3. (a) iii. (4)
set.seed(5701)
beta.star=c(1,1,1,beta4)
sigma.star=1
K=5
p=4
# create Sigma 
Sigma=covar.mat(n=p-1, sigma=sigma.star, rho=rho)
# create the design matrix X
Xtilde=mvrnorm(n=n, mu=rep(0, p-1), Sigma = Sigma)
X.mat=cbind(1, Xtilde)
# create y
y=X.mat%*%beta.star+sigma.star*rnorm(n)
count=0
for(r in 1:reps){
  if(olscv(X.mat,y,K)<olscv(X.mat[,1:3],y,K)){
    count=count+1
  }
  result=count/reps
}
p4=result
final=cbind(beta4,rho,n,p1,p2,p3,p4)
colnames(final)=c("beta4","rho","n","Type I Error","P(AIC makes incorrect choice)","P(BIC makes incorrect choice)","P(5-fold CV makes incorrect choice)")
final

### 3. (b)
# H0: beta_2=beta_3=beta_4, H1: otherwise
power.calculate=function(n, beta1, beta2, beta3, beta4, rho, alpha, reps) {
  beta.star=c(beta1=1, beta2, beta3, beta4)
  sigma.star=1
  p=4
  d=2
  
  # create Sigma 
  Sigma=covar.mat(n=p-1, sigma=sigma.star, rho=rho)
  # create the design matrix X
  Xtilde=mvrnorm(n=n, mu=rep(0, p-1), Sigma = Sigma)
  X=cbind(1, Xtilde)
  
  C=rbind(c(0,1,-1,0),c(0,0,-1,1))
  C=matrix(C, nrow=2,ncol=4)
  
  # mat=C'{C(X'X)^-1C'}^-1C
  mat=t(C)%*%solve(C%*%solve(t(X)%*%X)%*%t(C))%*%C
  qrX=qr(x=X)
  
  # to save p-values
  pval.list=numeric(reps)
  for(r in 1:reps){
    # generate y
    y=X%*%beta.star+sigma.star*rnorm(n)
    
    beta.hat=qr.coef(qrX,y=y)
    
    # rss1(under H_1)
    rss1=sum((y-X%*%beta.hat)^2)
    
    # calc F statistic and store current p-value
    f = as.numeric(t(beta.hat)%*%mat%*%beta.hat)*(n-p)/(d*rss1)
    pval.list[r]=1-pf(f, d, n-p)
  }
  return(mean(pval.list < alpha))
}


# parameter settings
alpha=0.01
reps=1e3
beta1=beta2=beta3=1
beta4.vec=c(1.5, 2.5)
rho.vec=c(0.3, 0.9)
n.vec=seq(from=10, by=10, length.out=100)

res.mat = matrix(0, nrow=length(beta4.vec)*length(rho.vec), ncol=4)
colnames(res.mat) = c("beta4", "rho", "n", "power")

set.seed(5701)
k=1
for(beta4 in beta4.vec){
  for(rho in rho.vec){
    dum=2
    final.n=0
    final.power=0
    for(n in n.vec){
      power=power.calculate(n, beta1, beta2, beta3, beta4, rho, alpha, reps)
      if(abs(power-0.85)<dum){
        final.n=n
        final.power=power
        dum=abs(power-0.85)
      }
    }
    # result matrix
    res.mat[k,] = c(beta4, rho, final.n, final.power)
    k=k+1
  }
}
res.mat