# 1. (b)

# A simulation to compute a simulated estimate of the coverage probability of a 100(1-alpha)% random confidence interval for the variance, based on an iid sample from N(mu, sigma^2)

# function that generates CI for sigma^2 with alpha/2
make.CI.for.sigma.1=function(x, alpha){
  n=length(x)
  variance=var(x)
  left.pt=(n-1)*variance/qchisq(1-alpha/2, n-1)
  right.pt=(n-1)*variance/qchisq(alpha/2, n-1)
  return(c(left.pt,right.pt))
}

# function that generates CI for sigma^2 with alpha/3
make.CI.for.sigma.2=function(x, alpha){
  n=length(x)
  variance=var(x)
  left.pt=(n-1)*variance/qchisq(1-alpha/3, n-1)
  right.pt=(n-1)*variance/qchisq(2*alpha/3, n-1)
  return(c(left.pt,right.pt))
}

# make score approximate confidence interval for the coverage probability of the random confidence interval for sigma^2 from X_1,,,.,X_n from Normal distribution with alpha/2
score.ci.for.coverage.1=function(reps, n, mu, sigma, alpha, score.c.int){
  dum=rnorm(n=n*reps, mean=mu, sd=sigma)
  x.mat=matrix(dum, nrow=reps, ncol=n)
  ind=apply(x.mat, 1, make.CI.for.sigma.1, alpha=alpha)
  indx=apply(ind, 2, function(x){1*(x[1]<sigma^2)&&(x[2]>sigma^2)})
  s.ci=prop.test(sum(indx), n=reps, p=1-alpha, conf.level=score.c.int)
  return(s.ci)
}
# make score approximate confidence interval for the coverage probability of the random confidence interval for sigma^2 with alpha/3
score.ci.for.coverage.2=function(reps, n, mu, sigma, alpha, score.c.int){
  dum=rnorm(n=n*reps, mean=mu, sd=sigma)
  x.mat=matrix(dum, nrow=reps, ncol=n)
  ind=apply(x.mat, 1, make.CI.for.sigma.2, alpha=alpha)
  indx=apply(ind, 2, function(x){1*(x[1]<sigma^2)&&(x[2]>sigma^2)})
  s.ci=prop.test(sum(indx), n=reps, p=1-alpha, conf.level=score.c.int)
  return(s.ci)
}

set.seed(5701)
# pick some values for the parameters
reps = 10000
n = 10
mu = 68
sigma = 3
alpha = 0.05
score.c.int = 0.99
# 99% score approximate confidence interval for the coverage probability of the random confidence interval for sigma^2 in equation (1)
score.ci.for.coverage.1(reps, n, mu, sigma, alpha, score.c.int)
# 99% score approximate confidence interval for the coverage probability of the random confidence interval for sigma^2 in equation (2)
score.ci.for.coverage.2(reps, n, mu, sigma, alpha, score.c.int)

# 1. (c) WRONG!! WORK IN PROGRESS
reps=10000
n.vec=c(10,50,500,1000,3000)
alpha.vec=c(0.01, 0.05)
score.c.int=0.99


set.seed(5701)
for(n in n.vec){
  for(alpha in alpha.vec){
    results1=score.ci.for.coverage.1(reps, n, mu, sigma, alpha, score.c.int)
    results2=score.ci.for.coverage.2(reps, n, mu, sigma, alpha, score.c.int)
    width1=results1$conf.int[2]-results1$conf.int[1]
    width2=results2$conf.int[2]-results2$conf.int[1]
    if(width1<width2){
      print(paste('With n=',n,'and alpha=',alpha,'(1) is narrower'))
    }
    else
      print(paste('With n=',n,'and alpha=',alpha,'(2) is narrower'))
  }
}

# 1. (d)
# make score approximate confidence interval for the coverage probability of the random confidence interval for sigma^2 from X_1,,,.,X_n from Normal distribution with alpha/2
score.ci.for.coverage.exp.1=function(reps, n, rate, alpha, score.c.int){
  dum=rexp(n=n*reps, rate=1)
  sigma=rate
  x.mat=matrix(dum, nrow=reps, ncol=n)
  ind=apply(x.mat, 1, make.CI.for.sigma.1, alpha=alpha)
  indx=apply(ind, 2, function(x){1*(x[1]<sigma^2)&&(x[2]>sigma^2)})
  s.ci=prop.test(sum(indx), n=reps, p=1-alpha, conf.level=score.c.int)
  return(s.ci)
}

# make score approximate confidence interval for the coverage probability of the random confidence interval for sigma^2 with alpha/3
score.ci.for.coverage.exp.2=function(reps, n, rate, alpha, score.c.int){
  dum=rexp(n=n*reps, rate=rate)
  sigma=rate
  x.mat=matrix(dum, nrow=reps, ncol=n)
  ind=apply(x.mat, 1, make.CI.for.sigma.2, alpha=alpha)
  indx=apply(ind, 2, function(x){1*(x[1]<sigma^2)&&(x[2]>sigma^2)})
  s.ci=prop.test(sum(indx), n=reps, p=1-alpha, conf.level=score.c.int)
  return(s.ci)
}

set.seed(5701)
# pick some values for the parameters
reps=10000
rate=1
n.vec=c(10,50,500)
alpha.vec=c(0.01, 0.05)
score.c.int=0.99
res.mat=matrix(0,length(n.vec)*length(alpha.vec),6)

set.seed(5701)
k=1
for(n in n.vec){
  for(alpha in alpha.vec){
    res.mat[k,1]=n
    res.mat[k,2]=1-alpha
    res.mat[k,3:4]=score.ci.for.coverage.exp.1(reps, n, rate, alpha, score.c.int)$conf.int
    res.mat[k,5:6]=score.ci.for.coverage.exp.2(reps, n, rate, alpha, score.c.int)$conf.int
    k=k+1
  }
}
colnames(res.mat)=c("n","1-alpha","CI with alpha/2 Lower","CI with alpha/2 Upper", "CI with alpha/3 Lower", "CI with alpha/3 Upper")
res.mat
# 1. (e)
library(MASS)
# make score approximate confidence interval for the coverage probability of the random confidence interval for sigma^2 from X_1,,,.,X_n from n-variate Normal distribution with alpha/2
score.ci.for.coverage.mvrnorm.1=function(reps, n, mu, sigma, rho, alpha, score.c.int){
  mean.vec=rep(mu,n)
  covar.mat=matrix(NA, nrow=n, ncol=n)
  for(i in 1:n){
    for(j in 1:n){
      covar.mat[i,j]=sigma^2*rho^abs(i-j)
    }
  }
  x.mat=mvrnorm(reps, mu=mean.vec, Sigma=covar.mat)
  ind=apply(x.mat, 1, make.CI.for.sigma.1, alpha=alpha)
  indx=apply(ind, 2, function(x){1*(x[1]<sigma^2)&&(x[2]>sigma^2)})
  s.ci=prop.test(sum(indx), n=reps, p=1-alpha, conf.level=score.c.int)
  return(s.ci)
}

# make score approximate confidence interval for the coverage probability of the random confidence interval for sigma^2 from X_1,,,.,X_n from n-variate Normal distribution with 2*alpha/3
score.ci.for.coverage.mvrnorm.2=function(reps, n, mu, sigma, rho, alpha, score.c.int){
  mean.vec=rep(mu,n)
  covar.mat=matrix(NA, nrow=n, ncol=n)
  for(i in 1:n){
    for(j in 1:n){
      covar.mat[i,j]=sigma^2*rho^abs(i-j)
    }
  }
  x.mat=mvrnorm(reps, mu=mean.vec, Sigma=covar.mat)
  ind=apply(x.mat, 1, make.CI.for.sigma.2, alpha=alpha)
  indx=apply(ind, 2, function(x){1*(x[1]<sigma^2)&&(x[2]>sigma^2)})
  s.ci=prop.test(sum(indx), n=reps, p=1-alpha, conf.level=score.c.int)
  return(s.ci)
}

set.seed(5701)
# pick some values for the parameters
reps=10000
n.vec=c(10,50,500)
alpha.vec=c(0.01, 0.05)
score.c.int=0.99
mu=68
sigma=3
rho=0.7
res.mat=matrix(0,length(n.vec)*length(alpha.vec),6)

set.seed(5701)
k=1
for(n in n.vec){
  for(alpha in alpha.vec){
    res.mat[k,1]=n
    res.mat[k,2]=1-alpha
    res.mat[k,3:4]=score.ci.for.coverage.mvrnorm.1(reps, n, mu, sigma, rho, alpha, score.c.int)$conf.int
    res.mat[k,5:6]=score.ci.for.coverage.mvrnorm.2(reps, n, mu, sigma, rho, alpha, score.c.int)$conf.int
    k=k+1
  }
}
colnames(res.mat)=c("n","1-alpha","CI with alpha/2 Lower","CI with alpha/2 Upper", "CI with alpha/3 Lower", "CI with alpha/3 Upper")
res.mat

# 2. (a)
# The following function produces a realization of reps independent copies of random p-value in a two-sample t-test of
# H_0: mu1=mu2
# These realizations could be used to estimate the distribution of the random p-value.
# Returns a vector of reps realization of the random p-value.
est.pval.dist=function(n1, n2, mu1, mu2, sigma1, sigma2, alpha, reps=1e4){
  # create a vector that will record the observed p-values
  pvalue.list=numeric(reps)
  for(r in 1:reps){
    # generate a realization of X_1,..,X_n iid N(mu1, sigma^2)
    x.list=mu1+sigma1*rnorm(n1)
    # compute the realization of the (random) sample mean
    xbar=mean(x.list)
    # compute the realization of the (random) sample standard dev.
    sigmax=sd(x.list)
    # generate a realization of Y_1,..,Y_n iid N(mu2, sigma^2)
    y.list=mu2+sigma2*rnorm(n2)
    # compute the realization of the (random) sample mean
    ybar=mean(y.list)
    # compute the realization of the (random) sample standard dev.
    sigmay=sd(y.list)
    
    # compute Sp
    pooled.var=((n1-1)*sigmax^2+(n2-1)*sigmay^2)/(n1+n2-2)
    # compute test statistic realization
    t=(xbar-ybar)/(sqrt(pooled.var)*sqrt(1/n1+1/n2))
    
    # compute the two-sided alternative observed p-value
    pvalue.list[r]=2*pt(-abs(t), n1+n2-2)
  }
  return(pvalue.list)
}
set.seed(5701)
pvals=est.pval.dist(n1=25, n2=25, mu1=68, mu2=68, sigma1=3, sigma2=3, alpha=0.05)

# make histogram of the observed p-values
hist(pvals)
probs=seq(0.01, 0.99, 0.01)
plot(probs, quantile(pvals, probs))
abline(0,1)

# 2. (b)
# The following function computes a simulated estimate of the power of a two-sample t-test of
# H_0: mu1=mu2
# performed with an iid sample of size n from N(mu, sigma^2). These realizations could be used to estimate the distribution of the random p-value.
# If mu1=mu2, then this function computes a simulated estimate of the Type I error prob.
# Returns the proportion of reps for which the null hypothesis was rejected.
fast.est.power.ttest=function(n1, n2, mu1, mu2, alpha, sigma1, sigma2, reps=1e3){
  x.mat=matrix(mu1+sigma1*rnorm(reps*n1), nrow=reps, ncol=n1)
  # compute the reps observed sample means
  xbar.list=apply(x.mat, 1, mean)
  # compute the reps observed sample sd's
  x.s.list=apply(x.mat, 1, sd)
  
  y.mat=matrix(mu2+sigma2*rnorm(reps*n2), nrow=reps, ncol=n2)
  # compute the reps observed sample means
  ybar.list=apply(y.mat, 1, mean)
  # compute the reps observed sample sd's
  y.s.list=apply(y.mat, 1, sd)
  
  # compute the reps observed realizations of T
  # compute Sp
  pooled.var=numeric(reps)
  t.list=numeric(reps)
  for(i in 1:reps){
    pooled.var[i]=((n1-1)*x.s.list[i]^2+(n2-1)*y.s.list[i]^2)/(n1+n2-2)
    # compute test statistic realization
    t.list[i]=(xbar.list[i]-ybar.list[i])/(sqrt(pooled.var[i])*sqrt(1/n1+1/n2))
  }
  #return(t.list)
  
  # compute the proportion of rejected null hypotheses
  prop.rejected=mean(abs(t.list)>qt(1-alpha/2, n1+n2-2))
  return(prop.rejected)
}

# make power curve for two-samaple t-test as mu2 increases with mu1=68, alpha=0.05, sigma=3
set.seed(5701)
numpts=100
mu2.seq=seq(from=68, to=71, length.out = numpts)
power.est=numeric(numpts)
for(i in 1:numpts){
  power.est[i]=fast.est.power.ttest(n1=50, n2=50, mu1=68, alpha=0.05, mu2=mu2.seq[i], sigma1=3, sigma2=3, reps=1e3)
}
plot(mu2.seq, power.est, type="o", xlab="true mean", ylab="Est. rejection probabilty")

# 2. (c)
mu1=68
mu2=68.5
sigma1=3
sigma2=3
alpha=0.05
reps=1e3
length.n=100
n.seq=seq(from=500, by=50, length.out=length.n)
results=numeric(length(n.seq))
for(i in 1:length(n.seq)){
  results[i]=fast.est.power.ttest(n1=n.seq[i], n2=n.seq[i], mu1=68, mu2=68.5, alpha=0.05, sigma1=3, sigma2=3, reps=1e3)
}
n.seq[which.min(abs(results - 0.95))]

# 2. (d)
pvals.1=est.pval.dist(n1=100, n2=100, mu1=68, mu2=68, sigma1=3, sigma2=6)
pvals.2=est.pval.dist(n1=20, n2=100, mu1=68, mu2=68, sigma1=3, sigma2=6)
pvals.3=est.pval.dist(n1=100, n2=20, mu1=68, mu2=68, sigma1=3, sigma2=6)

# make QQplot comparing the data percentiles of the observed p-values to the Unif(0,1) percentiles
drawQQplot=function(pval){
  probs=seq(0.01, 0.99, 0.01)
  plot(probs, quantile(pval, probs))
  abline(0,1)
}

par(mfrow=c(1,3))
drawQQplot(pvals.1)
drawQQplot(pvals.2)
drawQQplot(pvals.3)
quantile(pvals.1, c(0.01, 0.05))
quantile(pvals.2, c(0.01, 0.05))
quantile(pvals.3, c(0.01, 0.05))

# 2. (e)
library(MASS)
pval.mvrnorm=function(reps, n, mu1, mu2, sigma, rho){
  mu.vec=c(mu1, mu2)
  covar=matrix(0, nrow=2, ncol=2)
  covar[1,1]=covar[2,2]=sigma^2
  covar[1,2]=covar[2,1]=sigma^2*rho
  pvals=numeric(reps)
  for(r in 1:reps){
    xy=mvrnorm(n, mu.vec, covar)
    result=t.test(xy[,1], xy[,2], var.equal = TRUE)
    pvals[r]=result$p.value
  }
  pvals
}
set.seed(5701)
pvals.i=pval.mvrnorm(reps=1e3, n=20, mu1=68, mu2=68, sigma=3, rho=0.01)
pvals.ii=pval.mvrnorm(reps=1e3, n=20, mu1=68, mu2=68, sigma=3, rho=0.1)
par(mfrow=c(1,3))
drawQQplot(pvals.i)
drawQQplot(pvals.ii)
quantile(pvals.i, c(0.01, 0.05))
quantile(pvals.ii, c(0.01, 0.05))

# 3. (d)
# function to calculate simulated SE for shrinkage estimator
shrinkage.est=function(x,mu,a){
  (a*mean(x)-mu)^2
}
# function to calculate simulated SE for likelihood estimator
likelihood.est = function(x,mu){
  (mean(x)-mu)^2
}
# make function to return simulated SEs for likelihood and shrinkage estimates
get.se = function(x,mu,a){ 
  diff=likelihood.est(x,mu)-shrinkage.est(x,mu,a)
  return(c(shrinkage.est(x,mu,a), likelihood.est(x,mu), diff))
}

# make function that calculate simulated MSE for shrinkage est, likelihood est, and difference est with 1e4 reps.
get.mse=function(reps=1e4,mu,n){
  x.list=matrix(rexp(n*reps,1/mu), nrow=reps, ncol=n)
  ahat=n/(n+1)
  # squared errors
  ses=apply(x.list, 1, get.se, mu=mu, a=ahat)
  ci.1=t.test(x=ses[1,], conf.level=0.99)$conf.int
  ci.2=t.test(x=ses[2,], conf.level=0.99)$conf.int
  ci.3=t.test(x=ses[3,], conf.level=0.99)$conf.int
  # make mean squared error
  mses=apply(ses,1,mean)
  return(list("MSE"=mses,"CI.shrinkage"=ci.1,"CI.likelihood"=ci.2,"CI.differ"=ci.3))
}

k=1
n.vec=c(5,10,50)
mu.vec=c(0.5,1,10)
res.mat=matrix(0,length(n.vec)*length(mu.vec),12)
for(n in n.vec){
  for(mu in mu.vec){
    MSEs=get.mse(reps=1e3,mu=mu,n=n)
    res.mat[k,1:2]=c(n,mu)
    res.mat[k,3] = MSEs$MSE[1]
    res.mat[k,4] = (1/(n+1))*mu^2
    res.mat[k,5:6]= MSEs$CI.shrinkage
    res.mat[k,7] = MSEs$MSE[2]
    res.mat[k,8] = mu^2/n
    res.mat[k,9:10]= MSEs$CI.likelihood
    res.mat[k,11:12] = MSEs$CI.differ
    k=k+1
  }
}
colnames(res.mat)=c("n","mu","Sim. Shrinkage estimate","Formula estimate","Sim. shrinkage 99% Lower","Upper","Sim. Likelihood estimate","Formula estimate","Sim. Likelihood 99% Lower","Upper","Sim. Difference estimate 99% lower", "Upper")
res.mat
