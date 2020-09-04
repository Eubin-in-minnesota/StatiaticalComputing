#################### 1. (b)
rtri3=function(n){
  # generate a realization of U_1,...,U_n iid Uni(0,1)
  u.list=runif(n)
  # compute the corresponding realization of Z_1,...,Z_n iid g
  z.list=3*sqrt(u.list)
  return(z.list)
}

#################### 1. (c)
set.seed(5701)
probs=seq(0.01, 0.99, 0.01)
plot(3*sqrt(probs), quantile(rtri3(n=2000), probs), xlab="g distribution percentile", ylab="Data percentile")
abline(0,1)

#################### 2. (b)
# generate a realization of X_1,...,X_n iid f
# Arguments: n, the number of realizations
# Returns: x.list, the vector of n realizations
#          k.list, the vector of n iteration counts needed for each realization. 
#            
rquad3=function(n){
  # allocate the memory for the realization of the randon sample from f
  x.list=numeric(n)
  # allocate the memery for the number of iterations of the rejection sampling algorithm required for each of the n draws
  k.list=numeric(n)
  
  for(i in 1:n){
    # set the iteration counter to zero
    k=0
    accept=FALSE
    while(!accept) {
      k=k+1
      z=rtri3(1)
      u=runif(1)
      if (z>3*u)
        accept=TRUE
    }
    k.list[i]=k
    x.list[i]=z
  }
  return(list(x.list=x.list, k.list=k.list))
}

#################### 2. (c)
set.seed(5701)
n=1000
draw=rquad3(n=n)
hist(draw$x.list)
mean(draw$k.list)

#################### 3. (b)
set.seed(5701)
# make single realization from Normal distribution using box-muller
box.muller=function(mu, sigma){
  u=runif(n=2)
  tmp=sqrt(-2*log(u[1]))
  x1=tmp*cos(2*pi*u[2])
  x2=tmp*sin(2*pi*u[2])
  x=c(x1,x2) # two samples from standard normal
  a=sample(1:2, size=1, replace=T) # randomly choose one of them
  x.list=mu+sigma*x[a]
  return(x.list)
}

myrtnorm=function(n, mu, sigma, a, b){
  t.list=numeric(n)
  x=numeric(n)
  for(i in 1:n){
    accept=FALSE
    while(!accept){
      x=box.muller(mu=mu, sigma=sigma)
      if(a<x && x<b){
        t.list[i]=x
        accept=TRUE
      }
    }
  }
  return(t.list)
}

#################### 3. (d)
n=500
mu=0
sigma=1
a=-2
b=2
c=1/(pnorm(q=b, mean=mu, sd=sigma)-pnorm(q=a, mean=mu, sd=sigma))
c
myrtnorm(n, mu, sigma, a, b)
hist(myrtnorm(n, mu, sigma, a, b))

#################### 4. (a)
set.seed(5701)
# The following function generates a realization of X_1,...,X_n iid Exponential distribution with mu
# Argument: n, the sample size
#           mu, the user-specified mean of the exponential distribution
# Returns a vector of n entries with the generated realization of X_1,...,X_n
myrexp=function(n, mu){
  # generate a realization of U_1,..,U_n iid Uni(0,1)
  u.list=runif(n)
  # compute the corresponding realization of X_1,...,X_n iid Exp(mu)
  x.list=(-1)*mu*log(1-u.list)
  return(x.list)
}

#################### 4. (b)
set.seed(5701)
probs=seq(0.01,0.99,0.01)
mu=5
x.list=(-1)*mu*log(1-probs)
plot(x.list, quantile(myrexp(n=1000, mu=5), probs), xlab="Expected Exponential distribution percentile", ylab="Data percentile")
abline(0,1)

#################### 4. (c)
run.sim.exp=function(n, mu, reps){
  # ybar.list is the vector that stores realization of Y_1_bar,...,Y_n_bar
  ybar.list=numeric(reps)
  
  for(i in 1:reps){
    # generate a realization of Y_1_bar,...,Y_n_bar
    ybar.list[i]=mean(myrexp(n,mu))
  }
  
  probs=ppoints(reps)
  plot(qnorm(probs, mean=mu, sd=mu/sqrt(n)), quantile(ybar.list, probs), xlab="Expected Exponential distribution percentile", ylab="Data percentile")
  abline(0,1)
  return(ybar.list)
}

#################### 4. (d)
run.sim.exp(30,3.4,1e4)

#################### 5. (a)
mymvrnorm=function(n, mu, Sigma){
  # we generate an n by p matrix, with each row is drawn from the p-variate normal, mean mu covariance matrix sigma
  ## do an eigendecomposition of Sigma
  ei=eigen(Sigma, symmetric=TRUE)
  ## the square root of Sigma is
  Sigma.sqrt=ei$vectors %*% diag(ei$values^0.5)%*%t(ei$vectors)
  p=dim(Sigma)[1]
  mu=rep(mu, p)
  Z=matrix(rnorm(n*p), nrow=n, ncol=p)
  Y=rep(1,n) %*% t(mu) + Z %*% Sigma.sqrt
  return(Y)
}

#################### 5. (c)
n=10
mu=68
sigma=3
# covariance matrix for H_i
Sigma=matrix(NA, nrow=n, ncol=n)
for(i in 1:n){
  for(j in 1:n){
    if(i==j){
      Sigma[i,j]=sigma^2
    }
    else
      Sigma[i,j]=sigma^2*0.6
  }
}

reps=10000
h=mymvrnorm(reps, mu, Sigma)
hbar=apply(h,1,mean)
x=mymvrnorm(reps, mu, diag(rep(sigma^2), n))
xbar=apply(x,1,mean)

mean(hbar)
var(hbar)
mean(xbar)
var(xbar)
