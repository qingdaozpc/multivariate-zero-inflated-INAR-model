
# Simulation of trivariate INAR_MZIH
INAR_MZIH_sim <- function(n,time,p,pi0,pi,beta){
  data <- matrix(0,n*time,9)
  N.pre <- N0
  N.cur <- matrix(0,n,3)
  for (t in 1:time){
    lambda <- exp(X%*%beta)
    y <- 1+cbind(rpois(n,lambda[,1]),rpois(n,lambda[,2]),rpois(n,lambda[,3]))
    u <- cbind(rbinom(n,1,pi[1]),rbinom(n,1,pi[2]),rbinom(n,1,pi[3]))
    u0 <- (runif(n)<pi0)*1
    r <- u0*u*y
    for(i in 1:n){
      N.cur[i,] <- rbinom(3,N.pre[i,],p)+r[i,]
    }
    data[((t-1)*n+1):(t*n),] <- cbind(N.pre,N.cur,X)
    N.pre <- N.cur
    t <- t+1
  }
  return(data)
}

#The pmf of hurdle Poisson
HP_pmf <- function(z,pi,lambda){
  if(z==0)
    pmf <- 1-pi
  else
    pmf <- pi*dpois(z-1,lambda)
  return(pmf)
}

#The pmf of MZIH
MZIH_pmf <- function(z,pi0,pi,lambda){
  if(all(z==0))
    pmf <- 1-pi0 + pi0*prod(1-pi)
  else 
    pmf <- pi0*HP_pmf(z[1],pi[1],lambda[1])*HP_pmf(z[2],pi[2],lambda[2])*
      HP_pmf(z[3],pi[3],lambda[3])
  return(pmf)
}

# The pmf of INAR_MZIH
INAR_MZIH_pmf <- function(to,from,p,pi0,pi,lambda){
  pmf_array <- array(0,dim=c(to[1]+1,to[2]+1,to[3]+1))
  for(i in 0:(to[1])){
    for(j in 0:(to[2])){
      for(k in 0:(to[3])){
        pmf_array[i+1,j+1,k+1] <- MZIH_pmf(c(i,j,k),pi0,pi,lambda)*
          dbinom(to[1]-i,from[1],p[1])*
          dbinom(to[2]-j,from[2],p[2])*
          dbinom(to[3]-k,from[3],p[3])
      }
    }
  }
  return(sum(pmf_array))
}

# The pmf of INAR_MZIH with fixed margin1
INAR_MZIH1_pmf <- function(to,from,p,pi0,pi,lambda){
  pmf_matrix <- matrix(0,to[2]+1,to[3]+1)
  for(j in 0:(to[2])){
    for(k in 0:(to[3])){
      pmf_matrix[j+1,k+1] <- MZIH_pmf(c(0,j,k),pi0,pi,lambda)*
        dbinom(to[1],from[1],p[1])*
        dbinom(to[2]-j,from[2],p[2])*
        dbinom(to[3]-k,from[3],p[3])
    }
  }
  return(sum(pmf_matrix))
}

# The pmf of INAR_MZIH with fixed margin2
INAR_MZIH2_pmf <- function(to,from,p,pi0,pi,lambda){
  pmf_matrix <- matrix(0,to[1]+1,to[3]+1)
  for(i in 0:(to[1])){
    for(k in 0:(to[3])){
      pmf_matrix[i+1,k+1] <- MZIH_pmf(c(i,0,k),pi0,pi,lambda)*
        dbinom(to[1]-i,from[1],p[1])*
        dbinom(to[2],from[2],p[2])*
        dbinom(to[3]-k,from[3],p[3])
    }
  }
  return(sum(pmf_matrix))
}

# The pmf of INAR_MZIH with fixed margin3
INAR_MZIH3_pmf <- function(to,from,p,pi0,pi,lambda){
  pmf_matrix <- matrix(0,to[1]+1,to[2]+1)
  for(i in 0:(to[1])){
    for(j in 0:(to[2])){
      pmf_matrix[i+1,j+1] <- MZIH_pmf(c(i,j,0),pi0,pi,lambda)*
        dbinom(to[1]-i,from[1],p[1])*
        dbinom(to[2]-j,from[2],p[2])*
        dbinom(to[3],from[3],p[3])
    }
  }
  return(sum(pmf_matrix))
}

# loglik function of INAR_MZIH
INAR_MZIH_loglik <- function(data,p,pi0,pi,beta){
  N.pre <- data[,1:3]
  N.cur <- data[,4:6]
  X <- data[,7:9]
  lambda <- exp(X%*%beta)
  sum(sapply(1:nrow(data),function(k) 
    log(INAR_MZIH_pmf(N.cur[k,],N.pre[k,],p,pi0,pi,lambda[k,]))))
}

# EM for INAR_MZIH 
EM_INAR_MZIH <- function(data,par){
  N.pre <- data[,1:3]
  N.cur <- data[,4:6]
  X <- data[,7:9]
  nt <- nrow(data)
  
  # posterior expectation of y
  post_Ey <- function(p,pi0,pi,beta){
    indicator <- diag(3)
    y <- matrix(0,nt,3)
    lambda <- exp(X%*%beta)
    for(i in 1:nt){
      if(N.pre[i,1]==0|N.cur[i,1]==0){
        y1 <- 0
      }
      else{
        y1 <- p[1]*N.pre[i,1]*
          INAR_MZIH_pmf(N.cur[i,]-indicator[1,],N.pre[i,]-indicator[1,],p,pi0,pi,lambda[i,])/
          INAR_MZIH_pmf(N.cur[i,],N.pre[i,],p,pi0,pi,lambda[i,])
      }
      
      if(N.pre[i,2]==0|N.cur[i,2]==0){
        y2 <- 0
      }
      else{
        y2 <- p[2]*N.pre[i,2]*
          INAR_MZIH_pmf(N.cur[i,]-indicator[2,],N.pre[i,]-indicator[2,],p,pi0,pi,lambda[i,])/
          INAR_MZIH_pmf(N.cur[i,],N.pre[i,],p,pi0,pi,lambda[i,])
      }
      
      if(N.pre[i,3]==0|N.cur[i,3]==0){
        y3 <- 0
      }
      else{
        y3 <- p[3]*N.pre[i,3]*
          INAR_MZIH_pmf(N.cur[i,]-indicator[3,],N.pre[i,]-indicator[3,],p,pi0,pi,lambda[i,])/
          INAR_MZIH_pmf(N.cur[i,],N.pre[i,],p,pi0,pi,lambda[i,])
      }
      y[i,] <- c(y1,y2,y3)
    }
    return(y)
  }
  
  # posterior expectation of uv
  post_Euv <- function(p,pi0,pi,beta){
    lambda <- exp(X%*%beta)
    uv <- numeric(nt)
    for(i in 1:nt){
      uv[i] <- prod(dbinom(N.cur[i,],N.pre[i,],p))*(1-pi0)/
        INAR_MZIH_pmf(N.cur[i,],N.pre[i,],p,pi0,pi,lambda[i,])
    }
    return(uv)
  }
  
  # posterior expectation of tau
  post_Etau <- function(p,pi0,pi,beta){
    lambda <- exp(X %*% beta)
    tau <- matrix(0,nt,3)
    for(i in 1:nt){
      tau1 <- 1-INAR_MZIH1_pmf(N.cur[i,],N.pre[i,],p,pi0,pi,lambda[i,])/
        INAR_MZIH_pmf(N.cur[i,],N.pre[i,],p,pi0,pi,lambda[i,])
      tau2 <- 1-INAR_MZIH2_pmf(N.cur[i,],N.pre[i,],p,pi0,pi,lambda[i,])/
        INAR_MZIH_pmf(N.cur[i,],N.pre[i,],p,pi0,pi,lambda[i,])
      tau3 <- 1-INAR_MZIH3_pmf(N.cur[i,],N.pre[i,],p,pi0,pi,lambda[i,])/
        INAR_MZIH_pmf(N.cur[i,],N.pre[i,],p,pi0,pi,lambda[i,])
      tau[i,] <- c(tau1,tau2,tau3)
    }
    return(tau)
  }
  
  # initial values of parameters
  p <- par$p
  pi0 <- par$pi0
  pi <- par$pi
  beta <- par$beta
  name_par <- names(c(p=numeric(3),pi0=numeric(1),pi=numeric(3),beta1=numeric(3),
                      beta2=numeric(3),beta3=numeric(3)))
  # initial loglikelihood
  loglik <- INAR_MZIH_loglik(data,p,pi0,pi,beta)
  dif <- 1
  while(dif>1e-5){
    #E-step
    # posterior expectation of y
    y <- post_Ey(p,pi0,pi,beta)
    # posterior expectation of r
    r <- N.cur-y
    # posterior expectation of uv
    w <- post_Euv(p,pi0,pi,beta)
    # posterior expectation of tau
    tau <- post_Etau(p,pi0,pi,beta)
    
    #M-step
    # update the parameter p
    p <- apply(y,2,sum)/apply(N.pre,2,sum)
    # update the parameter pi0
    pi0 <- 1-sum(w)/nt
    # update the parameter pi
    pi <- apply(tau,2,sum)/ sum(1-w) 
    # update the parameter beta
    lambda <- exp(X%*%beta)
    s1 <- t(X)%*%(r[,1]-tau[,1]-tau[,1]*lambda[,1])
    h1 <- -t(X)%*%(tau[,1]*lambda[,1]*X)
    beta[,1] <- beta[,1]-solve(h1)%*%s1
    
    s2 <- t(X)%*%(r[,2]-tau[,2]-tau[,2]*lambda[,2])
    h2 <- -t(X)%*%(tau[,2]*lambda[,2]*X)
    beta[,2] <- beta[,2]-solve(h2)%*%s2
    
    s3 <- t(X)%*%(r[,3]-tau[,3]-tau[,3]*lambda[,3])
    h3 <- -t(X)%*%(tau[,3]*lambda[,3]*X)
    beta[,3] <- beta[,3]-solve(h3)%*%s3
    
    # update the loglikelihood
    loglik.new <- INAR_MZIH_loglik(data,p,pi0,pi,beta)
    dif <- loglik.new-loglik
    loglik <- loglik.new 
  }
  par <- c(p,pi0,pi,beta)
  names(par) <- name_par
  return(par)
}

#simulation
set.seed(123)
n=2000
time=5
X <- cbind(rep(1,n),rnorm(n),rbinom(n,1,0.5))
N0 <- cbind(rpois(n,0.1),rpois(n,0.2),rpois(n,0.3))

MC <- function(i){
  set.seed(i)
  data <-INAR_MZIH_sim(n=2000,time=5,p=c(0.1,0.2,0.3),pi0=0.5,pi=c(0.3,0.2,0.1),
                       beta=cbind(c(-3,-1,1),c(-2,-1,-1),c(-1,1,-1)))
  par <- EM_INAR_MZIH(data,par=list(p=rep(0.5,3),pi0=0.5,pi=rep(0.5,3),
                                    beta=matrix(0,3,3)))
  return(par)
}

library(doParallel)
cl <- makeCluster(20)
registerDoParallel(cl)

para.mat <- foreach(i = 1:100, .combine = rbind) %dopar% {MC(i)}
apply(para.mat,2,mean)
apply(para.mat,2,sd)/sqrt(100)

stopImplicitCluster()
stopCluster(cl)
