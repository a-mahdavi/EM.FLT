# R code for fitting finite mixtures (g components) of FLT distribution

mixFLT <- function(y, g=1, w=1, mu, s, del, nu, family="FLT", iter.max=100, tol=10^-6, CML=T, get.init = TRUE, group=F){  
  begin <- proc.time()[3]  ; y[which(y==0|y==1)]=10^(-30)
	dFLT<-function (x, mu = 0, sigma = 1, del, nu) 
		{
    ql <- qlogis(x)
          ifelse(x <= 0 | x >= 1, 0, dt(abs(ql-mu)/sigma-del,nu)/(2*x*(1-x)*sigma*pt(del,nu)) )
  	  }
 	dmixFLT <- function(y, w, mu, s, del, nu){
    	d <- 0 ; g <- length(w)
    	for ( j in 1:g)
      d <- d + w[j]*dFLT(y, mu[j], s[j], del[j], nu[j])
    return(d) }
	dFLN<-function (x, mu = 0, sigma = 1, del) 
	{
    ql <- qlogis(x)
          ifelse(x <= 0 | x >= 1, 0, dnorm(abs(ql-mu)/sigma-del)/(2*x*(1-x)*sigma*pnorm(del)) )
   	 }
  dmixFLN <- function(y, w, mu, s, del){
    d <- 0 ; g <- length(w)
    for ( j in 1:g)
      d <- d + w[j]*dFLN(y, mu[j], s[j], del[j])
    return(d) }
  dLN<-function (x, mu = 0, sigma = 1) 
  {
    ql <- qlogis(x)
    d <-  dnorm(ql, mu, sigma)/x/(1 - x)
	d[which(d==0)] <- 10^-30
	return(d)
  }
  dmixLN <- function(y, w, mu, s){
    d <- 0 ; g <- length(w)
    for ( j in 1:g)
      d <- d + w[j]*dLN(y, mu[j], s[j])
    return(d) }
  n <- length(y)   ;     dif <- 1;        count <- 0 
  if (get.init == TRUE) {
    init <- kmeans(y, g,  algorithm="Hartigan-Wong")
    w <- init$size/n ;mu <- s <- NULL
	for( j in 1:g){
	x <- y[init$cluster==j]
	mu[j] <- mean(log(x/(1-x)))
	s[j] <- 1/length(x)*sum((log(x/(1-x))-mu[j])^2)
		}
    del <- runif(g, -1,1)
 	 }
 if(family=="FLT"){
	nu <- runif(g,.5,5)
    LL <- 1 
    while ((dif > tol) && (count <= iter.max)) {
      z.hat <- matrix(0,n,g)
      # E step
      for (j in 1:g){
        z.hat[,j] <- w[j]*dFLT(y,mu[j],s[j], del[j], nu[j])/dmixFLT(y, w, mu, s, del, nu)
          # MCE steps
        w[j] <- sum(z.hat[,j])/n
	  tauh <- (nu[j]+1)/(nu[j]+(abs(qlogis(y)-mu[j])/s[j]-del[j])^2) 
        mu[j] <- sum(z.hat[,j]*tauh*(qlogis(y)-s[j]*del[j]*sign(qlogis(y)-mu[j])))/sum(tauh*z.hat[,j])
	a <- del[j]*sum(z.hat[,j]*tauh*abs(qlogis(y)-mu[j])) ; b <- sum(z.hat[,j]*tauh*(qlogis(y)-mu[j])^2)
	s1 <- (-a+sqrt(a^2+4*sum(z.hat[,j])*b))/(2*sum(z.hat[,j])) ; s2 <- (-a-sqrt(a^2+4*sum(z.hat[,j])*b))/(2*sum(z.hat[,j]))
		if ( s1 >0) 
		s[j] <- s1 else 
		s[j] <- s2
        del[j] <- (sum(z.hat[,j]*tauh*abs(qlogis(y)-mu[j])/s[j])-n*dt(del[j],nu[j])/pt(del[j],nu[j]))/sum(z.hat[,j]*tauh)
	nu[j] <- optim(nu[j],function(x){
            -sum(z.hat[,j]*log(dFLT(y,mu[j],s[j],del[j],x)))
          },method="L-BFGS-B",lower=0.1,upper=100)$par
       }
      LL.new <- sum(log(dmixFLT(y,w,mu,s,del,nu))) # log-likelihood function
      count <- count +1 
      dif <- abs(LL.new/LL-1)
      LL <- LL.new
      cat('iter =', count, '\tloglike =', LL.new, '\n')
    } 
    aic <- -2 * LL.new + 2 * (4*g+g-1)
    bic <- -2 * LL.new + log(n) * (4*g+g-1)
    edc <- -2 * LL.new + 0.2*sqrt(n) * (4*g+g-1)
    end <- proc.time()[3]
    time <- end-begin
    obj.out <- list(family=family,w=w, mu=mu, sigma=s , del=del, nu=nu, loglik=LL.new, aic=aic, bic=bic, edc=edc, iter=count,elapsed=as.numeric(time))
  }
  if(family=="FLN"){
    LL <- 1 
    while ((dif > tol) && (count <= iter.max)) {
      z.hat <- matrix(0,n,g)
      # E step
      for (j in 1:g){
        z.hat[,j] <- w[j]*dFLN(y,mu[j],s[j], del[j])/dmixFLN(y, w, mu, s, del)
          # MCE steps
        w[j] <- sum(z.hat[,j])/n
        mu[j] <- sum(z.hat[,j]*(qlogis(y)-s[j]*del[j]*sign(qlogis(y)-mu[j])))/sum(z.hat[,j])
	a <- del[j]*sum(z.hat[,j]*abs(qlogis(y)-mu[j])) ; b <- sum(z.hat[,j]*(qlogis(y)-mu[j])^2)
	s1 <- (-a+sqrt(a^2+4*sum(z.hat[,j])*b))/(2*sum(z.hat[,j])) ; s2 <- (-a-sqrt(a^2+4*sum(z.hat[,j])*b))/(2*sum(z.hat[,j]))
		if ( s1 >0) 
		s[j] <- s1 else 
		s[j] <- s2
	 del[j] <- sum(z.hat[,j]*abs(qlogis(y)-mu[j])/s[j])/sum(z.hat[,j])-dnorm(del[j])/pnorm(del[j])
       }
      LL.new <- sum(log(dmixFLN(y,w,mu,s,del))) # log-likelihood function
      count <- count +1 
      dif <- abs(LL.new/LL-1)
      LL <- LL.new
      cat('iter =', count, '\tloglike =', LL.new, '\n')
    } 
    aic <- -2 * LL.new + 2 * (3*g+g-1)
    bic <- -2 * LL.new + log(n) * (3*g+g-1)
    edc <- -2 * LL.new + 0.2*sqrt(n) * (3*g+g-1)
    end <- proc.time()[3]
    time <- end-begin
    obj.out <- list(family=family,w=w, mu=mu, sigma=s , del=del, loglik=LL.new, aic=aic, bic=bic, edc=edc, iter=count,elapsed=as.numeric(time))
  }
 if(family=="LT"){
	nu <- runif(g,.5,5); del<-rep(0,g)
    LL <- 1 
    while ((dif > tol) && (count <= iter.max)) {
      z.hat <- matrix(0,n,g)
      # E step
      for (j in 1:g){
        z.hat[,j] <- w[j]*dFLT(y,mu[j],s[j], 0, nu[j])/dmixFLT(y, w, mu, s, 0, nu)
          # MCE steps
        w[j] <- sum(z.hat[,j])/n
	  tauh <- (nu[j]+1)/(nu[j]+(abs(qlogis(y)-mu[j])/s[j]-del[j])^2) 
        mu[j] <- sum(z.hat[,j]*tauh*(qlogis(y)-s[j]*del[j]*sign(qlogis(y)-mu[j])))/sum(tauh*z.hat[,j])
	a <- del[j]*sum(z.hat[,j]*tauh*abs(qlogis(y)-mu[j])) ; b <- sum(z.hat[,j]*tauh*(qlogis(y)-mu[j])^2)
	s1 <- (-a+sqrt(a^2+4*sum(z.hat[,j])*b))/(2*sum(z.hat[,j])) ; s2 <- (-a-sqrt(a^2+4*sum(z.hat[,j])*b))/(2*sum(z.hat[,j]))
		if ( s1 >0) 
		s[j] <- s1 else 
		s[j] <- s2
	nu[j] <- optim(nu[j],function(x){
            -sum(z.hat[,j]*log(dFLT(y,mu[j],s[j],0,x)))
          },method="L-BFGS-B",lower=0.1,upper=100)$par
       }
      LL.new <- sum(log(dmixFLT(y,w,mu,s,del,nu))) # log-likelihood function
      count <- count +1 
      dif <- abs(LL.new/LL-1)
      LL <- LL.new
      cat('iter =', count, '\tloglike =', LL.new, '\n')
    } 
    aic <- -2 * LL.new + 2 * (3*g+g-1)
    bic <- -2 * LL.new + log(n) * (3*g+g-1)
    edc <- -2 * LL.new + 0.2*sqrt(n) * (3*g+g-1)
    end <- proc.time()[3]
    time <- end-begin
    obj.out <- list(family=family,w=w, mu=mu, sigma=s , nu=nu, loglik=LL.new, aic=aic, bic=bic, edc=edc, iter=count,elapsed=as.numeric(time))
  }
  if(family=="LN"){
    LL <- 1 
    while ((dif > tol) && (count <= iter.max)) {
      z.hat  <- matrix(0,n,g)
      # E step
      for (j in 1:g){
        z.hat[,j] <- w[j]*dLN(y,mu[j],s[j])/dmixLN(y, w, mu, s)
        # MCE steps
        w[j] <- sum(z.hat[,j])/n
        mu[j] <- sum(z.hat[,j]*qlogis(y))/sum(z.hat[,j])
        s[j] <- sqrt(sum(z.hat[,j]*(qlogis(y)-mu[j])^2)/sum(z.hat[,j]))
      }
      LL.new <- sum(log(dmixLN(y,w,mu,s))) # log-likelihood function
      count <- count +1 
      dif <- abs(LL.new/LL-1)
      LL <- LL.new
    cat('iter =', count, '\tloglike =', LL.new, '\n')
    } 
    aic <- -2 * LL.new + 2 * (2*g+g-1)
    bic <- -2 * LL.new + log(n) * (2*g+g-1)
    edc <- -2 * LL.new + 0.2*sqrt(n) * (2*g+g-1)
    end <- proc.time()[3]
    time <- end-begin
    obj.out <- list(family=family,w=w, mu=mu, sigma=s , loglik=LL.new, aic=aic, bic=bic, edc=edc, iter=count,elapsed=as.numeric(time))
  }
   if (group)
    obj.out$group <- apply(z.hat, 1, which.max)
  obj.out
}



