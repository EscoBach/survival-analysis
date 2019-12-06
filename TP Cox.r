## Main ##
mlcox <- function(n,N,a,beta){

	p <- length(beta)
	mle <- array(dim = c(N,p))

	upper <- rep(max(beta)+5,p)
	lower <- rep(min(beta)-5,p)

	for(i in 1:N){
		data <- gen(n,a,beta)

		optimr <- optim(beta,pl,data = data,control = list(fnscale = -1),upper = upper,lower = lower,method = 'L-BFGS-B')
		mle[i,] <- optimr$par
	}

	return(mle)
}
## Générateur de données simulées pour un modèle de Cox simple ##
gen <- function(n,a,beta){

	p <- length(beta)
	
	time <- vector(length = n)
	status <- vector(length = n)

	cnames <- sapply(1:p,function(t) sprintf('Z%g',t))
	z <- matrix(nrow = n, ncol = p)
	colnames(z) <- cnames

	for(i in 1:n){

		z[i,] <- runif(p)
		surv <- -log(1-runif(1))/(a*exp(sum(beta*z[i,])))
		cens <- rexp(1,a)

		time[i] <- ifelse(surv <= cens, surv, cens)
		status[i] <- ifelse(surv <= cens, 1, 0)
	}
	return(as.data.frame(cbind(time,status,z)))
}
## Estimateur de la fonction de risque commune ##
lambda0 <- function(data,t,beta){

	n <- dim(data)[1]
	p <- length(beta)

	v <- vector(length = n,mode='numeric')

	for(i in 1:n){
		v[i] <- (data$time[i] >= t)*exp(sum(beta*data[i,2+1:p]))
	}

	if(sum(v) > 0){
		return(sum(data$status == 1 & data$time == t)/sum(v))
	}else{
		return(0)
	}
}
## Log-vraisemblance partielle de Cox ##
pl <- function(beta,data){

	n <- dim(data)[1]
	p <- length(beta)

	lv <- array(dim = c(n,n))
	v <- vector(length = n,mode='numeric')

	for(t in 1:n){

		for(k in 1:n){
			v[k] <- (data$time[k] >= data$time[t])*exp(sum(beta*data[k,2+1:p]))
		}

		for(i in 1:n){

			lv[i,t] <- (data$time[i] == data$time[t] & data$status[i] == 1)*(sum(beta*data[i,2+1:p])-log(sum(v)))
		}
	}

	return(mean(lv))
}