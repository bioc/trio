poly4root <- function(a){
	if(length(a)!=5)
		stop("a must have length 5.")
	if(!is.numeric(a))
		stop("a must be numeric.")
	alpha <- -3/8 * (a[2]/a[1])^2 + a[3]/a[1]
	beta <- 0.125 * (a[2]/a[1])^3 - 0.5*a[2]*a[3]/a[1]^2 + a[4]/a[1]
	gamma <- -3/256 * (a[2]/a[1])^4 + a[3]*a[2]^2/(16*a[1]^3) - a[2]*a[4]/(4*a[1]^2) + a[5]/a[1]
	s <- c(-1, -1, 1, 1)
	r <- c(-1, 1, -1, 1)
	if(beta==0){
		x <- -a[2]/(4*a[1]) + s * sqrt((-alpha + r*sqrt(alpha^2 - 4*gamma)) / 2)
		return(x[!is.na(x)])
	}
	P <- -alpha^2/12 - gamma
	Q <- -alpha^3/108 + alpha*gamma/3 - beta^2/8
	if(P==0)
		y <- -5/6*alpha - Q^(1/3)
	else{
		d <- Q^2/4+P^3/27
		if(d<=0)
			return(casusIrred(alpha,beta, gamma, a[1], a[2]))
		U <- -(Q/2 + sqrt(Q^2/4 + P^3/27))^(1/3)
		y <- -5/6*alpha + U - P/(3*U)
	}
	w <- sqrt(alpha + 2*y)
	tmp <- -(3*alpha + 2*y + s*2*beta/w)
	if(any(abs(tmp) < 10^-8)){
		tmpids <- abs(tmp) < 10^-8
		tmp[tmpids] <- round(tmp[tmpids], 8)
	}
	x <- -a[2]/(4*a[1]) + 0.5*(s*w + r * sqrt(tmp))
	x[!is.na(x)]
}

casusIrred <- function(p4, q4, r4, a, b){
	r <- -2*p4
	s <- p4^2-4*r4
	t <- q4^2
	p <- s-r^2/3
	q <- 2*r^3/27 - r*s/3 + t
	R <- (q/2)^2 + (p/3)^3
	u <- sqrt(-(p/3)^3)
	tmp <- -q/(2*u)
	if(abs(abs(tmp)-1) < 10^-8)
		tmp <- round(tmp, 8)
	y <- (-1)^(0:2)*2*u^(1/3) * cos(acos(tmp)/3+(0:2)*pi/3)
	z <- sqrt(-(y-r/3))
	sigma <- -0.5 * q4/prod(z)
	mat <- matrix(c(1,1,-1,-1,1,-1,1,-1,1,-1,-1,1),4)
	yroot <- as.vector(sigma * mat %*% z)
	yroot - b/(4*a)
}



poly4rootMat <- function(amat){
	if(!is.matrix(amat))
		stop("amat must be a matrix.")
	if(ncol(amat) != 5)
		stop("amat must have 5 columns.")
	if(!is.numeric(amat))
		stop("amat must be a numeric matrix.")
	alpha <- -3/8 * (amat[,2]/amat[,1])^2 + amat[,3]/amat[,1]
	beta <- 0.125 * (amat[,2]/amat[,1])^3 - 0.5*amat[,2]*amat[,3]/amat[,1]^2 + amat[,4]/amat[,1]
	gamma <- -3/256 * (amat[,2]/amat[,1])^4 + amat[,3]*amat[,2]^2/(16*amat[,1]^3) - 
		amat[,2]*amat[,4]/(4*amat[,1]^2) + amat[,5]/amat[,1]
	r <- -2*alpha
	s <- alpha^2 - 4*gamma
	t <- beta^2
	p <- s - r^2/3
	q <- 2*r^3/27 - r*s/3 + t
	Rpq <- (q/2)^2 + (p/3)^3
	if(any(p==0))
		stop("p==0 currently not implemented.")
	idsBeta <- beta == 0
	idsR <- !idsBeta & (Rpq>0)
	mat.root <- matrix(NA, nrow(amat), 4)
	if(any(idsBeta))
		mat.root[idsBeta,] <- rootBeta(amat[idsBeta, 1], amat[idsBeta, 2], alpha[idsBeta], gamma[idsBeta], s[idsBeta])
	if(any(idsR))
		mat.root[idsR,] <- rootRleq0(amat[idsR,1], amat[idsR,2], alpha[idsR], beta[idsR], gamma[idsR]) 
	idsCI <- !idsBeta & !idsR
	mat.root[idsCI,] <- rootCasusIrr(amat[idsCI,1], amat[idsCI,2], beta[idsCI], p[idsCI], q[idsCI], r[idsCI]) 
	mat.root
}


rootBeta <- function(a, b, alpha, gamma, s){
	ssign <- c(-1, -1, 1, 1)
	rsign <- c(-1, 1, -1, 1)
	mat <- matrix(NA, length(a), 4)
	fac <- -b/(4*a)
	for(i in 1:length(a))
		mat[i,] <- fac[i] + ssign * sqrt((-alpha[i] + rsign*sqrt(s[i]))/2)
	mat
}

rootRleq0 <- function(a, b, alpha, beta, gamma){
	P <- -alpha^2/12 - gamma
	Q <- -alpha^3/108 + alpha*gamma/3 - beta^2/8
	qp <- Q^2/4 + P^3/27
	if(any(abs(qp) < 10^-8)){
		qpids <- abs(qp) < 10^-8
		qp[qpids] <- round(qp[qpids], 8)
	}
	u <- -(Q/2 + sqrt(qp))^(1/3)
	y <- -5/6 * alpha + u - P/(3*u)
	w <- sqrt(alpha + 2*y)
	rsign <- c(-1, -1, 1, 1)
	ssign <- c(-1, 1, -1, 1)
	mat <- matrix(NA, length(a), 4)
	fac <- -b/(4*a)
	for(i in 1:length(a)){
		tmp <- -(3*alpha[i] + 2*y[i] + ssign*2*beta[i]/w[i])
		if(any(abs(tmp) < 10^-8)){
			tmpids <- abs(tmp) < 10^-8
			tmp[tmpids] <- round(tmp[tmpids], 8)
		}
		mat[i,] <- fac[i] + 0.5 * (ssign * w[i] + rsign *sqrt(tmp))
	}
	mat
}

rootCasusIrr <- function(a, b, beta, P, Q, r){
	u <- sqrt(-(P/3)^3)
	fac1 <- 2*u^(1/3)
	tmp <- -Q/(2*u)
	if(any(abs(abs(tmp)-1) < 10^-8)){
		tmpids <- abs(abs(tmp)-1) < 10^-8
		tmp[tmpids] <- round(tmp[tmpids], 8)
	}
	fac2 <- acos(tmp)/3
	maty <- matrix(0, length(a), 3)
	maty[,1] <- fac1 * cos(fac2)
	maty[,2] <- -fac1 * cos(fac2 + pi/3)
	maty[,3] <- fac1 * cos(fac2 + 2*pi/3)
	matz <- sqrt(-(maty-r/3))
	zfac <- matz[,1] * matz[,2] * matz[,3]
	sigma <- -0.5 * beta / zfac
	mat <- matrix(c(1,1,-1,-1,1,-1,1,-1,1,-1,-1,1), 3, byrow=TRUE)
	mat.root <- sigma * matz %*% mat
	mat.root <- mat.root - b/(4*a)
	mat.root
}
		