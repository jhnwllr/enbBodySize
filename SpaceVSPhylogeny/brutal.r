 library(mvtnorm)
 library(MASS) 

# -----------------------------------------  Notes --------------------------------------------------
# In general these functions take a <formula> to fit to the data. <data> should contain the
# variables referenced in the model fomula. The code performs casewise deletion, so note that
# if <data> contains additonal variables with missing values, these will affect the analysis. 
# <phylomat> is the variance-covariance matrix for the phylogeny (e.g. use vcv.phylo) and
# for spatial analysis <spatialmat> is the spatial matrix. 
# IMPORTANT: <data> should have rownames that match the tip.labels of the phylogeny. An
# error will result if this is not the case. For generating a distance matrix lat and lon need to 
# be specified in addition to a list of rownames. 

# LIST OF FUNCTIONS:
# 				pglmSpatial			--- fits a spatial pglm
#				pglmSpatialFit		--- co-estimates lambda and phi whilst fitting a pglm
#				dist.mat 				--- generates a variance matrix based on distances using lats and lons
# See brutalDemo.r for examples based on the carnivore phylogeny.
# -----------------------------------------------------------------------------------------------------




# -----------------------------------------------
# - Simple function for calculating distances 
# -----------------------------------------------
geodetic <- function(l1, t1, l2, t2) {
	
	l1 <- l1 / 360.0 * 2.0 * pi
	l2 <- l2 / 360.0 * 2.0 * pi
	t1 <- t1 / 360.0 * 2.0 * pi
	t2 <- t2 / 360.0 * 2.0 * pi
	
	dist <- 6371.0 * acos( sin(t1) * sin(t2) + cos(t1) * cos(t2) * cos(l2 - l1) )

	return(dist)
	
	}
	
# -----------------------------------------------
# ------- Generates a distance matrix -------
# -----------------------------------------------
dist.mat <- function(lat, lon, rwnms) {
	n <- length( lat )
	mat <- matrix(0, n, n )
	
	for(i in 1:n) {
		for( j in 1:n) {
			mat[i, j] <- geodetic( lon[i], lat[i], lon[j], lat[j] )
			}
		}
		
	mdist <- geodetic(-90,-90, 90,90)
	mat <- mdist - mat
	diag(mat) <- mdist
	mat <- mat / mdist
	rownames(mat) <- rwnms
	return(mat)
	
	}

# -------------------------------------------------------------------------------------
# This fits a GLM, correcting for phylogeny and space, with the option
# of setting the value of lambda and rho, the indices of phylogenetic and spatial
# dependence (default is 1.0)
# -------------------------------------------------------------------------------------
pglmSpatial <- function(formula, data, phylomat, spatialmat, lambda = 1.0, phi = 0.0) {

	prune <- function(dat, Vmat, Dmat) {
	
# ------ Delete from here if you want to take the risk! ---------------------	
		# Makes sure data and matrix are in the same order
		nms <- row.names(dat)
		if(length(nms) == 0) stop("Need to supply row names for the data")
		idx <- sort(nms, index.return = TRUE)$ix
		dat <- dat[idx,]
		Vnms <- row.names(Vmat)
		if(length(Vnms) == 0) stop("Need to supply row names for the Variance matrix")
		idx <- sort(Vnms, index.return = TRUE)$ix
		Vmat <- Vmat[idx, idx]
		
		nms <- sort(nms)
		Vn <- sort(row.names(Vmat))
		idx <- which(nms != Vn)

		
		if(length(idx) > 0) stop("Error, taxon names do not match variance matrix")
		
		# Q&D to deall with spatial data
		Snms <- row.names( Dmat)
		if(length(Snms) == 0) stop("Need to supply row names for the Spatial matrix")
		idx <- sort(Snms, index.return = TRUE)$ix
		Dmat <- Dmat[idx, idx]
		
		
		
# ------ Delete to here if you want to take the risk! ---------------------
		
		complete <- complete.cases(dat)
		idx <- which(complete == TRUE)
		Vmat <- Vmat[idx, idx]
		Dmat <- Dmat[idx, idx]
		dat <- dat[idx,]

		return(list(dat = dat , Vmat = Vmat, Dmat = Dmat))
		}
	
	Dfun <- function(Cmat) {
		iCmat <- solve(Cmat,  tol = .Machine$double.eps)
		svdCmat <- La.svd(iCmat)
		D <- svdCmat$u %*% diag(sqrt( svdCmat$d )) %*% t(svdCmat$v)
		return( t(D) )
		}
			
	lamTrans <- function(Vmat, lambda) {
		V1 <- Vmat
		diag(V1) <- 0
		V2 <- diag( diag(Vmat), ncol = length(Vmat[1,]), nrow = length(Vmat[,1]))
		Vmat <- V1 * lambda + V2
		return(Vmat)
		}
			
	resVar <- function(y, Vmat, p) {
		iV <- solve(Vmat, tol = .Machine$double.eps)
		e <- y - p
		s2 <- crossprod(e, iV %*% e)
		if( s2 < 0) { cat("ERROR -- negative variance\n")
			s2 <- s2 * -1}
		n <- length(y) 
		return( s2 / (n- nx) )
		}
		
	# Estimates the GLS parameters for given data
	get.coeffs <- function(Y, V, X) {
		iV <- solve(V, tol = .Machine$double.eps)
		xVix <- crossprod(X, iV %*% X)
		xViy <- crossprod(X, iV %*% Y)
		mu <- solve(xVix, tol = .Machine$double.eps) %*% xViy 	#This is  a bad thing to do!!!!
		return(mu)
		}

	# Estimates the variance of a given trait (accounting for phylogeny)
	est.var <- function(y, V, x, mu ) {
		iV <- solve(V, tol = .Machine$double.eps)
		e <- y - x %*% mu
		s2 <- crossprod(e, iV %*% e)
		n <- length(y) 
		k <- length(x[1,])
		return( s2 / (n- k) )
		}
	
	# Full ML estimation for given x and V
	log.likelihood <- function(y, x, V) {
		mu <- get.coeffs(y, V, x)
		s2 <- est.var(y, V, x, mu)
		n <- length(x[,1])
		logDetV <- determinant(Vmat, logarithm = TRUE)$modulus[1]
		ll <- -n / 2.0 * log( 2 * pi) - n / 2.0 * log(s2) - logDetV / 2.0 - (n - 1)/2.0
		ypred <- x%*%mu	
		return( list(ll = ll, mu = mu, s2 = s2) )
		}
		
		
	null.var <- function(y, V) {
		X <- matrix(1, nrow = length(y))
		mu <- get.coeffs(y, V, X)
		return(est.var(y, V, X, mu))
		}

	# Scale matrices
	phylomat <- lamTrans(phylomat, lambda)
	spatialmat <- spatialmat / max(spatialmat) * max(phylomat)

	prune.dat <- prune(data, phylomat, spatialmat)
	phylomat <- prune.dat$Vmat
	spatialmat <- prune.dat$Dmat
	
	Vmat <- (1 - phi) * as.matrix(phylomat) + phi * as.matrix(spatialmat)
	

	
	data <- prune.dat$dat
	nm <- names(data)
	
	n <- length(data[,1])
	
	# Get the design matrix
	m <- model.frame(formula, data)
	y <- m[,1]
	x <- model.matrix(formula, m)
	k <- length(x[1,])
	
	namey <- names(m)[1]
	
	ll <- log.likelihood(y, x, Vmat)
	
	log.lik <- ll$ll	

	aic <- -2 * log.lik + 2 * k
	aicc <- -2 * log.lik + 2 * k + ((2*k*(k+1))/(n-k-1))
	
	coeffs <- ll$mu
	coeffs <- data.frame(t(coeffs))
	names(coeffs) <- colnames(x)
	varNames = names(m)

	
	pred <- x %*% ll$mu 
	
	res <- y - pred
	D <- Dfun(Vmat)
	pres <- D %*% res
	
	fm <- list(coef = coeffs, aic = aic, log.lik = log.lik)
	
	logDetV <- determinant(Vmat, logarithm = TRUE)$modulus[1]
 	
	logLikY <- -n / 2.0 * log( 2 * pi) - n / 2.0 * log( (n - k) * ll$s2 / n) - logDetV / 2.0  - n / 2.0
	
	RMS <- ll$s2
	RSSQ <- ll$s2 * (n - k)
	NMS <- RMS
	NSSQ <- RSSQ
	
	if(k > 0) {
		NMS <- null.var(y, Vmat)
		NSSQ <- NMS * (n - 1)
		}

	# Bits for parameter errors	
	errMat <- t(x)%*% solve(Vmat) %*% x  
	errMat <- solve(errMat) * RMS[1] 
	sterr <- diag(errMat)
	sterr <- sqrt(sterr)
	
	
	ret <- list(model = fm, formula = formula, logLikY = logLikY, RMS = RMS, NMS = NMS, NSSQ = NSSQ[1], RSSQ = RSSQ[1], 
	aic = aic, aicc = aicc, n = n, k = k, sterr = sterr, vcv = errMat, fitted = pred, residuals = res, phyres = pres, x = x, data = data,  
	varNames = varNames, y = y, V = Vmat, lambda = lambda, phi = phi, L0 = NULL, L1 = NULL, LamOptimised = FALSE, 
	phiOptimised = FALSE, namey = namey)
	
	class(ret) <- "pglm"
	return(ret)
	
	}
	
pglmSpatialFit <- function (formula, data, phylomat, spatialmat, plotit = FALSE) {
		 
		 
			# -----------------------------------------------
			# - Finds the joint ML values of lambda  & phi --
			# -----------------------------------------------
			optimizePhiLambda <- function( formula, data, phylomat, spatialmat) {
				optFun <- function(pars) -1*pglmSpatial( formula, data, phylomat, spatialmat, lambda = pars[1], phi = pars[2])$logLikY
				of <- optim( c(0.5, 0.5),  method = "L-BFGS-B", optFun, lower = c(0,0), upper = c(1,1) )
				return(of)
				}

			
			opt <- optimizePhiLambda( formula, data, phylomat, spatialmat)
			MLlambda <- opt$par[1]
			MLphi <- opt$par[2]
			ML <-  -1 * opt$value
			
			fm <- pglmSpatial( formula, data, phylomat, spatialmat, MLlambda, MLphi )
			fm$phiOptimised <- TRUE
			fm$LamOptimised <- TRUE
			return(fm)

	
	}













 



