library(ape)

# -----------------------------------------  Notes --------------------------------------------------
# In general these functions take a <formula> to fit to the data. <data> should contain the
# variables referenced in the model fomula. The code performs casewise deletion, so note that
# if <data> contains additonal variables with missing values, these will affect the analysis. 
# <phylo> is the phylogeny as a phylo object. 
# IMPORTANT: <data> should have rownames that match the tip.labels of the phylogeny. An
# error will result if this is not the case. In the case of the spatial models you need to also 
# include the latitude and longitude. These need to be lists of lat and lons, with rownames.

# LIST OF FUNCTIONS:
# 				cpglm 					--- fits a pglm by method of contrasts
#				cpglm.lambda		--- fits a pglm with additional <lambda> parameter
#				cpglm.estlambda	--- fits a pglm whilst estimating lambda parameter
#				cpglm.spatial			--- fits a pglm with both lambda and phi
# 				optimise.pglm		--- co-estimates lambda and phi whilst fitting a pglm
#
# See distance2.demo for examples based on the carnivore phylogeny.
# -----------------------------------------------------------------------------------------------------


# -----------------------------------------------
# --------PGLM based on contrasts -----------
# -----------------------------------------------
cpglm <- function( formula, data, phylo) {
	# Prune data & tree
	idx.retain <- which( complete.cases( data ) == TRUE)
	idx.drop <- which( complete.cases( data ) == FALSE)
	if(length(idx.drop) >0){
		phylo <- drop.tip( phylo, rownames(data)[ idx.drop ])}
	data <- data[idx.retain,]
	
	cpic <- function(x) pic(x, phylo)
	cdata <- apply( data, 2, cpic)
	cdata <- data.frame(cdata)
	colnames(cdata) <- colnames(data)
	fm <- lm(formula, cdata, y = TRUE)
	
	# Calculate the log-likelihood
	vars <- pic( data[,1], phylo, var = TRUE)[,2]
	V <- contrastData(data[,1], phylo)$V
	u <- residuals(fm) 
	
	n <- length(u) 
	sigma2 <- sum(u^2) / (n+1)
	logLik <- -0.5 *( (n + 1) * log( 2 * pi * sigma2) + sum(log( vars )) + sum(u^2) / sigma2  )
	logLik <- logLik - 0.5 * log(V)
	
	return(list( model = fm, logLik = logLik) )
	
	}

# -----------------------------------------------
# -- PGLM based on contrasts with lambda -
# -----------------------------------------------	
cpglm.lambda <- function( formula, data, phylo, lambda) {

# Prune data & tree
	idx.retain <- which( complete.cases( data ) == TRUE)
	idx.drop <- which( complete.cases( data ) == FALSE)
	if(length(idx.drop) >0){
		phylo <- drop.tip( phylo, rownames(data)[ idx.drop ])}
	data <- data[idx.retain,]
	
	phy <- lambda.trans( phylo, lambda )
	cpic <- function(x) pic(x, phy)
	cdata <- apply( data, 2, cpic)
	cdata <- data.frame(cdata)
	colnames(cdata) <- colnames(data)
	fm <- lm(formula, cdata, y = TRUE)
	
	# Calculate the log-likelihood
	vars <- pic( data[,1], phy, var = TRUE)[,2]
	V <- contrastData(data[,1], phy)$V	
	u <- residuals(fm) 
	n <- length(u)
	sigma2 <- sum(u^2) / (n+1)
	logLik <- -0.5 *( (n + 1) * log( 2 * pi * sigma2) + sum(log( vars )) + sum(u^2) / sigma2  )
	logLik <- logLik - 0.5 * log(V)
	
	return(list( model = fm, logLik = logLik) )
	
	}


# -----------------------------------------------
# -- PGLM based on contrasts with kappa -
# -----------------------------------------------	
cpglm.kappa <- function( formula, data, phylo, kappa) {

# Prune data & tree
	idx.retain <- which( complete.cases( data ) == TRUE)
	idx.drop <- which( complete.cases( data ) == FALSE)
	if(length(idx.drop) >0){
		phylo <- drop.tip( phylo, rownames(data)[ idx.drop ])}
	data <- data[idx.retain,]
	
	phy <- kappa.trans( phylo, kappa )
	cpic <- function(x) pic(x, phy)
	cdata <- apply( data, 2, cpic)
	cdata <- data.frame(cdata)
	colnames(cdata) <- colnames(data)
	fm <- lm(formula, cdata, y = TRUE)
	
	# Calculate the log-likelihood
	vars <- pic( data[,1], phy, var = TRUE)[,2]
	V <- contrastData(data[,1], phy)$V	
	u <- residuals(fm) 
	n <- length(u)
	sigma2 <- sum(u^2) / (n+1)
	logLik <- -0.5 *( (n + 1) * log( 2 * pi * sigma2) + sum(log( vars )) + sum(u^2) / sigma2  )
	logLik <- logLik - 0.5 * log(V)
	
	return(list( model = fm, logLik = logLik) )
	
	}
	
# -----------------------------------------------
# -- PGLM based on contrasts with kappa and lambda -
# -----------------------------------------------	
cpglm.lamkap <- function( formula, data, phylo, lambda, kappa) {
	
	# Prune data & tree
	idx.retain <- which( complete.cases( data ) == TRUE)
	idx.drop <- which( complete.cases( data ) == FALSE)
	if(length(idx.drop) >0){
		phylo <- drop.tip( phylo, rownames(data)[ idx.drop ])}
	data <- data[idx.retain,]
	
	phy <- kappa.trans( phylo, kappa )
	phy <- lambda.trans( phy, lambda)
	cpic <- function(x) pic(x, phy)
	cdata <- apply( data, 2, cpic)
	cdata <- data.frame(cdata)
	colnames(cdata) <- colnames(data)
	fm <- lm(formula, cdata, y = TRUE)
	
	# Calculate the log-likelihood
	vars <- pic( data[,1], phy, var = TRUE)[,2]
	V <- contrastData(data[,1], phy)$V	
	u <- residuals(fm) 
	n <- length(u)
	sigma2 <- sum(u^2) / (n+1)
	logLik <- -0.5 *( (n + 1) * log( 2 * pi * sigma2) + sum(log( vars )) + sum(u^2) / sigma2  )
	logLik <- logLik - 0.5 * log(V)
	
	return(list( model = fm, logLik = logLik) )
	
	}

# -----------------------------------------------
# -- PGLM based on contrasts for the -------
# ---- Estimation of lambda ------------------
cpglm.estlambda <- function( formula, data, phylo ) {
	
			optimizeLambda <- function( formula, data, phylo ) {
				optFun <- function(pars) -1 * cpglm.lambda( formula, data, phylo, lambda = pars[1])$logLik
				of <- optim( c(0.5),  method = "L-BFGS-B", optFun, lower = c(0), upper = c(1) )
				return(of)
				}
		
		# Prune data & tree
		idx.retain <- which( complete.cases( data ) == TRUE)
		idx.drop <- which( complete.cases( data ) == FALSE)
		if(length(idx.drop) >0){
			phylo <- drop.tip( phylo, rownames(data)[ idx.drop ])}
		data <- data[idx.retain,]
		
		opt <- optimizeLambda( formula, data, phylo )
		MLlambda <- opt$par[1]
		ML <-  -1 * opt$value
		
		fm <- cpglm.lambda( formula, data, phylo, MLlambda)$model
		
		return( list(model = fm, lambda = MLlambda, ML = ML)  )
	
	}



	
# -----------------------------------------------
# -- Spatial PGLM based on contrasts with -
# --- lambda and rho transformations ------
# -----------------------------------------------	
cpglm.spatial <- function( formula, data, phylo, lat, lon, lambda, phi) {


	# Prune data & tree
	idx.retain <- which( complete.cases( data ) == TRUE)
	idx.drop <- which( complete.cases( data ) == FALSE)
	if(length(idx.drop) >0){
				phylo <- drop.tip( phylo, rownames(data)[ idx.drop ])}
	data <- data[idx.retain,]
	phynames <- phylo$tip.label
	phyidx <- seq( 1:length(phynames) )
	datnames <- rownames(data)
	
	sortphy <- sort(phynames,index.return = TRUE)
	phyidx <- phyidx[sortphy$ix]
	
	sortdat <- sort(datnames, index.return = TRUE)
	data <- data[sortdat$ix,]
	lat <- lat[sortdat$ix]
	lon <- lon[sortdat$ix]
	
	data <- data[phyidx,]
	lat <- lat[phyidx]
	lon <- lon[phyidx]
	

	phy <- lambda.trans( phylo, lambda )
	cpic <- function(x) contrastData(x, phy)$con
	cdata <- apply( data, 2, cpic)
	
	names( lat ) <- rownames( data )
	names( lon) <- rownames( data)
	
	# get the nodal values & other contrasts information
	latCon <- contrastData(lat, phy)
	lonCon <- contrastData(lon, phy)
	
	# These are the nodal values
	lonDat1 <- lonCon$X1
	lonDat2 <- lonCon$X2
	latDat1 <- latCon$X1
	latDat2 <- latCon$X2
	
	# Get distances, contrasts and phylogenetic variances
	distances <- geodetic(lonDat1, latDat1, lonDat2, latDat2)
	phyVar <- latCon$vars
	
	# Scale these
	distances <- distances / max(distances) * max( phyVar )
	totalVar <-  phi * distances + (1 - phi) * phyVar
	
	for(i in 1:dim(cdata)[2]) { cdata[,i] <- cdata[,i] / sqrt (totalVar) }
	cdata <- data.frame(cdata)
	colnames(cdata) <- colnames(data)
	fm <- lm(formula, cdata) 
	
	V <- latCon$V	
	u <- residuals(fm)
	
	
	n <- length(u) 
	sigma2 <- sum(u^2) / (n + 1)
	logLik <- -0.5 *( (n + 1) * log( 2 * pi * sigma2) + sum(log( totalVar )) + sum(u^2) / sigma2  )
	logLik <- logLik - 0.5 * log(V)
	
	
	return(list(model = fm, logLik = logLik) )
	
	}


# -----------------------------------------------
# ------------- Fits phi & lambda -------------
# -----------------------------------------------	
optimise.cpglm <- function( formula, data, phylo, lat, lon) {
	
			optimizePhiLambda <- function( formula, data, phylo, lat, lon) {
				optFun <- function(pars) -1 * cpglm.spatial( formula, data, phylo, lat, lon, lambda = pars[1], phi = pars[2])$logLik
				of <- optim( c(0.5, 0.5), optFun, lower = c(0.01,0.01), upper = c(0.99,0.99) )
				return(of)
				}
			
			# Prune data & tree
			idx.retain <- which( complete.cases( data ) == TRUE)
			idx.drop <- which( complete.cases( data ) == FALSE)
			if(length(idx.drop) >0){
				phylo <- drop.tip( phylo, rownames(data)[ idx.drop ])}
			data <- data[idx.retain,]
			
			opt <- optimizePhiLambda( formula, data, phylo, lat, lon)
			MLlambda <- opt$par[1]
			MLphi <- opt$par[2]
			ML <-  -1 * opt$value
			
			fm <- cpglm.spatial( formula, data, phylo, lat, lon, MLlambda, MLphi )$model
			
			return( list(model = fm, lambda = MLlambda, phi = MLphi, ML = ML)  )
					
	}



# -----------------------------------------------
# ------------- Stolen from APE --------------
# -----------------------------------------------
tree.height <- function (phy, tol = .Machine$double.eps^0.5) 
{
    if (class(phy) != "phylo") 
        stop("object \"phy\" is not of class \"phylo\".")
    if (is.null(phy$edge.length)) 
        stop("the tree has no branch lengths.")
    n <- length(phy$tip.label)
    n.node <- phy$Nnode
    xx <- numeric(n + n.node)
    for (i in 1:dim(phy$edge)[1]) xx[phy$edge[i, 2]] <- xx[phy$edge[i, 
        1]] + phy$edge.length[i]
 	return(xx[1:n])
}


# -----------------------------------------------
# ------------- Lambda Transformation -----------
# -----------------------------------------------
lambda.trans <- function( phy, lambda, height = NULL ) {
    n <- length(phy$tip.label)
    if( is.null( height) == TRUE) height <- tree.height(phy)
    interns <- which( phy$edge[, 2] > n )
    externs <- which(  phy$edge[, 2] <= n)
    phy$edge.length[ interns ] <-  phy$edge.length[ interns ]  * lambda
    phy$edge.length[ externs ] <-  phy$edge.length[ externs ]  * lambda + (1 - lambda) * height
	return( phy )
	}


# -----------------------------------------------
# ------------- Kappa Transformation -----------
# -----------------------------------------------
kappa.trans <- function( phy, kappa) {
	idx <- which( phy$edge.length > 0)
	phy$edge.length[idx] <- ( phy$edge.length[idx] )^kappa
	return( phy )
	}



# -----------------------------------------------
# ----- My modification to Emmanuel's code 
# -----------------------------------------------
pic.Rob <- function (x, phy, scaled = TRUE, var.contrasts = TRUE) 
{
    if (class(phy) != "phylo") 
        stop("object 'phy' is not of class \"phylo\"")
        
    if (is.null(phy$edge.length)) 
        stop("your tree has no branch lengths: you may consider setting them equal to one, or using the function `compute.brlen'.")
       
    nb.tip <- length(phy$tip.label)
    nb.node <- phy$Nnode
    
    if (nb.node != nb.tip - 1) 
        stop("'phy' is not rooted and fully dichotomous")
        
    if (length(x) != nb.tip) 
        stop("length of phenotypic and of phylogenetic data do not match")
        
    if (any(is.na(x))) 
        stop("the present method cannot (yet) be used directly with missing data: you may consider removing the species with missing data from your tree with the function `drop.tip'.")
        
    phy <- reorder(phy, "pruningwise")
    
    phenotype <- numeric(nb.tip + nb.node)
    
    if (is.null(names(x))) {
        phenotype[1:nb.tip] <- x
    }
    else {
        if (all(names(x) %in% phy$tip.label)) 
            phenotype[1:nb.tip] <- x[phy$tip.label]
        else {
            phenotype[1:nb.tip] <- x
            warning("the names of argument \"x\" and the names of the tip labels did not match: the former were ignored in the analysis.")
        }
    }
    
    contr <- var.con <- numeric(nb.node)
    
    ans <- .C("pic", as.integer(nb.tip), as.integer(nb.node), 
        as.integer(phy$edge[, 1]), as.integer(phy$edge[, 2]), 
        as.double(phy$edge.length), as.double(phenotype), as.double(contr), 
        as.double(var.con), as.integer(var.contrasts), as.integer(scaled), 
        PACKAGE = "ape")
        
    contr <- ans[[7]]
    
    if (var.contrasts) {
        contr <- cbind(contr, ans[[8]])
        dimnames(contr) <- list(1:nb.node + nb.tip, c("contrasts", 
            "variance"))
    } else names(contr) <- 1:nb.node + nb.tip
        
   	idx <- which(ans[[3]] == (nb.tip+1) )
    root.v <- ans[[5]][idx]
    V <- root.v[1]*root.v[2]/(sum(root.v))
    
    return(list(contr = contr, root.v = root.v, V = V, ans = ans))
}


# -----------------------------------------------
# - Simple function for calculating distances ---
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
# ----- Explicitly extract data for contrasts ---
# -----------------------------------------------
contrastData <- function(x, phylo) {
	
	pR <-  pic.Rob(x, phylo)
	a <- pR$ans
	V <- pR$V
	
	
	nCon <- length( a[[7]] )
	
	vxy <- matrix(0, nrow = nCon)
	pos <- matrix(0, nrow = nCon)
	con <- matrix(0, nrow = nCon)
	X1 <- matrix(0, nrow = nCon)
	X2 <- matrix(0, nrow = nCon)
	
	for(i in 0:(nCon - 1)){
		
		xp <- i * 2 + 1
		yp <- i * 2 + 2
		
		vxy[i + 1] <- a[5][[1]][xp] + a[5][[1]][yp]
		pos[i + 1] <- a[[3]][xp] - (nCon + 1)
		
		xc <- a[[4]][xp]
		yc <- a[[4]][yp]
		
		con[i+1] <-  ( a[[6]][xc] - a[[6]][yc] ) 
		
		X1[i+1] <- a[[6]][xc]
		X2[i+1] <- a[[6]][yc] 
		
		}
	return( list( X1  = X1, X2 = X2, vars = vxy, contrasts = con, V = V, ans = a) )
	}




