

if(FALSE) { # simulate multivariate normal spatial data 

N = 100
lat = 1:N
lon = 1:N

library(tmvtnorm)
library(matrixcalc)
library(corpcor)
dm = as.matrix(dist(data.frame(lat,lon),upper = TRUE)) # distance matrix 
dm = make.positive.definite(dm) 

library(MASS)
library(ape)
library(geiger)

x = mvrnorm(n = 1,rep(0,N),dm) # lat and long

tree = rcoal(100)
y = sim.char(tree, 0.02, 1)[,1,]

e = rnorm(100,0,0.3)
y = 1*x + e

D = data.frame(y,x)
rownames(D) = tree$tip.label

source("C:/Users/Owner/Desktop/envBodySize/SpaceVsPhylogeny/brutal.r") # load library 


Vcarn = vcv.phylo(tree)
dat = read.table("C:/Users/Owner/Desktop/envBodySize/SpaceVsPhylogeny/carnDat.txt", header = TRUE)
# dat = order.D(dat)

# DM = dist.mat(lat,lon,rownames(D)) # distance matrix
DM = dm
V = vcv.phylo(tree)

fm1 = pglmSpatial(y ~ 1,D,V,DM)
fm2 = pglmSpatial(y ~ x,D,V,DM)

fm1$model$coef
fm2$model$coef

fm1$aic
fm2$aic
}


if(FALSE) { # simulate brownian motion characters 

library(geiger)
library(ape)

DL = list()
for(i in 1:1000) { 
tree = rcoal(100)
treeSignal = lambdaTree(tree, 0.1)

y = sim.char(tree, 0.02, 1)[,1,] # response variable simulation 
ySignal = sim.char(treeSignal, 0.02, 1)[,1,] # response variable simulation 

x = rnorm(100)
e = rnorm(100,0,0.3)
y = 0.1*x + e
ySignal = 0.1*x + e


D = data.frame(y,x,sp=tree$tip.label)

library(nlme)
V = corBrownian(1,phy=tree) 
rownames(D) = D$sp
fit = gls(y ~ x,correlation=V,data=D)
summary(fit)
pvalue1 = summary(fit)$tTable[2,4]

DSignal = data.frame(ySignal,x,sp=tree$tip.label)

V = corBrownian(1,phy=treeSignal) 
rownames(DSignal) = DSignal$sp
fit = gls(ySignal ~ x,correlation=V,data=DSignal)
pvalue2 = summary(fit)$tTable[2,4]

DL[[i]] = c(pvalue1,pvalue2)
}

D = do.call("rbind",DL)

jpeg("C:/Users/Owner/Desktop/plot.jpg")
boxplot(D)
dev.off()
}

# par(mfrow=c(2,1))
# plot(treeSignal)
# plot(tree)




# V = corPagel(1, phy = tree, form = ~1, fixed = FALSE)


# fit = gls(hwl ~ Temp + Prec + cover + div_brd + div_mam,correlation=V,data=D)




# library(nlme)
# V = corBrownian(1,phy=tree) 
# rownames(D) = D$GenusSpecies
# V = corPagel(1, phy = tree, form = ~1, fixed = FALSE)
# fit = gls(hwl ~ div_brd + div_mam + Temp + Prec,correlation=V,data=D)
# fit = gls(hwl ~ Temp + Prec + cover + div_brd + div_mam,correlation=V,data=D)



if(FALSE) {


library(ape)
source("C:/Users/Owner/Desktop/envBodySize/SpaceVsPhylogeny/brutal.r") # load library 


# list.files("C:/Users/Owner/Desktop/envBodySize/SpaceVsPhylogeny/brutal.r")

phy = read.tree("C:/Users/Owner/Desktop/envBodySize/SpaceVsPhylogeny/carn.phy")
Vcarn = vcv.phylo(phy)
dat = read.table("C:/Users/Owner/Desktop/envBodySize/SpaceVsPhylogeny/carnDat.txt", header = TRUE)
# dat = order.D(dat)

Dcarn = dist.mat( dat$AvgLat, dat$AvgLon, rownames(dat) )

# Vcarn <- order.V(Vcarn)
# Dcarn <- order.V(Dcarn)


# fm1 = pglmSpatial(Mass ~ 1, dat, Vcarn, Dcarn)
fm2 = pglmSpatial(Mass ~ AvTemp, dat, Vcarn, Dcarn)
fm2$model$coef

lm(Mass~AvTemp,dat)


# fm1$aic
# str(fm2)

# fm = pglm( Mass ~ AvTemp, dat, Vcarn)
# fm1 = pglmEstLambda( Mass ~ AvTemp, dat, Vcarn)

# fm2 = pglmSpatialFit( Mass ~ AvTemp, dat, Vcarn, Dcarn)

# Vcarn = Vcarn / max(Vcarn)
# Dcarn = Dcarn / max(Dcarn)
# fm3 <- regress( Mass ~ AvTemp, ~ Vcarn,identity = TRUE,  data =  dat )
}


