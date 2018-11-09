
# Try with carnivores

cphy <- read.tree("carn.phy")
dat <- read.table( "carnDat.txt", header = TRUE)

fm <- cpglm(Mass ~ AvTemp - 1, dat, cphy)
fm1 <- cpglm.lambda(Mass ~ AvTemp - 1, dat, cphy, 1)
fm2 <- pglm( Mass ~ AvTemp, dat, vcv.phylo(cphy) )
fm3 <-  cpglm.spatial( Mass ~ AvTemp, dat, cphy, dat$AvgLat, dat$AvgLon, 1, 0)

fit.temp <- optimise.cpglm( AvTemp ~  -1, dat, cphy, dat$AvgLat, dat$AvgLon)
fit.Mass <- optimise.cpglm( Mass ~  1, dat, cphy, dat$AvgLat, dat$AvgLon)
fit.Bergmann <- optimise.cpglm( Mass ~  AvTemp, dat, cphy, dat$AvgLat, dat$AvgLon)

