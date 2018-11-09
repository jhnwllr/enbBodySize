
# Manuscript on environmental determinates of body size in Odonates...
# envBodySize.r 

if(FALSE) { # Fig. 1B. # Temp vs Tropical Phylo Corrected

load(file = "C:/Users/Owner/Desktop/envBodySize/data/opdb/quickBodySize.rda") # load cleaned body size data
B = quickBodySize # B for body size
# B = B[c("GenusSpecies","Species","Genus","Family","SubOrder","tbl","hwl")] # drop fwl since often NA

B_hwl = B[c("GenusSpecies","Species","Genus","Family","SubOrder","hwl")] # drop another variable so that is easy to work with 
B_hwl = na.omit(B_hwl) # remove missing 

B_tbl = B[c("GenusSpecies","Species","Genus","Family","SubOrder","tbl")] # drop another variable so that is easy to work with 
B_tbl = na.omit(B_tbl) # remove missing 


load(file = "C:/Users/Owner/Desktop/envBodySize/data/gbif/climateLat.rda") # load un-cleaned climate data
C = climateLat # C for climate
C$GenusSpecies = gsub(" ","_",C$GenusSpecies) 

AL = aggregate(C$decimalLatitude ~ GenusSpecies, data = C, median) # average latitude
colnames(AL) = c("GenusSpecies","lat")
AL$climate = ifelse(AL$lat < 23 & AL$lat > -23,"trop","temp") # define tropical climates

D_hwl = merge(B_hwl,AL,id="GenusSpecies")
D_tbl = merge(B_tbl,AL,id="GenusSpecies")

library(ape) 
load("C:/Users/Owner/Desktop/envBodySize/data/trees/tree.rda") # load tree for pgls


# D_hwl = D_hwl[D_hwl$SubOrder == "Anisoptera",]
# D_tbl = D_tbl[D_tbl$SubOrder == "Anisoptera",]

treeD = function(tree,D_hwl,D_tbl) {  # function to process 

D_hwl = D_hwl[D_hwl$GenusSpecies %in% tree$tip.label,] # get only those with tip
D_tbl = D_tbl[D_tbl$GenusSpecies %in% tree$tip.label,] # get only those with tip

tip_hwl = tree$tip.label[!tree$tip.label %in% D_hwl$GenusSpecies] # drop tips without data
tip_tbl = tree$tip.label[!tree$tip.label %in% D_tbl$GenusSpecies] # drop tips without data

tree_hwl = drop.tip(tree, tip_hwl)
tree_tbl = drop.tip(tree, tip_tbl)

library(nlme)
rownames(D_hwl) = D_hwl$GenusSpecies
V_hwl = corBrownian(1,phy=tree_hwl) 
fit_hwl = gls(hwl ~ climate,correlation=V_hwl,data=D_hwl)

rownames(D_tbl) = D_tbl$GenusSpecies
V_tbl = corBrownian(1,phy=tree_tbl) 
fit_tbl = gls(tbl ~ climate,correlation=V_tbl,data=D_tbl)

tbl_se = summary(fit_tbl)$tTable[2,2]
hwl_se = summary(fit_hwl)$tTable[2,2]

pvalue_tbl = summary(fit_tbl)$tTable[2,4]
pvalue_hwl = summary(fit_hwl)$tTable[2,4]

print(summary(fit_tbl))
print(summary(fit_hwl))

mean_temp_tbl = coef(fit_tbl)[1]
mean_trop_tbl = coef(fit_tbl)[1] + coef(fit_tbl)[2]
tbl_means = c(mean_trop_tbl,mean_temp_tbl)

mean_temp_hwl = coef(fit_hwl)[1]
mean_trop_hwl = coef(fit_hwl)[1] + coef(fit_hwl)[2]
hwl_means = c(mean_trop_hwl,mean_temp_hwl)

# hwl_means = c(mean(D_hwl[D_hwl$climate == "trop",]$hwl),mean(D_hwl[D_hwl$climate == "temp",]$hwl))
# names(hwl_means) = c("trop","temp")
# tbl_means = c(mean(D_tbl[D_tbl$climate == "trop",]$tbl),mean(D_tbl[D_tbl$climate == "temp",]$tbl))
# names(tbl_means) = c("trop","temp")

# lower_tbl=2*tbl_se; upper_tbl=2*tbl_se # 95% ci
# lower_hwl=2*hwl_se; upper_hwl=2*hwl_se # 95% ci

climate = c("trop","temp","trop","temp")
morph = c("hwl","hwl","tbl","tbl")
means = c(hwl_means, tbl_means)

se = rep(c(hwl_se, tbl_se),each=2) # 95% ci

pvalue = rep(c(pvalue_hwl, pvalue_tbl),each=2)
 
return(
data.frame(
climate = climate,
morph = morph,
means = means,
se = se,
pvalue = pvalue
)
)

}



# Prepare data for ggplot

# Odonata 
D_hwl_in = D_hwl; D_tbl_in = D_tbl
d1 = treeD(tree,D_hwl_in,D_tbl_in) # process data
print("Odonata")
d1$group = "Odonata (both suborders)"


# Zygoptera 
D_hwl_in = D_hwl[D_hwl$SubOrder == "Zygoptera",]; D_tbl_in = D_tbl[D_tbl$SubOrder == "Zygoptera",]
d2 = treeD(tree,D_hwl_in,D_tbl_in) # process data
print("Zygoptera")
d2$group = "Zygoptera (damselflies)"

D_hwl_in = D_hwl[D_hwl$SubOrder == "Anisoptera",]; D_tbl_in = D_tbl[D_tbl$SubOrder == "Anisoptera",]
d3 = treeD(tree,D_hwl_in,D_tbl_in) # process data
print("Anisoptera")
d3$group = "Anisoptera (dragonflies)"

D = rbind(d1,d2,d3) # final plot data 
D$pvalueRES = round(D$pvalue,3)
D
D$group = factor(D$group, levels = c("Odonata (both suborders)", "Zygoptera (damselflies)", "Anisoptera (dragonflies)"))

group.colors = c(trop = "#e13d14", temp = "#0190b6")

library(ggplot2)
p = ggplot(D,aes(x=morph,y=means,fill=factor(climate))) +
geom_bar(stat="identity",position="dodge") + 
facet_wrap(~group) +
geom_errorbar(aes(ymin=means-se, ymax=means+se),width=.2,position=position_dodge(.9)) + 
theme_bw() + 
scale_fill_manual(values=group.colors) + 
guides(fill=guide_legend(title="Latitude")) + 
theme(strip.background = element_rect(fill="#f1e4d7")) +  
theme(legend.position=c(0.93,0.12)) + 
ylim(0,62)


ggsave("C:/Users/Owner/Desktop/plot.pdf")  

} 

if(FALSE) { # Fig. 1A. # Temp vs Tropical Raw Means 

load(file = "C:/Users/Owner/Desktop/envBodySize/data/opdb/quickBodySize.rda") # load cleaned body size data
B = quickBodySize # B for body size
# B = B[c("GenusSpecies","Species","Genus","Family","SubOrder","tbl","hwl")] # drop fwl since often NA

B_hwl = B[c("GenusSpecies","Species","Genus","Family","SubOrder","hwl")] # drop another variable so that is easy to work with 
B_hwl = na.omit(B_hwl) # remove missing 

B_tbl = B[c("GenusSpecies","Species","Genus","Family","SubOrder","tbl")] # drop another variable so that is easy to work with 
B_tbl = na.omit(B_tbl) # remove missing 


load(file = "C:/Users/Owner/Desktop/envBodySize/data/gbif/climateLat.rda") # load un-cleaned climate data
C = climateLat # C for climate
C$GenusSpecies = gsub(" ","_",C$GenusSpecies) 

AL = aggregate(C$decimalLatitude ~ GenusSpecies, data = C, median) # average latitude
colnames(AL) = c("GenusSpecies","lat")
AL$climate = ifelse(AL$lat < 23 & AL$lat > -23,"trop","temp") # define tropical climates

D_hwl = merge(B_hwl,AL,id="GenusSpecies")
D_tbl = merge(B_tbl,AL,id="GenusSpecies")

library(ape) 
load("C:/Users/Owner/Desktop/envBodySize/data/trees/tree.rda") # load tree for pgls

getD = function(D_hwl,D_tbl) { 

fit_hwl = t.test(hwl ~ climate,data=D_hwl)
fit_tbl = t.test(tbl ~ climate,data=D_tbl)


fitD_hwl = broom::tidy(fit_hwl) # fit data
fitD_tbl = broom::tidy(fit_tbl) # fit data

pvalue_tbl = fitD_tbl$p.value
pvalue_hwl = fitD_hwl$p.value

print(pvalue_tbl)
print(pvalue_hwl)

mean_temp_tbl = fitD_tbl$estimate1 
mean_trop_tbl = fitD_tbl$estimate2
tbl_means = c(mean_trop_tbl,mean_temp_tbl)

mean_temp_hwl = fitD_hwl$estimate1
mean_trop_hwl = fitD_hwl$estimate2
hwl_means = c(mean_trop_hwl,mean_temp_hwl)

hwl_se = sd(D_hwl$hwl)/sqrt(length(D_hwl$hwl))
tbl_se = sd(D_tbl$tbl)/sqrt(length(D_tbl$tbl))

se = rep(c(hwl_se, tbl_se),each=2) # 95% ci

climate = c("trop","temp","trop","temp")
morph = c("hwl","hwl","tbl","tbl")
means = c(hwl_means, tbl_means)

pvalue = rep(c(pvalue_hwl, pvalue_tbl),each=2)
 
return(
data.frame(
climate = climate,
morph = morph,
means = means,
se = se,
pvalue = pvalue
)
)

}

# Prepare data for ggplot

# Odonata 
D_hwl_in = D_hwl; D_tbl_in = D_tbl
d1 = getD(D_hwl_in,D_tbl_in) # process data
print("Odonata")
d1$group = "Odonata (both suborders)"

# Zygoptera 
D_hwl_in = D_hwl[D_hwl$SubOrder == "Zygoptera",]; D_tbl_in = D_tbl[D_tbl$SubOrder == "Zygoptera",]
d2 = getD(D_hwl_in,D_tbl_in) # process data
print("Zygoptera")
d2$group = "Zygoptera (damselflies)"

D_hwl_in = D_hwl[D_hwl$SubOrder == "Anisoptera",]; D_tbl_in = D_tbl[D_tbl$SubOrder == "Anisoptera",]
d3 = getD(D_hwl_in,D_tbl_in) # process data
print("Anisoptera")
d3$group = "Anisoptera (dragonflies)"


D = rbind(d1,d2,d3) # final plot data 
D$pvalueRES = round(D$pvalue,3)
D$group = factor(D$group, levels = c("Odonata (both suborders)", "Zygoptera (damselflies)", "Anisoptera (dragonflies)"))

group.colors = c(trop = "#e13d14", temp = "#0190b6")

library(ggplot2)
p = ggplot(D,aes(x=morph,y=means,fill=factor(climate))) +
geom_bar(stat="identity",position="dodge") + 
facet_wrap(~group) +
geom_errorbar(aes(ymin=means-se*2, ymax=means+se*2),width=.2,position=position_dodge(.9)) + 
theme_bw() + 
scale_fill_manual(values=group.colors) + 
guides(fill=guide_legend(title="Latitude")) + 
theme(strip.background = element_rect(fill="#f1e4d7")) +  
theme(legend.position=c(0.93,0.12)) + 
ylim(0,62)


ggsave("C:/Users/Owner/Desktop/plot.pdf")  

} 



if(FALSE) { # old 
pdf("C:/Users/Owner/Desktop/Fig1_TempVsTropical.pdf",width = 8.7, height = 8.7, useDingbats=FALSE)

par(mfrow = c(1,6))
lower_tbl=out_Z$lower_tbl; upper_tbl=out_Z$upper_tbl; tbl_means=out_Z$tbl_means
x = barplot(tbl_means,col=c("red","gray"),main="Odonates\ntotal body length",ylab="length (mm)",ylim=c(0,60))
arrows(x,tbl_means+upper_tbl,x,tbl_means-lower_tbl,angle=90,code=3,length=0.05)
lower_hwl=out_Z$lower_hwl; upper_hwl=out_Z$upper_hwl; hwl_means=out_Z$hwl_means
x = barplot(hwl_means,col=c("red","gray"),main="Odonates\nhind wing length",ylab="length (mm)",ylim=c(0,60))
arrows(x,hwl_means+upper_hwl,x,hwl_means-lower_hwl,angle=90,code=3,length=0.05)

#### END ODONATES 

D_hwl_in = D_hwl[D_hwl$SubOrder == "Zygoptera",]; D_tbl_in = D_tbl[D_tbl$SubOrder == "Zygoptera",]

out_Z = treeD(tree,D_hwl_in,D_tbl_in) # function to process 
print(out_Z)

lower_tbl=out_Z$lower_tbl; upper_tbl=out_Z$upper_tbl; tbl_means=out_Z$tbl_means
x = barplot(tbl_means,col=c("red","gray"),main="Damselflies\ntotal body length",ylab="length (mm)",ylim=c(0,60))
arrows(x,tbl_means+upper_tbl,x,tbl_means-lower_tbl,angle=90,code=3,length=0.05)
lower_hwl=out_Z$lower_hwl; upper_hwl=out_Z$upper_hwl; hwl_means=out_Z$hwl_means
x = barplot(hwl_means,col=c("red","gray"),main="Damselflies\nhind wing length",ylab="length (mm)",ylim=c(0,60))
arrows(x,hwl_means+upper_hwl,x,hwl_means-lower_hwl,angle=90,code=3,length=0.05)

D_hwl_in = D_hwl[D_hwl$SubOrder == "Anisoptera",]; D_tbl_in = D_tbl[D_tbl$SubOrder == "Anisoptera",]

out_Z = treeD(tree,D_hwl_in,D_tbl_in) # function to process 

lower_tbl=out_Z$lower_tbl; upper_tbl=out_Z$upper_tbl; tbl_means=out_Z$tbl_means
x = barplot(tbl_means,col=c("red","gray"),main="Dragonflies\ntotal body length",ylab="length (mm)",ylim=c(0,60))
arrows(x,tbl_means+upper_tbl,x,tbl_means-lower_tbl,angle=90,code=3,length=0.05)
lower_hwl=out_Z$lower_hwl; upper_hwl=out_Z$upper_hwl; hwl_means=out_Z$hwl_means
x = barplot(hwl_means,col=c("red","gray"),main="Dragonflies\nhind wing length",ylab="length (mm)",ylim=c(0,60))
arrows(x,hwl_means+upper_hwl,x,hwl_means-lower_hwl,angle=90,code=3,length=0.05)

dev.off()

}

if(FALSE) { # stats for fossils latitude trends
# load present day body sizes
load(file = "C:/Users/Owner/Desktop/envBodySize/data/opdb/quickBodySize.rda") # load cleaned body size data
B = quickBodySize # B for body size
B = B[c("GenusSpecies","Species","Genus","Family","SubOrder","tbl","hwl")] # drop fwl since often NA
B = B[c("GenusSpecies","Species","Genus","Family","SubOrder","hwl")] # drop another variable so that is easy to work with 
B = na.omit(B) # remove missing 

load(file = "C:/Users/Owner/Desktop/envBodySize/data/gbif/climateLat.rda") # load un-cleaned climate data
C = climateLat # C for climate
C$GenusSpecies = gsub(" ","_",C$GenusSpecies) 

AL = aggregate(decimalLatitude ~ GenusSpecies, data = C, mean) # average latitude

D = merge(B,AL,id="GenusSpecies") # final data.frame

# Run Some stats suggested by Erik

# make the data.frames match up for plotting
D$age_ma = 0
D$lat = D$decimalLatitude
D = D[c("GenusSpecies","Genus","Family","SubOrder","hwl","age_ma","lat")]
colnames(D) = c("GenusSpecies","Genus","Family","SubOrder","size_mm","age_ma","lat")

# D = D[dplyr::between(D$lat, 25, 50),] # get only in the range present in the fossil record 

load(file="C:/Users/Owner/Desktop/envBodySize/data/Fossils/FossilsLatTax.rda") # load fossil data
DF = FossilsLatTax # Data Fossils

DF$lat = DF$paleolat
DF = DF[c("GenusSpecies","Genus","Family","SubOrder","size_mm","age_ma","lat")]

D$type = "present"
DF$type = "fossil"

D = rbind(D,DF) # combine present and fossil datasets

D$age_cat = D$age_ma
D$age_cat = "80-0"
D$age_cat[D$age_ma > 80] = "150-80"
D$age_cat[D$age_ma > 150] = "210-150"
D$age_cat[D$type == "present"] = "0 extant"

D$age_cat = factor(D$age_cat,c("210-150","150-80","80-0","0 extant"))

# Stats for Fig 2. 
# str(D)
# D=D[D$SubOrder == "Zygoptera",]
D=D[D$SubOrder == "Anisoptera",]
D = D[D$age_cat == "0 extant",]
# D = D[D$age_cat == "150-80",]
# D = D[D$age_cat == "80-0",]
# D = D[D$age_cat == "210-150",]

str(D)
D$lat

fit=lm(size_mm~abs(lat),data=D)
summary(fit)

# [1] 0 extant 150-80   210-150  80-0
# fit=lm(size_mm~age_cat*abs(lat),data=D)
# summary(fit)
# summary(aov(fit))
}

if(FALSE) { # Plot present day with paleolat
library(ggplot2)
ggplot(D, aes(abs(lat),size_mm)) + 
geom_point(aes(colour=age_ma), size=2) + 
scale_colour_gradient(low="antiquewhite",high="tan4") + 
geom_point(shape = 1,size = 2,colour = "black") + 
geom_smooth(method="gam",se = FALSE,colour="gray31") + 
ggtitle("Fossils: wing length and latitude") + 
facet_grid(SubOrder~age_cat,scales="free_x") + 
xlab("Paleo-latitude") + 
theme_bw() + 
theme(legend.position=c(0.06, 0.25)) + 
ylab("Wing Length (mm)") 

ggsave("C:/Users/Owner/Desktop/PaleoLatPresent.pdf",width=10,height=5)

# fit = lme(size_mm ~ paleolat,random=~1|SubOrder/Family,data=D)
# summary(fit)

}

if(FALSE) { # Fig. 2. # polygons and maps

load(file = "C:/Users/Owner/Desktop/envBodySize/data/opdb/quickBodySize.rda") # load cleaned body size data
B = quickBodySize # B for body size
B = B[c("GenusSpecies","Species","Genus","Family","SubOrder","hwl")] # drop fwl since often NA
B = na.omit(B) # remove missing 

# load(file = "C:/Users/Owner/Desktop/envBodySize/data/gbif/climateLat.rda") # load un-cleaned climate data
# C = climateLat # C for climate
# C$GenusSpecies = gsub(" ","_",C$GenusSpecies) 
# C = C[c("GenusSpecies","decimalLatitude", "decimalLongitude","bio_max1","bio_max12")]
# colnames(C) = c("GenusSpecies","Lat","Lon","Temp","Prec")
# C$Temp = C$Temp/10
# D = merge(B,C,id="GenusSpecies",all.x=TRUE)
# D = na.omit(D)

load(file = "C:/Users/Owner/Desktop/envBodySize/data/gbif/climateLat.rda") # load un-cleaned climate data
C = climateLat # C for climate
C$GenusSpecies = gsub(" ","_",C$GenusSpecies) 
C = C[c("GenusSpecies","decimalLatitude", "decimalLongitude","bio_max1","bio_max12")]
colnames(C) = c("GenusSpecies","Lat","Lon","Temp","Prec")
C$Temp = C$Temp/10
C$Prec = C$Prec/10
D = merge(B,C,id="GenusSpecies",all.x=TRUE)
D = na.omit(D)


# load(file="C:/Users/Owner/Desktop/envBodySize/data/gbif/Birds.rda") # load Birds
# Brd = Birds # rename B for Birds
# Brd$GenusSpecies = gsub(" ","_",Brd$GenusSpecies) 
# str(Brd)
# Brd$diversity = Brd$all # select Bird variable to analyze 
# Brd = Brd[c("GenusSpecies", "Lat", "Lon", "diversity")]
# Brd = na.omit(Brd)
# D = merge(B,Brd,id="GenusSpecies")

# load(file="C:/Users/Owner/Desktop/envBodySize/data/gbif/Mammals.rda") # 
# M = Mammals # rename M for Mammals
# M$GenusSpecies = gsub(" ","_",M$GenusSpecies)
# M$mam_div = M$all_spp
# M = M[c("GenusSpecies", "Lat", "Lon", "mam_div")]
# M = na.omit(M)
# D = merge(B,M,id="GenusSpecies")

# load("C:/Users/Owner/Desktop/envBodySize/data/gbif/trees.rda")
# T = trees # T for Tree Cover
# T = T[!T$tree == 0,] # remove zero value
# T$GenusSpecies = gsub(" ","_",T$GenusSpecies)
# T$cover = T$trees # select Bird variable to analyze 
# T = T[c("GenusSpecies", "Lat", "Lon", "cover")]
# T = na.omit(T)
# D = merge(B,T,id="GenusSpecies")

# focalTrait = "diversity"
# focalTrait = "cover"
# focalTrait = "mam_div"
focalTrait = "Prec"


# D = D[D$SubOrder == "Zygoptera",] # only damselflies
# D = D[D$SubOrder == "Anisoptera",] # only damselflies


mapHexagons = function(D,crclrs=55) { # hexagons on map function
# data.frame with Lat, Lon, focalTrait, GenusSpecies
# start function 
D$focalTrait = D[,focalTrait]
D = D[c("GenusSpecies","Lat","Lon","focalTrait")]


ss = 5 # step size
binLt = seq(-90,90,ss)
binLn = seq(-180,180,ss)

# Create Bins
xLt = D$Lat
DL = list()
count = 1
for(i in 1:(length(binLt)-1)) {
d = D[xLt > binLt[i] & xLt < binLt[i+1],] # get latitude band

for(j in 1:(length(binLn)-1)) {
xLn = d$Lon

dd = d[xLn > binLn[j] & xLn < binLn[j+1],]
if(nrow(dd) != 0) {
dd$lat_point = binLt[i]
dd$lon_point = binLn[j]
}
DL[[count]] = dd

count = count + 1
}

}

dl = list() # output datalist
for(i in 1:length(DL)) {
d = DL[[i]] 

if(nrow(d) != 0) { 
y = d$lat_point[1]
x = d$lon_point[1]
n = length(d$GenusSpecies)

ft = mean(d[!duplicated(d$GenusSpecies),]$focalTrait,na.rm=TRUE)

} else {
x = NA
y = NA
n = NA
ft = NA
}

# dl[[i]] = data.frame(n,temp,hwl,x,y)
dl[[i]] = data.frame(n,ft,x,y)
}

D = do.call("rbind",dl)


# hexagon plot

library(RColorBrewer) #to use brewer.pal
library(fields) #to use designer.colors

Hexagon <- function (x, y, unitcell = 1, col = col) {
polygon(c(x, x, x + unitcell/2, x + unitcell, x + unitcell, 
x + unitcell/2), c(y + unitcell * 0.125, 
y + unitcell * 0.875, 
y + unitcell * 1.125, 
y + unitcell * 0.875, 
y + unitcell * 0.125, 
y - unitcell * 0.125), 
col = col, border=NA)
}#function

# Heatmap_Matrix = matrix(rnorm(1000),nrow=100)
# hwl = D$hwl
# hwl[is.na(hwl)] = 0 # replace with 0
# Heatmap_Matrix = t(matrix(hwl,ncol=(length(binLt)-1)))

ft = D$ft
ft[is.na(ft)] = 0 # replace with 0
Heatmap_Matrix = t(matrix(ft,ncol=(length(binLt)-1)))

x = as.vector(Heatmap_Matrix)

makeScale = function(x,crclrs,xlab) { # make the scale in separate pdf
colfunc = colorRampPalette(c("red","orange","yellow","green","lightsteelblue","blue","violet"))
ColRamp = rev(colfunc(crclrs))

ColorCode = rep("#FFFFFF", length(x)) #default is all white
Bins = seq(min(x, na.rm=T), max(x, na.rm=T), length=length(ColRamp))

pdf("C:/Users/Owner/Desktop/scale.pdf", height = 5, width = 5,useDingbats=FALSE)
plot(Bins,rep(0,length(Bins)),col=ColRamp,pch=19,xlab=xlab,bty="n",xaxt="n",yaxt="n",ylab="",ylim=c(0,1),cex=1.5)
axis(1, at = seq(min(x), max(x), by = 50), las=1)
dev.off()
}

makeScale(x,crclrs,xlab="num. sp")

SOM_Rows = dim(Heatmap_Matrix)[1]
SOM_Columns = dim(Heatmap_Matrix)[2]

# jpeg("C:/Users/Owner/Desktop/plot.jpg",units = "in",width=15,height=10,res=300)
pdf("C:/Users/Owner/Desktop/polygons.pdf",width=15,height=10)
plot(0, 0, type = "n", axes = FALSE, xlim=c(0, SOM_Columns), ylim=c(0, SOM_Rows), xlab="", ylab= "", asp=1)

# ColRamp <- rev(designer.colors(n=100, col=brewer.pal(9, "Spectral")))

colfunc = colorRampPalette(c("red","orange","yellow","green","lightsteelblue","blue","violet"))
# 55
ColRamp = rev(colfunc(crclrs))

ColorCode = rep("#FFFFFF", length(x)) #default is all white
Bins = seq(min(x, na.rm=T), max(x, na.rm=T), length=length(ColRamp))
for (i in 1:length(x)) if (!is.na(x[i])) ColorCode[i] = ColRamp[which.min(abs(Bins-x[i]))] 

background = names(rev(sort(table(ColorCode)))[1]) 
ColorCode[ColorCode == background] = "#00000000"

offset = 0.5 #offset for the hexagons when moving up a row
for (row in 1:SOM_Rows) {
for (column in 0:(SOM_Columns - 1)) 
Hexagon(column + offset, row - 1, col = ColorCode[row + SOM_Rows * column])
offset = ifelse(offset, 0, 0.5)
}

# image(Heatmap_Matrix)
dev.off()
}

mapHexagons(D)

# make the map that goes underneath
# pdf("C:/Users/Owner/Desktop/map.pdf",width=15,height=10)
# library(maps)
# map(col="gray")
# dev.off()

}

if(FALSE) { # Fig. 4. # plot shifts and climate on tree shifts on tree
load(file="C:/Users/Owner/Desktop/envBodySize/data/optima/optima.rda")
E = optima
E = E[E$tn == "tip",] # filter edge matrix

examples = c("Calopteryx_exul","Ictinogomphus_decoratus","Aeshna_eremita","Macromia_taeniolata","Megaloprepus_caerulatus","Mecistogaster_linearis","Ischnura_aurora") 
# examples = c("Calopteryx_exul","Ictinogomphus_decoratus","Aeshna_eremita","Macromia_taeniolata","Ischnura_aurora") 

example_optima = c(); for(i in 1:length(examples)) example_optima[i] = E[E$label == examples[i],]$optima_tip

E$shifts = ifelse(E$optima_tip %in% example_optima,"pos","neg")
E = E[c("label","shifts")]
colnames(E) = c("GenusSpecies","shifts") 

load(file = "C:/Users/Owner/Desktop/envBodySize/data/gbif/climateLat.rda") # load un-cleaned climate data
C = climateLat # C for climate
C$GenusSpecies = gsub(" ","_",C$GenusSpecies) 

AL = aggregate(C$decimalLatitude ~ GenusSpecies, data = C, mean) # average latitude
colnames(AL) = c("GenusSpecies","lat")

AL$climate = ifelse(AL$lat < 23 & AL$lat > -23,"trop","temp") # define tropical climates

D = merge(E,AL,id="GenusSpecies")

load(file = "C:/Users/Owner/Desktop/envBodySize/data/opdb/quickBodySize.rda")
BS = quickBodySize # body size 
BS = BS[c("GenusSpecies","tbl")]
D = merge(D,BS,id="GenusSpecies")

library(ape) 
load("C:/Users/Owner/Desktop/envBodySize/data/trees/tree.rda") # load tree for pgls
D = D[D$GenusSpecies %in% tree$tip.label,] # get only those with tip

tip = tree$tip.label[!tree$tip.label %in% D$GenusSpecies] # drop tips without data
tree = drop.tip(tree, tip)

shifts = setNames(D$shifts,D$GenusSpecies)
climate = setNames(D$climate,D$GenusSpecies)

library(phytools)

tree_map_shifts = make.simmap(tree,x=shifts,model="SYM",nsim=1)
cols_shifts = setNames(c("gray","#ff7373"),c("neg","pos"))
bar_cols_shifts = setNames(ifelse(D$shifts == "pos","#ff7373","gray"),D$GenusSpecies)
x_shifts = D$tbl
names(x_shifts) = D$GenusSpecies

tree_map_climate = make.simmap(tree,x=climate,model="SYM",nsim=1)
cols_climate = setNames(c("gray","#ff7373"),c("trop","temp"))
bar_cols_climate = setNames(ifelse(D$climate == "pos","#ff7373","gray"),D$GenusSpecies)
x_climate = rep(1,length(climate))
names(x_climate) = D$GenusSpecies


# jpeg("C:/Users/Owner/Desktop/plot.jpg")
pdf("C:/Users/Owner/Desktop/Fig4_mirrorTrees.pdf",height=10,width=10,useDingbats=FALSE)
par(mfrow=c(1,2))
plotTree.wBars(tree_map_shifts, x_shifts, scale=1, width=NULL, type="phylogram"
,method="plotSimmap", tip.labels=TRUE, col= bar_cols_shifts, border=NA,fsize=0.1,colors=cols_shifts)

plotTree.wBars(tree_map_climate, x_climate, scale=1, width=NULL, type="phylogram",method="plotSimmap", tip.labels=TRUE, col= bar_cols_climate, border=NA,fsize=0.1,colors=cols_climate,direction="leftwards")

dev.off()

# jpeg("C:/Users/Owner/Desktop/plot.jpg")
# plot(tree_map,cols,fsize=0.1)
# dev.off()

# library(phytools)
# fit = fitPagel(tree,x,y)
# round(fit$independent.Q,3)
# unclass(fit)
}

if(FALSE) { # Fig. 3 Path Analysis
library(ape)
load(file="C:/Users/Owner/Desktop/BodySize/Manuscript/Data/trees/tree.rda") # load supertree

load(file="C:/Users/Owner/Desktop/BodySize/Manuscript/Data/Development/Development.rda") # load development time
D = Development # rename

im = function(x) replace(x, is.na(x), mean(x, na.rm = TRUE))
D = plyr::ddply(D, ~ Genus, transform, x = im(body_lengths)) # replace missing taxa with mean of Genus
D = na.omit(D[c("GenusSpecies","dev","MeanLat","x")])
D$GenusSpecies = gsub(" ","_",D$GenusSpecies)

tip = tree$tip.label[!tree$tip.label %in% D$GenusSpecies] # drop tips without data
tree = drop.tip(tree, tip)

D = D[D$GenusSpecies %in% tree$tip.label,] # remove data without tip in tree

load(file = "C:/Users/Owner/Desktop/envBodySize/data/gbif/climateLat.rda") # load un-cleaned climate data
C = climateLat # C for climate
C$GenusSpecies = gsub(" ","_",C$GenusSpecies) 
AC = aggregate(cbind(bio_max1,bio_max12) ~ GenusSpecies, data = C, mean) # average 

D = merge(D,AC,id="GenusSpecies",all.x=TRUE) # 
colnames(D) = c("SP","dev","MeanLat","tbl","temp","prec")
rownames(D) = D$SP

library(phylopath)

models = list(
mod1 = DAG(tbl ~ dev,dev ~ temp),
mod2 = DAG(tbl ~ dev + temp, dev ~ temp),
mod3 = DAG(tbl ~ temp,dev ~ temp),
mod4 = DAG(tbl ~ dev, temp ~ temp),
mod5 = DAG(tbl ~ temp, dev ~ dev)
)

fit = phylo_path(models, data = D, tree = tree)
summary(fit)

# methods(class = "phylopath")

best_model = best(fit)
# best_model

# plot_model_set(models, labels = NULL, algorithm = "kk", text_size = 5,
# box_x = 12, box_y = 10, edge_width = 1, curvature = 0.05,
# rotation = 0, flip_x = FALSE, flip_y = FALSE, nrow = NULL,
# arrow = grid::arrow(type = "closed", 15, grid::unit(10, "points")))

# 12.628, 13.474, 21.902, 34.721, 43.994

# best_model
# jpeg("C:/Users/Owner/Desktop/plot.jpg")
pdf("C:/Users/Owner/Desktop/Fig3_coef.pdf",height=10,width=10,useDingbats=FALSE)
coef_plot(average(fit), reverse_order = FALSE)
# coef_plot(best_model, reverse_order = FALSE)
# plot_model_set(models,curvature=0)
# plot(best_model)
dev.off()
}


if(FALSE) { # Table 1 # Spatial Data, Mammals, Birds, Tree Cover, Temp, Prec

# list.files("C:/Users/Owner/Desktop/envBodySize/data/")
load("C:/Users/Owner/Desktop/envBodySize/data/gbif/trees.rda")
T = trees # T for Tree Cover
# T = T[!T$tree == 0,] # remove zero value
T$GenusSpecies = gsub(" ","_",T$GenusSpecies)
AT = aggregate(trees ~ GenusSpecies, data = T, mean) # AT average tree cover
colnames(AT) = c("GenusSpecies","cover")

load(file = "C:/Users/Owner/Desktop/envBodySize/data/gbif/climateLat.rda") # load un-cleaned climate data
C = climateLat # C for climate
C$GenusSpecies = gsub(" ","_",C$GenusSpecies) 

AC = aggregate(cbind(bio_max1,bio_max12) ~ GenusSpecies, data = C, mean) # average climate
colnames(AC) = c("GenusSpecies","Temp","Prec")

load(file = "C:/Users/Owner/Desktop/envBodySize/data/opdb/quickBodySize.rda") # load cleaned body size data
BS = quickBodySize # BS for body size DATA
BS = BS[c("GenusSpecies","Species","Genus","Family","SubOrder","tbl","hwl")] # drop fwl since often NA
BS = BS[c("GenusSpecies","Species","Genus","Family","SubOrder","hwl")] # drop another variable so that is easy to work with 
BS = na.omit(BS) # remove missing 

# load biotic grids 
load(file="C:/Users/Owner/Desktop/envBodySize/data/gbif/Mammals.rda") # 
A = Mammals # rename A for Amphibians

load(file="C:/Users/Owner/Desktop/envBodySize/data/gbif/Birds.rda") # load Birds
B = Birds # rename B for Birds

# str(B) # what variable are available
B$div_brd = B$all # select Amphibians variable to analyze 
A$div_mam = A$all_spp # select Bird variable to analyze 

AB = aggregate(div_brd ~ GenusSpecies,data=B,mean) # get average birds diversity 
AB$GenusSpecies = gsub(" ","_",AB$GenusSpecies)
D = merge(BS,AB,id="GenusSpecies")
AA = aggregate(div_mam ~ GenusSpecies,data=A,mean) # get average birds diversity 
AA$GenusSpecies = gsub(" ","_",AA$GenusSpecies)

D = merge(D,AA,id="GenusSpecies")
D = merge(D,AC,id="GenusSpecies") # merge diversity with average climates
D = merge(D,AT,id="GenusSpecies") # merge diversity with average climates


library(ape); load("C:/Users/Owner/Desktop/envBodySize/data/trees/tree.rda") # load tree and ape for pgls

# Scale Data
D$Temp = scale(D$Temp)
D$Prec = scale(D$Prec)
D$div_mam = scale(D$div_mam)
D$div_brd = scale(D$div_brd)
D$cover = scale(D$cover)

# naive analysis
fit1 = lm(hwl ~ Temp, data=D)
fit2 = lm(hwl ~ div_mam, data=D)
fit3 = lm(hwl ~ div_brd, data=D)
fit4 = lm(hwl ~ div_brd + div_mam, data=D)
fit5 = lm(hwl ~ div_brd + div_mam + Temp, data=D)
fit6 = lm(hwl ~ div_brd + div_mam + Temp + Prec, data=D)
fit7 = lm(hwl ~ div_brd + div_mam + Temp*Prec, data=D)
fit8 = lm(hwl ~ div_brd*div_mam + Temp*Prec, data=D)
fit9 = lm(hwl ~ div_brd*div_mam + Temp*Prec, data=D)
fit10 = lm(hwl ~ Temp*Prec, data=D)
fit11 = lm(hwl ~ Temp + Prec, data=D)
fit12 = lm(hwl ~ Prec, data=D)

MS = AIC(fit1,fit2,fit3,fit4,fit5,fit6,fit7,fit8,fit9,fit10,fit11,fit12)
# MS[order(MS$AIC),]


# mixed effect model 
library(nlme)
fit1 = lme(hwl ~ div_brd, random = ~ 1|SubOrder/Family,data=D)
fit2 = lme(hwl ~ div_mam, random = ~ 1|SubOrder/Family, data=D)
fit3 = lme(hwl ~ div_brd, random = ~ 1|SubOrder/Family, data=D)
fit4 = lme(hwl ~ div_brd + div_mam, random = ~ 1|SubOrder/Family, data=D)
fit5 = lme(hwl ~ div_brd + div_mam + Temp, random = ~ 1|SubOrder/Family, data=D)
fit6 = lme(hwl ~ div_brd + div_mam + Temp + Prec, random = ~ 1|SubOrder/Family, data=D)
fit7 = lme(hwl ~ div_brd + div_mam + Temp*Prec, random = ~ 1|SubOrder/Family, data=D)
fit8 = lme(hwl ~ div_brd*div_mam + Temp*Prec, random = ~ 1|SubOrder/Family, data=D)
fit9 = lme(hwl ~ div_brd*div_mam + Temp*Prec, random = ~ 1|SubOrder/Family, data=D)
fit10 = lme(hwl ~ Temp*Prec, random = ~ 1|SubOrder/Family, data=D)
fit11 = lme(hwl ~ Temp + Prec, random = ~ 1|SubOrder/Family, data=D)
fit12 = lme(hwl ~ Prec, random = ~ 1|SubOrder/Family, data=D)

MS = AIC(fit1,fit2,fit3,fit4,fit5,fit6,fit7,fit8,fit9,fit10,fit11,fit12)
# MS[order(MS$AIC),]

# remove those not in tree or data
D = D[D$GenusSpecies %in% tree$tip.label,] # get only those with tip
# D = D[!duplicated(D$GenusSpecies),] # not necessary but lets leave it in 
tip = tree$tip.label[!tree$tip.label %in% D$GenusSpecies] # drop tips without data
tree = drop.tip(tree, tip) # drop missing

# nrow(D)
# tree

library(nlme)
V = corBrownian(1,phy=tree) 
rownames(D) = D$GenusSpecies
# V = corPagel(1, phy = tree, form = ~1, fixed = FALSE)
# fit = gls(hwl ~ div_brd + div_mam + Temp + Prec,correlation=V,data=D)
fit = gls(hwl ~ Temp + Prec + cover + div_brd + div_mam,correlation=V,data=D)
# fit = gls(hwl ~ Temp,correlation=V,data=D)


# clean gls summary info
AIC = summary(fit)$AIC
N = nrow(D)
tt = as.data.frame(summary(fit)$tTable)

p = tt$"p-value"
p = p[2:nrow(tt)] # drop intercept
p = round(p,4) # round
# if(p == 0) p = "< 0.0001"
p = paste("P =", p)

ce = rownames(tt)[2:nrow(tt)] # coefficients

val = round(tt$Value[2:nrow(tt)],3) # drop intercept
val = paste(ce,"=",val,collapse=" ;")

se = round(tt$Std.Error[2:nrow(tt)],3)
se = paste("SE =",se)

N
val
se
p
AIC

}


if(FALSE) {

# list.files("C:/Users/Owner/Desktop/envBodySize/data/")
load("C:/Users/Owner/Desktop/envBodySize/data/gbif/trees.rda")
T = trees # T for Tree Cover
# T = T[!T$tree == 0,] # remove zero value
T$GenusSpecies = gsub(" ","_",T$GenusSpecies)
AT = aggregate(trees ~ GenusSpecies, data = T, mean) # AT average tree cover
colnames(AT) = c("GenusSpecies","cover")

load(file = "C:/Users/Owner/Desktop/envBodySize/data/gbif/climateLat.rda") # load un-cleaned climate data
C = climateLat # C for climate
C$GenusSpecies = gsub(" ","_",C$GenusSpecies) 

AC = aggregate(cbind(decimalLatitude,decimalLongitude,bio_max1,bio_max12) ~ GenusSpecies, data = C, mean) # average climate
colnames(AC) = c("GenusSpecies","lat","lon","Temp","Prec")

load(file = "C:/Users/Owner/Desktop/envBodySize/data/opdb/quickBodySize.rda") # load cleaned body size data
BS = quickBodySize # BS for body size DATA
BS = BS[c("GenusSpecies","Species","Genus","Family","SubOrder","tbl","hwl")] # drop fwl since often NA
BS = BS[c("GenusSpecies","Species","Genus","Family","SubOrder","hwl")] # drop another variable so that is easy to work with 
BS = na.omit(BS) # remove missing 

# load biotic grids 
load(file="C:/Users/Owner/Desktop/envBodySize/data/gbif/Mammals.rda") # 
A = Mammals # rename A for Amphibians

load(file="C:/Users/Owner/Desktop/envBodySize/data/gbif/Birds.rda") # load Birds
B = Birds # rename B for Birds

# str(B) # what variable are available
B$div_brd = B$all # select Amphibians variable to analyze 
A$div_mam = A$all_spp # select Bird variable to analyze 

AB = aggregate(div_brd ~ GenusSpecies,data=B,mean) # get average birds diversity 
AB$GenusSpecies = gsub(" ","_",AB$GenusSpecies)
D = merge(BS,AB,id="GenusSpecies")
AA = aggregate(div_mam ~ GenusSpecies,data=A,mean) # get average birds diversity 
AA$GenusSpecies = gsub(" ","_",AA$GenusSpecies)

D = merge(D,AA,id="GenusSpecies")
D = merge(D,AC,id="GenusSpecies") # merge diversity with average climates
D = merge(D,AT,id="GenusSpecies") # merge diversity with average climates


library(ape); load("C:/Users/Owner/Desktop/envBodySize/data/trees/tree.rda") # load tree and ape for pgls

# Scale Data
D$Temp = scale(D$Temp)
D$Prec = scale(D$Prec)
D$div_mam = scale(D$div_mam)
D$div_brd = scale(D$div_brd)
D$cover = scale(D$cover)

source("C:/Users/Owner/Desktop/envBodySize/SpaceVsPhylogeny/brutal.r") # load library 

# remove those not in tree or data
D = D[D$GenusSpecies %in% tree$tip.label,] # get only those with tip
# D = D[!duplicated(D$GenusSpecies),] # not necessary but lets leave it in 
tip = tree$tip.label[!tree$tip.label %in% D$GenusSpecies] # drop tips without data
tree = drop.tip(tree, tip) # drop missing

V = vcv.phylo(tree)
rownames(D) = D$GenusSpecies

DM = dist.mat(D$lat, D$lon, rownames(D)) # distance matrix

# phi = 0.01
# fit1 = pglmSpatialFit(hwl ~ Temp, D, V, DM)

fit1 = pglmSpatial(hwl ~ Temp, D, V, DM,phi=0.01)
fit2 = pglmSpatial(hwl ~ Temp + Prec, D, V, DM,phi=0.01)
fit3 = pglmSpatial(hwl ~ Temp + Prec + cover, D, V, DM,phi=0.01)
fit4 = pglmSpatial(hwl ~ div_brd + Temp + Prec + cover, D, V, DM,phi=0.01)
fit5 = pglmSpatial(hwl ~ div_mam + div_brd + Temp + Prec + cover, D, V, DM,phi=0.01)
fit6 = pglmSpatial(hwl ~ div_mam + Temp + Prec + cover, D, V, DM,phi=0.01)

fitList = list(fit1,fit2,fit3,fit4,fit5,fit6)

save(fitList,file="C:/Users/Owner/Desktop/fitList.rda")

# fit1 = pglmSpatial(hwl ~ div_brd, D, V, DM,phi=0.01)
# fit2 = pglmSpatial(hwl ~ div_mam, D, V, DM,phi=0.01)
# fit3 = pglmSpatial(hwl ~ div_brd + div_mam, D, V, DM,phi=0.01)
# fit4 = pglmSpatial(hwl ~ div_brd + div_mam + Temp, D, V, DM,phi=0.01)
# fit5 = pglmSpatial(hwl ~ div_brd + div_mam + Temp + Prec, D, V, DM,phi=0.01)
# fit6 = pglmSpatial(hwl ~ div_brd + div_mam + Temp*Prec, D, V, DM,phi=0.01)
# fit7 = pglmSpatial(hwl ~ div_brd*div_mam + Temp*Prec, D, V, DM,phi=0.01)
# fit8 = pglmSpatial(hwl ~ div_brd*div_mam + Temp*Prec, D, V, DM,phi=0.01)
# fit9 = pglmSpatial(hwl ~ Temp*Prec, D, V, DM,phi=0.01)
# fit10 = pglmSpatial(hwl ~ Temp + Prec, D, V, DM,phi=0.01)
# fit11 = pglmSpatial(hwl ~ Prec, D, V, DM,phi=0.01)
# fit12 = pglmSpatial(hwl ~ Temp, D, V, DM,phi=0.01)
# fit13 = pglmSpatial(hwl ~ Temp + Prec + div_brd, D, V, DM,phi=0.01)
# fit14 = pglmSpatial(hwl ~ Temp*Prec + div_brd, D, V, DM,phi=0.01)

# fitList = list(fit1,fit2,fit3,fit4,fit5,fit6,fit7,fit8,fit9,fit10,fit11,fit12,fit13,fit14)


save(fitList,file="C:/Users/Owner/Desktop/fitList.rda")
}


if(FALSE) { 
load(file="C:/Users/Owner/Desktop/fitList.rda")

aics = setNames(sapply(fitList,function(x) {x$model$aic}),paste("fit",1:6))

d = data.frame(aics)
d$dummy = 1
d = d[order(d$aics),]
d

jpeg("C:/Users/Owner/Desktop/plot.jpg")
plot(d$aics)
dev.off()

# fit 2  3982.600     1
# fit 12 3982.938     1
# fit 13 3983.002     1
# fit 1  3983.101     1
# fit 3  3983.101     1

# sort(aics)


setNames(sapply(fitList,function(x) {x$model$coef}),paste("fit",1:5))
# setNames(sapply(fitList,function(x) {x$phi}),paste("fit",1:15))
# setNames(sapply(fitList,function(x) {x$k}),paste("fit",1:15))

# create supplementary table 

}

if(FALSE) { # ignore bias produce simple species density 

load(file = "C:/Users/Owner/Desktop/envBodySize/data/gbif/climateLat.rda") # load un-cleaned climate data
C = climateLat # C for climate
C$GenusSpecies = gsub(" ","_",C$GenusSpecies) 
C = C[c("GenusSpecies", "decimalLatitude", "decimalLongitude")]
colnames(C) = c("GenusSpecies", "lat", "lon")

D = aggregate(cbind(lat,lon) ~ GenusSpecies, data = C, mean) # average 

if(FALSE) {
# D = D[1:500,]

# get simple density 
library(geosphere)
d = D
den = c()
for(j in 1:nrow(d)) {
p1 = d[c("lon","lat")][j,]
dis = c()
for(i in 1:nrow(d)) {
p2 = d[c("lon","lat")][i,]
dis[i] = distVincentySphere(p1, p2, r=6378137) # radius of earth
}
dis = dis[dis != 0]
den[j] = mean(dis)
print(den[j]) 
}
D$den = den
sppDensity = D
save(sppDensity,file="C:/Users/Owner/Desktop/sppDensity.rda")

}

}

if(FALSE) { # Extract numbers from Fig. 1 for results section 

load(file = "C:/Users/Owner/Desktop/envBodySize/data/opdb/quickBodySize.rda") # load cleaned body size data
B = quickBodySize # B for body size
# B = B[c("GenusSpecies","Species","Genus","Family","SubOrder","tbl","hwl")] # drop fwl since often NA

B_hwl = B[c("GenusSpecies","Species","Genus","Family","SubOrder","hwl")] # drop another variable so that is easy to work with 
B_hwl = na.omit(B_hwl) # remove missing 

B_tbl = B[c("GenusSpecies","Species","Genus","Family","SubOrder","tbl")] # drop another variable so that is easy to work with 
B_tbl = na.omit(B_tbl) # remove missing 


load(file = "C:/Users/Owner/Desktop/envBodySize/data/gbif/climateLat.rda") # load un-cleaned climate data
C = climateLat # C for climate
C$GenusSpecies = gsub(" ","_",C$GenusSpecies) 

AL = aggregate(C$decimalLatitude ~ GenusSpecies, data = C, median) # average latitude
colnames(AL) = c("GenusSpecies","lat")
AL$climate = ifelse(AL$lat < 23 & AL$lat > -23,"trop","temp") # define tropical climates

D_hwl = merge(B_hwl,AL,id="GenusSpecies")
D_tbl = merge(B_tbl,AL,id="GenusSpecies")

library(ape) 
load("C:/Users/Owner/Desktop/envBodySize/data/trees/tree.rda") # load tree for pgls


treeD = function(tree,D_hwl,D_tbl,group) {  # function to process 

D_hwl = D_hwl[D_hwl$GenusSpecies %in% tree$tip.label,] # get only those with tip
D_tbl = D_tbl[D_tbl$GenusSpecies %in% tree$tip.label,] # get only those with tip

tip_hwl = tree$tip.label[!tree$tip.label %in% D_hwl$GenusSpecies] # drop tips without data
tip_tbl = tree$tip.label[!tree$tip.label %in% D_tbl$GenusSpecies] # drop tips without data

tree_hwl = drop.tip(tree, tip_hwl)
tree_tbl = drop.tip(tree, tip_tbl)

library(nlme)
rownames(D_hwl) = D_hwl$GenusSpecies
V_hwl = corBrownian(1,phy=tree_hwl) 
fit_hwl = gls(hwl ~ climate,correlation=V_hwl,data=D_hwl)

rownames(D_tbl) = D_tbl$GenusSpecies
V_tbl = corBrownian(1,phy=tree_tbl) 
fit_tbl = gls(tbl ~ climate,correlation=V_tbl,data=D_tbl)

print(summary(fit_tbl))

hwl_means = c(mean(D_hwl[D_hwl$climate == "trop",]$hwl),mean(D_hwl[D_hwl$climate == "temp",]$hwl))
names(hwl_means) = c("trop","temp")
tbl_means = c(mean(D_tbl[D_tbl$climate == "trop",]$tbl),mean(D_tbl[D_tbl$climate == "temp",]$tbl))
names(tbl_means) = c("trop","temp")


cat(paste(group,": \n"))

cat("tbl;\n")
cat("trop:\n")
cat(paste(round(tbl_means[1],1),"\n"))
cat("temp:\n")
cat(paste(round(tbl_means[2],1),"\n"))
cat(paste("P = "))
cat(paste(round(summary(fit_tbl)$tTable[,4][2],4),"\n"))

cat("hwl;\n")
cat("trop:\n")
cat(paste(round(hwl_means[1],1),"\n"))
cat("temp:\n")
cat(paste(round(hwl_means[2],1),"\n"))
cat(paste("P = "))
cat(paste(round(summary(fit_hwl)$tTable[,4][2],4),"\n"))



}


D_hwl_in = D_hwl[D_hwl$SubOrder == "Zygoptera",]; D_tbl_in = D_tbl[D_tbl$SubOrder == "Zygoptera",]
out = treeD(tree,D_hwl_in,D_tbl_in,group="Zygoptera") # function to process 

D_hwl_in = D_hwl[D_hwl$SubOrder == "Anisoptera",]; D_tbl_in = D_tbl[D_tbl$SubOrder == "Anisoptera",]
out = treeD(tree,D_hwl_in,D_tbl_in,group="Anisoptera") # function to process 

}

if(FALSE) { # pagel and corDisc AIC test and positive shifts in body size optima

load(file="C:/Users/Owner/Desktop/envBodySize/data/optima/optima.rda")
E = optima
E = E[E$tn == "tip",] # filter edge matrix

examples = c("Calopteryx_exul","Ictinogomphus_decoratus","Aeshna_eremita","Macromia_taeniolata","Megaloprepus_caerulatus","Mecistogaster_linearis","Ischnura_aurora") 
# examples = c("Calopteryx_exul","Ictinogomphus_decoratus","Aeshna_eremita","Macromia_taeniolata","Ischnura_aurora") 

example_optima = c(); for(i in 1:length(examples)) example_optima[i] = E[E$label == examples[i],]$optima_tip

E$shifts = ifelse(E$optima_tip %in% example_optima,"pos","neg")
E = E[c("label","shifts")]
colnames(E) = c("GenusSpecies","shifts") 

load(file = "C:/Users/Owner/Desktop/envBodySize/data/gbif/climateLat.rda") # load un-cleaned climate data
C = climateLat # C for climate
C$GenusSpecies = gsub(" ","_",C$GenusSpecies) 

AL = aggregate(C$decimalLatitude ~ GenusSpecies, data = C, mean) # average latitude
colnames(AL) = c("GenusSpecies","lat")

AL$climate = ifelse(AL$lat < 23 & AL$lat > -23,"trop","temp") # define tropical climates

D = merge(E,AL,id="GenusSpecies")

library(ape) 
load("C:/Users/Owner/Desktop/envBodySize/data/trees/tree.rda") # load tree for pgls
D = D[D$GenusSpecies %in% tree$tip.label,] # get only those with tip

tip = tree$tip.label[!tree$tip.label %in% D$GenusSpecies] # drop tips without data
tree = drop.tip(tree, tip)

trait = data.frame(Genus_sp=D$GenusSpecies,T1=D$climate,T2=D$shifts)
trait$T1 = ifelse(trait$T1=="temp",1,0)
trait$T2 = ifelse(trait$T2=="pos",1,0)

library(corHMM)

# fit_ARD = corDISC(phy=tree,data=trait, ntraits=2, rate.mat=NULL, model="ARD", node.states="joint", lewis.asc.bias=FALSE, p=NULL, root.p=NULL, ip=1, lb=0, ub=100, diagn=FALSE)
# fit_ER = corDISC(phy=tree,data=trait, ntraits=2, rate.mat=NULL, model="ER", node.states="joint", lewis.asc.bias=FALSE, p=NULL, root.p=NULL, ip=1, lb=0, ub=100, diagn=FALSE)
# fit_SYM = corDISC(phy=tree,data=trait, ntraits=2, rate.mat=NULL, model="SYM", node.states="joint", lewis.asc.bias=FALSE, p=NULL, root.p=NULL, ip=1, lb=0, ub=100, diagn=FALSE)

# sort(setNames(c(fit_ARD$AIC,fit_ER$AIC,fit_SYM$AIC),c("ARD","ER","SYM"))) # look at AIC
     # ARD      SYM       ER 
# 662.6291 699.3973 804.0975


# (temp, pos)
# 0 = trop, 1 = pos
# 1 = temp, 1 pos
      # (0,0) (0,1) (1,0) (1,1)
# (0,0)    NA     3     5    NA
# (0,1)     1    NA    NA     7
# (1,0)     2    NA    NA     8
# (1,1)    NA     4     6    NA
            # (0,0)        (0,1)      (1,0)        (1,1)
# (0,0)          NA 0.0004337238 0.01026007           NA
# (0,1) 0.000465252           NA         NA 0.0238275184
# (1,0) 0.002754444           NA         NA 0.0004450216
# (1,1)          NA 0.0044976031 0.00000000           NA
# fit_ARD
# fit_ARD$index.mat
       #     (trop0,neg0)  (trop0,pos1)  (temp1,neg0) (temp1,pos1)
# (trop0,neg0)    NA         0.0004(3)        0.0102(5)          NA
# (trop0,pos1) 0.0004(1)          NA             NA         0.0238(7)
# (temp1,neg0) 0.0027(2)          NA             NA         0.0004(8)
# (temp1,pos1)   NA          0.0044(4)        0.0000(6)          NA

Q = rate.mat.maker(hrm=FALSE, ntraits=2, model="ARD")

fit1=corDISC(phy=tree,data=trait, ntraits=2, rate.mat=Q, model="ARD", node.states="joint", lewis.asc.bias=FALSE, p=NULL, root.p=NULL, ip=1, lb=0, ub=100, diagn=FALSE)

library(corrplot)
M = fit1$solution
M[is.na(M)] = 0
M = M*100
rownames(M) = c("trop & neutral", "trop & pos. shift", "temp & neutral", "temp & pos. shift")
colnames(M) = c("trop & neutral", "trop & pos. shift", "temp & neutral", "temp & pos. shift")

# jpeg("C:/Users/Owner/Desktop/plot.jpg")
pdf("C:/Users/Owner/Desktop/corPlot1.pdf",height=10,width=10,useDingbats=FALSE)
corrplot(M, method="circle",is.corr=FALSE)
dev.off() 


# Q2 = rate.par.drop(Q, c(7))

# fit2=corDISC(phy=tree,data=trait, ntraits=2, rate.mat=Q2, model="ARD", node.states="joint", lewis.asc.bias=FALSE, p=NULL, root.p=NULL, ip=1, lb=0, ub=100, diagn=FALSE)

# sort(setNames(c(fit1$AIC,fit2$AIC),c("fit1","fit2")))




# x = setNames(D$shifts,D$GenusSpecies)
# y = setNames(D$climate,D$GenusSpecies)

# Pagel alternative models 
# library(phytools)
# fit = fitPagel(tree,x,y)
# fit
# round(fit$independent.Q,3)

}

if(FALSE) { # create optima.r60r shifts from l1ou model 

load(file="C:/Users/Owner/Desktop/BodySize/Manuscript/Data/l1ou/fit_z.rda") # load supertree

Anisoptera = function() { # function to load tree and associated data
library(ape)
library(phytools)
library(mvMORPH)

load("C:/Users/Owner/Desktop/BodySize/Manuscript/Data/trees/Anisoptera_tree.rda") # load Zygoptera
load("C:/Users/Owner/Desktop/BodySize/Manuscript/Data/Species.rda") # load Species
D = read.table("C:/Users/Owner/Desktop/BodySize/Manuscript/Data/OdonateData.csv",header=TRUE,sep=";") # load OdonateData

D$GenusSpecies = paste(D$genus,D$species,sep="_") # GenusSpecies correct
Species$GenusSpecies = gsub(" ","_",Species$GenusSpecies) # GenusSpecies correct

D = merge(D,Species,by="GenusSpecies",all.x=TRUE)
D = D[c("SubOrder","Family","Genus","GenusSpecies","body_lengths")]

# im = function(x) replace(x, is.na(x), mean(x, na.rm = TRUE))
# D = plyr::ddply(D, ~ Genus, transform, x = im(body_lengths)) # impute mean place in x
D$x = D$body_lengths # add x
# D$x = rnorm(length(D$body_lengths)) # noise to test is some problem with tree

D = na.omit(D[c("SubOrder","Family","Genus","GenusSpecies","x")]) # remove missing
D = D[D$GenusSpecies %in% tree$tip.label,] # get only those with tip
D = D[!duplicated(D$GenusSpecies),] # remove duplicated 
D = D[D$SubOrder == "Anisoptera",] # get only dragonflies, remove outgroup
tip = tree$tip.label[!tree$tip.label %in% D$GenusSpecies] # drop tips without data
tree = drop.tip(tree, tip)

data = D[c("x")] 
rownames(data) = D$GenusSpecies

return(list(tree=tree,data=data)) 
}

Zygoptera = function() { # function to load tree and associated data
library(ape)
library(phytools)
library(mvMORPH)

load("C:/Users/Owner/Desktop/BodySize/Manuscript/Data/trees/Zygoptera_tree.rda") # load Zygoptera
load("C:/Users/Owner/Desktop/BodySize/Manuscript/Data/Species.rda") # load Species
D = read.table("C:/Users/Owner/Desktop/BodySize/Manuscript/Data/OdonateData.csv",header=TRUE,sep=";") # load OdonateData

D$GenusSpecies = paste(D$genus,D$species,sep="_") # GenusSpecies correct
Species$GenusSpecies = gsub(" ","_",Species$GenusSpecies) # GenusSpecies correct

D = merge(D,Species,by="GenusSpecies",all.x=TRUE)
D = D[c("SubOrder","Family","Genus","GenusSpecies","body_lengths")]

# im = function(x) replace(x, is.na(x), mean(x, na.rm = TRUE))
# D = plyr::ddply(D, ~ Genus, transform, x = im(body_lengths)) # impute mean place in x
D$x = D$body_lengths # add x
# D$x = rnorm(length(D$body_lengths)) # noise to test is some problem with tree

D = na.omit(D[c("SubOrder","Family","Genus","GenusSpecies","x")]) # remove missing
D = D[D$GenusSpecies %in% tree$tip.label,] # get only those with tip
D = D[!duplicated(D$GenusSpecies),] # remove duplicated 
D = D[D$SubOrder == "Zygoptera",] # get only dragonflies, remove outgroup
tip = tree$tip.label[!tree$tip.label %in% D$GenusSpecies] # drop tips without data
tree = drop.tip(tree, tip)

data = D[c("x")] # data for mvEB,mvBM,mvOU
rownames(data) = D$GenusSpecies

return(list(tree=tree,data=data)) 
}

library(phytools)
library(ape)
library(methods)

fEdge = function(fit,obj) { # function to create edge matrix data 

tree = obj$tree
tree_ou = fit$tree

Ef = function(tree) {
E = data.frame(num1=tree$edge[,1],num2=tree$edge[,2],l=tree$edge.length)
E$tn = ifelse(E$num2 >= length(tree$tip.label)+1,"node","tip")
E$label[E$tn == "tip"] = tree$tip.label
# E = E[order(E$num2),]
E$id = paste0(E$num1,"_",E$num2)
return(E)
}

E = Ef(tree)
E$ord = 1:nrow(E)
E = E[c("id","num1","num2","l","tn","ord","label")]

E_ou = Ef(tree_ou)
E_ou$optima = fit$edge.optima[,1]
E_ou$map = as.numeric(as.factor(E_ou$optima))
E_ou = E_ou[c("id","optima","map")]
E = merge(E,E_ou,by="id")
E = E[order(E$ord),]
E$x[E$tn == "tip"] = obj$data[,1]

E$optima_tip[E$tn == "tip"] = fit$optima

return(E)
}


load(file="C:/Users/Owner/Desktop/BodySize/Manuscript/Data/l1ou/fit_a.rda") # load fitted l1ou model for Anisoptera
obj_a = Anisoptera()
E_a = fEdge(fit,obj_a)

obj_z = Zygoptera()
load(file="C:/Users/Owner/Desktop/BodySize/Manuscript/Data/l1ou/fit_z.rda") # load fitted l1ou model for Zygoptera
E_z = fEdge(fit,obj_z)

tree_z = obj_z$tree
tree_a = obj_a$tree

max_ta = max(nodeHeights(tree_a)[,2])
max_tz = max(nodeHeights(tree_z)[,2])

# Odonata 237 - root in R
# Zygotera 132.9 x
root_age = 237

tree_arm = list(edge = matrix(c(3,3,1,2),2),edge.length=c(root_age-max_ta, root_age-max_tz),tip.label=c("ta","tz"),Nnode=1)
class(tree_arm)="phylo"

tree_arm$tip.label[tree_arm$tip.label == "ta"] = "NA"
tree_a$root.edge = 0
tree = paste.tree(tree_arm, tree_a)

tree$tip.label[tree$tip.label == "tz"] = "NA"
tree_z$root.edge = 0
tree = paste.tree(tree, tree_z)

arm_a = data.frame(id=NA,num1=NA,num2=NA,l=tree$edge.length[1],tn="node",ord=NA,label=NA,optima=mean(E_a$x,na.rm=TRUE),map=NA,x=NA,optima_tip=NA) # 729 - 730
E_a = rbind(arm_a,E_a)

arm_z = data.frame(id=NA,num1=NA,num2=NA,l=tree$edge.length[784],tn="node",ord=NA,label=NA,optima=mean(E_z$x,na.rm=TRUE),map=NA,x=NA,optima_tip=NA) 
E_z = rbind(arm_z,E_z)

E = rbind(E_a,E_z) # edge data

optima = E
# save(optima,file="C:/Users/Owner/Desktop/optima.rda")
}

if(FALSE) { # fit pagel example and intuition building

library(phytools)
tree = pbtree(n=300,scale=1)
Q = matrix(c(0,0.4,0.4,0,2,0,0,2,2,0,0,2,0,0.4,0.4,0),4,4,byrow=TRUE)

Q = rbind(
c(0,0.4,0.4,0),
c(2,0.0,0.0,2),
c(2,0.0,0.0,2),
c(0,0.4,0.4,0)
)

rownames(Q) = colnames(Q)<-c("red","blue","green","black")
diag(Q) <- -rowSums(Q)

Q
# aa = both red
# bb = both blue
# ab = left red, right blue
# ba = left blue, right red

tt = sim.history(tree,Q,anc = "green")

# "red","blue","green","black"
# "aa","bb","ab","ba"

tt$states = getStates(tt,"tips")
pdf("C:/Users/Owner/Desktop/plot.pdf")
plotSimmap(tt,setNames(c("red","blue","green","black"),c("red","blue","green","black")),lwd=1,ftype="off")
dev.off()

# t1 = mergeMappedStates(tt,c("aa","ab"),"a")
# t1 = mergeMappedStates(t1,c("ba","bb"),"b")
# t2 = mergeMappedStates(tt,c("aa","ba"),"a")
# t2 = mergeMappedStates(t2,c("ab","bb"),"b")
# t1$states = getStates(t1,"tips")
# t2$states = getStates(t2,"tips")

# pdf("C:/Users/Owner/Desktop/plot.pdf")
# par(mfrow=c(1,2))
# plotSimmap(t1,setNames(c("red","blue"),letters[1:2]),lwd=1,ftype="off")
# plotSimmap(t2,setNames(c("red","blue"),letters[1:2]),lwd=1,ftype="off",direction="leftwards")
# dev.off()

# x = getStates(t1,"tips")
# y = getStates(t2,"tips")

# fit = fitPagel(tree,x,y)
# fit
}


if(FALSE) { # mvMORPH; do temperate and tropical species have different body size optima? 

load(file = "C:/Users/Owner/Desktop/envBodySize/data/gbif/climateLat.rda") # load un-cleaned climate data
C = climateLat # C for climate
C$GenusSpecies = gsub(" ","_",C$GenusSpecies) 

AL = aggregate(C$decimalLatitude ~ GenusSpecies, data = C, mean) # average latitude
colnames(AL) = c("GenusSpecies","lat")

AL$climate = ifelse(AL$lat < 23 & AL$lat > -23,"tropical","temperate") # define tropical climates

load(file = "C:/Users/Owner/Desktop/envBodySize/data/opdb/quickBodySize.rda") # load cleaned body size data
B = quickBodySize # B for body size
B = B[c("GenusSpecies","Species","Genus","Family","SubOrder","tbl","hwl")] # drop fwl since often NA
# B = B[c("GenusSpecies","Species","Genus","Family","SubOrder","hwl")] # drop another variable so that is easy to work with 
B = na.omit(B) # remove missing 

D = merge(B,AL,id="GenusSpecies")

library(ape) 
load("C:/Users/Owner/Desktop/envBodySize/data/trees/tree.rda") # load tree for pgls

D = D[D$GenusSpecies %in% tree$tip.label,] # get only those with tip

tip = tree$tip.label[!tree$tip.label %in% D$GenusSpecies] # drop tips without data
tree = drop.tip(tree, tip)

state = setNames(D$climate,D$GenusSpecies) # set up climate variable for acestral state reconstruction

library(phytools)
tree = make.simmap(tree, state, model="ER", nsim=1)

data = D[c("tbl","hwl")]
rownames(data) = D$GenusSpecies

library(mvMORPH)
OUM = mvOU(tree, data, model="OUM")
OU1 = mvOU(tree, data, model="OU1")

AIC(OUM)
AIC(OU1)

# library(mvMORPH)
# trait1_OU1 = mvOU(tree, data[,1], model="OU1", diagnostic=FALSE, echo=FALSE)
# trait1_OUM = mvOU(tree, data[,1], model="OUM", diagnostic=FALSE, echo=FALSE)
# AIC(trait1_OU1)
# AIC(trait1_OUM)

# test for different optima 

# test if temperate are bigger than tropical; yes they are by around 10 mm
# library(nlme)
# V = corBrownian(1,phy=tree) 
# rownames(D) = D$GenusSpecies
# V = corPagel(1, phy = tree, form = ~1, fixed = FALSE)
# fit = gls(hwl ~ climate,correlation=V,data=D)
# summary(fit)

# rownames(D) = D$GenusSpecies
# library(phylolm)
# fit = phyloglm(climate ~ hwl,phy=tree,data=D,boot=1000)
# summary(fit)$coefficients
# round(summary(fit)$coefficients[,4],3)



# jpeg("C:/Users/Owner/Desktop/plot.jpg")
# plot(AL$lat,col=as.factor(AL$climate))
# dev.off()

# colnames(C)
# jpeg("C:/Users/Owner/Desktop/plot.jpg")
# hist(C$decimalLatitude)
# dev.off()

}

if(FALSE) { # mvMORPH multiple optima for different environments example

library(mvMORPH)

tree = pbtree(n=100)
state = as.vector(c(rep("Forest",60),rep("Savannah",40))); names(state)<-tree$tip.label
state
tree = make.simmap(tree, state, model="ER", nsim=1)

# col<-c("blue","orange"); names(col)<-c("Forest","Savannah")
# plotSimmap(tree,col, fsize=0.6, node.numbers=FALSE, lwd=3, pts=FALSE)

set.seed(101)
alpha = matrix(c(1.1,-0.9,-0.9,1),2)
sigma = matrix(c(0.35,0.06,0.06,0.35),2)
theta = c(5.5,5.1,1.2,1.4)
data = mvSIM(tree, param=list(sigma=sigma, alpha=alpha, ntraits=2, theta=theta,
names_traits = c("limb.length","limb.width")), model="OUM", nsim=1)

data

# Fitting the Ornstein Uhlenbeck on the whole tree
trait1_OU1 = mvOU(tree, data[,1], model="OU1", diagnostic=FALSE, echo=FALSE)
trait2_OU1 = mvOU(tree, data[,2], model="OU1", diagnostic=FALSE, echo=FALSE)
# Fitting the Ornstein Uhlenbeck with multiple optimums
trait1_OUM = mvOU(tree, data[,1], model="OUM", diagnostic=FALSE, echo=FALSE)
trait2_OUM = mvOU(tree, data[,2], model="OUM", diagnostic=FALSE, echo=FALSE)
# Compare the AIC values between models fit

# AIC(trait1_OUM)
# AIC(trait1_OU1)
# AIC(trait2_OUM)
# AIC(trait2_OU1)

OUM = mvOU(tree, data, model="OUM")
OU1 = mvOU(tree, data, model="OU1")

AIC(OUM)
AIC(OU1)
}

if(FALSE) { # logistic regression of shifts vs temperature

load(file="C:/Users/Owner/Desktop/BodySize/Manuscript/Data/l1ou/fit_z.rda") # load supertree

Anisoptera = function() { # function to load tree and associated data
library(ape)
library(phytools)
library(mvMORPH)

load("C:/Users/Owner/Desktop/BodySize/Manuscript/Data/trees/Anisoptera_tree.rda") # load Zygoptera
load("C:/Users/Owner/Desktop/BodySize/Manuscript/Data/Species.rda") # load Species
D = read.table("C:/Users/Owner/Desktop/BodySize/Manuscript/Data/OdonateData.csv",header=TRUE,sep=";") # load OdonateData

D$GenusSpecies = paste(D$genus,D$species,sep="_") # GenusSpecies correct
Species$GenusSpecies = gsub(" ","_",Species$GenusSpecies) # GenusSpecies correct

D = merge(D,Species,by="GenusSpecies",all.x=TRUE)
D = D[c("SubOrder","Family","Genus","GenusSpecies","body_lengths")]

# im = function(x) replace(x, is.na(x), mean(x, na.rm = TRUE))
# D = plyr::ddply(D, ~ Genus, transform, x = im(body_lengths)) # impute mean place in x
D$x = D$body_lengths # add x
# D$x = rnorm(length(D$body_lengths)) # noise to test is some problem with tree

D = na.omit(D[c("SubOrder","Family","Genus","GenusSpecies","x")]) # remove missing
D = D[D$GenusSpecies %in% tree$tip.label,] # get only those with tip
D = D[!duplicated(D$GenusSpecies),] # remove duplicated 
D = D[D$SubOrder == "Anisoptera",] # get only dragonflies, remove outgroup
tip = tree$tip.label[!tree$tip.label %in% D$GenusSpecies] # drop tips without data
tree = drop.tip(tree, tip)

data = D[c("x")] 
rownames(data) = D$GenusSpecies

return(list(tree=tree,data=data)) 
}

Zygoptera = function() { # function to load tree and associated data
library(ape)
library(phytools)
library(mvMORPH)

load("C:/Users/Owner/Desktop/BodySize/Manuscript/Data/trees/Zygoptera_tree.rda") # load Zygoptera
load("C:/Users/Owner/Desktop/BodySize/Manuscript/Data/Species.rda") # load Species
D = read.table("C:/Users/Owner/Desktop/BodySize/Manuscript/Data/OdonateData.csv",header=TRUE,sep=";") # load OdonateData

D$GenusSpecies = paste(D$genus,D$species,sep="_") # GenusSpecies correct
Species$GenusSpecies = gsub(" ","_",Species$GenusSpecies) # GenusSpecies correct

D = merge(D,Species,by="GenusSpecies",all.x=TRUE)
D = D[c("SubOrder","Family","Genus","GenusSpecies","body_lengths")]

# im = function(x) replace(x, is.na(x), mean(x, na.rm = TRUE))
# D = plyr::ddply(D, ~ Genus, transform, x = im(body_lengths)) # impute mean place in x
D$x = D$body_lengths # add x
# D$x = rnorm(length(D$body_lengths)) # noise to test is some problem with tree

D = na.omit(D[c("SubOrder","Family","Genus","GenusSpecies","x")]) # remove missing
D = D[D$GenusSpecies %in% tree$tip.label,] # get only those with tip
D = D[!duplicated(D$GenusSpecies),] # remove duplicated 
D = D[D$SubOrder == "Zygoptera",] # get only dragonflies, remove outgroup
tip = tree$tip.label[!tree$tip.label %in% D$GenusSpecies] # drop tips without data
tree = drop.tip(tree, tip)

data = D[c("x")] # data for mvEB,mvBM,mvOU
rownames(data) = D$GenusSpecies

return(list(tree=tree,data=data)) 
}

library(phytools)
library(ape)
library(methods)

fEdge = function(fit,obj) { # function to create edge matrix data 

tree = obj$tree
tree_ou = fit$tree

Ef = function(tree) {
E = data.frame(num1=tree$edge[,1],num2=tree$edge[,2],l=tree$edge.length)
E$tn = ifelse(E$num2 >= length(tree$tip.label)+1,"node","tip")
E$label[E$tn == "tip"] = tree$tip.label
# E = E[order(E$num2),]
E$id = paste0(E$num1,"_",E$num2)
return(E)
}

E = Ef(tree)
E$ord = 1:nrow(E)
E = E[c("id","num1","num2","l","tn","ord","label")]

E_ou = Ef(tree_ou)
E_ou$optima = fit$edge.optima[,1]
E_ou$map = as.numeric(as.factor(E_ou$optima))
E_ou = E_ou[c("id","optima","map")]
E = merge(E,E_ou,by="id")
E = E[order(E$ord),]
E$x[E$tn == "tip"] = obj$data[,1]

E$optima_tip[E$tn == "tip"] = fit$optima

return(E)
}


load(file="C:/Users/Owner/Desktop/BodySize/Manuscript/Data/l1ou/fit_a.rda") # load fitted l1ou model for Anisoptera
obj_a = Anisoptera()
E_a = fEdge(fit,obj_a)

obj_z = Zygoptera()
load(file="C:/Users/Owner/Desktop/BodySize/Manuscript/Data/l1ou/fit_z.rda") # load fitted l1ou model for Zygoptera
E_z = fEdge(fit,obj_z)

tree_z = obj_z$tree
tree_a = obj_a$tree

max_ta = max(nodeHeights(tree_a)[,2])
max_tz = max(nodeHeights(tree_z)[,2])

# Odonata 237 - root in R
# Zygotera 132.9 x
root_age = 237

tree_arm = list(edge = matrix(c(3,3,1,2),2),edge.length=c(root_age-max_ta, root_age-max_tz),tip.label=c("ta","tz"),Nnode=1)
class(tree_arm)="phylo"

tree_arm$tip.label[tree_arm$tip.label == "ta"] = "NA"
tree_a$root.edge = 0
tree = paste.tree(tree_arm, tree_a)

tree$tip.label[tree$tip.label == "tz"] = "NA"
tree_z$root.edge = 0
tree = paste.tree(tree, tree_z)

arm_a = data.frame(id=NA,num1=NA,num2=NA,l=tree$edge.length[1],tn="node",ord=NA,label=NA,optima=mean(E_a$x,na.rm=TRUE),map=NA,x=NA,optima_tip=NA) # 729 - 730
E_a = rbind(arm_a,E_a)

arm_z = data.frame(id=NA,num1=NA,num2=NA,l=tree$edge.length[784],tn="node",ord=NA,label=NA,optima=mean(E_z$x,na.rm=TRUE),map=NA,x=NA,optima_tip=NA) 
E_z = rbind(arm_z,E_z)

E = rbind(E_a,E_z) # edge data
E = E[E$tn == "tip",] # filter edge matrix
# "Calopteryx_exul" = 1

# E$shifts = NA
examples = c("Calopteryx_exul","Ictinogomphus_decoratus","Aeshna_eremita","Macromia_taeniolata","Megaloprepus_caerulatus","Mecistogaster_linearis","Ischnura_aurora") 
examples = c("Calopteryx_exul","Ictinogomphus_decoratus","Aeshna_eremita","Macromia_taeniolata","Ischnura_aurora") 

example_optima = c(); for(i in 1:length(examples)) example_optima[i] = E[E$label == examples[i],]$optima_tip

E$shifts = ifelse(E$optima_tip %in% example_optima,1,0)

D = E[c("label","x","shifts")]
colnames(D) = c("GenusSpecies","x","shifts") 

load(file = "C:/Users/Owner/Desktop/envBodySize/data/gbif/climateLat.rda") # load un-cleaned climate data
C = climateLat # C for climate
C$GenusSpecies = gsub(" ","_",C$GenusSpecies) 
AC = aggregate(cbind(bio_max1,bio_max12) ~ GenusSpecies, data = C, mean) # average 

D = merge(D,AC,id="GenusSpecies",all.x=TRUE) # 
colnames(D) = c("GenusSpecies","x","shifts","temp","prec")

load(file="C:/Users/Owner/Desktop/envBodySize/data/gbif/Birds.rda") # load Birds
B = Birds # rename B for Birds
# str(B) # what variable are available
B$diversity = B$all # select Bird variable to analyze 
AB = aggregate(diversity ~ GenusSpecies,data=B,mean) # get average birds diversity 
AB$GenusSpecies = gsub(" ","_",AB$GenusSpecies)

D = merge(D,AB,id="GenusSpecies",all.x=TRUE) # 
D = na.omit(D)
rownames(D) = D$GenusSpecies

tip = tree$tip.label[!tree$tip.label %in% D$GenusSpecies] # drop tips without data
tree = drop.tip(tree, tip)

# naive analysis
# fit = glm(shifts ~ diversity,data=D,family=binomial)
# summary(fit)
  
# normalize climate
# D$temp = scale(D$temp)  
# D$prec = scale(D$prec)  
# D$diversity = scale(D$diversity)  

library(phylolm)
fit = phyloglm(shifts ~ diversity,phy=tree,data=D,boot=1000)
summary(fit)$coefficients
round(summary(fit)$coefficients[,4],3)

    
# library(nlme)
# V = corBrownian(1,phy=tree) 
# rownames(D) = D$GenusSpecies
# V = corPagel(1, phy = tree, form = ~1, fixed = FALSE)
# fit = gls(shifts ~ temp,correlation=V,data=D)
# summary(fit)
} 
 
if(FALSE) { # phylo path analysis temp, dev, tbl 
library(ape)
load(file="C:/Users/Owner/Desktop/BodySize/Manuscript/Data/trees/tree.rda") # load supertree

load(file="C:/Users/Owner/Desktop/BodySize/Manuscript/Data/Development/Development.rda") # load development time
D = Development # rename

im = function(x) replace(x, is.na(x), mean(x, na.rm = TRUE))
D = plyr::ddply(D, ~ Genus, transform, x = im(body_lengths)) # replace missing taxa with mean of Genus
D = na.omit(D[c("GenusSpecies","dev","MeanLat","x")])
D$GenusSpecies = gsub(" ","_",D$GenusSpecies)

tip = tree$tip.label[!tree$tip.label %in% D$GenusSpecies] # drop tips without data
tree = drop.tip(tree, tip)

D = D[D$GenusSpecies %in% tree$tip.label,] # remove data without tip in tree

load(file = "C:/Users/Owner/Desktop/envBodySize/data/gbif/climateLat.rda") # load un-cleaned climate data
C = climateLat # C for climate
C$GenusSpecies = gsub(" ","_",C$GenusSpecies) 
AC = aggregate(cbind(bio_max1,bio_max12) ~ GenusSpecies, data = C, mean) # average 

D = merge(D,AC,id="GenusSpecies",all.x=TRUE) # 
colnames(D) = c("SP","dev","MeanLat","tbl","temp","prec")
rownames(D) = D$SP

library(phylopath)

models = list(
mod1 = DAG(tbl ~ temp + dev,temp ~ dev),
mod2 = DAG(tbl ~ dev,dev ~ temp),
mod3 = DAG(tbl ~ temp,dev ~ temp),
mod4 = DAG(tbl ~ dev, temp ~ temp),
mod5 = DAG(tbl ~ dev + temp, dev ~ temp),
mod6 = DAG(tbl ~ temp, dev ~ dev)
)

fit = phylo_path(models, data = D, tree = tree)
best_model = best(fit)
best_model
# best_model
jpeg("C:/Users/Owner/Desktop/plot.jpg")
# plot_model_set(models)
plot(best_model)
dev.off()
}

if(FALSE) { # phylo path analysis temp, dev, tbl, birds 
library(ape)
load(file="C:/Users/Owner/Desktop/BodySize/Manuscript/Data/trees/tree.rda") # load supertree

load(file="C:/Users/Owner/Desktop/BodySize/Manuscript/Data/Development/Development.rda") # load development time
D = Development # rename

im = function(x) replace(x, is.na(x), mean(x, na.rm = TRUE))
D = plyr::ddply(D, ~ Genus, transform, x = im(body_lengths)) # replace missing taxa with mean of Genus
D = na.omit(D[c("GenusSpecies","dev","MeanLat","x")])
D$GenusSpecies = gsub(" ","_",D$GenusSpecies)

tip = tree$tip.label[!tree$tip.label %in% D$GenusSpecies] # drop tips without data
tree = drop.tip(tree, tip)

D = D[D$GenusSpecies %in% tree$tip.label,] # remove data without tip in tree

load(file = "C:/Users/Owner/Desktop/envBodySize/data/gbif/climateLat.rda") # load un-cleaned climate data
C = climateLat # C for climate
C$GenusSpecies = gsub(" ","_",C$GenusSpecies) 
AC = aggregate(cbind(bio_max1,bio_max12) ~ GenusSpecies, data = C, mean) # average 
D = merge(D,AC,id="GenusSpecies",all.x=TRUE) # 

# add birds to path analysis 
load(file="C:/Users/Owner/Desktop/envBodySize/data/gbif/Birds.rda") # load Birds
B = Birds # rename B for Birds
# str(B) # what variable are available
B$diversity = B$all # select Bird variable to analyze 
AB = aggregate(diversity ~ GenusSpecies,data=B,mean) # get average birds diversity 
AB$GenusSpecies = gsub(" ","_",AB$GenusSpecies)

D = merge(D,AB,id="GenusSpecies",all.x=TRUE) # 
colnames(D) = c("SP","dev","MeanLat","tbl","temp","prec","birds")
rownames(D) = D$SP

library(phylopath)

models = list(
mod1 = DAG(tbl ~ birds + temp + dev,dev ~ temp,birds ~ temp),
mod2 = DAG(tbl ~ birds + dev,dev ~ temp,birds ~ temp),
mod3 = DAG(tbl ~ birds + dev,dev ~ temp),
mod4 = DAG(tbl ~ birds + dev + temp),
mod5 = DAG(tbl ~ birds + temp,dev ~ dev),
mod6 = DAG(tbl ~ birds + temp,dev ~ temp)
)

fit = phylo_path(models, data = D, tree = tree)
best_model = best(fit)

jpeg("C:/Users/Owner/Desktop/plot.jpg")
# plot_model_set(models)
plot(best_model)
dev.off()

}

if(FALSE) { # phylopath, temp, tbl, birds, no development

library(ape)
load(file="C:/Users/Owner/Desktop/BodySize/Manuscript/Data/trees/tree.rda") # load supertree

load(file = "C:/Users/Owner/Desktop/envBodySize/data/gbif/climateLat.rda") # load un-cleaned climate data
C = climateLat # C for climate
C$GenusSpecies = gsub(" ","_",C$GenusSpecies) 
AC = aggregate(cbind(bio_max1,bio_max12) ~ GenusSpecies, data = C, mean) # average 

# add birds to path analysis 
load(file="C:/Users/Owner/Desktop/envBodySize/data/gbif/Birds.rda") # load Birds
B = Birds # rename B for Birds
# str(B) # what variable are available
B$diversity = B$all # select Bird variable to analyze 
AB = aggregate(diversity ~ GenusSpecies,data=B,mean) # get average birds diversity 
AB$GenusSpecies = gsub(" ","_",AB$GenusSpecies)

D = merge(AC,AB,id="GenusSpecies") # merge climage and birds

load(file = "C:/Users/Owner/Desktop/envBodySize/data/opdb/quickBodySize.rda")
BS = quickBodySize # body size 
BS = BS[c("GenusSpecies","hwl")]
D = merge(D,BS,id="GenusSpecies") # merge body size and climate and birds
D = na.omit(D) # remove missing

D = D[D$GenusSpecies %in% tree$tip.label,] # remove data without tip in tree
tip = tree$tip.label[!tree$tip.label %in% D$GenusSpecies] # drop tips without data
tree = drop.tip(tree, tip)

colnames(D) = c("SP","temp","prec","birds","hwl")
rownames(D) = D$SP

nrow(D)

library(phylopath)

models = list(
mod1 = DAG(hwl ~ birds + temp,birds ~ temp),
mod2 = DAG(hwl ~ birds + temp)
)

fit = phylo_path(models, data = D, tree = tree)
best_model = best(fit)

jpeg("C:/Users/Owner/Desktop/plot.jpg")
# plot_model_set(models)
plot(best_model)
dev.off()
}

if(FALSE) { # simulate phylo path
library(phylopath)
library(ape)
N = 50
tree = rcoal(N)

SP = tree$tip.label
temp = rnorm(N)
dev = temp + rnorm(N)
tbl = dev + rnorm(N)

D = data.frame(SP = SP, tbl = tbl, dev = dev, temp = temp)
rownames(D) = SP

models = list(
mod1 = DAG(tbl ~ temp + dev,temp ~ dev),
mod2 = DAG(tbl ~ dev,dev ~ temp),
mod3 = DAG(tbl ~ temp,dev ~ temp),
mod4 = DAG(tbl ~ dev, temp ~ temp)
)

# rhino

fit = phylo_path(models, data = D, tree = tree)
best_model = best(fit)

# best_model
jpeg("C:/Users/Owner/Desktop/plot.jpg")
# plot_model_set(models)
plot(best_model)
dev.off()


# models = list(
# one   = DAG(LS ~ BM, NL ~ BM, DD ~ NL, RS ~ DD),
# two   = DAG(LS ~ BM, NL ~ BM, DD ~ NL, RS ~ LS + DD),
# three = DAG(LS ~ BM, NL ~ BM, DD ~ NL, RS ~ NL),
# four  = DAG(LS ~ BM, NL ~ BM, DD ~ NL, RS ~ BM + NL),
# five  = DAG(LS ~ BM, NL ~ BM, DD ~ NL, RS ~ BM + NL + DD),
# six   = DAG(LS ~ BM, NL ~ BM + RS, DD ~ NL, RS ~ BM),
# seven = DAG(LS ~ BM, NL ~ BM + RS, DD ~ NL, RS ~ LS + BM),
# eight = DAG(LS ~ BM, NL ~ BM + RS, DD ~ NL),
# nine  = DAG(LS ~ BM, NL ~ BM + RS, DD ~ NL, RS ~ LS)
# )

# SP 

# result = phylo_path(models, data = rhino, tree = rhino_tree, order = c('BM', 'NL', 'DD', 'LS', 'RS'))

# summary(result)

# best_model = best(result)
# best_model
# jpeg("C:/Users/Owner/Desktop/plot.jpg")
# plot_model_set(models)
# plot(best_model)
# dev.off()

}

if(FALSE) { # BIOCLIM definitions
# BIO1 = Annual Mean Temperature
# BIO2 = Mean Diurnal Range (Mean of monthly (max temp - min temp))
# BIO3 = Isothermality (BIO2/BIO7) (* 100)
# BIO4 = Temperature Seasonality (standard deviation *100)
# BIO5 = Max Temperature of Warmest Month
# BIO6 = Min Temperature of Coldest Month
# BIO7 = Temperature Annual Range (BIO5-BIO6)
# BIO8 = Mean Temperature of Wettest Quarter
# BIO9 = Mean Temperature of Driest Quarter
# BIO10 = Mean Temperature of Warmest Quarter
# BIO11 = Mean Temperature of Coldest Quarter
# BIO12 = Annual Precipitation
# BIO13 = Precipitation of Wettest Month
# BIO14 = Precipitation of Driest Month
# BIO15 = Precipitation Seasonality (Coefficient of Variation)
# BIO16 = Precipitation of Wettest Quarter
# BIO17 = Precipitation of Driest Quarter
# BIO18 = Precipitation of Warmest Quarter
# BIO19 = Precipitation of Coldest Quarter
}

if(FALSE) { # Create quickBodySize.rda; cleaned body size and wing data. 

D = read.table(file = "C:/Users/Owner/Desktop/envBodySize/data/opdb/OdonateData.csv",sep=";",header=TRUE)
load(file = "C:/Users/Owner/Desktop/envBodySize/data/Species.rda") # add Species data
Species$GenusSpecies = gsub(" ","_",Species$GenusSpecies)

D = D[D$sex == "male",] # get only males for quick body size
D$GenusSpecies = paste0(D$genus,"_",D$species) 
D = D[c("GenusSpecies","body_lengths","forewing_lengths","hindwing_lengths")]
D = merge(D,Species,id ="GenusSpecies") # new merged data
D = D[c("GenusSpecies","Species","Genus","Family","SubOrder","body_lengths","forewing_lengths","hindwing_lengths")]
D = D[!is.na(D$body_lengths) | !is.na(D$forewing_lengths)| !is.na(D$hindwing_lengths),] # remove rows with totally missing data
colnames(D) = c("GenusSpecies","Species","Genus","Family","SubOrder","tbl","fwl","hwl")
quickBodySize = D # rename

# save(quickBodySize, file = "C:/Users/Owner/Desktop/quickBodySize.rda") # save to desktop to avoid overwriting
}

if(FALSE) { # Is there a relationship between body size and temperature? Yes negative

load(file = "C:/Users/Owner/Desktop/envBodySize/data/opdb/quickBodySize.rda") # load cleaned body size data
B = quickBodySize # B for body size
B = B[c("GenusSpecies","Species","Genus","Family","SubOrder","tbl","hwl")] # drop fwl since often NA
B = B[c("GenusSpecies","Species","Genus","Family","SubOrder","hwl")] # drop another variable so that is easy to work with 
B = na.omit(B) # remove missing 

load(file = "C:/Users/Owner/Desktop/envBodySize/data/gbif/climateLat.rda") # load un-cleaned climate data
C = climateLat # C for climate
C$GenusSpecies = gsub(" ","_",C$GenusSpecies) 

# str(C)

library(ape) 
load("C:/Users/Owner/Desktop/envBodySize/data/trees/tree.rda") # load tree for pgls

AC = aggregate(bio_max1/10 ~ GenusSpecies, data = C, mean) # average climate
colnames(AC) = c("GenusSpecies","temp")
D = merge(B,AC,id="GenusSpecies") # final data.frame

# remove those not in tree or data
D = D[D$GenusSpecies %in% tree$tip.label,] # get only those with tip
D = D[!duplicated(D$GenusSpecies),] # not necessary but lets leave it in 
tip = tree$tip.label[!tree$tip.label %in% D$GenusSpecies] # drop tips without data
tree = drop.tip(tree, tip) # drop missing

# naive analysis 

library(nlme)
V = corBrownian(1,phy=tree) 
rownames(D) = D$GenusSpecies

# V = corPagel(1, phy = tree, form = ~1, fixed = FALSE)
fit = gls(hwl ~ temp,correlation=V,data=D)

}


if(FALSE) { # Is there a relationship between body size and precipitation?
load(file = "C:/Users/Owner/Desktop/envBodySize/data/opdb/quickBodySize.rda") # load cleaned body size data
B = quickBodySize # B for body size
B = B[c("GenusSpecies","Species","Genus","Family","SubOrder","tbl","hwl")] # drop fwl since often NA
B = B[c("GenusSpecies","Species","Genus","Family","SubOrder","hwl")] # drop another variable so that is easy to work with 
B = na.omit(B) # remove missing 

load(file = "C:/Users/Owner/Desktop/envBodySize/data/gbif/climateLat.rda") # load un-cleaned climate data
C = climateLat # C for climate
C$GenusSpecies = gsub(" ","_",C$GenusSpecies) 

library(ape) 
load("C:/Users/Owner/Desktop/envBodySize/data/trees/tree.rda") # load tree for pgls

# colnames(C)

AC = aggregate(cbind(bio_max1,bio_max12) ~ GenusSpecies, data = C, mean) # average climate
# colnames(AC) = c("GenusSpecies","bio_max1","bio_max12") # not necessary for formula aggregate
D = merge(B,AC,id="GenusSpecies") # final data.frame

# remove those not in tree or data
D = D[D$GenusSpecies %in% tree$tip.label,] # get only those with tip
D = D[!duplicated(D$GenusSpecies),] # not necessary but lets leave it in 
tip = tree$tip.label[!tree$tip.label %in% D$GenusSpecies] # drop tips without data
tree = drop.tip(tree, tip) # drop missing

# naive analysis 
D = na.omit(D[c("GenusSpecies","Family","Genus","SubOrder","hwl","bio_max1","bio_max12")])
D$bio_max1 = scale(D$bio_max1) # standardize
D$bio_max12 = scale(D$bio_max12)
fit = lm(hwl ~ bio_max1*bio_max12,data=D) # larger living in colder regions
# summary(fit)

library(nlme)
V = corBrownian(1,phy=tree) 
# rownames(D) = D$GenusSpecies
# V = corPagel(1, phy = tree, form = ~1, fixed = FALSE)
fit = gls(hwl ~ bio_max1*bio_max12,correlation=V,data=D)
# summary(fit)
}

if(FALSE) { # Bird diversity plots
Files = list.files("C:/Users/Owner/Desktop/envBodySize/data/gbif/biodiversitymapping_TIFFs/Birds/") # get only GeoTiffs
Files = Files[!grepl("aux|ovr|xml|dbf",Files)]
Files = Files[grepl(".tif",Files)]

library(raster)
library(rgdal)

Paths = paste0("C:/Users/Owner/Desktop/envBodySize/data/gbif/biodiversitymapping_TIFFs/Birds/",Files)

DL = list() 
for(i in 1:length(Paths)){
DL[[i]] = raster(Paths[i])
}

Names = gsub("richness_10km_","",Files) # Clean names
Names = gsub("_raster.tif","",Names)

i = 11
r = DL[[i]]

newproj = "+proj=longlat +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_def" # set projection coordinates
r2 = projectRaster(r, crs=newproj) # data projected to new coordinates

jpeg("C:/Users/Owner/Desktop/plot.jpg",units = "in",width=15,height=10,res=300)
plot(r2,main=Names[i])
dev.off()
}


if(FALSE) { # Process Bird data to produce Birds.rda
Files = list.files("C:/Users/Owner/Desktop/envBodySize/data/gbif/biodiversitymapping_TIFFs/Birds/") # get only GeoTiffs
Files = Files[!grepl("aux|ovr|xml|dbf",Files)]
Files = Files[grepl(".tif",Files)]

library(raster)
library(rgdal)

Paths = paste0("C:/Users/Owner/Desktop/envBodySize/data/gbif/biodiversitymapping_TIFFs/Birds/",Files)

DL = list() 
for(i in 1:length(Paths)){
DL[[i]] = raster(Paths[i])
}

Names = gsub("richness_10km_","",Files) # Clean names
Names = gsub("_raster.tif","",Names)
Names

load(file = "C:/Users/Owner/Desktop/envBodySize/data/gbif/climateLat.rda") # load un-cleaned climate data
C = climateLat[c("GenusSpecies","decimalLatitude","decimalLongitude")]
colnames(C) = c("GenusSpecies","Lat","Lon")

# list.files("C:/Users/Owner/Desktop/envBodySize/data/gbif/biodiversitymapping_TIFFs/Birds/richness_10km_all_raster.tif")
# r = raster("C:/Users/Owner/Desktop/envBodySize/data/gbif/biodiversitymapping_TIFFs/Birds/richness_10km_all_raster.tif",values=TRUE)
# library(raster)


extractBirds = function(r,C) {
# newproj = "+proj=eck4 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_def"
newproj = "+proj=longlat +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_def" # set projection coordinates
r2 = projectRaster(r, crs=newproj) # data projected to new coordinates
# rasterToPoints(r2) 

# point = cbind(-100,50) # Lat Long
point = C[c("Lon","Lat")]
b = extract(r2, point) # B for birds
return(b)
}

L = list()
for(i in 1:length(DL)) {
L[[i]] = extractBirds(DL[[i]],C)
}

B = do.call("cbind",L) # B for birds
colnames(B) = Names
Birds = cbind(C,B)

save(Birds,file = "C:/Users/Owner/Desktop/Birds.rda")
}

if(FALSE) { # What is the effect of bird diversity on body size? 
load(file = "C:/Users/Owner/Desktop/envBodySize/data/opdb/quickBodySize.rda") # load cleaned body size data
BS = quickBodySize # BS for body size DATA
BS = BS[c("GenusSpecies","Species","Genus","Family","SubOrder","tbl","hwl")] # drop fwl since often NA
BS = BS[c("GenusSpecies","Species","Genus","Family","SubOrder","hwl")] # drop another variable so that is easy to work with 
BS = na.omit(BS) # remove missing 

load(file="C:/Users/Owner/Desktop/envBodySize/data/gbif/Birds.rda") # load Birds
B = Birds # rename B for Birds

str(B) # what variable are available
B$diversity = B$all # select Bird variable to analyze 

AB = aggregate(diversity ~ GenusSpecies,data=B,mean) # get average birds diversity 
AB$GenusSpecies = gsub(" ","_",AB$GenusSpecies)
D = merge(BS,AB,id="GenusSpecies")

library(ape); load("C:/Users/Owner/Desktop/envBodySize/data/trees/tree.rda") # load tree and ape for pgls

# remove those not in tree or data
D = D[D$GenusSpecies %in% tree$tip.label,] # get only those with tip
D = D[!duplicated(D$GenusSpecies),] # not necessary but lets leave it in 
tip = tree$tip.label[!tree$tip.label %in% D$GenusSpecies] # drop tips without data
tree = drop.tip(tree, tip) # drop missing

# naive analysis
# fit = lm(hwl ~ diversity, data=D)
# summary(fit)

library(nlme)
V = corBrownian(1,phy=tree) 
# rownames(D) = D$GenusSpecies
# V = corPagel(1, phy = tree, form = ~1, fixed = FALSE)
fit = gls(hwl ~ diversity,correlation=V,data=D)
summary(fit)
}
 
if(FALSE) { # What is the effect of bird diversity after we control for climate?

load(file = "C:/Users/Owner/Desktop/envBodySize/data/gbif/climateLat.rda") # load un-cleaned climate data
C = climateLat # C for climate
C$GenusSpecies = gsub(" ","_",C$GenusSpecies) 

AC = aggregate(cbind(bio_max1,bio_max12) ~ GenusSpecies, data = C, mean) # average climate
colnames(AC) = c("GenusSpecies","Temp","Prec")

load(file = "C:/Users/Owner/Desktop/envBodySize/data/opdb/quickBodySize.rda") # load cleaned body size data
BS = quickBodySize # BS for body size DATA
BS = BS[c("GenusSpecies","Species","Genus","Family","SubOrder","tbl","hwl")] # drop fwl since often NA
BS = BS[c("GenusSpecies","Species","Genus","Family","SubOrder","hwl")] # drop another variable so that is easy to work with 
BS = na.omit(BS) # remove missing 

load(file="C:/Users/Owner/Desktop/envBodySize/data/gbif/Birds.rda") # load Birds
B = Birds # rename B for Birds

# str(B) # what variable are available
B$diversity = B$all # select Bird variable to analyze 

AB = aggregate(diversity ~ GenusSpecies,data=B,mean) # get average birds diversity 
AB$GenusSpecies = gsub(" ","_",AB$GenusSpecies)
D = merge(BS,AB,id="GenusSpecies")
D = merge(D,AC) 

library(ape); load("C:/Users/Owner/Desktop/envBodySize/data/trees/tree.rda") # load tree and ape for pgls

# Scale Data
D$Temp = scale(D$Temp)
D$Prec = scale(D$Prec)
D$diversity = scale(D$diversity)

# naive analysis
fit1 = lm(hwl ~ Temp, data=D)
fit2 = lm(hwl ~ diversity, data=D)
fit3 = lm(hwl ~ Prec, data=D)
fit4 = lm(hwl ~ Prec + Temp, data=D)
fit5 = lm(hwl ~ diversity + Temp, data=D)
fit6 = lm(hwl ~ diversity + Temp + Prec, data=D)
fit7 = lm(hwl ~ diversity + Temp*Prec, data=D)
fit8 = lm(hwl ~ diversity*Temp*Prec, data=D)
fit9 = lm(hwl ~ Temp*Prec, data=D)
fit10 = lm(hwl ~ diversity*Temp, data=D)
fit11 = lm(hwl ~ diversity*Prec, data=D)

MS = AIC(fit1,fit2,fit3,fit4,fit5,fit6,fit7,fit8,fit9,fit10,fit11)
MS[order(MS$AIC),]
summary(fit10)

# mixed effect model 
library(nlme)
fit1 = lme(hwl ~ Temp, random = ~ 1|SubOrder/Family,data=D)
fit2 = lme(hwl ~ diversity, random = ~ 1|SubOrder/Family,data=D)
fit3 = lme(hwl ~ Prec,random = ~ 1|SubOrder/Family,data=D)
fit4 = lme(hwl ~ Prec + Temp,random = ~ 1|SubOrder/Family,data=D)
fit5 = lme(hwl ~ diversity + Temp,random = ~ 1|SubOrder/Family,data=D)
fit6 = lme(hwl ~ diversity + Temp + Prec,random = ~ 1|SubOrder/Family,data=D)
fit7 = lme(hwl ~ diversity + Temp*Prec,random = ~ 1|SubOrder/Family,data=D)
fit8 = lme(hwl ~ diversity*Temp*Prec, random = ~ 1|SubOrder/Family,data=D)
fit9 = lme(hwl ~ Temp*Prec, random = ~ 1|SubOrder/Family,data=D)
fit10 = lme(hwl ~ diversity*Temp,random = ~ 1|SubOrder/Family,data=D)
fit11 = lme(hwl ~ diversity*Prec,random = ~ 1|SubOrder/Family,data=D)

MS = AIC(fit1,fit2,fit3,fit4,fit5,fit6,fit7,fit8,fit9,fit10,fit11)
MS[order(MS$AIC),]
summary(fit6)
summary(fit5)
summary(fit10)

# remove those not in tree or data
D = D[D$GenusSpecies %in% tree$tip.label,] # get only those with tip
# D = D[!duplicated(D$GenusSpecies),] # not necessary but lets leave it in 
tip = tree$tip.label[!tree$tip.label %in% D$GenusSpecies] # drop tips without data
tree = drop.tip(tree, tip) # drop missing

nrow(D)
tree


library(nlme)
V = corBrownian(1,phy=tree) 
rownames(D) = D$GenusSpecies
# V = corPagel(1, phy = tree, form = ~1, fixed = FALSE)
fit = gls(hwl ~ diversity + Temp + Prec,correlation=V,data=D)
summary(fit)
fit = gls(tbl ~ diversity + Temp + Prec,correlation=V,data=D)
summary(fit)
}


if(FALSE) { # Process Amphibians data to produce Amphibians.rda

Files = list.files("C:/Users/Owner/Desktop/envBodySize/data/gbif/biodiversitymapping_TIFFs/Amphibians/") # get only GeoTiffs
Files = Files[!grepl("aux|ovr|xml|dbf",Files)]
Files = Files[grepl(".tif",Files)]

library(raster)
library(rgdal)

Paths = paste0("C:/Users/Owner/Desktop/envBodySize/data/gbif/biodiversitymapping_TIFFs/Amphibians/",Files)

DL = list() 
for(i in 1:length(Paths)){
DL[[i]] = raster(Paths[i])
}

Names = gsub("richness_10km_","",Files) # Clean names
Names = gsub("_raster.tif","",Names)
Names

load(file = "C:/Users/Owner/Desktop/envBodySize/data/gbif/climateLat.rda") # load un-cleaned climate data
C = climateLat[c("GenusSpecies","decimalLatitude","decimalLongitude")]
colnames(C) = c("GenusSpecies","Lat","Lon")

# list.files("C:/Users/Owner/Desktop/envBodySize/data/gbif/biodiversitymapping_TIFFs/Birds/richness_10km_all_raster.tif")
# r = raster("C:/Users/Owner/Desktop/envBodySize/data/gbif/biodiversitymapping_TIFFs/Birds/richness_10km_all_raster.tif",values=TRUE)
# library(raster)


extractBirds = function(r,C) { # or extract Amphibians
# newproj = "+proj=eck4 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_def"
newproj = "+proj=longlat +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_def" # set projection coordinates
r2 = projectRaster(r, crs=newproj) # data projected to new coordinates
# rasterToPoints(r2) 

# point = cbind(-100,50) # Lat Long
point = C[c("Lon","Lat")]
b = extract(r2, point) # B for birds
return(b)
}

L = list()
for(i in 1:length(DL)) {
L[[i]] = extractBirds(DL[[i]],C)
}

B = do.call("cbind",L) # B for birds or aphiBians :)
colnames(B) = Names
Amphibians = cbind(C,B)

# save(Amphibians,file = "C:/Users/Owner/Desktop/Amphibians.rda")
 
}

if(FALSE) { # What is the effect of Amphibian diversity after we control for climate?

load(file = "C:/Users/Owner/Desktop/envBodySize/data/gbif/climateLat.rda") # load un-cleaned climate data
C = climateLat # C for climate
C$GenusSpecies = gsub(" ","_",C$GenusSpecies) 

AC = aggregate(cbind(bio_max1,bio_max12) ~ GenusSpecies, data = C, mean) # average climate
colnames(AC) = c("GenusSpecies","Temp","Prec")

load(file = "C:/Users/Owner/Desktop/envBodySize/data/opdb/quickBodySize.rda") # load cleaned body size data
BS = quickBodySize # BS for body size DATA
BS = BS[c("GenusSpecies","Species","Genus","Family","SubOrder","tbl","hwl")] # drop fwl since often NA
BS = BS[c("GenusSpecies","Species","Genus","Family","SubOrder","hwl")] # drop another variable so that is easy to work with 
BS = na.omit(BS) # remove missing 

load(file="C:/Users/Owner/Desktop/envBodySize/data/gbif/Amphibians.rda") # load Birds
B = Amphibians # rename B for Birds

# str(B) # what variable are available
B$diversity = B$all_spp # select Bird variable to analyze 

AB = aggregate(diversity ~ GenusSpecies,data=B,mean) # get average birds diversity 
AB$GenusSpecies = gsub(" ","_",AB$GenusSpecies)
D = merge(BS,AB,id="GenusSpecies")
D = merge(D,AC) 

library(ape); load("C:/Users/Owner/Desktop/envBodySize/data/trees/tree.rda") # load tree and ape for pgls

# Scale Data
D$Temp = scale(D$Temp)
D$Prec = scale(D$Prec)
D$diversity = scale(D$diversity)

# naive analysis
fit1 = lm(hwl ~ Temp, data=D)
fit2 = lm(hwl ~ diversity, data=D)
fit3 = lm(hwl ~ Prec, data=D)
fit4 = lm(hwl ~ Prec + Temp, data=D)
fit5 = lm(hwl ~ diversity + Temp, data=D)
fit6 = lm(hwl ~ diversity + Temp + Prec, data=D)
fit7 = lm(hwl ~ diversity + Temp*Prec, data=D)
fit8 = lm(hwl ~ diversity*Temp*Prec, data=D)
fit9 = lm(hwl ~ Temp*Prec, data=D)
fit10 = lm(hwl ~ diversity*Temp, data=D)
fit11 = lm(hwl ~ diversity*Prec, data=D)

MS = AIC(fit1,fit2,fit3,fit4,fit5,fit6,fit7,fit8,fit9,fit10,fit11)
# MS[order(MS$AIC),]
# summary(fit10)

# mixed effect model 
library(nlme)
fit1 = lme(hwl ~ Temp, random = ~ 1|SubOrder/Family,data=D)
fit2 = lme(hwl ~ diversity, random = ~ 1|SubOrder/Family,data=D)
fit3 = lme(hwl ~ Prec,random = ~ 1|SubOrder/Family,data=D)
fit4 = lme(hwl ~ Prec + Temp,random = ~ 1|SubOrder/Family,data=D)
fit5 = lme(hwl ~ diversity + Temp,random = ~ 1|SubOrder/Family,data=D)
fit6 = lme(hwl ~ diversity + Temp + Prec,random = ~ 1|SubOrder/Family,data=D)
fit7 = lme(hwl ~ diversity + Temp*Prec,random = ~ 1|SubOrder/Family,data=D)
fit8 = lme(hwl ~ diversity*Temp*Prec, random = ~ 1|SubOrder/Family,data=D)
fit9 = lme(hwl ~ Temp*Prec, random = ~ 1|SubOrder/Family,data=D)
fit10 = lme(hwl ~ diversity*Temp,random = ~ 1|SubOrder/Family,data=D)
fit11 = lme(hwl ~ diversity*Prec,random = ~ 1|SubOrder/Family,data=D)

MS = AIC(fit1,fit2,fit3,fit4,fit5,fit6,fit7,fit8,fit9,fit10,fit11)
# MS[order(MS$AIC),]
# summary(fit6)
# summary(fit5)
# summary(fit10)

# remove those not in tree or data
D = D[D$GenusSpecies %in% tree$tip.label,] # get only those with tip
# D = D[!duplicated(D$GenusSpecies),] # not necessary but lets leave it in 
tip = tree$tip.label[!tree$tip.label %in% D$GenusSpecies] # drop tips without data
tree = drop.tip(tree, tip) # drop missing

# nrow(D)
# tree

library(nlme)
V = corBrownian(1,phy=tree) 
rownames(D) = D$GenusSpecies
# V = corPagel(1, phy = tree, form = ~1, fixed = FALSE)
fit = gls(hwl ~ diversity + Temp + Prec,correlation=V,data=D)
AIC(fit)
summary(fit)

# fit = gls(tbl ~ diversity + Temp + Prec,correlation=V,data=D)
# summary(fit)
}

if(FALSE) { # What is the effect of bird diversity after we control for Amphibian diversity? It stays strong
load(file = "C:/Users/Owner/Desktop/envBodySize/data/gbif/climateLat.rda") # load un-cleaned climate data
C = climateLat # C for climate
C$GenusSpecies = gsub(" ","_",C$GenusSpecies) 

AC = aggregate(cbind(bio_max1,bio_max12) ~ GenusSpecies, data = C, mean) # average climate
colnames(AC) = c("GenusSpecies","Temp","Prec")

load(file = "C:/Users/Owner/Desktop/envBodySize/data/opdb/quickBodySize.rda") # load cleaned body size data
BS = quickBodySize # BS for body size DATA
BS = BS[c("GenusSpecies","Species","Genus","Family","SubOrder","tbl","hwl")] # drop fwl since often NA
BS = BS[c("GenusSpecies","Species","Genus","Family","SubOrder","hwl")] # drop another variable so that is easy to work with 
BS = na.omit(BS) # remove missing 

# load biotic grids 
load(file="C:/Users/Owner/Desktop/envBodySize/data/gbif/Amphibians.rda") # 
A = Amphibians # rename A for Amphibians

load(file="C:/Users/Owner/Desktop/envBodySize/data/gbif/Birds.rda") # load Birds
B = Birds # rename B for Birds

# str(B) # what variable are available
B$div_brd = B$all # select Amphibians variable to analyze 
A$div_amp = A$all_spp # select Bird variable to analyze 

AB = aggregate(div_brd ~ GenusSpecies,data=B,mean) # get average birds diversity 
AB$GenusSpecies = gsub(" ","_",AB$GenusSpecies)
D = merge(BS,AB,id="GenusSpecies")
AA = aggregate(div_amp ~ GenusSpecies,data=A,mean) # get average birds diversity 
AA$GenusSpecies = gsub(" ","_",AA$GenusSpecies)
D = merge(D,AA,id="GenusSpecies")

D = merge(D,AC,id="GenusSpecies") # merge diversity with average climates

library(ape); load("C:/Users/Owner/Desktop/envBodySize/data/trees/tree.rda") # load tree and ape for pgls

# Scale Data
D$Temp = scale(D$Temp)
D$Prec = scale(D$Prec)
D$div_amp = scale(D$div_amp)
D$div_brd = scale(D$div_brd)

# naive analysis
fit1 = lm(hwl ~ Temp, data=D)
fit2 = lm(hwl ~ div_amp, data=D)
fit3 = lm(hwl ~ div_brd, data=D)
fit4 = lm(hwl ~ div_brd + div_amp, data=D)
fit5 = lm(hwl ~ div_brd + div_amp + Temp, data=D)
fit6 = lm(hwl ~ div_brd + div_amp + Temp + Prec, data=D)
fit7 = lm(hwl ~ div_brd + div_amp + Temp*Prec, data=D)
fit8 = lm(hwl ~ div_brd*div_amp + Temp*Prec, data=D)
fit9 = lm(hwl ~ div_brd*div_amp + Temp*Prec, data=D)
fit10 = lm(hwl ~ Temp*Prec, data=D)
fit11 = lm(hwl ~ Temp + Prec, data=D)
fit12 = lm(hwl ~ Prec, data=D)

MS = AIC(fit1,fit2,fit3,fit4,fit5,fit6,fit7,fit8,fit9,fit10,fit11,fit12)
# MS[order(MS$AIC),]


# mixed effect model 
library(nlme)
fit1 = lme(hwl ~ div_brd, random = ~ 1|SubOrder/Family,data=D)
fit2 = lme(hwl ~ div_amp, random = ~ 1|SubOrder/Family, data=D)
fit3 = lme(hwl ~ div_brd, random = ~ 1|SubOrder/Family, data=D)
fit4 = lme(hwl ~ div_brd + div_amp, random = ~ 1|SubOrder/Family, data=D)
fit5 = lme(hwl ~ div_brd + div_amp + Temp, random = ~ 1|SubOrder/Family, data=D)
fit6 = lme(hwl ~ div_brd + div_amp + Temp + Prec, random = ~ 1|SubOrder/Family, data=D)
fit7 = lme(hwl ~ div_brd + div_amp + Temp*Prec, random = ~ 1|SubOrder/Family, data=D)
fit8 = lme(hwl ~ div_brd*div_amp + Temp*Prec, random = ~ 1|SubOrder/Family, data=D)
fit9 = lme(hwl ~ div_brd*div_amp + Temp*Prec, random = ~ 1|SubOrder/Family, data=D)
fit10 = lme(hwl ~ Temp*Prec, random = ~ 1|SubOrder/Family, data=D)
fit11 = lme(hwl ~ Temp + Prec, random = ~ 1|SubOrder/Family, data=D)
fit12 = lme(hwl ~ Prec, random = ~ 1|SubOrder/Family, data=D)

MS = AIC(fit1,fit2,fit3,fit4,fit5,fit6,fit7,fit8,fit9,fit10,fit11,fit12)
MS[order(MS$AIC),]

# remove those not in tree or data
D = D[D$GenusSpecies %in% tree$tip.label,] # get only those with tip
# D = D[!duplicated(D$GenusSpecies),] # not necessary but lets leave it in 
tip = tree$tip.label[!tree$tip.label %in% D$GenusSpecies] # drop tips without data
tree = drop.tip(tree, tip) # drop missing

# nrow(D)
# tree

library(nlme)
V = corBrownian(1,phy=tree) 
rownames(D) = D$GenusSpecies
# V = corPagel(1, phy = tree, form = ~1, fixed = FALSE)
fit = gls(hwl ~ div_brd + div_amp + Temp + Prec,correlation=V,data=D)
summary(fit)
}


if(FALSE) { # Process Mammals data to produce Mammals.rda

Files = list.files("C:/Users/Owner/Desktop/envBodySize/data/gbif/biodiversitymapping_TIFFs/Mammals/") # get only GeoTiffs
Files = Files[!grepl("aux|ovr|xml|dbf",Files)]
Files = Files[grepl(".tif",Files)]

library(raster)
library(rgdal)

Paths = paste0("C:/Users/Owner/Desktop/envBodySize/data/gbif/biodiversitymapping_TIFFs/Mammals/",Files)

DL = list() 
for(i in 1:length(Paths)){
DL[[i]] = raster(Paths[i])
}

Names = gsub("richness_10km_","",Files) # Clean names
Names = gsub("_raster.tif","",Names)
Names

load(file = "C:/Users/Owner/Desktop/envBodySize/data/gbif/climateLat.rda") # load un-cleaned climate data
C = climateLat[c("GenusSpecies","decimalLatitude","decimalLongitude")]
colnames(C) = c("GenusSpecies","Lat","Lon")

# list.files("C:/Users/Owner/Desktop/envBodySize/data/gbif/biodiversitymapping_TIFFs/Birds/richness_10km_all_raster.tif")
# r = raster("C:/Users/Owner/Desktop/envBodySize/data/gbif/biodiversitymapping_TIFFs/Birds/richness_10km_all_raster.tif",values=TRUE)
# library(raster)


extractBirds = function(r,C) { # or extract Amphibians
# newproj = "+proj=eck4 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_def"
newproj = "+proj=longlat +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_def" # set projection coordinates
r2 = projectRaster(r, crs=newproj) # data projected to new coordinates
# rasterToPoints(r2) 

# point = cbind(-100,50) # Lat Long
point = C[c("Lon","Lat")]
b = extract(r2, point) # B for birds
return(b)
}

L = list()
for(i in 1:length(DL)) {
L[[i]] = extractBirds(DL[[i]],C)
}

B = do.call("cbind",L) # B for birds or aphiBians :)
colnames(B) = Names
Mammals = cbind(C,B)

# save(Mammals,file = "C:/Users/Owner/Desktop/Mammals.rda")
 
}

if(FALSE) { # What is the effect of Mammals diversity after we control for climate?

load(file = "C:/Users/Owner/Desktop/envBodySize/data/gbif/climateLat.rda") # load un-cleaned climate data
C = climateLat # C for climate
C$GenusSpecies = gsub(" ","_",C$GenusSpecies) 

AC = aggregate(cbind(bio_max1,bio_max12) ~ GenusSpecies, data = C, mean) # average climate
colnames(AC) = c("GenusSpecies","Temp","Prec")

load(file = "C:/Users/Owner/Desktop/envBodySize/data/opdb/quickBodySize.rda") # load cleaned body size data
BS = quickBodySize # BS for body size DATA
BS = BS[c("GenusSpecies","Species","Genus","Family","SubOrder","tbl","hwl")] # drop fwl since often NA
BS = BS[c("GenusSpecies","Species","Genus","Family","SubOrder","hwl")] # drop another variable so that is easy to work with 
BS = na.omit(BS) # remove missing 

load(file="C:/Users/Owner/Desktop/envBodySize/data/gbif/Mammals.rda") # load Birds
B = Mammals # rename B for Birds

# str(B) # what variable are available
B$diversity = B$all_spp # select Bird variable to analyze 

AB = aggregate(diversity ~ GenusSpecies,data=B,mean) # get average birds diversity 
AB$GenusSpecies = gsub(" ","_",AB$GenusSpecies)
D = merge(BS,AB,id="GenusSpecies")
D = merge(D,AC) 

library(ape); load("C:/Users/Owner/Desktop/envBodySize/data/trees/tree.rda") # load tree and ape for pgls

# Scale Data
D$Temp = scale(D$Temp)
D$Prec = scale(D$Prec)
D$diversity = scale(D$diversity)

# naive analysis
fit1 = lm(hwl ~ Temp, data=D)
fit2 = lm(hwl ~ diversity, data=D)
fit3 = lm(hwl ~ Prec, data=D)
fit4 = lm(hwl ~ Prec + Temp, data=D)
fit5 = lm(hwl ~ diversity + Temp, data=D)
fit6 = lm(hwl ~ diversity + Temp + Prec, data=D)
fit7 = lm(hwl ~ diversity + Temp*Prec, data=D)
fit8 = lm(hwl ~ diversity*Temp*Prec, data=D)
fit9 = lm(hwl ~ Temp*Prec, data=D)
fit10 = lm(hwl ~ diversity*Temp, data=D)
fit11 = lm(hwl ~ diversity*Prec, data=D)

MS = AIC(fit1,fit2,fit3,fit4,fit5,fit6,fit7,fit8,fit9,fit10,fit11)
# MS[order(MS$AIC),]
# summary(fit10)

# mixed effect model 
library(nlme)
fit1 = lme(hwl ~ Temp, random = ~ 1|SubOrder/Family,data=D)
fit2 = lme(hwl ~ diversity, random = ~ 1|SubOrder/Family,data=D)
fit3 = lme(hwl ~ Prec,random = ~ 1|SubOrder/Family,data=D)
fit4 = lme(hwl ~ Prec + Temp,random = ~ 1|SubOrder/Family,data=D)
fit5 = lme(hwl ~ diversity + Temp,random = ~ 1|SubOrder/Family,data=D)
fit6 = lme(hwl ~ diversity + Temp + Prec,random = ~ 1|SubOrder/Family,data=D)
fit7 = lme(hwl ~ diversity + Temp*Prec,random = ~ 1|SubOrder/Family,data=D)
fit8 = lme(hwl ~ diversity*Temp*Prec, random = ~ 1|SubOrder/Family,data=D)
fit9 = lme(hwl ~ Temp*Prec, random = ~ 1|SubOrder/Family,data=D)
fit10 = lme(hwl ~ diversity*Temp,random = ~ 1|SubOrder/Family,data=D)
fit11 = lme(hwl ~ diversity*Prec,random = ~ 1|SubOrder/Family,data=D)

MS = AIC(fit1,fit2,fit3,fit4,fit5,fit6,fit7,fit8,fit9,fit10,fit11)
# MS[order(MS$AIC),]
# summary(fit6)
# summary(fit5)
# summary(fit10)

# remove those not in tree or data
D = D[D$GenusSpecies %in% tree$tip.label,] # get only those with tip
# D = D[!duplicated(D$GenusSpecies),] # not necessary but lets leave it in 
tip = tree$tip.label[!tree$tip.label %in% D$GenusSpecies] # drop tips without data
tree = drop.tip(tree, tip) # drop missing

# nrow(D)
# tree

library(nlme)
V = corBrownian(1,phy=tree) 
rownames(D) = D$GenusSpecies
# V = corPagel(1, phy = tree, form = ~1, fixed = FALSE)
fit = gls(hwl ~ diversity + Temp + Prec,correlation=V,data=D)
AIC(fit)
summary(fit)

# fit = gls(tbl ~ diversity + Temp + Prec,correlation=V,data=D)
# summary(fit)
}

if(FALSE) { # What is the effect of bird diversity after we control for Mammals diversity? It stays strong

# list.files("C:/Users/Owner/Desktop/envBodySize/data/")
load("C:/Users/Owner/Desktop/envBodySize/data/gbif/trees.rda")
T = trees # T for Tree Cover
# T = T[!T$tree == 0,] # remove zero value
T$GenusSpecies = gsub(" ","_",T$GenusSpecies)
AT = aggregate(trees ~ GenusSpecies, data = T, mean) # AT average tree cover
colnames(AT) = c("GenusSpecies","cover")

load(file = "C:/Users/Owner/Desktop/envBodySize/data/gbif/climateLat.rda") # load un-cleaned climate data
C = climateLat # C for climate
C$GenusSpecies = gsub(" ","_",C$GenusSpecies) 

AC = aggregate(cbind(bio_max1,bio_max12) ~ GenusSpecies, data = C, mean) # average climate
colnames(AC) = c("GenusSpecies","Temp","Prec")

load(file = "C:/Users/Owner/Desktop/envBodySize/data/opdb/quickBodySize.rda") # load cleaned body size data
BS = quickBodySize # BS for body size DATA
BS = BS[c("GenusSpecies","Species","Genus","Family","SubOrder","tbl","hwl")] # drop fwl since often NA
BS = BS[c("GenusSpecies","Species","Genus","Family","SubOrder","hwl")] # drop another variable so that is easy to work with 
BS = na.omit(BS) # remove missing 

# load biotic grids 
load(file="C:/Users/Owner/Desktop/envBodySize/data/gbif/Mammals.rda") # 
A = Mammals # rename A for Amphibians

load(file="C:/Users/Owner/Desktop/envBodySize/data/gbif/Birds.rda") # load Birds
B = Birds # rename B for Birds

# str(B) # what variable are available
B$div_brd = B$all # select Amphibians variable to analyze 
A$div_mam = A$all_spp # select Bird variable to analyze 

AB = aggregate(div_brd ~ GenusSpecies,data=B,mean) # get average birds diversity 
AB$GenusSpecies = gsub(" ","_",AB$GenusSpecies)
D = merge(BS,AB,id="GenusSpecies")
AA = aggregate(div_mam ~ GenusSpecies,data=A,mean) # get average birds diversity 
AA$GenusSpecies = gsub(" ","_",AA$GenusSpecies)

D = merge(D,AA,id="GenusSpecies")
D = merge(D,AC,id="GenusSpecies") # merge diversity with average climates
D = merge(D,AT,id="GenusSpecies") # merge diversity with average climates


library(ape); load("C:/Users/Owner/Desktop/envBodySize/data/trees/tree.rda") # load tree and ape for pgls

# Scale Data
D$Temp = scale(D$Temp)
D$Prec = scale(D$Prec)
D$div_mam = scale(D$div_mam)
D$div_brd = scale(D$div_brd)
D$cover = scale(D$cover)

# naive analysis
fit1 = lm(hwl ~ Temp, data=D)
fit2 = lm(hwl ~ div_mam, data=D)
fit3 = lm(hwl ~ div_brd, data=D)
fit4 = lm(hwl ~ div_brd + div_mam, data=D)
fit5 = lm(hwl ~ div_brd + div_mam + Temp, data=D)
fit6 = lm(hwl ~ div_brd + div_mam + Temp + Prec, data=D)
fit7 = lm(hwl ~ div_brd + div_mam + Temp*Prec, data=D)
fit8 = lm(hwl ~ div_brd*div_mam + Temp*Prec, data=D)
fit9 = lm(hwl ~ div_brd*div_mam + Temp*Prec, data=D)
fit10 = lm(hwl ~ Temp*Prec, data=D)
fit11 = lm(hwl ~ Temp + Prec, data=D)
fit12 = lm(hwl ~ Prec, data=D)

MS = AIC(fit1,fit2,fit3,fit4,fit5,fit6,fit7,fit8,fit9,fit10,fit11,fit12)
# MS[order(MS$AIC),]


# mixed effect model 
library(nlme)
fit1 = lme(hwl ~ div_brd, random = ~ 1|SubOrder/Family,data=D)
fit2 = lme(hwl ~ div_mam, random = ~ 1|SubOrder/Family, data=D)
fit3 = lme(hwl ~ div_brd, random = ~ 1|SubOrder/Family, data=D)
fit4 = lme(hwl ~ div_brd + div_mam, random = ~ 1|SubOrder/Family, data=D)
fit5 = lme(hwl ~ div_brd + div_mam + Temp, random = ~ 1|SubOrder/Family, data=D)
fit6 = lme(hwl ~ div_brd + div_mam + Temp + Prec, random = ~ 1|SubOrder/Family, data=D)
fit7 = lme(hwl ~ div_brd + div_mam + Temp*Prec, random = ~ 1|SubOrder/Family, data=D)
fit8 = lme(hwl ~ div_brd*div_mam + Temp*Prec, random = ~ 1|SubOrder/Family, data=D)
fit9 = lme(hwl ~ div_brd*div_mam + Temp*Prec, random = ~ 1|SubOrder/Family, data=D)
fit10 = lme(hwl ~ Temp*Prec, random = ~ 1|SubOrder/Family, data=D)
fit11 = lme(hwl ~ Temp + Prec, random = ~ 1|SubOrder/Family, data=D)
fit12 = lme(hwl ~ Prec, random = ~ 1|SubOrder/Family, data=D)

MS = AIC(fit1,fit2,fit3,fit4,fit5,fit6,fit7,fit8,fit9,fit10,fit11,fit12)
# MS[order(MS$AIC),]

# remove those not in tree or data
D = D[D$GenusSpecies %in% tree$tip.label,] # get only those with tip
# D = D[!duplicated(D$GenusSpecies),] # not necessary but lets leave it in 
tip = tree$tip.label[!tree$tip.label %in% D$GenusSpecies] # drop tips without data
tree = drop.tip(tree, tip) # drop missing

# nrow(D)
# tree

library(nlme)
V = corBrownian(1,phy=tree) 
rownames(D) = D$GenusSpecies
# V = corPagel(1, phy = tree, form = ~1, fixed = FALSE)
# fit = gls(hwl ~ div_brd + div_mam + Temp + Prec,correlation=V,data=D)
fit = gls(hwl ~ Temp + Prec + cover + div_brd + div_mam,correlation=V,data=D)
# fit = gls(hwl ~ Temp,correlation=V,data=D)


# clean gls summary info
AIC = summary(fit)$AIC
N = nrow(D)
tt = as.data.frame(summary(fit)$tTable)

p = tt$"p-value"
p = p[2:nrow(tt)] # drop intercept
p = round(p,4) # round
# if(p == 0) p = "< 0.0001"
p = paste("P =", p)

ce = rownames(tt)[2:nrow(tt)] # coefficients

val = round(tt$Value[2:nrow(tt)],3) # drop intercept
val = paste(ce,"=",val,collapse=" ;")

se = round(tt$Std.Error[2:nrow(tt)],3)
se = paste("SE =",se)

N
val
se
p
AIC

}

if(FALSE) { # Process Binary Tree cover to matrix
library(caTools)
X = read.ENVI("C:/Users/Owner/Desktop/envBodySize/data/gbif/treeCover/tree.img")# make sure to have tree.img.hdr file in same directory!!
# save(X,file = "C:/Users/Owner/Desktop/X.rda")
}

if(FALSE) { # convert tree cover X matrix to r raster
load("C:/Users/Owner/Desktop/envBodySize/data/gbif/treeCover/X.rda")
library(raster)
r = raster(X)
# save(r,file = "C:/Users/Owner/Desktop/r.rda")
}

if(FALSE) { # assign spatial extent
library(raster)

load("C:/Users/Owner/Desktop/envBodySize/data/gbif/treeCover/r.rda")
e = extent(-180, 180, -90, 90)
extent(r) = e
save(r, file = "C:/Users/Owner/Desktop/r.rda") # save here to prevent overwriting
}

if(FALSE) { # extract tree cover data
library(raster)
load("C:/Users/Owner/Desktop/envBodySize/data/gbif/treeCover/r.rda")

load(file = "C:/Users/Owner/Desktop/envBodySize/data/gbif/climateLat.rda") # load un-cleaned climate data
C = climateLat[c("GenusSpecies","decimalLatitude","decimalLongitude")]
colnames(C) = c("GenusSpecies","Lat","Lon")

# point = cbind(55,13)
# point = cbind(-100,50) # Lon Lat
# point = cbind(-100,50) # Lat Long

# point = cbind(-52,-9.24) # Lon Lat # is in reverse order thanks hjmajiz!
# point = cbind(14.88,59.946) # Lon Lat # is in reverse order thanks hjmajiz!
# point = cbind(14.88,59.946) # Lon Lat # is in reverse order thanks hjmajiz!

point = C[c("Lon","Lat")]
# point = cbind(-165.684746, 18.152153)

# jpeg("C:/Users/Owner/Desktop/plot.jpg",units = "in",width=15,height=10,res=600)
# plot(r,main="Tree Cover")
# points(C$Lon,C$Lat,cex=0.1,pch=19,col="red")
# dev.off()
trees = extract(r, point) # B for birds
trees = cbind(C,trees)

# Value	Label
# 0-100	Percent of pixel area covered by tree cover
# 200	Water
# 210	Cloud
# 211	Shadow
# 220	Filled Value

trees = trees[!trees$trees == 200,] # clean up non-data values
trees = trees[!trees$trees == 210,]
trees = trees[!trees$trees == 210,]
trees = trees[!trees$trees == 220,]
trees = trees[!trees$trees == 211,]

# save(trees,file="C:/Users/Owner/Desktop/trees.rda") # save in neutral location
}

if(FALSE) { # How related is bird diversity to tree cover? 

# list.files("C:/Users/Owner/Desktop/envBodySize/data/")
load("C:/Users/Owner/Desktop/envBodySize/data/gbif/trees.rda")
T = trees # T for Tree Cover

load(file="C:/Users/Owner/Desktop/envBodySize/data/gbif/Birds.rda") # load Birds
B = Birds # rename B for Birds

D = merge(B,T,id="GenusSpecies")
D = aggregate(cbind(all,trees) ~ GenusSpecies,data=D,mean)

cor.test(D$all,D$trees)
summary(lm(all ~ trees,data=D))

# jpeg("C:/Users/Owner/Desktop/plot.jpg",units = "in",width=15,height=10,res=300)
# par(mfrow=c(1,2))
# hist(D$trees)
# hist(D$all)
# plot(D$trees,D$all,main="Trees Vs Birds")
# points(C$Lon,C$Lat,cex=0.1,pch=19,col="red")
# dev.off()

# load(file = "C:/Users/Owner/Desktop/envBodySize/data/gbif/climateLat.rda") # load un-cleaned climate data
# C = climateLat # C for climate
# C$GenusSpecies = gsub(" ","_",C$GenusSpecies) 
}

if(FALSE) { # Is Bird diversity related to body size after controlling for tree density? 

# list.files("C:/Users/Owner/Desktop/envBodySize/data/")
load("C:/Users/Owner/Desktop/envBodySize/data/gbif/trees.rda")
T = trees # T for Tree Cover
# T = T[!T$tree == 0,] # remove zero value
T$GenusSpecies = gsub(" ","_",T$GenusSpecies)
AT = aggregate(trees ~ GenusSpecies, data = T, mean) # AT average tree cover
colnames(AT) = c("GenusSpecies","cover")

load(file="C:/Users/Owner/Desktop/envBodySize/data/gbif/Birds.rda") # load Birds
B = Birds # rename B for Birds
B$GenusSpecies = gsub(" ","_",B$GenusSpecies)
AB = aggregate(all ~ GenusSpecies, data = B, mean) # AB average birds
colnames(AB) = c("GenusSpecies","diversity")

load(file = "C:/Users/Owner/Desktop/envBodySize/data/gbif/climateLat.rda") # load un-cleaned climate data
C = climateLat # C for climate
C$GenusSpecies = gsub(" ","_",C$GenusSpecies) 
AC = aggregate(cbind(bio_max1,bio_max12) ~ GenusSpecies, data = C, mean) # average climate
colnames(AC) = c("GenusSpecies","Temp","Prec")

load(file = "C:/Users/Owner/Desktop/envBodySize/data/opdb/quickBodySize.rda") # load cleaned body size data
BS = quickBodySize # BS for body size DATA
BS = BS[c("GenusSpecies","Species","Genus","Family","SubOrder","tbl","hwl")] # drop fwl since often NA
BS = BS[c("GenusSpecies","Species","Genus","Family","SubOrder","hwl")] # drop another variable so that is easy to work with 
BS = na.omit(BS) # remove missing 

D = merge(BS,AC,id="GenusSpecies")
D = merge(D,AT,id="GenusSpecies")
D = merge(D,AB,id="GenusSpecies")

library(ape); load("C:/Users/Owner/Desktop/envBodySize/data/trees/tree.rda") # load tree and ape for pgls

# Scale Data
D$Temp = scale(D$Temp)
D$Prec = scale(D$Prec)
D$diversity = scale(D$diversity)
D$cover = scale(D$cover)

# naive analysis
fit1 = lm(hwl ~ cover, data=D)
fit2 = lm(hwl ~ diversity, data=D)
fit3 = lm(hwl ~ Temp, data=D)
fit4 = lm(hwl ~ Prec, data=D)
fit5 = lm(hwl ~ diversity + cover, data=D)
fit6 = lm(hwl ~ diversity*cover, data=D)
fit7 = lm(hwl ~ diversity + cover + Temp, data=D)
fit8 = lm(hwl ~ diversity + cover + Temp + Prec, data=D)
fit9 = lm(hwl ~ diversity + cover + Temp*Prec, data=D)
fit10 = lm(hwl ~ diversity + cover*Temp*Prec, data=D)
fit11 = lm(hwl ~ diversity*cover*Temp*Prec, data=D)
fit12 = lm(hwl ~ diversity*cover + Temp*Prec, data=D)
fit13 = lm(hwl ~ diversity + Temp*Prec, data=D)
fit14 = lm(hwl ~ diversity*Temp*Prec, data=D)
fit15 = lm(hwl ~ diversity*Temp + Prec, data=D)
fit15 = lm(hwl ~ diversity*cover + Temp, data=D)
fit16 = lm(hwl ~ diversity*cover*Temp, data=D)
fit17 = lm(hwl ~ SubOrder, data=D)
fit18 = lm(hwl ~ SubOrder*diversity, data=D)

MS = AIC(fit1,fit2,fit3,fit4,fit5,fit6,fit7,fit8,fit9,fit10,fit11,fit12,fit13,fit14,fit15,fit16,fit17,fit18)
MS[order(MS$AIC),]
summary(fit18)


# mixed effect model 
library(nlme)
fit1 = lme(hwl ~ cover,random = ~ 1|SubOrder/Family,data=D)
fit2 = lme(hwl ~ diversity,random = ~ 1|SubOrder/Family,data=D)
fit3 = lme(hwl ~ Temp,random = ~ 1|SubOrder/Family,data=D)
fit4 = lme(hwl ~ Prec,random = ~ 1|SubOrder/Family,data=D)
fit5 = lme(hwl ~ diversity + cover,random = ~ 1|SubOrder/Family,data=D)
fit6 = lme(hwl ~ diversity*cover,random = ~ 1|SubOrder/Family,data=D)
fit7 = lme(hwl ~ cover + Temp,random = ~ 1|SubOrder/Family,data=D)
fit8 = lme(hwl ~ diversity + cover + Temp + Prec,random = ~ 1|SubOrder/Family,data=D)
fit9 = lme(hwl ~ diversity + cover + Temp*Prec,random = ~ 1|SubOrder/Family,data=D)
fit10 = lme(hwl ~ diversity + cover*Temp*Prec,random = ~ 1|SubOrder/Family,data=D)
fit11 = lme(hwl ~ diversity*cover*Temp*Prec,random = ~ 1|SubOrder/Family,data=D)
fit12 = lme(hwl ~ diversity*cover + Temp*Prec,random = ~ 1|SubOrder/Family,data=D)
fit13 = lme(hwl ~ diversity + Temp*Prec,random = ~ 1|SubOrder/Family,data=D)
fit14 = lme(hwl ~ diversity*Temp*Prec,random = ~ 1|SubOrder/Family,data=D)
fit15 = lme(hwl ~ diversity*Temp + Prec,random = ~ 1|SubOrder/Family,data=D)
fit15 = lme(hwl ~ diversity*cover + Temp,random = ~ 1|SubOrder/Family,data=D)
fit16 = lme(hwl ~ diversity*cover*Temp,random = ~ 1|SubOrder/Family,data=D)

MS = AIC(fit1,fit2,fit3,fit4,fit5,fit6,fit7,fit8,fit9,fit10,fit11,fit12,fit13,fit14,fit15,fit16)
MS[order(MS$AIC),]
summary(fit8)

# remove those not in tree or data
D = D[D$GenusSpecies %in% tree$tip.label,] # get only those with tip
# D = D[!duplicated(D$GenusSpecies),] # not necessary but lets leave it in 
tip = tree$tip.label[!tree$tip.label %in% D$GenusSpecies] # drop tips without data
tree = drop.tip(tree, tip) # drop missing

library(nlme)
V = corBrownian(1,phy=tree) 
rownames(D) = D$GenusSpecies
# V = corPagel(1, phy = tree, form = ~1, fixed = FALSE)
# fit = gls(hwl ~ diversity + Temp + cover,correlation=V,data=D)
fit = gls(hwl ~ diversity + cover + Temp + Prec,correlation=V,data=D)
# fit = gls(hwl ~ Temp + cover,correlation=V,data=D)
summary(fit)

}


# Consider using mammal diversity as a control for biodiversity hotspots


if(FALSE) { # Combine Body Size, Taxonomic Opinions,and Paleo-latitude --> FossilsLatTax.rda

D = read.table("C:/Users/Owner/Desktop/envBodySize/data/Fossils/pbdb_occs.csv",sep=",",header=TRUE)
# D = D[c("taxon_name","collection_no","early_age","late_age","paleolat")]

cn = unique(D$collection_no)
DL = list()
for(i in 1:length(cn)) {
paleolat = unique(D[D$collection_no == cn[i],]$paleolat)
collection_no = cn[i]
DL[[i]] = data.frame(collection_no,paleolat)
}

D = do.call("rbind",DL)

load("C:/Users/Owner/Desktop/envBodySize/data/Fossils/Fossils.rda")
F = Fossils
F = F[!is.na(F$collection_no),]

D = merge(F,D,id="collection_no")


fg = unique(D[D$type2 == "fossil",]$Genus)

O = read.table("C:/Users/Owner/Desktop/envBodySize/data/Fossils/opinions.csv",header=TRUE,sep=",",fill=TRUE) # O for opinions

DL = list()
for(i in 1:length(fg)) {
ffg = fg[i] # iterated fossil genus
d = O[O$child_name == ffg,] # sub data.frame

# print(d)
if(nrow(d) != 0) {
n_opinions = nrow(d)
d = d[nrow(d),]
child_name = d$child_name
parent_name = d$parent_name
status = d$status
} else {
child_name = NA
parent_name = NA
status = NA
n_opinions = NA
}

DL[[i]] = data.frame(fg = ffg,child_name=child_name,parent_name=parent_name,status=status,n_opinions=n_opinions)
}

FT = do.call("rbind",DL) # Fossil Taxonomy
colnames(FT) = c("fg","Genus","Family","status","n_opinions")
FT = FT[c("fg","Family","Genus")]

if(TRUE) { # Corrections...
FT = FT[order(FT$Family),]
FT$Family = as.character(FT$Family)
FT$Genus = as.character(FT$Genus)
FT$Genus[FT$fg == "Heterophlebia"] = "Heterophlebia"
FT$Family[FT$Genus == "Heterophlebia"] = "Heterophlebiidae"

FT = FT[c("Family","Genus")]
FT$SubOrder = NA
FT = FT[c("SubOrder","Family","Genus")]
FT$Family[FT$Family == "Odonata"] = NA
FT$SubOrder[FT$Family == "Zygoptera"] = "Zygoptera"
FT$Family[FT$Family == "Zygoptera"] = NA

FT$Family[FT$Family == "Libellula (Heterophlebia)"] = "Libellulidae"
FT$SubOrder[FT$Family == "Aeschnidiidae"] = "Anisoptera"
FT$Family[FT$Family == "Aeshninae"] = "Aeshnidae"
FT$SubOrder[FT$Family == "Aeshnidae"] = "Anisoptera"
FT$SubOrder[FT$Family == "Agrionidae"] = "Zygoptera"
FT$SubOrder[FT$Family == "Aktassiinae"] = "Anisoptera"
FT$SubOrder[FT$Family == "Anisoptera"] = "Anisoptera"
FT$SubOrder[FT$Family == "Anisozygoptera"] = "Anisozygoptera"
FT$Family[FT$Family == "Anisozygoptera"] = NA
FT$SubOrder[FT$Family == "Araripelibellulidae"] = "Anisoptera"
FT$SubOrder[FT$Family == "Architheminae"] = "Anisozygoptera"
FT$SubOrder[FT$Family == "Bolcacorduliidae"] = "Anisoptera"
FT$SubOrder[FT$Family == "Bolcathoridae"] = "Zygoptera"
FT$Family[FT$Family == "Parabrachydiplax"] = "Libellulidae" # http://fossilworks.org/cgi-bin/bridge.pl?a=taxonInfo&taxon_no=178765
FT$SubOrder[FT$Family == "Asiopteridae"] = "Anisoptera" # http://fossilworks.org/bridge.pl?a=taxonInfo&taxon_no=188917
FT$SubOrder[FT$Family == "Brachydiplacini"] = "Anisoptera" # http://fossilworks.org/cgi-bin/bridge.pl?a=referenceInfo&reference_no=35360
FT$SubOrder[FT$Family == "Campterophlebiidae"] = "Anisoptera" # http://www.tandfonline.com/doi/full/10.1080/03115518.2017.1318170
FT$SubOrder[FT$Family == "Cordulagomphinae"] = "Anisoptera" # http://fossilworks.org/?a=taxonInfo&taxon_no=175158
FT$SubOrder[FT$Family == "Corduliidae"] = "Anisoptera" # present
FT$SubOrder[FT$Family == "Corduliinae"] = "Anisoptera" # http://fossilworks.org/bridge.pl?a=taxonInfo&taxon_no=180523
FT$SubOrder[FT$Family == "Corduliinae"] = "Anisoptera" # http://fossilworks.org/bridge.pl?a=taxonInfo&taxon_no=180523
FT$SubOrder[FT$Family == "Cretacoenagrionidae"] = "Zygoptera" # http://fossilworks.org/bridge.pl?a=taxonInfo&taxon_no=181187
FT$SubOrder[FT$Family == "Cretapetaluridae"] = "Anisoptera" # http://fossilworks.org/bridge.pl?a=taxonInfo&taxon_no=221181
FT$SubOrder[FT$Family == "Cymatophlebiina"] = "Anisoptera" # http://fossilworks.org/bridge.pl?a=taxonInfo&taxon_no=178397
FT$SubOrder[FT$Family == "Enigmaeshnidae"] = "Anisoptera" # http://fossilworks.org/bridge.pl?a=taxonInfo&taxon_no=177573
FT$SubOrder[FT$Family == "Eodichromini"] = "Zygoptera" # http://fossilworks.org/bridge.pl?a=taxonInfo&taxon_no=177573
FT$SubOrder[FT$Family == "Epallaginae"] = "Zygoptera" # 
FT$SubOrder[FT$Family == "Eumorbaeschnidae"] = "Anisoptera" # 
FT$SubOrder[FT$Family == "Euzygoptera"] = "Zygoptera" # 
FT$SubOrder[FT$Family == "Gallophlebiidae"] = "Anisoptera" # 
FT$SubOrder[FT$Family == "Gomphaeschnaoidini"] = "Anisoptera" # 
FT$SubOrder[FT$Family == "Gomphaeschnidae"] = "Anisoptera" # 
FT$SubOrder[FT$Family == "Gomphaeschninae"] = "Anisoptera" # 
FT$SubOrder[FT$Family == "Gomphidae"] = "Anisoptera" # 
FT$SubOrder[FT$Family == "Gomphini"] = "Anisoptera" # 
FT$SubOrder[FT$Family == "Heterophlebiidae"] = "Anisoptera" # 
FT$SubOrder[FT$Family == "Valdicorduliidae"] = "Anisoptera" # 
FT$SubOrder[FT$Family == "Valdaeshninae"] = "Anisoptera" # 
FT$SubOrder[FT$Family == "Triassolestidae"] = "Anisoptera" # 
FT$SubOrder[FT$Family == "Trameini"] = "Anisoptera" # 
FT$SubOrder[FT$Family == "Trameini"] = "Anisoptera" # 
FT$SubOrder[FT$Family == "Stenophlebiinae"] = "Anisoptera" # 
FT$SubOrder[FT$Family == "Steleopteridae"] = "Zygoptera" # 
FT$SubOrder[FT$Family == "Trameinae"] = "Anisoptera" # 
FT$SubOrder[FT$Family == "Sphenophlebiidae"] = "Anisoptera" # 
FT$SubOrder[FT$Family == "Sieblosiidae"] = "Zygoptera" # 
FT$SubOrder[FT$Family == "Rhyothemistinae"] = "Anisoptera" # 
FT$SubOrder[FT$Family == "Selenothemidae"] = "Anisoptera" # 
FT$SubOrder[FT$Family == "Pseudocymatophlebiinae"] = "Anisoptera" # 
FT$SubOrder[FT$Family == "Protoneuridae"] = "Zygoptera" # 
FT$SubOrder[FT$Family == "Protolindeniina"] = "Anisoptera" # 
FT$SubOrder[FT$Family == "Proterogomphidae"] = "Anisoptera" # 
FT$SubOrder[FT$Family == "Prostenophlebiidae"] = "Anisoptera" # 
FT$SubOrder[FT$Family == "Prohemeroscopidae"] = "Anisoptera" # 
FT$SubOrder[FT$Family == "Priscalestidae"] = "Zygoptera" # 
FT$SubOrder[FT$Family == "Parastenophlebiidae"] = "Anisoptera" # 
FT$SubOrder[FT$Family == "Palaeomacromiidae"] = "Anisoptera" # 
FT$SubOrder[FT$Genus == "Oryctothemis"] = "Anisoptera" # 
FT$SubOrder[FT$Genus == "Syrrhoe"] = "Anisoptera" # 
FT$SubOrder[FT$Genus == "Magnasupplephlebia"] = "Anisoptera" # 
FT$SubOrder[FT$Family == "Nannogomphidae"] = "Anisoptera" # 
FT$SubOrder[FT$Family == "Myopophlebiidae"] = "Anisoptera" # 
FT$SubOrder[FT$Family == "Mesuropetalidae"] = "Anisoptera" # 
FT$SubOrder[FT$Family == "Mesochlorogomphidae"] = "Anisoptera" # 
FT$SubOrder[FT$Family == "Megapodagrionidae"] = "Zygoptera" # 
FT$SubOrder[FT$Family == "Liupanshaniidae"] = "Anisoptera" # 
FT$SubOrder[FT$Family == "Libellulidae"] = "Anisoptera" # 
FT$SubOrder[FT$Family == "Liassostenophlebiidae"] = "Anisoptera" # 
FT$SubOrder[FT$Family == "Liassophlebiidae"] = "Anisoptera" # 
FT$SubOrder[FT$Family == "Liassogomphidae"] = "Anisoptera" # 
FT$SubOrder[FT$Family == "Lestoidea"] = "Zygoptera" # 
FT$SubOrder[FT$Family == "Juracorduliidae"] = "Anisoptera" # 
FT$SubOrder[FT$Family == "Isophlebiidae"] = "Anisoptera" # 
FT$SubOrder[FT$Family == "Idionychidae"] = "Anisoptera" # 
FT$SubOrder[FT$Family == "Hypolestinae"] = "Zygoptera" # 
FT$SubOrder[FT$Family == "Henrotayiidae"] = "Anisoptera" # 
FT$SubOrder[FT$Family == "Hemiphlebiidae"] = "Zygoptera" # 

FT$SubOrder[FT$SubOrder == "Anisozygoptera"] = "Anisoptera"
}

D = merge(FT,D,id="Genus")
D$GenusSpecies = D$species
D = D[c("SubOrder","Family","Genus","GenusSpecies","size_mm","body_part","age_ma","paleolat","specimen_no","collection_name","collection_no")]

FossilsLatTax = D
# save(FossilsLatTax,file="C:/Users/Owner/Desktop/FossilsLatTax.rda")
}

if(FALSE) { # Fossils vs PaleoLat Plot
load(file="C:/Users/Owner/Desktop/envBodySize/data/Fossils/FossilsLatTax.rda")
D = FossilsLatTax

D = D[!duplicated(D$GenusSpecies),]


library(nlme)

# D = D[!is.na(D$Family),]

D$age_cat = D$age_ma
D$age_cat = "80-0"
D$age_cat[D$age_ma > 80] = "150-80"
D$age_cat[D$age_ma > 150] = "210-150"

# levels(D$age_cat) = 
D$age_cat = factor(D$age_cat,c("210-150","150-80","80-0"))


library(ggplot2)
ggplot(D, aes(paleolat,size_mm)) + 
geom_point(aes(colour=SubOrder), size=2) + 
geom_point(shape = 1,size = 2,colour = "black") + 
geom_smooth(method="gam",se = FALSE,colour="gray31") + 
ggtitle("Fossils: wing length and latitude") + 
facet_wrap(~age_cat) + 
xlab("Paleo-latitude") + 
ylab("Wing Length (mm)") + 
theme_bw()

# scale_colour_gradient(low="antiquewhite",high="tan4") + 
ggsave("C:/Users/Owner/Desktop/PaleoLat.pdf",width=10,height=5)
}

if(FALSE) { # Plot present day with paleolat

# load present day body sizes
load(file = "C:/Users/Owner/Desktop/envBodySize/data/opdb/quickBodySize.rda") # load cleaned body size data
B = quickBodySize # B for body size
B = B[c("GenusSpecies","Species","Genus","Family","SubOrder","tbl","hwl")] # drop fwl since often NA
B = B[c("GenusSpecies","Species","Genus","Family","SubOrder","hwl")] # drop another variable so that is easy to work with 
B = na.omit(B) # remove missing 

load(file = "C:/Users/Owner/Desktop/envBodySize/data/gbif/climateLat.rda") # load un-cleaned climate data
C = climateLat # C for climate
C$GenusSpecies = gsub(" ","_",C$GenusSpecies) 

AL = aggregate(decimalLatitude ~ GenusSpecies, data = C, mean) # average latitude

D = merge(B,AL,id="GenusSpecies") # final data.frame

# make the data.frames match up for plotting
D$age_ma = 0
D$lat = D$decimalLatitude
D = D[c("GenusSpecies","Genus","Family","SubOrder","hwl","age_ma","lat")]
colnames(D) = c("GenusSpecies","Genus","Family","SubOrder","size_mm","age_ma","lat")

# D = D[dplyr::between(D$lat, 25, 50),] # get only in the range present in the fossil record 

load(file="C:/Users/Owner/Desktop/envBodySize/data/Fossils/FossilsLatTax.rda") # load fossil data
DF = FossilsLatTax # Data Fossils

DF$lat = DF$paleolat
DF = DF[c("GenusSpecies","Genus","Family","SubOrder","size_mm","age_ma","lat")]

D$type = "present"
DF$type = "fossil"

D = rbind(D,DF) # combine present and fossil datasets

D$age_cat = D$age_ma
D$age_cat = "80-0"
D$age_cat[D$age_ma > 80] = "150-80"
D$age_cat[D$age_ma > 150] = "210-150"
D$age_cat[D$type == "present"] = "0 extant"

D$age_cat = factor(D$age_cat,c("210-150","150-80","80-0","0 extant"))


library(ggplot2)
ggplot(D, aes(abs(lat),size_mm)) + 
geom_point(aes(colour=age_ma), size=2) + 
scale_colour_gradient(low="antiquewhite",high="tan4") + 
geom_point(shape = 1,size = 2,colour = "black") + 
geom_smooth(method="gam",se = FALSE,colour="gray31") + 
ggtitle("Fossils: wing length and latitude") + 
facet_grid(SubOrder~age_cat,scales="free_x") + 
xlab("Paleo-latitude") + 
theme_bw() + 
theme(legend.position=c(0.06, 0.25)) + 
ylab("Wing Length (mm)") 

ggsave("C:/Users/Owner/Desktop/PaleoLatPresent.pdf",width=10,height=5)

# fit = lme(size_mm ~ paleolat,random=~1|SubOrder/Family,data=D)
# summary(fit)

}

if(FALSE) { # altitude we don't have data for this yet
load(file = "C:/Users/Owner/Desktop/envBodySize/data/gbif/climateLat.rda") # load un-cleaned climate data
C = climateLat # C for climate
C$GenusSpecies = gsub(" ","_",C$GenusSpecies) 
colnames(C)
} 


if(FALSE) { # plot our two oxygen datasets...

O1 = read.table("C:/Users/Owner/Desktop/envBodySize/data/oxygen/oxygen.txt",header=TRUE)
O2 = read.table("C:/Users/Owner/Desktop/envBodySize/data/oxygen/oxygenCharcoal.txt",header=TRUE,sep=";")

jpeg("C:/Users/Owner/Desktop/plot.jpg")
plot(-1*O1$age,O1$oxygen,type="l",col="red")
points(-1*O1$age,O1$oxygen,col="red")
points(-1*O2$age,O2$oxygen,type="l")
points(-1*O2$age,O2$oxygen)

dev.off()
}

if(FALSE) { # Simulated data OU,EB,ENV,kappa,trend...

library(RPANDA)
# data(Cetacea)
data(InfTemp)
str(InfTemp)

# set.seed(1)
tree = rbdtree(birth=0.07, death=0.01,Tmax = max(InfTemp$Age))

tree

# first parameter is sigma2
trait = sim_t_env(tree, param=c(3,0.2), env_data=InfTemp, model="EnvExp",root.value=0, step=0.001, plot=TRUE)

ENV = fit_t_env(tree, trait, env_data=InfTemp,model = "EnvExp", scale=TRUE)

ENV

library(geiger)
dat = data.frame(trait)
phy = tree
BM = fitContinuous(phy, dat, model=c("BM"))
OU = fitContinuous(phy, dat, model=c("OU"))
EB = fitContinuous(phy, dat, model=c("EB"))
trend = fitContinuous(phy, dat, model=c("trend"))
white = fitContinuous(phy, dat, model=c("white"))
Kappa = fitContinuous(phy, dat, model=c("kappa")) # number of speciation events between, Pagel 199

AIC = setNames(c(BM$opt$aic,OU$opt$aic,EB$opt$aic,trend$opt$aic,Kappa$opt$aic, white$opt$aic,ENV$aic),c("BM","OU","EB","trend","Kappa","white","ENV"))



plotENV = function(x, steps = 100,...) {

if (is.function(x$model)) {
fun_temp <- function(x, temp, model, param) {
	rate_fun <- function(x) {
		model(x, temp, param)
	}
	rate <- rate_fun(x)
	return(rate)
}
}
else if (x$model == "EnvExp") {
fun_temp <- function(x, temp, model, param) {
	sig <- param[1]
	beta <- param[2]
	rate <- (sig * exp(beta * temp(x)))
	return(rate)
}
}
else if (x$model == "EnvLin") {
fun_temp <- function(x, temp, model, param) {
	sig <- param[1]
	beta <- param[2]
	rate <- sig + (beta - sig) * temp(x)
	return(rate)
}
}


t <- seq(0, x$tot_time, length.out = steps)
t = InfTemp$Age
rates <- fun_temp(x = t, temp = x$env_func, model = x$model, 
param = x$param)
plot(-t, rates, type = "l", xlab = "Times", ylab = bquote(paste("Evolutionary rates ", sigma)), ...)
results <- list(time_steps = t, rates = rates)
return(results)
}

ENV
x = plotENV(ENV,steps = 17632)

# str(x)

D = data.frame(Temp = InfTemp$Temperature,Rate = x$rates,Age = InfTemp$Age,Age_x = x$time_steps)

library(ggtree)
p = ggtree(tree)
d = data.frame(id = names(trait), value = trait)

facet_plot(p, panel='bar', data=d, geom=geom_segment, aes(x=0, xend=value, y=y, yend=y), size=1, color='steelblue') + 
theme_tree2()
ggsave("C:/Users/Owner/Desktop/plot2.pdf")


library(phytools)
pdf("C:/Users/Owner/Desktop/plot.pdf",width=10,height=15)

par(mfrow=c(4,1))
phenogram(tree,trait,lwd=1)
plot(-1*D$Age,D$Temp)
plot(-1*D$Age,D$Rate)
plot(D$Temp,D$Rate)

# plot(ENV)
# plotENV(ENV,steps = 17632)
# plot(InfTemp$Age,-1*InfTemp$Temperature)
# plot(tree)
# axisPhylo()
dev.off()

sort(AIC)
}

if(FALSE) { # start simple fit various models of phenotypic evolution to extant species...

# load present day body sizes
load(file = "C:/Users/Owner/Desktop/envBodySize/data/opdb/quickBodySize.rda") # load cleaned body size data
B = quickBodySize # B for body size
B = B[c("GenusSpecies","Species","Genus","Family","SubOrder","tbl","hwl")] # drop fwl since often NA
B = B[c("GenusSpecies","Species","Genus","Family","SubOrder","hwl")] # drop another variable so that is easy to work with 
B = na.omit(B) # remove missing

D = B # rename; this is our main data

library(ape); load("C:/Users/Owner/Desktop/envBodySize/data/trees/tree.rda") # load tree

# remove those not in tree or data
D = D[D$GenusSpecies %in% tree$tip.label,] # get only those with tip
# D = D[!duplicated(D$GenusSpecies),] # not necessary but lets leave it in 
tip = tree$tip.label[!tree$tip.label %in% D$GenusSpecies] # drop tips without data
tree = drop.tip(tree, tip) # drop missing


library(geiger)
dat = data.frame(setNames(D$hwl,D$GenusSpecies))
phy = tree
BM = fitContinuous(phy, dat, model=c("BM"))
OU = fitContinuous(phy, dat, model=c("OU"))
EB = fitContinuous(phy, dat, model=c("EB"))
white = fitContinuous(phy, dat, model=c("white"))
# trend = fitContinuous(phy, dat, model=c("trend"))
# Kappa = fitContinuous(phy, dat, model=c("kappa")) # number of speciation events between, Pagel 199

# AIC = setNames(c(BM$opt$aic,OU$opt$aic,EB$opt$aic,trend$opt$aic,Kappa$opt$aic, white$opt$aic),c("BM","OU","EB","trend","Kappa","white"))


print("OU")
paste("alpha = ", OU$opt$alpha)
paste("sg2 = ", OU$opt$sigsq)
paste("lnLik = ", round(OU$opt$lnL,3))
paste("n. param. = ",OU$opt$k)
paste("AIC = ", round(OU$opt$aic,3))


print("BM")
paste("sg2 = ", BM$opt$sigsq)
paste("lnLik = ", round(BM$opt$lnL,3))
paste("n. param. = ",BM$opt$k)
paste("AIC = ", round(BM$opt$aic,3))

print("EB")
paste("a = ", EB$opt$a)
paste("sg2 = ", EB$opt$sigsq)
paste("lnLik = ", round(EB$opt$lnL,3))
paste("n. param. = ",EB$opt$k)
paste("AIC = ", round(EB$opt$aic,3))

}



if(FALSE) { # Look at oxygen with our extant data

# load present day body sizes
load(file = "C:/Users/Owner/Desktop/envBodySize/data/opdb/quickBodySize.rda") # load cleaned body size data
B = quickBodySize # B for body size
B = B[c("GenusSpecies","Species","Genus","Family","SubOrder","tbl","hwl")] # drop fwl since often NA
B = B[c("GenusSpecies","Species","Genus","Family","SubOrder","hwl")] # drop another variable so that is easy to work with 
B = na.omit(B) # remove missing
B$hwl = log(B$hwl) # log body size

D = B # rename; this is our main data

library(ape); load("C:/Users/Owner/Desktop/envBodySize/data/trees/tree.rda") # load tree

# remove those not in tree or data
D = D[D$GenusSpecies %in% tree$tip.label,] # get only those with tip
# D = D[!duplicated(D$GenusSpecies),] # not necessary but lets leave it in 
tip = tree$tip.label[!tree$tip.label %in% D$GenusSpecies] # drop tips without data
tree = drop.tip(tree, tip) # drop missing

O = read.table("C:/Users/Owner/Desktop/envBodySize/data/oxygen/oxygen.txt",header=TRUE) # O for oxygen
O = O[O$age < 237,] # 237 # age of tree

# env_data Environmental data, given as a time continuous function (see, e.g. splinefun) or
# a data frame with two columns. The first column is time, the second column is
# the environmental data (temperature for instance).

env_data = O[c("age","oxygen")]

library(RPANDA)
trait = setNames(D$hwl,D$GenusSpecies)

# simulate brownian motion with our dataset
# trait = sim_t_env(tree, param=c(3,0), env_data=env_data, model="EnvExp",root.value=0, step=0.001, plot=TRUE)

ENV = fit_t_env(tree, trait, env_data=env_data,model = "EnvExp", scale=TRUE)

cat("other models...\n")
library(geiger)
dat = data.frame(setNames(D$hwl,D$GenusSpecies))
phy = tree
BM = fitContinuous(phy, dat, model=c("BM"))
OU = fitContinuous(phy, dat, model=c("OU"))
EB = fitContinuous(phy, dat, model=c("EB"))
trend = fitContinuous(phy, dat, model=c("trend"))
white = fitContinuous(phy, dat, model=c("white"))
Kappa = fitContinuous(phy, dat, model=c("kappa")) # number of speciation events between, Pagel 199

AIC = setNames(c(BM$opt$aic,OU$opt$aic,EB$opt$aic,trend$opt$aic,Kappa$opt$aic, white$opt$aic,ENV$aic),c("BM","OU","EB","trend","Kappa","white","ENV"))

sort(AIC)

}

if(FALSE) { # Compare simulation with our tree and oxygen

# load present day body sizes
load(file = "C:/Users/Owner/Desktop/envBodySize/data/opdb/quickBodySize.rda") # load cleaned body size data
B = quickBodySize # B for body size
B = B[c("GenusSpecies","Species","Genus","Family","SubOrder","tbl","hwl")] # drop fwl since often NA
B = B[c("GenusSpecies","Species","Genus","Family","SubOrder","hwl")] # drop another variable so that is easy to work with 
B = na.omit(B) # remove missing
B$hwl = log(B$hwl) # log body size

D = B # rename; this is our main data

library(ape); load("C:/Users/Owner/Desktop/envBodySize/data/trees/tree.rda") # load tree

# remove those not in tree or data
D = D[D$GenusSpecies %in% tree$tip.label,] # get only those with tip
# D = D[!duplicated(D$GenusSpecies),] # not necessary but lets leave it in 
tip = tree$tip.label[!tree$tip.label %in% D$GenusSpecies] # drop tips without data
tree = drop.tip(tree, tip) # drop missing

O = read.table("C:/Users/Owner/Desktop/envBodySize/data/oxygen/oxygen.txt",header=TRUE) # O for oxygen
O = O[O$age < 237,] # 237 # age of tree

# env_data Environmental data, given as a time continuous function (see, e.g. splinefun) or
# a data frame with two columns. The first column is time, the second column is
# the environmental data (temperature for instance).

env_data = O[c("age","oxygen")]

library(RPANDA)
trait = setNames(D$hwl,D$GenusSpecies)

# simulate brownian motion with our dataset
trait = sim_t_env(tree, param=c(10,0), env_data=env_data, model="EnvExp",root.value=0, step=0.001, plot=TRUE)
print("taco")
ENV = fit_t_env(tree, trait, env_data=env_data,model = "EnvExp", scale=TRUE)

cat("other models...\n")
library(geiger)
# dat = data.frame(setNames(D$hwl,D$GenusSpecies))
dat = data.frame(trait) # for simulation 
phy = tree
BM = fitContinuous(phy, dat, model=c("BM"))
OU = fitContinuous(phy, dat, model=c("OU"))
EB = fitContinuous(phy, dat, model=c("EB"))
trend = fitContinuous(phy, dat, model=c("trend"))
white = fitContinuous(phy, dat, model=c("white"))
Kappa = fitContinuous(phy, dat, model=c("kappa")) # number of speciation events between, Pagel 199

AIC = setNames(c(BM$opt$aic,OU$opt$aic,EB$opt$aic,trend$opt$aic,Kappa$opt$aic, white$opt$aic,ENV$aic),c("BM","OU","EB","trend","Kappa","white","ENV"))

sort(AIC)
}


if(FALSE) { # simulate 
library(TreeSim)
library(ape)


simModels = function(N) { 

L = list()
for(i in 1:N) {

tree = NA
class(tree) = "" 
while(class(tree) != "phylo") {
trees = sim.bd.age(age=237, numbsim =1, lambda=0.02, mu=0.001, frac = 1, mrca = FALSE,complete = TRUE, K = 0)
tree = trees[[1]]
print(tree)
}

O = read.table("C:/Users/Owner/Desktop/envBodySize/data/oxygen/oxygen.txt",header=TRUE) # O for oxygen
O = O[O$age < 237,] # 237 # age of tree

env_data = O[c("age","oxygen")]

library(RPANDA)
# trait = setNames(D$hwl,D$GenusSpecies)

# simulate brownian motion with our dataset
trait = sim_t_env(tree, param=c(10,0), env_data=env_data, model="EnvExp",root.value=0, step=0.001, plot=TRUE)
print(i)
ENV = fit_t_env(tree, trait, env_data=env_data,model = "EnvExp", scale=TRUE)
cat("other models...\n")

library(geiger)
# dat = data.frame(setNames(D$hwl,D$GenusSpecies))
dat = data.frame(trait) # for simulation 
phy = tree
BM = fitContinuous(phy, dat, model=c("BM"))
OU = fitContinuous(phy, dat, model=c("OU"))
EB = fitContinuous(phy, dat, model=c("EB"))
trend = fitContinuous(phy, dat, model=c("trend"))
white = fitContinuous(phy, dat, model=c("white"))
Kappa = fitContinuous(phy, dat, model=c("kappa")) # number of speciation events between, Pagel 199

L[[i]] = list(dat = dat, phy = phy, BM = BM, OU = OU, EB = EB, trend = trend, white = white, Kappa = Kappa, ENV = ENV)

# AIC = setNames(c(BM$opt$aic,OU$opt$aic,EB$opt$aic,trend$opt$aic,Kappa$opt$aic, white$opt$aic,ENV$aic),c("BM","OU","EB","trend","Kappa","white","ENV"))

# L[[i]] = sort(AIC)
}

return(L)
}

library(foreach)
library(doParallel)
library(random)
library(muscle)

cl = makePSOCKcluster(6)             
clusterSetRNGStream(cl)
registerDoParallel(cl)

X = foreach(i=1:6, .packages=c('TreeSim','geiger'), .multicombine=TRUE) %dopar% {  
L = simModels(50)
return(L)
}
stopCluster(cl)

save(X, file = "C:/Users/Owner/Desktop/X.rda")
}


if(FALSE) { # analyze data from simulations
load(file = "C:/Users/Owner/Desktop/X.rda")

for(j in 1:3) {
for(i in 1:length(X)) {
x = X[[i]][[j]]
OU = x$OU; BM = x$BM; EB = x$EB; trend = x$trend; Kappa = x$Kappa; white = x$white; ENV = x$ENV
AIC = setNames(c(BM$opt$aic,OU$opt$aic,EB$opt$aic,trend$opt$aic,Kappa$opt$aic,white$opt$aic,ENV$aic),c("BM","OU","EB","trend","Kappa","white","ENV"))
print(sort(AIC))
}
}

}


# list.files("C:/Users/Owner/Desktop/envBodySize/data/Fossils/")
# load("C:/Users/Owner/Desktop/envBodySize/data/Fossils/FossilsLatTax.rda")
# load(file = "C:/Users/Owner/Desktop/envBodySize/data/trees/tree.rda")
# ls()
# tree$tip.label
# str(FossilsLatTax)


if(FALSE) { # plotting

plotENV = function(x, steps = 100,...) {

if (is.function(x$model)) {
fun_temp <- function(x, temp, model, param) {
	rate_fun <- function(x) {
		model(x, temp, param)
	}
	rate <- rate_fun(x)
	return(rate)
}
}
else if (x$model == "EnvExp") {
fun_temp <- function(x, temp, model, param) {
	sig <- param[1]
	beta <- param[2]
	rate <- (sig * exp(beta * temp(x)))
	return(rate)
}
}
else if (x$model == "EnvLin") {
fun_temp <- function(x, temp, model, param) {
	sig <- param[1]
	beta <- param[2]
	rate <- sig + (beta - sig) * temp(x)
	return(rate)
}
}


t <- seq(0, x$tot_time, length.out = steps)
t = env_data$age
rates <- fun_temp(x = t, temp = x$env_func, model = x$model, 
param = x$param)
plot(-t, rates, type = "l", xlab = "Times", ylab = bquote(paste("Evolutionary rates ", sigma)), ...)
results <- list(time_steps = t, rates = rates)
return(results)
}

x = plotENV(ENV,steps = nrow(env_data))

# str(x)

D = data.frame(Oxygen = env_data$oxygen,Rate = x$rates,Age = env_data$age,Age_x = x$time_steps)

library(ggtree)
p = ggtree(tree)
d = data.frame(id = names(trait), value = trait)

facet_plot(p, panel='bar', data=d, geom=geom_segment, aes(x=0, xend=value, y=y, yend=y), size=1, color='steelblue') + 
theme_tree2()
ggsave("C:/Users/Owner/Desktop/plot2.pdf")


library(phytools)
pdf("C:/Users/Owner/Desktop/plot.pdf",width=10,height=15)

par(mfrow=c(4,1))
phenogram(tree,trait,lwd=1,spread.labels=FALSE, fsize=0)
plot(-1*D$Age,D$Oxygen)
plot(-1*D$Age,D$Rate)
plot(D$Oxygen,D$Rate)
dev.off()

# Environmental determinants of body size? 
# Helen Morlon rate of evolution stuff??? 
}

if(FALSE) {

load(file = "C:/Users/Owner/Desktop/envBodySize/data/gbif/climateLat.rda") # load un-cleaned climate data
C = climateLat # C for climate
C$GenusSpecies = gsub(" ","_",C$GenusSpecies) 

AC = aggregate(cbind(bio_max1,bio_max12) ~ GenusSpecies, data = C, mean) # average climate
colnames(AC) = c("GenusSpecies","Temp","Prec")

load(file = "C:/Users/Owner/Desktop/envBodySize/data/opdb/quickBodySize.rda") # load cleaned body size data
BS = quickBodySize # BS for body size DATA
BS = BS[c("GenusSpecies","Species","Genus","Family","SubOrder","tbl","hwl")] # drop fwl since often NA
BS = BS[c("GenusSpecies","Species","Genus","Family","SubOrder","hwl")] # drop another variable so that is easy to work with 
BS = na.omit(BS) # remove missing 

load(file="C:/Users/Owner/Desktop/envBodySize/data/gbif/Birds.rda") # load Birds
B = Birds # rename B for Birds

# str(B) # what variable are available
B$diversity = B$all # select Bird variable to analyze 

AB = aggregate(diversity ~ GenusSpecies,data=B,mean) # get average birds diversity 
AB$GenusSpecies = gsub(" ","_",AB$GenusSpecies)
D = merge(BS,AB,id="GenusSpecies")
D = merge(D,AC) 

library(ape); load("C:/Users/Owner/Desktop/envBodySize/data/trees/tree.rda") # load tree and ape for pgls

}


if(FALSE) { # Simple Heat map example with points


load(file = "C:/Users/Owner/Desktop/envBodySize/data/gbif/climateLat.rda") # load un-cleaned climate data
C = climateLat # C for climate
C = C[c("GenusSpecies", "decimalLatitude", "decimalLongitude")]
colnames(C) = c("GenusSpecies","Lat","Lon")
C$GenusSpecies = gsub(" ","_",C$GenusSpecies)

load(file = "C:/Users/Owner/Desktop/envBodySize/data/opdb/quickBodySize.rda") # load cleaned body size data
B = quickBodySize # B for body size
B = B[c("GenusSpecies","Species","Genus","Family","SubOrder","hwl")] # drop fwl since often NA
B = na.omit(B) # remove missing 

D = merge(B,C,id="GenusSpecies")
# D = D[D$SubOrder == "Zygoptera",] # only damselflies
D = D[D$SubOrder == "Anisoptera",] # only damselflies


ss = 5 # step size
binLt = seq(-90,90,ss)
binLn = seq(-180,180,ss)


xLt = D$Lat

DL = list()
count = 1
for(i in 1:(length(binLt)-1)) {

d = D[xLt > binLt[i] & xLt < binLt[i+1],] # get latitude band


for(j in 1:(length(binLn)-1)) {
xLn = d$Lon

dd = d[xLn > binLn[j] & xLn < binLn[j+1],]
if(nrow(dd) != 0) {
dd$lat_point = binLt[i]
dd$lon_point = binLn[j]
}
DL[[count]] = dd

count = count + 1
}

}

dl = list() # output datalist
for(i in 1:length(DL)) {
d = DL[[i]] 

if(nrow(d) != 0) { 
y = d$lat_point[1]
x = d$lon_point[1]
n = length(d$GenusSpecies)
hwl = mean(d[!duplicated(d$GenusSpecies),]$hwl)

} else {
x = NA
y = NA
n = NA
hwl = NA
}

dl[[i]] = data.frame(n,hwl,x,y)
}

D = na.omit(do.call("rbind",dl))

colfunc = colorRampPalette(c("red","orange","yellow","green","lightsteelblue","blue","darkblue","violet"))
Colors = rev(colfunc(round(max(D$hwl))))


cols = Colors[round(D$hwl)]
cols
jpeg("C:/Users/Owner/Desktop/plot.jpg",units = "in",width=15,height=10,res=300)
plot(D$x,D$y,col = cols,pch=19,cex=3,ylim = c(-80,80))
abline(h = 0)
# image(t(t(X)))
dev.off()

}

if(FALSE) { # Heat Map with Hexagons
load(file = "C:/Users/Owner/Desktop/envBodySize/data/gbif/climateLat.rda") # load un-cleaned climate data
C = climateLat # C for climate
C = C[c("GenusSpecies", "decimalLatitude", "decimalLongitude")]
colnames(C) = c("GenusSpecies","Lat","Lon")
C$GenusSpecies = gsub(" ","_",C$GenusSpecies)

load(file = "C:/Users/Owner/Desktop/envBodySize/data/opdb/quickBodySize.rda") # load cleaned body size data
B = quickBodySize # B for body size
B = B[c("GenusSpecies","Species","Genus","Family","SubOrder","hwl")] # drop fwl since often NA
B = na.omit(B) # remove missing 

D = merge(B,C,id="GenusSpecies")
# D = D[D$SubOrder == "Zygoptera",] # only damselflies
D = D[D$SubOrder == "Anisoptera",] # only damselflies


ss = 5 # step size
binLt = seq(-90,90,ss)
binLn = seq(-180,180,ss)

# COL <- adjustcolor(c("red", "blue", "darkgreen")[iris$Species], alpha.f = 0.5)
xLt = D$Lat

DL = list()
count = 1
for(i in 1:(length(binLt)-1)) {

d = D[xLt > binLt[i] & xLt < binLt[i+1],] # get latitude band


for(j in 1:(length(binLn)-1)) {
xLn = d$Lon

dd = d[xLn > binLn[j] & xLn < binLn[j+1],]
if(nrow(dd) != 0) {
dd$lat_point = binLt[i]
dd$lon_point = binLn[j]
}
DL[[count]] = dd

count = count + 1
}

}

dl = list() # output datalist
for(i in 1:length(DL)) {
d = DL[[i]] 

if(nrow(d) != 0) { 
y = d$lat_point[1]
x = d$lon_point[1]
n = length(d$GenusSpecies)
hwl = mean(d[!duplicated(d$GenusSpecies),]$hwl)

} else {
x = NA
y = NA
n = NA
hwl = NA
}

dl[[i]] = data.frame(n,hwl,x,y)
}

D = do.call("rbind",dl)



# hexagon plot

library(RColorBrewer) #to use brewer.pal
library(fields) #to use designer.colors

Hexagon <- function (x, y, unitcell = 1, col = col) {
  polygon(c(x, x, x + unitcell/2, x + unitcell, x + unitcell, 
            x + unitcell/2), c(y + unitcell * 0.125, 
                               y + unitcell * 0.875, 
                               y + unitcell * 1.125, 
                               y + unitcell * 0.875, 
                               y + unitcell * 0.125, 
                               y - unitcell * 0.125), 
          col = col, border=NA)
}#function


# Heatmap_Matrix = matrix(rnorm(1000),nrow=100)

hwl = D$hwl
hwl[is.na(hwl)] = 0 # replace with 0
# Heatmap_Matrix = matrix(hwl,ncol=(length(binLt)-1))
Heatmap_Matrix = t(matrix(hwl,ncol=(length(binLt)-1)))

x <- as.vector(Heatmap_Matrix)
SOM_Rows <- dim(Heatmap_Matrix)[1]
SOM_Columns <- dim(Heatmap_Matrix)[2]

jpeg("C:/Users/Owner/Desktop/plot.jpg",units = "in",width=15,height=10,res=300)
# par(mfrow = c(1,2))
plot(0, 0, type = "n", axes = FALSE, xlim=c(0, SOM_Columns), ylim=c(0, SOM_Rows), xlab="", ylab= "", asp=1)


ColRamp <- rev(designer.colors(n=20, col=brewer.pal(9, "Spectral")))

ColorCode <- rep("#FFFFFF", length(x)) #default is all white
Bins <- seq(min(x, na.rm=T), max(x, na.rm=T), length=length(ColRamp))
for (i in 1:length(x))
  if (!is.na(x[i])) ColorCode[i] <- ColRamp[which.min(abs(Bins-x[i]))] 


offset <- 0.5 #offset for the hexagons when moving up a row
for (row in 1:SOM_Rows) {
  for (column in 0:(SOM_Columns - 1)) 
    Hexagon(column + offset, row - 1, col = ColorCode[row + SOM_Rows * column])
  offset <- ifelse(offset, 0, 0.5)
}

# image(Heatmap_Matrix)
dev.off()
}

if(FALSE) { # process Continents into temperate and tropical


L = as.character(D$Continents)

L = lapply(L, function(x) {
if(x != "")
substring(x, seq(1,nchar(x),2), seq(2,nchar(x),2))
else ""
})

L = lapply(L, function(x) gsub("AS","As",x))
L = lapply(L, function(x) gsub("NA","Na",x))
L = lapply(L, function(x) gsub("EU","Eu",x))
L = lapply(L, function(x) gsub("SA","Sa",x))
L = lapply(L,function(x) gsub("As|Eu|Na","temp",x)) # temperate
L = lapply(L,function(x) gsub("Af|Au|Sa","trop",x)) # tropical
L = lapply(L, unique)
L = lapply(L, function(x) if(length(x) > 1) "both" else x)
lat = do.call("c",L)

tbl = D$Male.average.body.length..mm.
tbl = as.numeric(gsub(",",".",tbl))
tbl[tbl == 0] = NA 

D$tbl = tbl # add tbl back into data
D$lat = lat # add back to data

D = D[!D$lat == "both",]
D = D[!D$lat == "",]
D = D[c("lat","tbl","GenusSpecies","Genus","Family","SubOrder")]
D$lat = as.factor(D$lat)
D = na.omit(D) # remove missing


# Naive
fit = lm(tbl ~ lat, data = D)
summary(fit)

# mixed model
library(nlme)
fit = lme(tbl ~ lat, random = ~ 1|SubOrder/Family/Genus, data = D)
summary(fit)

# PGLS 
library(ape)
load(file = "C:/Users/Owner/Desktop/BodySize/workingtree.rda")

D = D[D$GenusSpecies %in% tree$tip.label,] # drop data not in tree

tip = tree$tip.label[!tree$tip.label %in% D$GenusSpecies] # drop outgroup
tree = drop.tip(tree, tip)

library(nlme)
V = corBrownian(1,phy=tree)
fit = gls(tbl~lat,correlation=V,data=D)
summary(fit)

}

if(FALSE) { # output simple table for miguel
# list.files("C:/Users/Owner/Desktop/envBodySize/data/")
load("C:/Users/Owner/Desktop/John/envBodySize/data/gbif/trees.rda")
T = trees # T for Tree Cover
# T = T[!T$tree == 0,] # remove zero value
T$GenusSpecies = gsub(" ","_",T$GenusSpecies)
AT = aggregate(trees ~ GenusSpecies, data = T, mean) # AT average tree cover
colnames(AT) = c("GenusSpecies","cover")

load(file = "C:/Users/Owner/Desktop/John/envBodySize/data/gbif/climateLat.rda") # load un-cleaned climate data
C = climateLat # C for climate
C$GenusSpecies = gsub(" ","_",C$GenusSpecies) 

AC = aggregate(cbind(decimalLatitude,decimalLongitude,bio_max1,bio_max12) ~ GenusSpecies, data = C, mean) # average climate
colnames(AC) = c("GenusSpecies","lat","lon","Temp","Prec")

load(file = "C:/Users/Owner/Desktop/John/envBodySize/data/opdb/quickBodySize.rda") # load cleaned body size data
BS = quickBodySize # BS for body size DATA
BS = BS[c("GenusSpecies","Species","Genus","Family","SubOrder","tbl","hwl")] # drop fwl since often NA
BS = BS[c("GenusSpecies","Species","Genus","Family","SubOrder","hwl")] # drop another variable so that is easy to work with 
BS = na.omit(BS) # remove missing 

# load biotic grids 
load(file="C:/Users/Owner/Desktop/John/envBodySize/data/gbif/Mammals.rda") # 
A = Mammals # rename A for Amphibians

load(file="C:/Users/Owner/Desktop/John/envBodySize/data/gbif/Birds.rda") # load Birds
B = Birds # rename B for Birds

# str(B) # what variable are available
B$div_brd = B$all # select Amphibians variable to analyze 
A$div_mam = A$all_spp # select Bird variable to analyze 

AB = aggregate(div_brd ~ GenusSpecies,data=B,mean) # get average birds diversity 
AB$GenusSpecies = gsub(" ","_",AB$GenusSpecies)
D = merge(BS,AB,id="GenusSpecies")
AA = aggregate(div_mam ~ GenusSpecies,data=A,mean) # get average birds diversity 
AA$GenusSpecies = gsub(" ","_",AA$GenusSpecies)

D = merge(D,AA,id="GenusSpecies")
D = merge(D,AC,id="GenusSpecies") # merge diversity with average climates
D = merge(D,AT,id="GenusSpecies") # merge diversity with average climates

str(D)

write.table(D,"C:/USers/Owner/Desktop/John/envBodySize/spatialData.csv",sep=",")
}


# GBIF.org (09 November 2018) GBIF Occurrence Download https://doi.org/10.15468/dl.ynwpkx
# make odonate diversity grids 

library(dplyr)

D = data.table::fread("C:/Users/Owner/Desktop/OdonateOcc.csv")

D = D %>% select(decimalLatitude,decimalLongitude,speciesKey)


str(D)

pryr::mem_used()









