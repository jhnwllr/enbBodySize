
library(RColorBrewer) #to use brewer.pal
library(fields) #to use designer.colors

length(18:52)

# ColRamp <- rev(designer.colors(n=20, col=brewer.pal(9, "Spectral")))
# ColRamp

load(file = "C:/Users/Owner/Desktop/x.rda")


colfunc = colorRampPalette(c("red","orange","yellow","green","lightsteelblue","blue","violet"))
crclrs = 55
ColRamp = rev(colfunc(crclrs))

ColorCode = rep("#FFFFFF", length(x)) #default is all white
Bins = seq(min(x, na.rm=T), max(x, na.rm=T), length=length(ColRamp))

pdf("C:/Users/Owner/Desktop/scale.pdf", height = 5, width = 5,useDingbats=FALSE)
plot(Bins,rep(0,length(Bins)),col=ColRamp,pch=19,xlab="Temperature",bty="n",xaxt="n",yaxt="n",ylab="",ylim=c(0,1))
axis(1, at = seq(-5, 35, by = 5), las=1)
dev.off()

# for (i in 1:length(x)) if (!is.na(x[i])) ColorCode[i] = ColRamp[which.min(abs(Bins-x[i]))] 

# ColorCode
#FF0000

# Make scale 
# Color Gradient 
# colfunc = colorRampPalette(c("red","orange","yellow","green","lightsteelblue","blue","violet"))
# Colors = rev(colfunc(55))
# Colors
# pdf("C:/Users/Owner/Desktop/scale.pdf", height = 5, width = 5,useDingbats=FALSE)
# plot(x,rep(0,length(-5:35)),col=Colors[0:40],pch=19,yaxt = "n",ylab="",xlab="mm",ylim=c(0,1),bty="n",xaxt="n",xlim=c(-5,35),cex=1.5)
# axis(1, at = seq(-5, 35, by = 5), las=1)

# axis(labels=c("30"),side = 1, at = c(30,40,60,80,100,120))
# dev.off()


