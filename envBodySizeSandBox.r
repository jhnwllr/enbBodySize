

# if(FALSE) {


load("C:/Users/Owner/Desktop/envBodySize/data/Fossils/FossilsLatTax.rda")
load(file = "C:/Users/Owner/Desktop/envBodySize/data/trees/tree.rda")
load(file = "C:/Users/Owner/Desktop/envBodySize/data/Species.rda")

F = FossilsLatTax # fossils
F = F[!duplicated(F$Genus),]
S = Species # species 
F = F[!F$Family %in% S$Family,]

FL = split(F,F$Family) # fossil list 

FL = FL[sapply(FL, function(x) nrow(x) > 1)]

# length(FL)
f = FL[[1]]
f$Family = as.factor(f$Family)
f$Genus = as.factor(f$Genus)
f$GenusSpecies = as.factor(f$GenusSpecies)

library(ape)
tree = as.phylo(~GenusSpecies, data=f)
rownames(f) = f$GenusSpecies 
f = f[tree$tip.label,] # reorder data.frame according to tip.label 
tree$edge.length = (max(f$age_ma) - f$age_ma)
# max(f$age_ma) - f[f$GenusSpecies == "Aeschnidium bubas",]$age_ma # check 

# set up simple donar tree for testing
set.seed(1)
tree_z = rcoal(4) 
# tree_z$root.edge = 237
tree_z$root.edge = 2


# tree2 = bind.tree(tree_z, tree,position=max(f$age_ma))

# tree_xx = rcoal(2)
# tree2 = bind.tree(tree_z, tree_xx,position=0.1)

# tree2 = tree_z

# tree2$edge
# tree2$edge.length
# xxx = node.depth.edgelength(tree2)[tree2$edge[,2] <= tree2$Nnode - 1]

# length(xxx)
# xxx = node.depth.edgelength(tree2)
# xxx[1]

# tree2$root.edge = NULL
# 237 age of tree
jpeg("C:/Users/Owner/Desktop/plot.jpg")
# par(mfrow = c(1,2))
plot(tree_z)
axisPhylo()
nodelabels()
abline(v=0.5) # check ages
# plot(tree2)
# nodelabels()
# axisPhylo()
dev.off()
# }


if(FALSE) {
# > library(phytools)
# > # first, let's say we have a phylogeny with 20 tips
# > tree<-pbtree(n=20)
# > # plot the tree with node labels
# > plotTree(tree,node.numbers=T,pts=F)

library(phytools)
# ls("package:phytools")
set.seed(3)
tree = rcoal(5)

tips = c("t2","t3","t8")


node = findMRCA(tree, tips=tips, type=c("node","height"))
age = max(node.depth.edgelength(tree))
fossil_age = 0.3 # 1 mya
node_age = node.depth.edgelength(tree)[node]

node.depth.edgelength(tree)
el = (age-fossil_age) - node_age # edge length
if(el <= 0) el = 0 # check if fossil is older than group
tree$root.edge = 2
tree2 = bind.tip(tree, tip.label ="fossil", edge.length= el,where=5,position=0.5)



jpeg("C:/Users/Owner/Desktop/plot.jpg")
par(mfrow=c(1,2))
plot(tree,root.edge=TRUE)
axisPhylo()
nodelabels()
# abline(v=1.221)
# abline(v=0.95)
plot(tree2,root.edge=TRUE)
nodelabels()
axisPhylo()
dev.off()

}

# getNode = function(genus) {
# nl = list() # node list
# for(i in 1:length(genus)) {
# node = fastMRCA(tree, genus[1], genus[i])
# node_height = node.height(tree)[node]
# nl[[i]] = data.frame(node,node_height)
# }
# nld = do.call("rbind",nl) # node list data
# node = nld[order(nld$node_height),]$node[1]
# return(node)
# }




# library(TreeSim)

# set.seed(1)
# trees = sim.bd.age(age=10, numbsim =1, lambda=0.2, mu=0, frac = 1, mrca = FALSE,complete = TRUE, K = 0)
# tree = trees[[1]]

# tree$node.label = "x"

# tree = unroot(tree)
# print(tree)

# tree = rcoal(10)
# unclass(tree)
# nodeHeights(tree)
# node.height(tree)[11]
# unclass(tree)

# genus = c("t4","t6","t2")

# node = getNode(genus)
# node

# tip=list(edge=matrix(c(2L,1L),1,2),tip.label="added", edge.length=0,Nnode=1L)
# class(tip)="phylo"

# tree2 = bind.tree(tree, tip, where = 12, position = 0, interactive = FALSE)

# tree2 = bind.tip(tree, tip.label ="added", edge.length=0,where=15)
# tree2

# bind.tree


# bind.tip


