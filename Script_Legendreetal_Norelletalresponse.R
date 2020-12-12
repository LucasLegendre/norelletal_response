# Script compiled under R 4.0.3 (2020-10-10)

# Analyses performed by Lucas Legendre for the manuscript:
# Legendre LJ, Rubilar-Rogers D, Vargas AO, Clarke JA. The first dinosaur egg remains a mystery.

# All original data and trees from:
# Legendre, L.J., Rubilar-Rogers, D., Musser, G.M., Davis, S.N., Otero, R.A., Vargas, A.O.,
#   Clarke, J.A., 2020. A giant soft-shelled egg from the Late Cretaceous of Antarctica.
#   Nature 583, 411–414. https://doi.org/10.1038/s41586-020-2377-7
# Norell, M.A., Wiemann, J., Fabbri, M., Yu, C., Marsicano, C.A., Moore-Nall, A.,
#   Varricchio, D.J., Pol, D., Zelenitsky, D.K., 2020. The first dinosaur egg was soft.
#   Nature 583, 406–410. https://doi.org/10.1038/s41586-020-2412-8

# WARNING: set your working directory to your preferred folder before executing this script

library(ape)
library(phytools)
library(evobiR)
library(ggtree)
library(castor)

# data and tree
norelltree<-read.nexus("soft.trees.nex")
norelldata<-read.table("Norelldata.txt", header=T, row.names="Species")
norelldata2<-read.table("Norelldata.txt", header=T)

# data as factor / remove NAs / add branch lengths
norelldata<-na.omit(norelldata); norelldata2<-na.omit(norelldata2)

norelltree<-drop.tip(norelltree,setdiff(norelltree$tip.label,norelldata2$Species))
norelltree$edge.length<-rep(1,182)
plotTree(norelltree)


## Maximum likelihood (ML) ancestral states (in ape)
fitER<-ace(norelldata$Eggshell,norelltree,model="ER",type="discrete")
round(fitER$lik.anc,3)

# To visualize node numbers
ggtree(norelltree) + geom_text2(aes(subset=!isTip, label=node), hjust=-.3) + geom_tiplab()

# Order data as tips in the tree
norelldata<-ReorderData(norelltree, norelldata, taxa.names="row names")

# plot ancestral states
cols<-setNames(c("royalblue","green2","red3"),sort(unique(norelldata$Eggshell)))
plotTree(norelltree,fsize=0.7,ftype="i",lwd=1)
nodelabels(node=1:norelltree$Nnode+Ntip(norelltree),pie=fitER$lik.anc,piecol=cols,cex=0.5)
tiplabels(pie=to.matrix(norelldata$Eggshell,sort(unique(norelldata$Eggshell))),piecol=cols,cex=0.3)
add.simmap.legend(colors=cols,prompt=FALSE,x=25,y=30,fsize=0.8)

## ML marginal ancestral states, with rerootingMethod (phytools)
norellRR<-norelldata[,1]; names(norellRR)<-rownames(norelldata)
marganc<-rerootingMethod(norelltree,norellRR,model="ER")
plot(marganc)

# Same with ultrametric tree
ultranorell<-force.ultrametric(norelltree)
ultranorell$edge.length[ultranorell$edge.length<0.001]<-1; ultranorell<-force.ultrametric(ultranorell)
# repeat the above line three times so that no branch length is equal or too close to 0

fitultra<-ace(norelldata$Eggshell,ultranorell,model="ER",type="discrete")
plotTree(ultranorell,fsize=0.7,ftype="i",lwd=1)
nodelabels(node=1:ultranorell$Nnode+Ntip(ultranorell),pie=fitultra$lik.anc,piecol=cols,cex=0.5)
tiplabels(pie=to.matrix(norelldata$Eggshell,sort(unique(norelldata$Eggshell))),piecol=cols,cex=0.3)
add.simmap.legend(colors=cols,x=25,y=30,fsize=0.8)
# recovers hard-shelled as the most likely ancestral state for archosaurs

marganc<-rerootingMethod(ultranorell,norellRR,model="ER")
plot(marganc)
# same; branch length information greatly influences ancestral state reconstruction


# ====================== Change values of Mussaurus and Lufengosaurus =================
norelldata[c(51,52),]<-2
fitNew<-ace(norelldata$Eggshell,norelltree,model="ER",type="discrete", marginal=T)

# original tree
plotTree(norelltree,fsize=0.7,ftype="i",lwd=1)
nodelabels(node=1:norelltree$Nnode+Ntip(norelltree),pie=fitNew$lik.anc,piecol=cols,cex=0.5)
tiplabels(pie=to.matrix(norelldata$Eggshell,sort(unique(norelldata$Eggshell))),piecol=cols,cex=0.3)
add.simmap.legend(colors=cols,prompt=FALSE,x=25,y=30,fsize=0.8)

norellRR<-norelldata[,1]; names(norellRR)<-rownames(norelldata)
marganc<-rerootingMethod(norelltree,norellRR,model="ER")
plot(marganc)

# ultrametric tree
plotTree(ultranorell,fsize=0.7,ftype="i",lwd=1)
nodelabels(node=1:ultranorell$Nnode+Ntip(ultranorell),pie=fitNew$lik.anc,piecol=cols,cex=0.5)
tiplabels(pie=to.matrix(norelldata$Eggshell,sort(unique(norelldata$Eggshell))),piecol=cols,cex=0.3)
add.simmap.legend(colors=cols,prompt=FALSE,x=25,y=30,fsize=0.8)

norellRR<-norelldata[,1]; names(norellRR)<-rownames(norelldata)
marganc<-rerootingMethod(ultranorell,norellRR,model="ER")
plot(marganc)

# ============ Ancestral state reconstruction (ASR) using SIMMAP in phytools ==============

norelldata[c(51,52),]<-0
norellRR<-norelldata[,1]; names(norellRR)<-rownames(norelldata)

# Using make.simmap in phytools
mtrees<-make.simmap(norelltree, norellRR, model="ER", nsim=1000)
par(mfrow=c(10,10))
null<-sapply(mtrees,plotSimmap,colors=cols,lwd=1,ftype="off")
pd<-summary(mtrees); pd
par(mfrow=c(1,1))
pdf("SIMMAP-NorelltreeB12.pdf")
plot(pd,fsize=0.6,ftype="i",colors=cols,ylim=c(-2,Ntip(norelltree)))
add.simmap.legend(colors=cols,prompt=FALSE,x=20,y=0,vertical=FALSE)
dev.off()

## With corrected values for dibamid and sauropodomorphs
norelldata[c(51,52),]<-2
norelltree<-drop.tip(norelltree, "Anelytropsis_papillosus")
ultranorell<-drop.tip(ultranorell, "Anelytropsis_papillosus")
norelldata[c(29),]<-NA; norelldata<-na.omit(norelldata)
norellRR<-norelldata[,1]; names(norellRR)<-rownames(norelldata)

# original tree
mtrees<-make.simmap(norelltree, norellRR, model="ER", nsim=1000)
par(mfrow=c(10,10))
null<-sapply(mtrees,plotSimmap,colors=cols,lwd=1,ftype="off")
pd<-summary(mtrees); pd
par(mfrow=c(1,1))
plot(pd,fsize=0.6,ftype="i",colors=cols,ylim=c(-2,Ntip(norelltree)))
add.simmap.legend(colors=cols,prompt=FALSE,x=20,y=0,vertical=FALSE)

# ultrametric tree
mtrees<-make.simmap(ultranorell, norellRR, model="ER", nsim=1000)
par(mfrow=c(10,10))
null<-sapply(mtrees,plotSimmap,colors=cols,lwd=1,ftype="off")
pd<-summary(mtrees); pd
par(mfrow=c(1,1))
pdf("SIMMAP-Norelltreeultrametric.pdf")
plot(pd,fsize=0.6,ftype="i",colors=cols,ylim=c(-2,Ntip(ultranorell)))
add.simmap.legend(colors=cols,prompt=FALSE,x=20,y=0,vertical=FALSE)
dev.off()

# ==================== Maximum parsimony ASR, using castor ========================

norelldata[c(51,52),]<-1
norelldataMP<-norelldata+1
MP<-asr_max_parsimony(norelltree, norelldataMP[,1], Nstates=3)
pdf("MP-Norelltree.pdf")
plotTree(norelltree,fsize=0.4,lwd=1,ftype="i",colors=cols)
nodelabels(node=1:norelltree$Nnode+Ntip(norelltree),pie=MP$ancestral_likelihoods,piecol=cols,cex=0.5)
tiplabels(pie=to.matrix(norelldata$Eggshell,sort(unique(norelldata$Eggshell))),piecol=cols,cex=0.3)
add.simmap.legend(colors=cols,prompt=FALSE,x=25,y=30,fsize=0.8)
dev.off()

# ============= ASR with prismatic layer absent/present, using our sample ===============

wholesample<-read.table("Datawhole.txt", header=T)
wholetree<-read.nexus("treewholenew.trees.nex")
wholetree<-drop.tip(wholetree, setdiff(wholetree$tip.label,wholesample$Taxon))
wholesample<-read.table("Datawhole.txt", header=T, row.names="Taxon")
newRR<-wholesample[,10]; names(newRR)<-rownames(wholesample)
cols<-setNames(c("royalblue","red3"),levels(newRR))

# SIMMAP
mtrees<-make.simmap(wholetree, newRR, model="ER", nsim=1000)
pd<-summary(mtrees); pd
pdf("SIMMAP-Legendreneweggs.pdf")
plot(pd,fsize=0.4,ftype="i",colors=cols,ylim=c(-2,Ntip(wholetree)))
add.simmap.legend(colors=cols,prompt=FALSE,x=20,y=0,vertical=FALSE)
dev.off()

# SIMMAP with all branch lengths of 1
wholetree1<-wholetree; wholetree1$edge.length<-rep(1,315)
mtrees<-make.simmap(wholetree1, newRR, model="ER", nsim=1000)
pd<-summary(mtrees); pd
pdf("SIMMAP-Legendreneweggs-uncalibrated.pdf")
plot(pd,fsize=0.4,ftype="i",colors=cols,ylim=c(-2,Ntip(wholetree1)))
add.simmap.legend(colors=cols,prompt=FALSE,x=20,y=0,vertical=FALSE)
dev.off()

# =========================================================================================================
# New figures with new calibrated tree for our sample, with more soft-shelled eggs from Norell added

newdata<-read.table("Datawhole2.txt", header=T, stringsAsFactors = T); newdata<-newdata[-179,]
treenew<-read.nexus("treewholenew2.DGE.trees.nex")
newdata<-ReorderData(treenew, newdata, taxa.names="Species")
treenew2<-drop.tip(treenew, setdiff(treenew$tip.label, newdata$Taxon))
newdata<-read.table("Datawhole2.txt", header=T, row.names="Taxon")
newRR<-newdata[,10]; names(newRR)<-rownames(newdata)
cols<-setNames(c("royalblue","red3"),c("Absent","Present"))

# SIMMAP
mtrees<-make.simmap(treenew2, newRR, model="ER", nsim=1000)
pd<-summary(mtrees); pd
pdf("SIMMAP-Legendreneweggs-Norell.pdf")
plot(pd,fsize=0.5,lwd=1,type="fan",ftype="i",colors=cols,ylim=c(-2,Ntip(treenew2)),part=0.95)
add.simmap.legend(colors=cols,prompt=FALSE,x=20,y=0,vertical=FALSE)
dev.off()

# Same with Lufengosaurus and Massospondylus as hard-shelled
newdata[c(116,118),10]<-"Present"
newRR<-newdata[,10]; names(newRR)<-rownames(newdata)
mtrees<-make.simmap(treenew2, newRR, model="ER", nsim=1000)
pd<-summary(mtrees); pd
pdf("SIMMAP-Legendreneweggs-Norell-2hard.pdf",paper="a4")
plot(pd,type="fan",fsize=0.5,lwd=1,ftype="i",colors=cols,ylim=c(-2,Ntip(treenew2)),part=0.95)
add.simmap.legend(colors=cols,prompt=FALSE,x=20,y=0,vertical=FALSE)
dev.off()

# posterior probabilities of ancestral states for each node
describe.simmap(mtrees)$ace
ggtree(treenew2) + geom_text2(aes(subset=!isTip, label=node), hjust=-.3) + geom_tiplab()

# SIMMAP with all branch lengths of 1
wholetree1<-treenew2; wholetree1$edge.length<-rep(1,321)
mtrees<-make.simmap(wholetree1, newRR, model="ER", nsim=1000)
pd<-summary(mtrees); pd
pdf("SIMMAP-Legendreneweggs-uncalibrated.pdf")
plot(pd,fsize=0.4,ftype="i",colors=cols,ylim=c(-2,Ntip(wholetree1)))
add.simmap.legend(colors=cols,prompt=FALSE,x=20,y=0,vertical=FALSE)
dev.off()

