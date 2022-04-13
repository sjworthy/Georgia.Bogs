setwd("~/Documents/Bogs/Melanie.Bogs")

library(ape)
library(picante)
library(geiger)
library(vegan)
library(seqinr)
library(phangorn)
library(phytools)
library(abind)
library(vegan)

#### Tree with only species that have rbcL and matK ####

# Building the phylogenies

# export .fasta files from Geneious
# sequences must be in the same order (namewise) for both gene regions so that alignments are made in that order
# and so that cbind() correctly matches the sequences for the concatenations. It is easiest to just sort sequences 
# in alphabetical order in Geneious before exporting.

# run MAFFT alignment in Terminal
# place .fasta files into home/username folder
# Terminal code = mafft bog.rbcL.fasta > rbcl.aln.fasta
# Terminal code = mafft bog.matK.fasta > matk.aln.fasta

# Alignment Strategy Used = FFT-NS-2

#### Read in alignments and concatenate them ####

rbcl.aln=read.dna("rbcl.aln.fasta", format="fasta")
matk.aln=read.dna("matk.aln.fasta", format="fasta")

write.csv(rbcl.aln, file = "rbcL.aln.test.csv")
write.csv(matk.aln, file = "matk.aln.test.csv")

new.concat=cbind(rbcl.aln, matk.aln, fill.with.gaps=TRUE, check.names=FALSE)
write.dna(new.concat, "new.concat.fasta", format="fasta")
write.csv(new.concat, file="new.concat.test.csv")

both.seq.concat=read.phyDat("new.concat.fasta", format = "fasta", type="DNA")

# Modeltest - test different models to see which one is best for building the phylogeny.
# Model with the lowest AIC score is the best

modelTest(both.seq.concat) # GTR + G + I is the best

#### Maximum likelihood tree ####

dist.both.seq=dist.logDet(both.seq.concat)
both.seq.nj.tree=NJ(dist.both.seq)
both.seq.ml.model=pml(both.seq.nj.tree, both.seq.concat)
both.seq.ml.tree=optim.pml(both.seq.ml.model, model="GTR", optGamma = TRUE, optInv = TRUE)
both.seq.optim.ml.tree=optim.pml(both.seq.ml.model, model="GTR", optGamma = TRUE, optInv = TRUE, optNni = TRUE)
both.seq.optim.ml.tree.rooted=root(both.seq.optim.ml.tree$tree, outgroup="Pinus_taeda Pinus taeda  PNCA265 ribulose-1,5-bisphosphate carboxylase/oxygenase large subunit (rbcL) gene KC156886.1", resolve.root=T)
plot(both.seq.optim.ml.tree.rooted, cex=0.3)

#### Bootstrap trees ####

both.seq.boot=bootstrap.pml(both.seq.optim.ml.tree, bs=1000, optNni=TRUE)
both.seq.bootsrooted <- lapply(both.seq.boot, function(x) root(x, resolve.root=TRUE, outgroup="Pinus_taeda Pinus taeda  PNCA265 ribulose-1,5-bisphosphate carboxylase/oxygenase large subunit (rbcL) gene KC156886.1"))
class(both.seq.bootsrooted) <- "multiPhylo"
both.seqtree=plotBS(both.seq.optim.ml.tree$tree, both.seq.bootsrooted, p=75, type="phylo")
both.seqtree.rooted=root(both.seqtree, outgroup="Pinus_taeda Pinus taeda  PNCA265 ribulose-1,5-bisphosphate carboxylase/oxygenase large subunit (rbcL) gene KC156886.1", resolve.root=T)
both.seqtree.rooted$tip.label=all.tips
write.tree(both.seqtree.rooted, file="both.tree.bs..old.names.2.txt")

plot(both.seqtree.rooted, cex=0.3)

#### Phylogenetic Analyses ####

# add underscores (_) into names of tips of phylogeny to match cdm
all.phylo=both.seqtree.rooted
all.phylo$tip.label=all.tips.underscore

# drop the outgroup tip label for analyses

all.phylo=drop.tip(all.phylo, "Pinus_taeda")

#### Making the community data matrix (species as columns, bogs as rows) ####

bog.cdm=matrix(data=NA, nrow=3, ncol=161)
row.names(bog.cdm)=c("Colquitt", "Wheeler", "Charlton")
colnames(bog.cdm)=all.phylo$tip.label

write.csv(bog.cdm, file="bog.cdm.csv")

bog.cdm=read.csv("bog.cdm.csv", header=T, row.names=1)

#### Alpha Diversity Metrics ####
#### Faith's Index (PD) ####

pd.all.root=ses.pd(bog.cdm, all.phylo, null.model="taxa.labels", runs=999, iterations=1000)
write.csv(pd.all.root, file = "Faith.PD.csv")

#### Mean pairwise distance ####

mpd.all=ses.mpd(bog.cdm, cophenetic(all.phylo), null.model="taxa.labels", abundance.weighted=FALSE,
                runs=999, iterations=1000)
write.csv(mpd.all, file = "MPD.csv")

#### Mean nearest neighbor distance ####

mntd.all=ses.mntd(bog.cdm, cophenetic(all.phylo), null.model="taxa.labels", abundance.weighted=FALSE,
                  runs=999, iterations=1000)
write.csv(mntd.all, file = "MNTD.csv")

#### Beta Diversity Metrics ####
#### Phylogenetic Sorenson's Index ####
# tells you the phylogenentic similarity of communities

# getting observed phyloSor values
obs.phyloSor=phylosor(bog.cdm,all.phylo)

# getting random phyloSor values
phylo.Sor=phylosor.rnd(bog.cdm, all.phylo, null.model = "taxa.labels", runs = 999, iterations = 1000)

# https://pedrohbraga.github.io/CommunityPhylogenetics-Workshop/CommunityPhylogenetics-Workshop.html#phylogenetic-patterns-as-proxies-of-community-assembly-mechanisms
# create function to calculate the standardized effect size of a phylogenetic beta-diversity metric. In this case, we are using it for PhyloSor.

ses.PBD <- function(obs, rand){
  
  # first, we make sure our observed PhyloSor values are numeric
  pbd.obs <- as.numeric(obs) 
  
  # then, we take the mean of the 999 null expectations we generated
  rand <- t(as.data.frame(lapply(rand, as.vector)))
  pbd.mean <- apply(rand, MARGIN = 2, 
                    FUN = mean, 
                    na.rm = TRUE)
  # as well as their standard deviation
  pbd.sd <- apply(rand, 
                  MARGIN = 2, 
                  FUN = sd, 
                  na.rm = TRUE)
  
  # now, we can calculate the standardized effect size (SES)!
  pbd.ses <- (pbd.obs - pbd.mean)/pbd.sd
  
  # rank observed PhyloSor (we use this to calculate p-values for SES)
  pbd.obs.rank <- apply(X = rbind(pbd.obs, rand), 
                        MARGIN = 2, 
                        FUN = rank)[1, ]
  pbd.obs.rank <- ifelse(is.na(pbd.mean), NA, pbd.obs.rank)
  
  # return results in a neat dataframe
  data.frame(pbd.obs, 
             pbd.mean, 
             pbd.sd, 
             pbd.obs.rank, 
             pbd.ses, 
             pbd.obs.p = pbd.obs.rank/(dim(rand)[1] + 1))
}

phylo.Sor.all=ses.PBD(obs = obs.phyloSor, rand = phylo.Sor)
write.csv(phylo.Sor.all, file="Phylo.Sor.csv")

#### Pairwise Phylogenetic Dissimilairty (Dpw) ####

# randomizes community data matrix one time 
comdist.shuff <- function(x){
  as.matrix(comdist(bog.cdm,cophenetic(tipShuffle(x)),abundance.weighted = F))
}

# generate null distribution

nulls.Dpw=replicate(999, comdist.shuff(all.phylo))

# calculate mean and standard deviation

null.Dpw.mean=apply(nulls.Dpw,c(1:2),mean,na.rm=T)
null.Dpw.sd=apply(nulls.Dpw,c(1:2),sd,na.rm=T)

# calculate observed values

obs.Dpw=as.matrix(comdist(bog.cdm,cophenetic(all.phylo),abundance.weighted = F))
write.csv(obs.Dpw, file="obs.Dpw.csv")

# Calculate standardized effect sizes (SES)

ses.Dpw=(obs.Dpw-null.Dpw.mean)/null.Dpw.sd
write.csv(ses.Dpw, file = "SES.Dpw.csv")

# combine observed data with nulls

obs.nulls.Dpw=abind(obs.Dpw,nulls.Dpw)

# Calculate the rank of the values
# Rank values would need to be greater than 975 or less than 25 to be significant

temp.rank.Dpw=array(dim=dim(obs.nulls.Dpw),t(apply(apply(obs.nulls.Dpw,c(1,2),rank),3,t)))
rank.Dpw=temp.rank.Dpw[,,1]
colnames(rank.Dpw)=colnames(obs.Dpw)
row.names(rank.Dpw)=colnames(obs.Dpw)
write.csv(rank.Dpw, file="rank.Dpw.csv")

#### Nearest Neighbor Phylogenetic Dissimilarity ####

# randomizes community data matrix one time 
comdist.nt.shuff <- function(x){
  as.matrix(comdistnt(bog.cdm,cophenetic(tipShuffle(x)),abundance.weighted = F))
}

# generate null distribution

nulls.Dnn=replicate(999, comdist.nt.shuff(all.phylo))

# calculate mean and standard deviation

null.Dnn.mean=apply(nulls.Dnn,c(1:2),mean,na.rm=T)
null.Dnn.sd=apply(nulls.Dnn,c(1:2),sd,na.rm=T)

# calculate observed values

obs.Dnn=as.matrix(comdistnt(bog.cdm,cophenetic(all.phylo),abundance.weighted = F))
write.csv(obs.Dnn, file="obs.Dnn.csv")

# Calculate standardized effect sizes (SES)

ses.Dnn=(obs.Dnn-null.Dnn.mean)/null.Dnn.sd
write.csv(ses.Dnn, file = "SES.Dnn.csv")

# combine observed data with nulls

obs.nulls.Dnn=abind(obs.Dnn,nulls.Dnn)

# Calculate the rank of the values
# Rank values would need to be greater than 975 or less than 25 to be significant

temp.rank.Dnn=array(dim=dim(obs.nulls.Dnn),t(apply(apply(obs.nulls.Dnn,c(1,2),rank),3,t)))
rank.Dnn=temp.rank.Dnn[,,1]
colnames(rank.Dnn)=colnames(obs.Dnn)
row.names(rank.Dnn)=colnames(obs.Dnn)
write.csv(rank.Dnn, file="rank.Dnn.csv")

#### Taxonomic Indicies ####

# Sorenson's Dissimilarity Index
# Sorenson of species used in phylogenetic analyses
sorenson.output=vegdist(bog.cdm, binary = TRUE)

# Sorenson of all species in all bogs
taxa.cdm=read.csv("taxonomic.cdm.csv", header=T, row.names = 1)
sorenson.all.taxa=vegdist(taxa.cdm, binary=TRUE)

# Sorenson's similarity Index
designdist(bog.cdm, method ="(2*J)/(A+B)",terms = c("binary"))

### Plotting ####

setwd("~/Documents/Bogs/Melanie.Bogs/Results")

# Figure 3a
beta.results=read.csv("Sorenson.Value.csv", header=T, row.names=1)
par(mar=c(5.1, 4.1, 4.1, 5))
plot(beta.results[,1]~row.names(beta.results), pch=19, type="b",  cex.axis=1.5, cex=2, col="black", cex.lab=1.5,
     xaxt="n", ylab="Dissimilarity", xlab="Bogs", ylim=c(0.5, 0.8))
par(new=TRUE)
plot(beta.results[,2]~row.names(beta.results), pch=19, type="b", cex.axis=1.5, cex=2, col="gray", cex.lab=1.5,
     xaxt="n", ylab="", xlab="",yaxt="n",ylim=c(0.5, 0.8))
xtick=seq(300,500,by=100)
axis(side=1, at=xtick, cex.axis=1.5, labels=c("Colquitt-Wheeler", "Colquitt-Charlton", "Wheeler-Charlton"))

# Figure 3b
par(mar=c(5.1, 4.1, 4.1, 5))
phy.Sor=read.csv("Phylo.Sor.csv", header=T, row.names = 1)
plot(phy.Sor[,1]~row.names(phy.Sor), pch=19, type="b", cex=2, col="black", xaxt="n", 
     ylab="Similarity", xlab="Bogs", cex.axis=1.5, cex.lab=1.5,ylim=c(0.5, 0.8))
par(new=TRUE)
xtick=seq(1,3, by=1)
axis(side=1, at=xtick, cex.axis=1.5,labels=c("Colquitt-Wheeler", "Colquitt-Charlton", "Wheeler-Charlton"))

# Figure 5a
richness=read.csv("SR.table.csv", header=T, row.names = 1)
plot(richness[,1]~row.names(richness), cex.lab=1.5, cex.axis=1.5,pch=19, type="b", ylim=c(50,100), cex=2, col="black", xaxt="n", ylab="Species Richness", xlab="Site")
par(new=TRUE)
plot(richness[,2]~row.names(richness), pch=19, type="b", ylim=c(50,100), yaxt="n", cex=2, col="gray", xaxt="n", ylab="", xlab="")
par(new=TRUE)
xtick=seq(1,3, by=1)
axis(side=1, at=xtick, cex.axis=1.5,labels=c("Charlton", "Colquitt", "Wheeler"))

# Figure 5b
plot(richness[,3]~row.names(richness), pch=19, type="b", cex=2, cex.lab=1.5, cex.axis=1.5, col="black", xaxt="n", ylab="PD", xlab="Site")
par(new=TRUE)
xtick=seq(1,3, by=1)
axis(side=1, at=xtick, cex.axis=1.5, labels=c("Charlton", "Colquitt", "Wheeler"))



