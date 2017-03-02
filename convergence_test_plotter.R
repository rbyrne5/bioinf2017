###########################
#Note this is not a self sufficient script, run each line of code on objects created by Lawson's script
#First generate two meancoincidence files using Lawson's plotting script.
#Also need a tdend object from it
##Plotting both matrixes on one diagonal ordered by a single tree
##requires that you have built tdend as in lawson script and loaded two mean coincidence matrixes
fullorder<-labels(tdend) # the order according to the tree, choose one.
mcmcmatrixraw<-read.csv(meancoincidencefile,header = TRUE, row.names = 1, stringsAsFactors = FALSE)# read in the pairwise coincidence file we created earlier
mcmcmatrixraw2<-read.csv(meancoincidencefile2,header = TRUE, row.names = 1, stringsAsFactors = FALSE)#create this file first using second MCMC file

colnames(mcmcmatrixraw) = rownames(mcmcmatrixraw)
colnames(mcmcmatrixraw2) = rownames(mcmcmatrixraw2)

mcmcmatrix = as.matrix(mcmcmatrixraw[fullorder,fullorder])
mcmcmatrix_tri = mcmcmatrix
mcmcmatrix2 = as.matrix(mcmcmatrixraw2[fullorder,fullorder])
matrix_tri2 = mcmcmatrix2
mcmcmatrix_tri[lower.tri(mcmcmatrix_tri)] <- matrix_tri2[lower.tri(matrix_tri2)] 

#########################
## PLOT: Pairwise coincidence, both pairwise coincidences ordered by the same map

source("~/R_scripts/FinestructureLibrary.R")
pdf(file="convergence_test.pdf",height=12,width=12)
plotFinestructure(mcmcmatrix_tri,dimnames(row.names(mcmcmatrix_tri)),dend=tdend,cex.axis=0.6,edgePar=list(p.lwd=0,t.srt=90,t.off=-0.1,t.cex=0.8),main = "Independent runs")

dev.off()
 
