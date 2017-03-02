##Hacks to be used in conjunction with Lawson's example script. If any don't work just ask me, they are by no means tested on all 
#datasets
   
###### 1.) function for mapstate of length N based on our tdend object######
## Supply with n to return a mapstatlist like object (list of n vectors each containing names of individuals in those "populations")
## usage: sub_tree = mapstate_n(n)
 
mapstate_n = function(n){
  x = cutree(tdend, n)
  ##make an orderd list of all such names
  list = list(length = n)
  
  
  for(i in 1:n){
    list[[i]] = names(x[x == i])
  }
  list
}


####### 2.) Plotting loop: Colours PCA by each branch split sequentially for 2 - N branches (uses mapstate_n)########

pcares<-mypca(dataraw)##Start by generating the pcares dataframe
rcols<-rainbow(max(21)) ##Set your colouramp with the maximum N value 

pdf("your_pca_split.pdf",height=16,width=12) ##open a pdf
par(mfrow=c(4,3))
##Here choose a suitable range ie for(i in 2:N)
for(i in 2:21){
  mapstate_t = mapstate_n(i)
  pcapops<-getPopIndices(rownames(dataraw),mapstate_t)
  pcanames<-rownames(dataraw)
  plot(pcares$vectors[,1],pcares$vectors[,2],col=rcols[pcapops],xlab=paste("PC",1),ylab=paste("PC",2),main=paste(i," : " ,"PC",1,"vs",2),pch = rcols[pcapops])
  text(pcares$vectors[,1],pcares$vectors[,2],labels=as.numeric(pcapops),col=rcols[pcapops],cex=0.5,pos=1)
  
}
dev.off()


####### 3.) Colouring the lables of your Dendogram by order of your inferred clusters ###########
##Making coloured 
mapstate<-extractValue(treexml,"Pop") # map state as a finestructure clustering
mapstatelist<-popAsList(mapstate)
library(dendextend)
label = labels(tdend)
pcapops<-getPopIndices(label,mapstatelist) ### Note mapstatelist may be replaced by mapstate_n(n)
#labels(tdend) = pcapops ##optional changes your lables to cluster number
colours = colors()
set.seed(73)

N = 22    ##Choose your N to reflect the number of clusters (ie length(mapstatelist))
col_ramp = sample(colours, N, replace = FALSE) 
labels_colors(tdend) = col_ramp[as.numeric(pcapops)]

pdf(file="your_dendogram_name.pdf",height=6,width=14)
par(mar=c(6,0,2,0),mfrow=c(1,1))
plot.dendrogram(tdend,horiz=FALSE,nodePar=list(cex=0,lab.cex=0.6),edgePar=list(p.lwd=0,t.srt=90,t.off=-0.5),axes=F)
dev.off() 


######## 4.) Colouring the lables of your popdend dendogram by inferred clusters ##########
tree = unique(as.numeric(pcapops))
labels_colors(popdend) = col_ramp[tree]  ##colours the labels of the tree by population

pdf(file="Your_coloured_poptree.pdf",height=6,width=12)
par(mar=c(8,2,2,2),cex=0.8)
plot.dendrogram(popdend,horiz=FALSE,nodePar=list(cex=0,lab.cex=1.2,las=2),edgePar=list(p.lwd=0,t.srt=0,t.off=0.3),yaxt="n",height=0.5,dLeaf=0.2)
dev.off()
