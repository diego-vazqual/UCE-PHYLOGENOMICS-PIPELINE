###########################################################################################################
# Authors: Tauana Cunha and Bruno de Medeiros
# Date: 20-Apr-2017
# This script takes a 0/1 occupancy matrix and creates a nice image
# Columns are species names, plus column named 'filename' (gene alignment files)
# Rows are the genes/OGs
# Dimensions are ordered so higher occupancy is on top left corner
# You should edit color, dimensions of pdf image, x axis ticks as desired
# Currently has 3 colors: black for presence, and 2 for empty parts of two matrices of different sizes
###########################################################################################################

```{r}
setwd("/Users/Diego/OneDrive/Escritorio/TFM/occupancy_matrix")
matriz <- read.table("occ.matrix.tsv", header=T, sep="\t")
rownames(matriz) <- matriz$filename # Gives OG names as rownames
matriz <- matriz[,colnames(matriz)!="filename"] # Removes column of OG file names
```

```{r}
### Simple, unordered matrix
#pdf(file="/Users/Diego/OneDrive/Escritorio/TFM/occupancy_matrix/matrix_unordered.pdf", width=4.5, height=2.25) # Image size in inches
par(mar=c(2,1,1,1)) # Margin without species names
image(x=1:dim(matriz)[1],y=1:dim(matriz)[2],as.matrix(matriz),col=c("orange2", "black"),yaxt='n',xaxt='n')
axis(1,at=c(997,1634), tick=T, las=1, cex.axis=1) # X axis ticks
dev.off()
```

```{r}
###### Ordered by occupancy
# Order matrix to show higher occupancy on top left corner
occ <- matriz[,order(colSums(matriz), decreasing=F)]
occ <- occ[order(rowSums(occ), decreasing=T),]
# Changes 0 for 2 in first X (e.g. 88) genes, for different color of smaller matrix
occ[1:797,][occ[1:797,]=="0"] <- as.integer("2")
occ[798:1634,][occ[798:1634,]=="0"] <- as.integer("0")
```

```{r}
#pdf(file="/Users/Diego/OneDrive/Escritorio/TFM/occupancy_matrix/matrix_occ.pdf", width=4.3, height=2.25) # Image size in inches
par(mar=c(2,1,1,1)) # Margin without species names
image(x=1:dim(occ)[1],y=1:dim(occ)[2],as.matrix(occ),col=c("coral", "black", "bisque"),yaxt='n',xaxt='n',ylab="")
axis(1,at=c(797,1634), tick=T, las=1, cex.axis=0.7, mgp=c(2,0.4,0)) # X axis ticks
mtext(side=3,text="50% occupancy: 1634 loci", adj=1, cex=0.4)
mtext(side=3,text="75% occupancy: 797 loci", adj=0, line=0, cex=0.4)
dev.off() #mtext(side=3,text="149 genes", adj=0, line=1, cex=0.75)

# With species names
pdf(file="/Users/Diego/OneDrive/Escritorio/TFM/occupancy_matrix/matrix_occ_spnames8.pdf", width=6.55, height=5.25)
par(mar=c(2,2.3,1,1)) # Margin with species names
image(x=1:dim(occ)[1],y=1:dim(occ)[2],as.matrix(occ),col=c("coral", "black", "bisque"),yaxt='n',xaxt='n',ylab="", xlab="Genes", mgp=c(1,0.1,0), cex.lab=0.7)
axis(2,at=seq(1,dim(occ)[2],1), colnames(occ), tick=F, las=1, cex.axis=0.2, mgp=c(3,0.1,0), font=3) # Species names
#mgp=c(2,0.3,0) # For distance between axis titles and labels to axis
axis(1,at=c(798,1634), tick=T, las=1, cex.axis=0.7, mgp=c(2,0.3,0)) # X axis ticks
mtext(side=3,text="50% occupancy: 1648 loci", adj=1, cex=0.4)
mtext(side=3,text="75% occupancy: 797 loci", adj=0, line=0, cex=0.4)
dev.off()
```
