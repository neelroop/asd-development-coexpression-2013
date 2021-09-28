## Parikshak et al. code for network analysis
rm(list = ls()) ## Clear workspace
library(WGCNA) ## Load the WGCNA package in R
options(stringsAsFactors=FALSE)

## 0) Load all data and make directories in current folder for outputting the analyses
load(file="./AllSupData-11-1-2013.Rdata")

dir.create("./data")
dir.create("./figures")
dir.create("./tables")

## 1) Check the normalized data, metadata, gene model information from GENCODEV10, and gene lists for enrichment analysis
options(stringsAsFactors=FALSE) ## Set this so R uses the proper format for certain variables
print(dim(datExpr)); print(dim(datMeta));
##[1] 22084 transripts, 331 measured characteristics
##[1] 331 samples, 10 measured characteristics

## Check gene annotation data - it corresponds to GENECODEV10 or ENSEMBL65
ENSGExpr <- as.character(geneDat[match(rownames(datExpr),geneDat[,1]),1])
symbolExpr <- toupper(as.character(geneDat[match(rownames(datExpr),geneDat[,1]),2]))
dim(geneDat)
head(geneDat)
table(substr(geneDat[,3],nchar(geneDat[,3])-13,nchar(geneDat[,3]))=="protein_coding") # Number of protein coding transcripts
length(unique(geneDat[substr(geneDat[,3],nchar(geneDat[,3])-13,nchar(geneDat[,3]))=="protein_coding",2])) # 15585 unique symbols

## From the tab-delimitted version of Table S1, extract set memberships for 17 gene sets - note that the RDNV sets all reflect the "combined set" from the analysis
colorMat <- t(listMembership[,c(39:54,56)]) #transfer list membership to a matrix that's marked for membership
colnames(colorMat) <- listMembership[,"Gene.Symbol"]

## 2) Check the distribution of the data
apply(datExpr,2,quantile) # view distribution. Note that these values are from GC and gene length normalizaiton via conditional quantile normalization (package cqn in R), which leads to some negative values, which were set to 0. The normalized RPKM values were then log2(RPKM+1) transformed and corrected for batch effects relating to site of processing using ComBat.

## 4) Pick data that has RIN >= 9, age < 45 months (3 years post-natal age)
print(colnames(datMeta))
bestrin <- datMeta[,"RIN"]>=9
print(sort(datMeta[,"MonthsAge"]))
table(datMeta[,"MonthsAge"]<=45)  ## 9mo. + 12*3 mo. = 3 years of age

keep.for.dev <- datMeta[,8]<=45

print(bicorAndPvalue(datMeta[keep.for.dev,"RIN"],datMeta[keep.for.dev,"MonthsAge"])) # < 45 months with RIN >8 has a notable RIN correlation - p = 1.3e-6, r = -0.33

print(bicorAndPvalue(datMeta[keep.for.dev&bestrin,"RIN"],datMeta[keep.for.dev&bestrin,"MonthsAge"])) # < 45 months with RIN > 9 has low RIN correlation - p = 0.24, r = -0.10 - we use these to ensure network structure is not heavily influenced by RIN, and that the co-expression we observe is not driven by RNA quality

## 5) Pick samples for trajectory time points
devExpr <- vector(mode = "list", length = 2)
setLabels <- c("Prenatal2Child","Child2Adult")

for ( i in 1:length(devExpr) ) {
  devExpr[[i]]$name <- setLabels[i]
}
devExpr[[1]]$datExpr <- datExpr[,keep.for.dev&bestrin]
devExpr[[2]]$datExpr <- datExpr[,!keep.for.dev&bestrin]

devExpr[[1]]$datMeta <- datMeta[keep.for.dev&bestrin,]
devExpr[[2]]$datMeta <- datMeta[!keep.for.dev&bestrin,]

## 6) Remove outliers in each set using a sample network and connectivity-based metrics - if a sample is 5 standard deviations away from the mean connectivity to other samples (in a sample-sample network), then we remove it
sdout <- 5
for (i in 1:2) {
  ## Remove outliers
  ##Calculate signed, weighted biweight midcorrelation
  normadj <- (0.5+0.5*bicor(devExpr[[i]]$datExpr))^2

  ## Calculate connectivity
  netsummary <- fundamentalNetworkConcepts(normadj)
  ku <- netsummary$Connectivity
  z.ku <- ku-(mean(ku))/sqrt(var(ku))
  ## Declare as outliers those samples which are more than sdout sd above the mean connectivity based on the chosen measure
  outliers <- (z.ku > mean(z.ku)+sdout*sd(z.ku))|(z.ku < mean(z.ku)-sdout*sd(z.ku))
  print(paste("There are ",sum(outliers)," outliers samples based on a bicor distance sample network connectivity standard deviation above ",sdout,sep=""))
  print(colnames(devExpr[[i]]$datExpr)[outliers])
  print(table(outliers))
  devExpr[[i]]$datExpr <- devExpr[[i]]$datExpr[,!outliers]
  devExpr[[i]]$datMeta <- devExpr[[i]]$datMeta[!outliers,]
}

## This rmoved one outlier in the adult data, now check to ensure there's not a significant age-RIN correlation in the data
for (i in 1:length(setLabels)) {
  print(setLabels[i])
  print(dim(devExpr[[i]]$datExpr))
  print(bicorAndPvalue(devExpr[[i]]$datMeta[,"RIN"],devExpr[[i]]$datMeta[,"MonthsAge"]))
  print(table(goodGenes(t(devExpr[[i]]$datExpr))))##There are some genes that are perhaps all 0 at some time windows
}

############### Now, with high quality expression data, the steps for network construction

## 7) Pick Soft Threshold - this step is slow if we use the entire set and the blockwise approach saves time, so use blcoks of size 5000 - see WGCNA documentation for more info
bsize=5000

if (!file.exists(paste("./data/signed.sft.thresh.RData",sep=""))) {
  for (n in 1:2) {
    devExpr[[n]]$softSignedThresh <- pickSoftThreshold(data=t(devExpr[[n]]$datExpr),networkType="signed",corFnc="bicor",verbose=3,powerVector=c(seq(8, 30, by = 2)),blockSize=bsize)

    ##[1] "0to45moRIN9":
    ##     Power SFT.R.sq slope truncated.R.sq mean.k. median.k. max.k.
    #1      8    0.650 -1.85          0.959   537.0    463.00   1690
    #2     10    0.723 -1.85          0.961   304.0    240.00   1220
    #3     12    0.759 -1.89          0.958   185.0    132.00    924
    #4     14    0.780 -1.93          0.956   119.0     76.00    719
    #5     16    0.757 -2.01          0.933    79.4     45.50    573
    #6     18    0.773 -2.02          0.939    55.1     28.10    466
    #7     20    0.791 -2.01          0.948    39.4     17.80    385
    #8     22    0.784 -2.04          0.942    28.9     11.70    323
    #9     24    0.799 -2.01          0.949    21.6      7.73    274
    #10    26    0.808 -2.00          0.954    16.5      5.22    235 <- Scale-free fit >0.8 at power of 26, use this.
    #11    28    0.818 -1.98          0.958    12.8      3.58    203
    #12    30    0.829 -1.95          0.964    10.0      2.48    176

    print(setLabels[n])
    ##  print(devExpr[[n]]$softSignedThresh)
    ##    Power SFT.R.sq  slope truncated.R.sq mean.k. median.k. max.k.
    #1      8    0.177  0.963          0.974  403.00   351.000 1160.0
    #2     10    0.623 -1.540          0.967  211.00   168.000  806.0
    #3     12    0.790 -1.800          0.980  120.00    86.100  588.0
    #4     14    0.836 -1.860          0.982   73.20    46.400  445.0 <- best fit for adult network
    #5     16    0.851 -1.870          0.981   47.10    26.200  347.0
    #6     18    0.855 -1.850          0.978   31.70    15.200  276.0
    #7     20    0.857 -1.830          0.976   22.10     9.130  223.0
    #8     22    0.863 -1.790          0.975   15.90     5.600  183.0
    #9     24    0.876 -1.730          0.979   11.70     3.510  152.0
    #10    26    0.881 -1.690          0.980    8.87     2.270  128.0
    #11    28    0.901 -1.630          0.989    6.83     1.480  108.0
    #12    30    0.894 -1.640          0.988    5.34     0.989   95.8
  }
    save(devExpr,file=paste("./data/signed.sft.thresh.RData",sep=""))
} else {
    load(paste("./data/signed.sft.thresh.RData",sep=""))
  }


## 8) TO matrix calculation - this step takes a lot of memory, and it is preferable to run the entire TO calculation in one block. Requires at least 16G RAM - see WGCNA documentation for more info
INDS <- c(1,2)
softPower=c(26,14)

for (n in INDS) {
  cat(paste("Plotting ",devExpr[[n]]$name,sep=""))

  tmpMulti <- vector(mode = "list", length = 1) ## Store data in temporary structure
    tmpMulti[[1]]$data <- t(devExpr[[n]]$datExpr)

    if (!file.exists(paste("./data/signed",devExpr[[n]]$name,"-block.1.RData",sep=""))) {
      
      ## The cutting parameters do not matter - we cut the tree using cutreeHybrid in the next step
      ## Key parameters are corType="bicor", networkType="signed"
      devExpr[[n]]$netData <- blockwiseModules(datExpr=tmpMulti[[1]]$data,maxBlockSize=25000,networkType="signed",
                                               power=softPower[n],mergeCutHeight=0.15,nThreads=10,
                                               saveTOMFileBase = paste("./data/signed",devExpr[[n]]$name,sep=""),saveTOMs=TRUE,
                                               corType="bicor",minModuleSize=20,pamStage=FALSE, ## The module size and other parameters do not matter right now - we re-cut the tree in the next step
                                               reassignThreshold = 1e-10,verbose=3,deepSplit=2)
    }
}

## 9) Plot TOMs and cut clusters - Load TOMs from these analyses, plot with different deep split values - see WGCNA documentation for more info

for (n in INDS) {
  TOMfile <- paste("./data/signed",devExpr[[n]]$name,"-block.1.RData",sep="")
  load(TOMfile)
  geneTree= hclust(1-TOM,method="average");
  
  pdf(file=paste("./figures/",devExpr[[n]]$name,"_TOMcuts.pdf",sep=""),height=8,width=11)
  
  
  devExpr[[n]]$datExpr <- devExpr[[n]]$datExpr[goodGenes(t(devExpr[[n]]$datExpr)),]
  mColorh <- mLabelh <- colorLabels <- NULL
  minModSize <- 200 # Modules are at least 200 genes large
  dthresh <- 0.15 # MEs are no more than 0.85 correlated, if they are then the modules are merged and the ME is re-calculated
  ds <- 2 # deep split parameter to determine how finely to cut the tree
  tree = cutreeHybrid(dendro = geneTree, pamStage=FALSE,
    minClusterSize = minModSize, cutHeight = 0.99999, ## Consider high maximum joining heights for dendrogram
    deepSplit = ds, distM = as.matrix(1-TOM))
  
  merged <- mergeCloseModules(exprData = t(devExpr[[n]]$datExpr),colors = tree$labels,
                              cutHeight = dthresh)
  mColorh <- labels2colors(merged$colors)
  mLabelh <- "DS=2,MMS=200,DCOR=0.15"
  
  save(geneTree,mColorh,mLabelh,file=paste("./data/",devExpr[[n]]$name,"_TOMcuts.Rdata",sep=""))
  
  plotDendroAndColors(geneTree, mColorh, groupLabels = mLabelh,addGuide=TRUE,dendroLabels=FALSE,main="Dendrogram With Defined Modules")
  dev.off()
}

## 10) Over-representation analysis function
## Odds-ratio estimator
OR <- function(q,k,m,t) {
  ## 2 x 2 table:
  ##         inTest   !inTest
  ## inRef     q        k
  ## !inRef    m        t
  
  fisher.out <- fisher.test(matrix(c(q, k-q, m-q, t-m-k+q), 2, 2),conf.int=TRUE)
  OR <- fisher.out$estimate
  pval <- fisher.out$p.value
  upCI <- fisher.out$conf.int[1]
  downCI <- fisher.out$conf.int[2]
  
  output <- c(OR,pval,upCI,downCI)
  names(output) <- c("OR","Fisher p","-95%CI","+95%CI")
  return(output)
}

## count overlaps and run the analysis
ORA <- function(testpath,refpath,testbackground,refbackground) {
  q <- length(intersect(testpath,refpath)) ## overlapped pathway size
  k <- length(intersect(refpath,testbackground))  ## input gene set
  m <- length(intersect(testpath,refbackground)) ## input module
  t <- length(intersect(testbackground,refbackground)) ## Total assessed background (intersect reference and test backgrounds)
  
  empvals <- OR(q,k,m,t)
  
  tmpnames <- names(empvals)
  empvals <- as.character(c(empvals,q,k,m,t,100*signif(q/k,3)))
  names(empvals) <- c(tmpnames,"Overlap","Reference List","Input List","Background","% List Overlap")
  return(empvals)
}


## Add annotation of protein-coding genes (also used as background for most sets)
pcList <- as.character(geneDat[substr(geneDat[,3],(nchar(geneDat[,3])+1)-nchar("protein_coding"),nchar(geneDat[,3]))=="protein_coding",2]) # extract GC10 protein-coding genes only
pcColors <- rep("white",length=ncol(colorMat))
names(pcColors) <- colnames(colorMat)
pcColors[!is.na(match(names(pcColors),pcList))] <- "blue"
colorMat <- rbind(colorMat,pcColors)
colorMat[colorMat=="N"] <- "white"
colorMat[colorMat=="Y"] <- "black"
rownames(colorMat)[nrow(colorMat)] <- "Protein Coding"

## 11) getting gene-wise statistics

for (net in INDS) {
  print(devExpr[[net]]$name)
  
  ## Open files and load appropriate data
  pdf(file=paste("./figures/",devExpr[[net]]$name,"_TOMenrich.pdf",sep=""),height=16,width=18)
  dim(devExpr[[net]]$datExpr)
  load(file=paste("./data/",devExpr[[net]]$name,"_TOMcuts.Rdata",sep=""))
  thisMeta <- devExpr[[net]]$datMeta
  keep <- goodGenes(t(devExpr[[net]]$datExpr))
  thisExpr <- devExpr[[net]]$datExpr[keep,]
  
  ## These are the metadata to save for further analysis
  devExpr[[net]]$traitmat <- cbind(as.factor(thisMeta[,"Brain.code"]),as.factor(thisMeta[,"Region"]),
                                   as.numeric(thisMeta[,"MonthsAge"]),as.numeric(thisMeta[,"RIN"]))
  
  rownames(devExpr[[net]]$traitmat) <- rownames(thisMeta)
  colnames(devExpr[[net]]$traitmat) <- c("ID","Region","AgeInMo.","RIN")
  geneSigs <- matrix(NA,nrow=3,ncol=nrow(thisExpr))
  
  ## Find adjusted multiple R^2 for each gene with each categorical variable, bicor with quantitative variables
  for (i in 1:ncol(geneSigs)) {
    exprvec <- as.numeric(thisExpr[i,])
    regionr <- sqrt(max(summary(lm(exprvec ~ as.factor(devExpr[[net]]$traitmat[,"Region"])))$adj.r.squared,0))
    ager <- bicor(devExpr[[net]]$traitmat[,"AgeInMo."],exprvec)
    rinr <- bicor(devExpr[[net]]$traitmat[,"RIN"],exprvec)
    geneSigs[,i] <- c(regionr,ager,rinr)
  }
  devExpr[[net]]$genesigs <- geneSigs ## Save information into this data structure
  
  ## Convert numbers to colors to plot a dendrogram with these correlations on the bottom
  geneSigs[1,] <- numbers2colors(as.numeric(geneSigs[1,]),signed=FALSE,centered=FALSE,blueWhiteRed(100)[51:100],lim=c(0,1))
  geneSigs[2,] <- numbers2colors(as.numeric(geneSigs[2,]),blueWhiteRed(100),signed=TRUE,centered=TRUE,lim=c(-1,1))
  geneSigs[3,] <- numbers2colors(as.numeric(geneSigs[3,]),blueWhiteRed(100),signed=TRUE,centered=TRUE,lim=c(-1,1))
  rownames(geneSigs) <- c("Region","AgeInMo.","RIN")
  devExpr[[net]]$genecols <- geneSigs
  
  ## Save some of the network data
  devExpr[[net]]$netData$netName <- c(paste("Signed bicor network at power of 10",ncol(devExpr[[net]]$datExpr))," samples")
  devExpr[[net]]$netData$TOMdendrogram <- geneTree
  devExpr[[net]]$netData$moduleColors <- mColorh
  devExpr[[net]]$netData$cutParameters <- mLabelh
  
  keptSymbols <- as.character(geneDat[match(rownames(devExpr[[net]]$datExpr),geneDat[,1]),2])
  
  ## Plot these on a dendrogram to visualize how network analysis grouped some of the other correlational relationships to genes in the data
  thesecolors <- colorMat[,match(keptSymbols,colnames(colorMat))]
  mColorh1 <- cbind(mColorh,t(thesecolors[,keep]),t(geneSigs))
  mLabelh1 <- c(mLabelh,rownames(colorMat),rownames(geneSigs))
  plotDendroAndColors(geneTree, mColorh1, groupLabels = mLabelh1,addGuide=TRUE,dendroLabels=FALSE,main=devExpr[[net]]$netData$netName)
  multiData <- devExpr[[net]]
  save(multiData,file=paste("./tables/",devExpr[[net]]$name,sep=""))
  dev.off()
}

## 12) Now for module-level statistics
matchExpr <- vector(mode = "list", length = 2)

ref <- 1 ## Anchor colors on the developmental colors if making the later aging network
load(file=paste("./tables/",devExpr[[ref]]$name,sep=""))
refExpr <- multiData
keepref <- goodGenes(t(refExpr$datExpr))

for (net in INDS) {
  load(file=paste("./tables/",setLabels[net],sep=""))
  keep <- goodGenes(t(multiData$datExpr))
  
  ## Match up the genes present in both datasets if recoloring the dendrogram
  goodgenes <- intersect(rownames(refExpr$datExpr[keepref,]),rownames(multiData$datExpr[keep,]))
  keepmap <- match(goodgenes,rownames(multiData$datExpr)[keep])
  keepmapref <- match(goodgenes,rownames(refExpr$datExpr)[keepref])
  
  multiData$netData$matchedColors <- as.vector(matchLabels(source=multiData$netData$moduleColors[keepmap],reference=refExpr$netData$moduleColors[keepmapref],pThreshold=0.05))
  names(multiData$netData$matchedColors) <- rownames(multiData$datExpr)[keep][keepmap]
  multiData$netData$finalColors <- rep("grey",nrow(multiData$datExpr[keep,]))
  names(multiData$netData$finalColors) <- rownames(multiData$datExpr)[keep]
  
  multiData$netData$finalColors[!is.na((match(names(multiData$netData$finalColors),names(multiData$netData$matchedColors))))] <- multiData$netData$matchedColors[na.omit(match(names(multiData$netData$finalColors),names(multiData$netData$matchedColors)))]
  
  
  ## Compute the eigengenes
  multiData$netData$MEs <- moduleEigengenes(expr=t(multiData$datExpr[keep,]),colors=multiData$netData$finalColors)
  MEmat <- t(multiData$netData$MEs$eigengenes)
  colnames(MEmat) <- colnames(multiData$datExpr)

  ## Compute the ranked connectivity of each gene by correlation to the ME - kMEs
  ## Go through all the modules and match the labels, output kME tables containing:
  ## Gene Symbol, ENSEMBL ID, Module Color, Module Label, Trait Measures, kMEs
  multiData$netData$kMEtable <- signedKME(datExpr=t(multiData$datExpr[keep,]),datME=multiData$netData$MEs$eigengenes,corFnc="bicor")
  allGenes <- substr(rownames(multiData$netData$kMEtable),17,100)
  
  rownames(multiData$genesigs) <- rownames(multiData$genecols)
  keptSymbols <- geneDat[match(rownames(multiData$datExpr),geneDat[,1]),2]
  thesecolors <- colorMat[,match(keptSymbols,colnames(colorMat))]
  geneTable <- cbind(substr(rownames(multiData$netData$kMEtable),1,15),keptSymbols[keep],allGenes,multiData$netData$finalColors,t(multiData$genesigs),t(thesecolors[,keep]),multiData$netData$kMEtable)
  
  ## Eigengene correlations with variables
  MESigs <- matrix(NA,nrow=6,ncol=nrow(MEmat))
  ## Find adjusted multiple R^2 for each gene withe ach categorical variable
  for (i in 1:ncol(MESigs)) {
    exprvec <- as.numeric(MEmat[i,])
    lmsum <- summary(lm(exprvec ~ as.factor(multiData$traitmat[,"Region"])))
    aovsum <- kruskal.test(exprvec ~ as.factor(multiData$traitmat[,"Region"]))
    regionr <- sqrt(max(lmsum$adj.r.squared,0))
    regionp <- aovsum$p.value
    agesum <- bicorAndPvalue(multiData$traitmat[,"AgeInMo."],exprvec)
    ager <- agesum$bicor
    agep <- agesum$p
    rinsum <- bicorAndPvalue(multiData$traitmat[,"RIN"],exprvec)
    rinr <- rinsum$bicor
    rinp <- rinsum$p
    MESigs[,i] <- c(regionr,regionp,ager,agep,rinr,rinp)
  }
  rownames(MESigs) <- c("Region Adjusted r","Region KW p","Age bicor","Age p","RIN bicor","RIN p")
  
  ## Include some specific background sets
  refbackground <- pcList ## Reference background is protein-coding genes only
  refbackground <- matrix(NA,ncol=nrow(colorMat),nrow=length(pcList))
  colnames(refbackground) <- paste(rownames(colorMat),"background",sep="_")
  for (n in c(1:6,9:16)) { ## All but Voineagu et al and FMRP sets get the protein coding gene background
    refbackground[,n] <- pcList
  }
  
  for (n in c(7,8)) {
    refbackground[1:nrow(asdexpr),n] <- asdexpr[,1]
  }
  
  for (n in c(17)) {
    refbackground[1:nrow(mouseorthologs),n] <- mouseorthologs[,6]
  }
  refbackMat <- refbackground
  
  ## Compute enrichments in each module
  enrichLists <- colorMat[,keep]
  modAll <- modEP <- modEZ <- matrix(NA,nrow=nrow(enrichLists),ncol=nrow(MEmat))
  
  if (net == 1) { i.inds <- c(2,3,4,8,9,11,13,14,15,16,17,18) } ## Only the ones that passed 2/3 independent replication criteria
  if (net == 2) { i.inds <- c(1:nrow(MEmat)) } ## Keeping all of the modules in adult, for exploratory purposes
  ## Find adjusted multiple R^2 for each gene withe ach categorical variable
  for (i in i.inds) {
    for (j in c(1:17)) {
      testpath <- as.character(geneTable[geneTable[,"multiData$netData$finalColors"]==substr(rownames(MEmat),3,100)[i],"keptSymbols[keep]"]) ## Genes in module
      refpath <- colnames(enrichLists)[enrichLists[j,]=="black"] ## Candidate gene set
      testbackground <- as.character(geneTable[,"keptSymbols[keep]"])
      refbackground <- unique(na.omit(refbackMat[,j]))
      
      oraout <- ORA(testpath,refpath,testbackground,refbackground)
      modEP[j,i] <- oraout[2]
      modEZ[j,i] <- oraout[1]
      modAll[j,i] <- paste(oraout,collapse="\t")
    }
  }
  rownames(modAll) <- rownames(modEP) <- rownames(modEZ) <- rownames(enrichLists)
  rownames(modAll) <- paste(paste(rownames(modAll),paste(names(oraout),collapse="\t")))
  colnames(modEP) <- paste(substr(rownames(MEmat),3,100),sep="")
  
  rownames(MESigs) <- c("Region Adjusted r","Region KW p","Age bicor","Age p","RIN bicor","RIN p")
  
  ## Eigengene table containing
  ## Module Color, Module Label, p-value/grouping for each trait, enrichment statistics for lists, and eigengene values
  MEtable <- cbind(rownames(MEmat),t(MESigs),t(modEP),t(modEZ),t(modAll),MEmat)
  multiData$netData$MEtable <- MEtable

  ## Replot the dendrogram and interesting relationships with these colors
  pdf(paste("./figures/MatchedTOMPlot",multiData$name,".pdf",sep=""),width=6.8,height=3.4)
  keepind <- c(1,2,9,10,7,8,11,12,15)
  mColorh1 <- cbind(multiData$netData$finalColors,
                    t(thesecolors[keepind,keep]),t(multiData$genecols))
  mLabelh1 <- c("Module",rownames(colorMat)[keep],rownames(multiData$genecols))
  plotDendroAndColors(multiData$netData$TOMdendrogram, mColorh1, groupLabels = mLabelh1,addGuide=TRUE,dendroLabels=FALSE,main=multiData$name)
  dev.off()
  
  ## Replot the dendrogram with only protein-coding and non-grey genes
  
  keep1 <- !(multiData$netData$finalColors=="grey")#&(colorMat["Protein Coding",]=="blue")
  
  pdf(paste("./figures/MatchedTOMPlotPCNoGrey",multiData$name,".pdf",sep=""),width=6.8,height=3.4)
  keepind <- c(9,10,7,8,11,12,1,2,5,6,15)
  mColorh1 <- cbind(multiData$netData$finalColors[keep1],
                    t(thesecolors[keepind,keep][,keep1]),t(multiData$genecols[,keep1]))
  mLabelh1 <- c("Module",rownames(colorMat)[keep],rownames(multiData$genecols))
  
  TOMfile <- paste("./data/signed",multiData$name,"-block.1.RData",sep="")
  load(TOMfile)
  PCgeneTree= hclust(1-as.dist(as.matrix(TOM)[keep1,keep1]),method="average");
  
  multiData$netData$TOMdendroPC <- PCgeneTree
  multiData$netData$PCcolors <- mColorh1
  multiData$netData$PCcolorsname <- mLabelh1
  
  plotDendroAndColors(PCgeneTree, mColorh1, groupLabels = mLabelh1,addGuide=TRUE,dendroLabels=FALSE,main=multiData$name,cex.colorLabels=0.5)
  dev.off()
  
  multiData$netData$TOMdendroPC <- PCgeneTree
  multiData$netData$PCcolors <- mColorh1
  multiData$netData$PCcolorsname <- mLabelh1
  
  save(multiData,file=paste("./tables/DataAggregate",multiData$name,".Rdata",sep=""))
  write.table(geneTable,file=paste("./tables/kMEtable_",multiData$name,".txt",sep=""),row.names=TRUE,col.names=TRUE,quote=TRUE,sep=",")
  write.table(MEtable,file=paste("./tables/MEtable_",multiData$name,".txt",sep=""),row.names=TRUE,col.names=TRUE,quote=TRUE,sep=",")
}

