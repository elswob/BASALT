### Ben Elsworth July 2015

# source("http://bioconductor.org/biocLite.R")
# biocIfnot <- function (packageName){
#   if (!(require(packageName, character.only=TRUE))) biocLite(packageName)
# }
# installifnot <- function (packageName){
#   if (!(require(packageName)))
#     install.packages(packageName, dependencies=TRUE, repos="http://cran.rstudio.com/")
# }
#installifnot("ggplot2")
#installifnot("gplots")
#installifnot("heatmap.plus")
#installifnot("GMD")
#installifnot("reshape")
#biocIfnot("limma")
#biocIfnot("org.Hs.eg.db")
#biocIfnot("genefu")
#biocIfnot("ctc")

library(ggplot2)
library(limma)
library(gplots)
library(org.Hs.eg.db)
library(genefu)
library(heatmap.plus)
library(reshape)
#library(GMD)

####################################################################
#' Set up the analysis
setup=function(){
  #load PAM50 genes
  p50<<-c('ACTR3B','ANLN','BAG1','BCL2','BIRC5','BLVRA','CCNB1','CCNE1','CDC20','CDC6','CDCA1','CDH3','CENPF','CEP55','CXXC5','EGFR','ERBB2','ESR1','EXO1','FGFR4','FOXA1','FOXC1','GPR160','GRB7','KIF2C','KNTC2','KRT14','KRT17','KRT5','MAPT','MDM2','MELK','MIA','MKI67','MLPH','MMP11','MYBL2','MYC','NAT1','ORC6L','PGR','PHGDH','PTTG1','RRM2','SFRP1','SLC39A6','TMEM45B','TYMS','UBE2C','UBE2T')

  #PAM50 gene groups
  ER_Signalling<<-c('BAG1','BCL2','BLVRA','CXXC5','ESR1','FOXA1','GPR160','MAPT','MLPH','NAT1','PGR1','SLC39A6')
  Growth_Factor_Signalling<<-c('EGFR','ERBB2','FGFR4','GRB7')
  Proliferation<<-c('ANLN','BIRC5','CCNB1','CCNE1','CDC20','CDC6','CDCA1','CENPF','CEP55','EX01','KIF2C','KNTC2','MELK','MKI67','MYBL2','ORC6','PTTG1','RRM2','TYMS','UBE2C','UBE2T')
  Invasion<<-c('MMP11')
  Miscellaneous<<-c('ACTR3B','MDM2','TMEM4B')
  Basal_pehnotype<<-c('CDH3','FOXC1','KRT14','KRT17','KRT5','MIA','MYC','PHGDH','SFRP1')
  Total<<-c('BAG1','BCL2','BLVRA','CXXC5','ESR1','FOXA1','GPR160','MAPT','MLPH','NAT1','PGR1','SLC39A6','EGFR','ERBB2','FGFR4','GRB7','ANLN','BIRC5','CCNB1','CCNE1','CDC20','CDC6','CDCA1','CENPF','CEP55','EX01','KIF2C','KNTC2','MELK','MKI67','MYBL2','ORC6','PTTG1','RRM2','TYMS','UBE2C','UBE2T','MMP11','ACTR3B','MDM2','TMEM4B','CDH3','FOXC1','KRT14','KRT17','KRT5','MIA','MYC','PHGDH','SFRP1')
  cat("Total in PAM50 groups = ",length(Total),"\n")

  #set location of original inputDir
  oldInputDir<<-inputDir

  #set percentage of samples required for a gene to be kept
  sampleNum<<-5

  #set pam50 confidence cutoff
  pamCC<<-0.75

  #read in the data
  p_file=paste(inputDir,"/",inputFile,sep="")
  cat("Reading ",p_file,"\n")

  m<-read.delim(p_file,header=T, sep = "\t")
  print(head(m[,0:5]))
  print(dim(m))

  #re-order to match pam50 gene groups
  m=m[match(Total, m[,1]),]
  print("Ordered by pam50 groups")
  #remove any rows not matched
  m=m[!is.na(m[,1]),]
  print(m[0:nrow(m),0:5])

  #print(head(m[,0:5]))

  #check for redundant names and create rownames from first column
  m=m[!duplicated(m[,1]),]
  print(dim(m))
  rownames(m)=m[,1]
  m[1]=NULL
  print(head(m[,0:5]))

  #identify missing PAM50 genes
  #cat(rownames(m))
  missing<-setdiff(sort(p50),sort(rownames(m)))
  cat("\nThe PAM50 genes",missing,"are missing from the data set\n")

  #remove genes with no values in less than x% of samples
  cat("Total dimensions = ",dim(m),"\n")
  mo=m
  print("Filtering genes based on percentage:")
  print(sampleNum)
  m=(m[rowSums(m != 0)>(ncol(m)/100)*(100-sampleNum),])
  #m=(m[rowSums(!is.na(m))>(ncol(m)/100)*(100-sampleNum),])
  cat("NA genes removed = ",dim(m),"\n")
  missing<-setdiff(sort(rownames(mo)),sort(rownames(m)))
  if(length(missing)>0){
    cat("The genes",missing,"have been removed from the data set\n")
  }else{
    cat("No genes were removed for being under represented\n")
  }

  #cat("Zero genes removed = ",dim(m),"\n")
  #m=(m[colSums(m != 0)>0,])
  #cat("Zero cells removed = ",dim(m),"\n")

  #run on a subset or all
  #m<-m[,0:10]
  #random number
  #m=m[,sort(sample(ncol(m),100))]

  m<<-m[order(rownames(m)),]

  #create master df
  master_df<<-data.frame(Sample=colnames(m))
}
#' Define the row colours for the heatmap.
#'
#' @param m A data frame with gene symbols as row names.
heatmap_row_colours=function(m){
  #get pam50 gene group counts
  g1=rownames(m) %in% ER_Signalling
  g2=rownames(m) %in% Growth_Factor_Signalling
  g3=rownames(m) %in% Proliferation
  g4=rownames(m) %in% Invasion
  g5=rownames(m) %in% Miscellaneous
  g6=rownames(m) %in% Basal_pehnotype
  gNames = c('ER Signalling','Growth Factor Signalling','Proliferation','Invasion','Miscellaneous','Basal pehnotype')

  gCol=c('gray', 'blue', 'black','red','orange','purple')

  #initialise vector
  gVec=c(rep("black",nrow(m)))

  gVec = replace(gVec,which(g1),"gray")
  gVec = replace(gVec,which(g2),"blue")
  gVec = replace(gVec,which(g3),"black")
  gVec = replace(gVec,which(g4),"red")
  gVec = replace(gVec,which(g5),"orange")
  gVec = replace(gVec,which(g6),"purple")
  #print(gVec)

  gLegText=character()
  gLegCol=character()
  gCount=0
  for (g in gCol){
    gCount=gCount+1
    if(g %in% gVec){
      gLegText = c(gLegText,gNames[gCount])
      gLegCol = c(gLegCol,g)
    }
  }
  #print(gLegText)
  #print(gLegCol)
  output<-list(gVec,gLegText,gLegCol)
  return(output)
}

#' Define the column colours for the heatmap.
#'
#' @param m A data frame with gene symbols as row names.
heatmap_col_colours=function(m){
  gCol=c('chocolate2','coral1','cyan','darkgoldenrod1','chartreuse4','chocolate','darkgray','darkmagenta','darkred','aquamarine','aquamarine4','beige','bisque4','black','blue','blueviolet','brown','burlywood','cadetblue1','chartreuse')  #print("# of col colours")
  #print(length(colData))
  c=colnames(m)
  cs=unique(sub("_.*","",c))
  csn=as.numeric(as.factor(sub("_.*","",c)))
  gVec = gCol[csn]
  gLegText=unique(as.factor(sub("_.*","",c)))
  gLegCol=unique(gVec)
  output<-list(gVec,gLegText,gLegCol)
  return(output)
}

#heatmap function
#' Generate the heatmaps.
#'
#' @param m A data frame with gene symbols as row names.
#' @param title The title for the heatmap
#' @param pam.res The PAM50 results file
makeHeatmap=function(m,title,pam.res){

  rowData=heatmap_row_colours(m)
  rgVec=rowData[[1]]
  rgLegText=rowData[[2]]
  rgLegCol=rowData[[3]]

  colData=heatmap_col_colours(m)
  cgVec=colData[[1]]
  cgLegText=colData[[2]]
  cgLegCol=colData[[3]]

  pScores=pam.res$Confidence

  print(head(m[,0:5]))

  #pdf(paste(inputDir,"/heatmap_plus.pdf",sep=""))

  subtypeCols=c("green","blue","red","orange","purple")
  g=grey.colors(10,start=0.9,end=0)
  gCols=g[pScores*10]

  pCalls=as.numeric(as.factor(pam.res$Call))
  pCols=subtypeCols[pCalls]
  pText=unique(as.factor(pam.res$Call))
  pCol=unique(pCols)

  cMat = cbind(cgVec,pCols,gCols)
  colnames(cMat)=c("Sample","Call","Confidence")
  rMat = as.matrix(t(rgVec))
  #rownames(rMat)=c("PAM50 group")

  print(head(cMat))
  print(dim(m))
  print(nrow(rMat))
  print(dim(as.matrix(cMat)))
  print(dim(as.matrix(rgVec)))

  #Set a working directory for output files
  #setwd("/Users/ogriffit/git/biostar-tutorials/Heatmaps")

  #Define custom dist and hclust functions for use with heatmaps
  mydist=function(c) {dist(c,method="euclidian")}
  myclust=function(c) {hclust(c,method="average")}

  #Create heatmap using custom heatmap.3 source code loaded above
  main_title="PAM50 subtype classification"
  par(cex.main=1)
  #col.pal <- brewer.pal(9,"YlOrRd")
  #col.pal <- brewer.pal(9,"PuRd")
  #col.pal=colorRampPalette(c("white", "orange", "red"),bias=1)
  col.pal=colorRampPalette(c("green", "white", "red"),bias=1.5)
  heatmap.3(
    m,
    col=col.pal,
    hclustfun=myclust,
    distfun=mydist,
    na.rm = TRUE,
    scale="none",
    dendrogram="both",
    margins=c(6,12),
    Rowv=TRUE,
    Colv=TRUE,
    ColSideColors=as.matrix(cMat),
    RowSideColors=rMat,
    symbreaks=FALSE,
    key=TRUE,
    symkey=FALSE,
    density.info="none",
    trace="none",
    main=main_title,
    labCol=FALSE,
    cexRow=1,
    ColSideColorsSize=6,
    RowSideColorsSize=1,
    KeyValueName="log2 Expression"
  )

  legend(
    0.9, 0.93,
    #"topright",      # location of the legend on the heatmap plot
    legend = cgLegText, # category labels
    col = cgLegCol,  # color key
    lty= 1,             # line style
    lwd = 5,            # line width
    cex = 0.4,
    #inset = 0.02
    title="Sample"

  )
  legend(
    0.1, 0.8,
    legend = pText, # category labels
    col = pCol,  # color key
    lty= 1,             # line style
    lwd = 5,            # line width
    cex = 0.4,
    inset = 0.1,
    title="Subtype Call"
  )
  legend(
    0.01, 0.35,
    legend = rgLegText, # category labels
    col = rgLegCol,  # color key
    lty= 1,             # line style
    lwd = 5,            # line width
    cex = 0.4,
    title="PAM50 gene type"
  )
}


########### PAM50 standard #################

#' Run the standard PAM50 analysis.
run_p50=function(){

  cat("\n ---- Running standard PAM50 analysis ----\n")
  #print(head(m))

  outDir=paste(inputDir,"/PAM50",sep="")
  dir.create(outDir,showWarnings = FALSE)

  #get pam50 set
  m<-m[rownames(m) %in% p50,]
  print(head(m[0:5]))

  #m[m==0]<-NA

  write.table(m,file=paste(outDir,"/pam50_un-normalised.txt",sep=""),sep="\t",quote=F,row.names=T,col.names=NA)

  #don't need to median normalise as does this by default but will do it for the heatmap
  x<-apply(m,1,median,na.rm=T)
  m_n<-m-x

  print(head(m)[,0:5])

  #write to file
  final_p_file=paste(outDir,"/pam50_median_normalised.txt",sep="")
  write.table(m_n,file=final_p_file,sep="\t",quote=F,row.names=T,col.names=NA)

  #inputFile="pam50_median_normalised_10_each.txt"
  inputFile<<-"pam50_median_normalised.txt"
  short<<-"res"
  inputDir <<- paste(inputDir,"/PAM50",sep="")

  ### Run PAM50
  print("Running PAM50")
  bio_sub()
  pam50_wrapper(inputDir,inputFile,short)

  ### Plot the output
  pam.result.file <- paste(outDir,"/",short,"_pam50scores.txt",sep="")
  pam.res <<- read.delim(pam.result.file, stringsAsFactors=F, row.names=1, check=F)

  #find calls with low confidence
  levels(pam.res$Call) <- c(levels(pam.res$Call), "LC")
  pam.res$Call[pam.res$Confidence<pamCC]="LC"

  #heatmaps
  pdf(paste(outDir,"/heatmap_raw_data.pdf",sep=""))
  h=makeHeatmap(m,"PAM50 raw heatmap",pam.res)
  dev.off()
  png(paste(outDir,"/heatmap_raw_data.png",sep=""))
  h=makeHeatmap(m,"PAM50 raw heatmap",pam.res)
  dev.off()

  pdf(paste(outDir,"/heatmap_median_centered.pdf",sep=""))
  h=makeHeatmap(m_n,"PAM50 median centered heatmap",pam.res)
  dev.off()
  png(paste(outDir,"/heatmap_median_centered.png",sep=""))
  h=makeHeatmap(m_n,"PAM50 median centered heatmap",pam.res)
  dev.off()

  #add to master file
  master_df<<-merge(master_df,pam.res[6],by.x="Sample",by.y="row.names")
  names(master_df)[names(master_df)=="Call"]<<-"PAM50"
  print(head(master_df))

  x=pam.res$Confidence; names(x) <- rownames(pam.res)
  qplot(x, main="Classification into subtypes", xlab="Confidence")

  pdf(paste(outDir,"/classification_plot_grouped.pdf",sep=""))
  g<-ggplot(data = pam.res, aes(x = sub("_.*","",rownames(pam.res)), fill = Call)) + geom_bar(position="fill") + labs(title = "PAM50 classification counts", y = "Classification Percentage", x = "Cell type", fill = "PAM50 Subtype") + theme(text = element_text(size=10), axis.text.x = element_text(angle = 45, hjust = 1))
  print(g)
  dev.off()

  qplot(Confidence, Call, data=pam.res, col=Call, main="Classification into subtypes") + geom_vline(xintercept=0.9, col="red", lty=2)
  g <- ggplot(pam.res, aes(factor(rownames(pam.res)), Confidence)) + geom_bar(stat = "identity",aes(fill = Call)) + labs(title = "PAM50 classification", y = "Classification Confidence", x = "Cell type", fill = "PAM50 Subtype") + theme(text = element_text(size=10), axis.text.x = element_text(angle = 45, hjust = 1))
  print(g)
  pdf(paste(outDir,"/classification_plot_ungrouped.pdf",sep=""))
  print(g)
  dev.off()

  #plot the correlation coefficients
  cor_plot(pam.res,".",short)
}

################# SCMGENE #######################

#' Run SCMGENE functions
run_scmgene=function(){
  cat("\n ---- Running SCMGENE analysis ----\n")

  outDir=paste(inputDir,"/scmgene",sep="")
  dir.create(outDir,showWarnings = FALSE)

  #create annotation df of gene symbol and EntrezGene.ID
  anno_df<<-select(org.Hs.eg.db, rownames(m), c("ENTREZID"), "ALIAS")
  #remove duplicates
  anno_df=anno_df[!duplicated(anno_df[,1]),]
  print(head(anno_df))
  #rename columns
  colnames(anno_df)[1]="probe"
  colnames(anno_df)[2]="EntrezGene.ID"
  rownames(anno_df)=anno_df$probe
  anno_m<<-as.matrix(anno_df)
  print(head(anno_m))

  #create model
  load("/Users/ben/Software/PAM50/SCMGENE/data/EXPO.RData")
  modgene <- lapply(scmod1.robust$mod, function(x) { return(x[1, , drop=FALSE]) })
  print(modgene)
  pdf(paste(outDir,"/scmgene_fit_EXPO.pdf",sep=""), width=7, height=7)
  tt <- subtype.cluster(module.ESR1=modgene$ESR1, module.ERBB2=modgene$ERBB2,module.AURKA=modgene$AURKA, data=data, annot=annot, do.mapping=FALSE, do.scale=TRUE,rescale.q=0.05, plot=TRUE, filen=paste(outDir,"/",sprintf("scmgene_model_EXPO"),sep=""))
  dev.off()
  scmgene.expo <<- tt$model

  #flip rows and columns of expression data
  m_flip=t(m)

  #rorS
  #print("#### rorS ####")
  #print(m[0:5,0:5])
  #r=rorS(m_flip, verbose=F, annot=anno_m,do.mapping = T)
  #print(r)

  # subtype clustering with SCMGENE
  cat("--- SCMGENE ---\n")
  pdf(paste(outDir,"/scmgene_classif.pdf", sep=""), width=7, height=7)
  #sc.out<<-subtype.cluster.predict(sbt.model=scmgene.expo, data=m_flip,annot=anno_m,do.mapping=TRUE, plot=TRUE, verbose=TRUE, logged2=TRUE)
  sc.out <<- subtype.cluster.predict(sbt.model=scmgene.robust, data=m_flip, annot=anno_m, do.mapping=TRUE, verbose=TRUE, plot=TRUE)
  dev.off()
  print(table(sc.out$subtype2))
  scm_df=as.data.frame(sc.out$subtype2)
  write.table(scm_df,paste(outDir,"/scmgene_out.txt",sep=""),sep="\t",quote=F,row.names=T,col.names=NA)
  #add to master file
  master_df<<-merge(master_df,scm_df,by.x="Sample",by.y="row.names")
  names(master_df)[names(master_df)=="sc.out$subtype2"]<<-"SCMGENE"
  print(head(master_df))

  # subtype clustering with SCMOD2
  cat("--- SCMOD2 ---\n")
  pdf(paste(outDir,"/scmgene_scmod2.pdf", sep=""), width=7, height=7)
  sc.mod2 <<- subtype.cluster.predict(sbt.model=scmod2.robust, data=m_flip, annot=anno_m, do.mapping=TRUE, verbose=TRUE, plot=TRUE)
  dev.off()
  #add to master file
  scm__mod2_df=as.data.frame(sc.mod2$subtype)
  master_df<<-merge(master_df,scm__mod2_df,by.x="Sample",by.y="row.names")
  names(master_df)[names(master_df)=="sc.mod2$subtype"]<<-"SCMGENE_SCMOD2"
  print(head(master_df))

  # subtype clustering with SCMOD1
  cat("--- SCMOD1 ---\n")
  pdf(paste(outDir,"/scmgene_scmod1.pdf", sep=""), width=7, height=7)
  sc.mod1 <<- subtype.cluster.predict(sbt.model=scmod1.robust, data=m_flip, annot=anno_m, do.mapping=TRUE, verbose=TRUE, plot=TRUE)
  dev.off()
  #add to master file
  scm__mod1_df=as.data.frame(sc.mod1$subtype)
  master_df<<-merge(master_df,scm__mod1_df,by.x="Sample",by.y="row.names")
  names(master_df)[names(master_df)=="sc.mod1$subtype"]<<-"SCMGENE_SCMOD1"
  print(head(master_df))

  #pam50 scale
  cat("--- SCMGENE PAM50 scale---\n")
  mysbt.pam50 <<- intrinsic.cluster.predict(sbt.model=pam50.scale, data=m_flip, annot=anno_m, do.mapping=TRUE, verbose=TRUE)
  print(table(mysbt.pam50$subtype))
  scm_p50=as.data.frame(mysbt.pam50$subtype)
  write.table(scm_p50,paste(outDir,"/scmgene_pam50_scale_out.txt",sep=""),sep="\t",quote=F,row.names=T,col.names=NA)
  #add to master file
  master_df<<-merge(master_df,scm_p50,by.x="Sample",by.y="row.names")
  names(master_df)[names(master_df)=="mysbt.pam50$subtype"]<<-"SCMGENE_P50_scale"
  print(head(master_df))

  #pam50 robust
  cat("--- SCMGENE PAM50 robust---\n")
  mysbt.pam50 <<- intrinsic.cluster.predict(sbt.model=pam50.robust, data=m_flip, annot=anno_m, do.mapping=TRUE, verbose=TRUE)
  print(table(mysbt.pam50$subtype))
  scm_p50=as.data.frame(mysbt.pam50$subtype)
  write.table(scm_p50,paste(outDir,"/scmgene_pam50_robust_out.txt",sep=""),sep="\t",quote=F,row.names=T,col.names=NA)
  #add to master file
  master_df<<-merge(master_df,scm_p50,by.x="Sample",by.y="row.names")
  names(master_df)[names(master_df)=="mysbt.pam50$subtype"]<<-"SCMGENE_P50_robust"
  print(head(master_df))

  # ssp2006 clustering
  cat("--- SCMGENE SSP2006---\n")
  mysbt.ssp2006 <<- intrinsic.cluster.predict(sbt.model=ssp2006.robust, data=m_flip, annot=anno_m, do.mapping=TRUE, verbose=TRUE)
  print(table(mysbt.ssp2006$subtype))
  scm_ssp2006=as.data.frame(mysbt.ssp2006$subtype)
  write.table(scm_ssp2006,paste(outDir,"/scmgene_ssp2006_out.txt",sep=""),sep="\t",quote=F,row.names=T,col.names=NA)
  #add to master file
  master_df<<-merge(master_df,scm_ssp2006,by.x="Sample",by.y="row.names")
  names(master_df)[names(master_df)=="mysbt.ssp2006$subtype"]<<-"SCMGENE_SSP2006"
  print(head(master_df))

  # ssp2003 clustering
  cat("--- SCMGENE SSP2003---\n")
  mysbt.ssp2003 <- intrinsic.cluster.predict(sbt.model=ssp2003.robust, data=m_flip, annot=anno_m, do.mapping=TRUE, verbose=TRUE)
  print(table(mysbt.ssp2003$subtype))
  scm_ssp2003=as.data.frame(mysbt.ssp2003$subtype)
  write.table(scm_ssp2003,paste(outDir,"/scmgene_ssp2003_out.txt",sep=""),sep="\t",quote=F,row.names=T,col.names=NA)
  #add to master file
  master_df<<-merge(master_df,scm_ssp2003,by.x="Sample",by.y="row.names")
  names(master_df)[names(master_df)=="mysbt.ssp2003$subtype"]<<-"SCMGENE_SSP2003"
  print(head(master_df))

  #########################################################

  #plot the pam50 robust data
  pam.result.file <- paste(outDir,"/scmgene_pam50_robust_out.txt",sep="")
  pam.res <<- read.delim(pam.result.file, stringsAsFactors=F, row.names=1, check=F)
  colnames(pam.res)[1]="Call"

  pdf(paste(outDir,"/pam50_robust_classification_plot_grouped.pdf",sep=""))
  g<-ggplot(data = pam.res, aes(x = sub("_.*","",rownames(pam.res)), fill = Call)) + geom_bar(position="fill") + labs(title = "SCMGENE PAM50 robust classification counts", y = "Classification Percentage", x = "Cell type", fill = "PAM50 Subtype") + theme(text = element_text(size=10), axis.text.x = element_text(angle = 45, hjust = 1))
  print(g)
  dev.off()

  #plot the scmgene data
  pam.result.file <- paste(outDir,"/scmgene_out.txt",sep="")
  pam.res <<- read.delim(pam.result.file, stringsAsFactors=F, row.names=1, check=F)
  colnames(pam.res)[1]="Call"

  pdf(paste(outDir,"/scmgene_plot_grouped.pdf",sep=""))
  g<-ggplot(data = pam.res, aes(x = sub("_.*","",rownames(pam.res)), fill = Call)) + geom_bar(position="fill") + labs(title = "SCMGENE classification counts", y = "Classification Percentage", x = "Cell type", fill = "PAM50 Subtype") + theme(text = element_text(size=10), axis.text.x = element_text(angle = 45, hjust = 1))
  print(g)
  dev.off()

}

#run the functions, make sure to run p50 last as it resets the inputDir variable
#' Run the subtype analysis
#'
#' @param inputDir The directory containing the input file
#' @param inputFile A dataframe of expression data with row names as gene symbols and column names as unique sample IDs
#' @param short A short name for the analysis
run_analysis=function(d,f,s){
  inputDir<<-d
  inputFile<<-f
  short<<-s
  setup()
  run_scmgene()
  run_p50()
  #print master dataframe to file
  write.table(master_df,paste(oldInputDir,"/subtype_summary.tsv",sep=""),sep="\t",quote=F,row.names=F)
}
