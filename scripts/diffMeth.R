#!/usr/bin/env Rscript

#Author: Johanna Klughammer
#Date: 26.07.2015

#options
options(error=quote({message(paste0("Error. Possibly due to missing replicates. Continuing to the end anyways."))})) #prevents RefFreeDMA from running on in case of error.

#libraries
suppressPackageStartupMessages(library(VennDiagram))
suppressPackageStartupMessages(library(hexbin))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(simpleCache))
suppressPackageStartupMessages(library(limma))
suppressPackageStartupMessages(library(MASS))
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(reshape2))


#get parameters
args=commandArgs(trailingOnly = TRUE)
wd=args[1]
RRBSdir=args[2]
species=args[3]
dedRef_id=args[4]
sampleAnnotation=args[5]
comp_col=args[6]
groups_col=args[7]
nTopDiffMeth=args[8]
scripts=args[9]
motif=args[10]
uc_tag=args[11]

source(paste0(scripts,"/ggbiplot_custom.R"))

#initialize
setwd(wd)
theme_set(theme_bw())
options(stringsAsFactors=FALSE)

system(paste0("mkdir -p ",getwd(),"/RCache/"))
setCacheDir(paste0(getwd(),"/RCache/"))

out_dir=paste0(RRBSdir,"/diffMeth_",motif)
system(paste0("mkdir -p ",out_dir))

stats=data.frame()


#--------------------------------
# functions		
#--------------------------------

# collect methylation tables
collectMethData=function(dir,motif){
  RRBSfiles=system(paste0("ls ",dir, "/*/*/*",motif,"Methylation*.bed"), intern=TRUE)
  RRBSfiles=RRBSfiles[!grepl(uc_tag,RRBSfiles)]
  Meth_bed=NULL
  for (file in RRBSfiles){
    if (grepl("__",file)==TRUE){
	type=unlist(strsplit(file,"__|\\.bed"))[3]}else{
	type=unlist(strsplit(file,"Methylation_|\\.bed"))[2]}
	
	print(type)
    input=fread(file)
    input=unique(input)
    input[,V4:=gsub("'","",V4),]
    spl=unlist(strsplit(input$V4,"/"))
    cov.idx=seq(1,length(spl),by=2)
    input[,V7:=as.numeric(spl[cov.idx+1]),]
    input[,V8:=as.numeric(spl[cov.idx]),]
    #input[,V5:=V5/10,] #calculate new from ration because of bug in chh and chg methylation calling in biseqMehcalling (not for CpG meth though)
    input[,V5:=V8/V7*100,]
    input[,V4:=NULL,]
    setnames(input,names(input),c("chr","start","end",paste0(type,".meth"),"strand",paste0(type,".cov"),paste0(type,".reads_meth")))
    if (is.null(Meth_bed)){Meth_bed=input}
    else {Meth_bed=merge(Meth_bed,input,by=c("chr","start","end","strand"),all=TRUE)}    
  }
  message("Done.")
  return(Meth_bed)
}

#functions for differential methylation calling from RnBeads
#limmaP
limmaP<-function (X, inds.g1, inds.g2 = -inds.g1, adjustment.table = NULL, 
                  fun.conversion = rnb.beta2mval, paired = FALSE) 
{
  if (is.logical(inds.g1)) 
    inds.g1 <- which(inds.g1)
  if (is.logical(inds.g2)) 
    inds.g2 <- which(inds.g2)
  n.g1 <- length(inds.g1)
  n.g2 <- length(inds.g2)
  n <- n.g1 + n.g2
  if (!is.null(adjustment.table)) {
    if (!(is.data.frame(adjustment.table) && nrow(adjustment.table) == 
            n && (!any(is.na(adjustment.table))))) {
      stop("invalid value for adjustment.table")
    }
    m <- ncol(adjustment.table)
    if (m == 0) {
      adjustment.table <- NULL
    }
    else {
      colnames(adjustment.table) <- paste0("x", 1:m, "x")
    }
  }
  ind.vec <- c(inds.g1, inds.g2)
  if (length(ind.vec) < 2) 
    stop("need at least two samples indices to compare")
  X.m <- fun.conversion(X[, ind.vec, drop = FALSE])
  
  df <- data.frame(xg = factor(rep(c("group1", "group2"), c(n.g1, 
                                                            n.g2)), levels = c("group1", "group2")))
  if (!is.null(adjustment.table)) {
    df <- cbind(df, adjustment.table)
  }
  if (paired) {
    if (n.g1 != n.g2) {
      stop("Could not conduct paired limma analysis: unequal groupsizes")
    }
    df$xp <- as.factor(rep(1:n.g1, 2))
  }
  formula.text <- paste0(c("~0", colnames(df)), collapse = "+")
  design.m <- model.matrix(as.formula(formula.text), data = df)
  colnames(design.m) <- make.names(colnames(design.m), unique = TRUE)
  colnames(design.m)[1:2] <- c("group1", "group2")
  fit <- limma::lmFit(X.m, design.m)
  contrasts.m <- makeContrasts(group1vs2 = group1 - group2, 
                               levels = design.m)
  fit <- limma::contrasts.fit(fit, contrasts.m)
  fit <- limma::eBayes(fit)
  return(fit$p.value[, "group1vs2"])
}

#rnb.beta2mval
rnb.beta2mval<-function (betas, epsilon = 1e-05) 
{
  if (!is.numeric(betas)) {
    stop("invalid value for betas")
  }
  if (!(is.numeric(epsilon) && length(epsilon) == 1 && (!is.na(epsilon)))) {
    stop("invalid value for epsilon")
  }
  if (epsilon < 0 || epsilon > 0.5) {
    stop("invalid value for epsilon; expected 0 <= epsilon <= 0.5")
  }
  betas[betas < epsilon] <- epsilon
  betas[betas > (1 - epsilon)] <- 1 - epsilon
  return(log2(betas/(1 - betas)))
}

#combine p values
combineTestPvalsMeth<-function (pvalues, testWeights = NULL, correlated = FALSE, methExpectedTestCorrelation = 0.8) 
{
  if (is.null(pvalues)) {
    return(NA)
  }
  if (!is.numeric(pvalues)) {
    #  logger.warning(c("Non numeric value for pvalues in combination:", pvalues))
    return(NA)
  }
  if (!is.null(testWeights)) {
    if (length(pvalues) != length(testWeights)) 
      stop("Number of items in <pvalues> and in <testWeights> must be identical if weights are to be used")
    if (sum(is.na(testWeights)) > 0) 
      stop("NA values are not permitted for the test weights")
    if (sum(testWeights < 0) > 0) 
      stop("Weights must be positive")
    testWeights = testWeights/sum(testWeights)
  }
  else {
    testWeights = rep(1/length(pvalues), length(pvalues))
  }
  pvalues[is.na(pvalues)] = 1
  if (length(pvalues) < 1) {
    return(NA)
  }
  else if (length(pvalues) == 1) {
    return(pvalues[1])
  }
  else if (length(pvalues) > sqrt(.Machine$integer.max)) {
    #  logger.info(c("Too many p-values to combine --> using subsampling"))
    nn <- trunc(sqrt(.Machine$integer.max))
    ss <- sample(length(pvalues), nn)
    pvalues <- pvalues[ss]
    testWeights <- testWeights[ss]
    testWeights <- testWeights/sum(testWeights)
  }
  if (correlated == FALSE & is.null(testWeights)) {
    tcombined = sum(-2 * log(pvalues))
    return(pchisq(tcombined, 2 * length(pvalues), lower.tail = FALSE))
  }
  else {
    r = ifelse(correlated, methExpectedTestCorrelation, 0)
    m = length(pvalues)
    M.Fm = sum(-2 * log(pvalues) * testWeights)
    var.M.Fm = 4 * sum(testWeights^2)
    ij.pairs <- expand.grid(1:m, 1:m)
    ij.pairs <- ij.pairs[ij.pairs[, 1] != ij.pairs[, 2], 
                         ]
    tw.i <- testWeights[ij.pairs[, 1]]
    tw.j <- testWeights[ij.pairs[, 2]]
    vv <- tw.i * tw.j * (3.25 * r + 0.75 * r^2)
    var.M.Fm <- var.M.Fm + sum(vv)
    nu = 8/var.M.Fm
    tcombined = nu * M.Fm/2
    return(pchisq(tcombined, nu, lower.tail = FALSE))
  }
}


#---------------------------------
#load data
#---------------------------------

#read sample annptation
sample_annot=fread(sampleAnnotation)
sample_annot[sample_annot==""]=NA
sample_annot=sample_annot[!grepl(uc_tag,Sample_Name),]

#combine meth bed files by caching (Nathan)
simpleCache(recreate=TRUE,paste0(motif,"_combinedMeth_",species),instruction="collectMethData(RRBSdir,motif)")
Meth_bed=get(paste0(motif,"_combinedMeth_",species))
#deduced geneome fragment .bed file
refBed_DT=fread(paste0("reduced/consensus/",RRBSdir,".bed"))
setnames(refBed_DT,names(refBed_DT),c("seqnames","start","end","meta","score","consN","mean_diffNts","max_diffNts","min_diffNts"))


#-----------------------------------------
#Visualize on individual motif level
#-----------------------------------------
#plot MDS on single motifs
meth.dist=dist(t(as.data.frame(Meth_bed)[,seq(5,ncol(Meth_bed),by=3)]),method="eucl")
meth.dist[meth.dist<=0]=0.00001
meth.mds=isoMDS(meth.dist)
MDS=data.table(Sample_Name=as.character(lapply(row.names(meth.mds$points),function(x){unlist(strsplit(x,"\\."))[1]})),MDS1=meth.mds$points[,1],MDS2=meth.mds$points[,2])
MDS_annot=merge(MDS,sample_annot,by="Sample_Name")
pdf(paste0(out_dir,"/",species,"_MDS_motif.pdf"),height=6,width=7.5)
ggplot(MDS_annot,aes(x=MDS1,y=MDS2,label=Sample_Name))+geom_point(aes(colour=factor(get(groups_col))))+geom_text(aes(colour=factor(get(groups_col))),size=3.5)+theme(legend.title=element_blank())+scale_x_continuous(expand=c(0.5,0.5))+scale_y_continuous(expand=c(0.5,0.5))
dev.off()

#plot PCA on single motifs
pca=prcomp(t(na.omit(as.data.frame(Meth_bed)[,seq(5,ncol(Meth_bed),by=3)])))
PC=data.table(Sample_Name=as.character(lapply(row.names(pca$x),function(x){unlist(strsplit(x,"\\."))[1]})),PC1=pca$x[,"PC1"],PC2=pca$x[,"PC2"])
PC_annot=merge(PC,sample_annot,by="Sample_Name")
pdf(paste0(out_dir,"/",species,"_PCA_motif.pdf"),height=6,width=7.5)
ggpl= ggplot(PC_annot,aes(x=PC1,y=PC2,label=Sample_Name))+geom_point(aes(colour=factor(get(groups_col))))+geom_text(aes(colour=factor(get(groups_col))),size=3.5)+theme(legend.title=element_blank())+scale_x_continuous(expand=c(0.5,0.5))+scale_y_continuous(expand=c(0.5,0.5))
print(ggpl) 
dev.off()

#-------------------------------------------------------
#Combine methylation scores by deduced genome fragment
#-------------------------------------------------------
Meth_bed[,meth.ID:=1:nrow(Meth_bed),]
meth_gr=with(Meth_bed,GRanges(chr, IRanges(start,end),strand=strand,meth.ID=meth.ID))

refBed_DT[,ref.ID:=1:nrow(refBed_DT),]
ref_gr <- with(refBed_DT,GRanges(seqnames,IRanges(start,end),ref.ID=ref.ID))

overlap_meth=data.table(as.data.frame(findOverlaps(meth_gr,ref_gr,type="within",select="all")))

overlap_meth[,meth.ID:=meth_gr$meth.ID[queryHits],]
overlap_meth[,ref.ID:=ref_gr$ref.ID[subjectHits],]

merge1_meth=merge(overlap_meth,refBed_DT,by="ref.ID",all.x=TRUE)
merge2_meth=merge(merge1_meth,Meth_bed,by="meth.ID",all.x=TRUE)


# combine methylation values uning a weighted mean
cols_dt=data.table(meth=grep("\\.meth",colnames(Meth_bed),value=TRUE),cov=grep("\\.cov",colnames(Meth_bed),value=TRUE))

cols_dt[,command:=paste0(meth,"=weighted.mean(",meth,",",cov,",na.rm=TRUE)",",",cov,"=mean(",cov,",na.rm=TRUE)"),]
command_weighted_mean=paste0("list(",paste0(cols_dt$command,collapse=","),")")

mean_meth=merge2_meth[,eval(parse(text=command_weighted_mean)),by=c("meta","score","consN","mean_diffNts","max_diffNts","min_diffNts")]
write.table(mean_meth,paste0(out_dir,"/",species,"_mean_meth.tsv"),quote=FALSE,sep="\t",row.names=FALSE)
#mean_meth1=fread(paste0(out_dir,"/",species,"_mean_meth.tsv"),sep="\t")

#-----------------------------------------
#Visualize on individual motif level
#-----------------------------------------
#PCA on fragments
noNA=as.data.table(na.omit(mean_meth[,c(1,seq(7,ncol(mean_meth),by=2)),with=FALSE]))
colNames=noNA$meta
t=t(noNA[,-1,with=FALSE])
colnames(t)=colNames

pca=prcomp(t)
PC=data.table(Sample_Name=as.character(lapply(row.names(pca$x),function(x){unlist(strsplit(x,"\\."))[1]})),PC1=pca$x[,"PC1"],PC2=pca$x[,"PC2"])
PC_annot=merge(PC,sample_annot,by="Sample_Name")
pdf(paste0(out_dir,"/",species,"_PCA_frag.pdf"),height=6,width=7.5)
ggpl=ggplot(PC_annot,aes(x=PC1,y=PC2,label=Sample_Name))+geom_point(aes(colour=factor(get(groups_col))))+geom_text(aes(colour=factor(get(groups_col))),size=3.5)+theme(legend.title=element_blank())+scale_x_continuous(expand=c(0.5,0.5))+scale_y_continuous(expand=c(0.5,0.5))
print(ggpl)
dev.off()

#biplot on fragments (with labels)
pcaSub=pca
pcaSub$x=pcaSub$x[order(row.names(pcaSub$x)),]
PC_annot=PC_annot[c(order(Sample_Name)),]

#PC1 and PC2
pcaSub$rotation=pca$rotation[order(sqrt((pca$rotation[,"PC1"])^2+(pca$rotation[,"PC2"])^2),decreasing=TRUE),][1:40,]
pdf(paste0(out_dir,"/",species,"_biplot_1-2_frag.pdf"),height=6,width=7)
ggbiplot(pcaSub,c(1,2),labels=PC_annot$Sample_Name,group=unlist(PC_annot[,groups_col,with=FALSE]),varname.size=2.5)
dev.off()
#PC3 and PC4
pcaSub$rotation=pca$rotation[order(sqrt((pca$rotation[,"PC3"])^2+(pca$rotation[,"PC4"])^2),decreasing=TRUE),][1:40,]
pdf(paste0(out_dir,"/",species,"_biplot_3-4_frag.pdf"),height=6,width=7)
ggbiplot(pcaSub,c(3,4),labels=PC_annot$Sample_Name,group=unlist(PC_annot[,groups_col,with=FALSE]),varname.size=2.5)
dev.off()

#biplot on fragments (without labels)
#PC1 and PC2
pcaSub$rotation=pca$rotation[order(sqrt((pca$rotation[,"PC1"])^2+(pca$rotation[,"PC2"])^2),decreasing=TRUE),][1:40,]
pdf(paste0(out_dir,"/",species,"_biplot_1-2_frag_noLab.pdf"),height=6,width=7)
ggbiplot(pcaSub,c(1,2),group=unlist(PC_annot[,groups_col,with=FALSE]),size=3.5)
dev.off()
#PC3 and PC4
pcaSub$rotation=pca$rotation[order(sqrt((pca$rotation[,"PC3"])^2+(pca$rotation[,"PC4"])^2),decreasing=TRUE),][1:40,]
pdf(paste0(out_dir,"/",species,"_biplot_3-4_frag_noLab.pdf"),height=6,width=7)
ggbiplot(pcaSub,c(3,4),group=unlist(PC_annot[,groups_col,with=FALSE]),size=3.5)
dev.off()

#MDS on fragments
meth.dist=dist(t(as.data.frame(mean_meth)[,seq(7,ncol(mean_meth),by=2)]),method="eucl")

meth.dist[meth.dist<=0]=0.00001
meth.mds=isoMDS(meth.dist)
MDS=data.table(Sample_Name=as.character(lapply(row.names(meth.mds$points),function(x){unlist(strsplit(x,"\\."))[1]})),MDS1=meth.mds$points[,1],MDS2=meth.mds$points[,2])
MDS_annot=merge(MDS,sample_annot,by="Sample_Name")
pdf(paste0(out_dir,"/",species,"_MDS_frag.pdf"),height=6,width=7.5)
ggpl=ggplot(MDS_annot,aes(x=MDS1,y=MDS2,label=Sample_Name))+geom_point(aes(colour=factor(get(groups_col))))+geom_text(aes(colour=factor(get(groups_col))),size=3.5)+theme(legend.title=element_blank())+scale_x_continuous(expand=c(0.5,0.5))+scale_y_continuous(expand=c(0.5,0.5))
print(ggpl)
dev.off()


#-------------------------------------------------------------------------
#prepare the replicate groupings for differential methylation analysis
#-------------------------------------------------------------------------

present_samples=unlist(lapply(cols_dt$meth,function(x){unlist(strsplit(x,".meth"))[1]}))
comp_matrix=sample_annot[!is.na(get(comp_col)),matrix(Sample_Name),by=comp_col][V1%in%present_samples]

setnames(comp_matrix,names(comp_matrix),c("group","sample"))
comp_matrix[,meth.cols:=paste0(sample,".meth"),]
comp_matrix[,cov.cols:=paste0(sample,".cov"),]
merge.cols=names(merge2_meth)
comp_matrix[,meth.idx:=unlist(lapply(meth.cols,function(x){which(merge.cols==x)})),]
comp_matrix[,cov.idx:=unlist(lapply(cov.cols,function(x){which(merge.cols==x)})),]
comp_matrix=split(comp_matrix,comp_matrix$group)

#-------------------------------------------------
#Visualize replicate concordance by Venn-diagrams
#-------------------------------------------------

for (i in 1:length(comp_matrix)){
	for (j in 1:floor(nrow(comp_matrix[[i]])-1)){

		scale=10^floor(log10(nrow(Meth_bed)))

		area1=round(nrow(na.omit(Meth_bed[,comp_matrix[[i]]$meth.cols[j],with=FALSE]))/scale,2)
		area2=round(nrow(na.omit(Meth_bed[,comp_matrix[[i]]$meth.cols[j+1],with=FALSE]))/scale,2)
		area.cross=round(nrow(na.omit(Meth_bed[,c(comp_matrix[[i]]$meth.cols[j],comp_matrix[[i]]$meth.cols[j+1]),with=FALSE]))/scale,2)

		pdf(paste0(out_dir,"/Venn_motif_",comp_matrix[[i]]$sample[j],"-",comp_matrix[[i]]$sample[j+1],"_x",scale,".pdf"),height=3,width=5.5)
		grid.newpage()
		draw.pairwise.venn(area1,area2,area.cross,category=c(comp_matrix[[i]]$sample[j],comp_matrix[[i]]$sample[j+1]),fill=c("blue","orange"),alpha=c(0.4,0.4),cat.dist=c(0.1,0.1),fontface=c("plain","plain","plain"),fontfamily=c("sans","sans","sans"),cat.fontface=c("plain","plain"),cat.fontfamily=c("sans","sans"),cat.pos=c(330,30),cat.col=c("blue","orange"),cat.cex=1.5,mar=c(0.05,0.05,0.05,0.05),cex=c(1.3,1.3,1.3))
		dev.off()

		stats=rbind(stats,data.frame(assay="diffMeth",condition=paste0(comp_matrix[[i]]$sample[j],"-",comp_matrix[[i]]$sample[j+1]),percent=area.cross/(area1+area2-area.cross),number=area.cross*scale,correlation=round(cor(na.omit(Meth_bed[,c(comp_matrix[[i]]$meth.cols[j],comp_matrix[[i]]$meth.cols[j+1]),with=FALSE]))[1,2],3)))
	}
}

#-----------------------------------------------------------------
#Calculate concordance in motif coverage between different samples
#-----------------------------------------------------------------
sampleCombies=as.data.table(t(combn(sample_annot$Sample_Name,2)))
sampleCombies[,overlap:=nrow(na.omit(Meth_bed[,c(paste0(V1,".meth"),paste0(V2,".meth")),with=FALSE])),by=1:nrow(sampleCombies)]

sampleCombies_extended=rbindlist(list(sampleCombies,with(sampleCombies,data.table(V1=V2,V2=V1,overlap=overlap))))
overlap_stats=sampleCombies_extended[,list(mean_overlap=mean(overlap),min_overlap=min(overlap),max_overlap=max(overlap),min_sample=V2[which.min(overlap)],max_sample=V2[which.max(overlap)]),by="V1"]
setnames(overlap_stats,"V1","Sample_Name")

write.table(overlap_stats,paste0(out_dir,"/",species,"_overlap_stats.tsv"),quote=FALSE,sep="\t",row.names=FALSE)

sampleCombies_extended[,V1:=factor(V1,levels=overlap_stats[order(mean_overlap)]$Sample_Name),]
sampleCombies_extended[,V2:=factor(V2,levels=overlap_stats[order(mean_overlap)]$Sample_Name),]
pdf(paste0(out_dir,"/",species,"_overlap_stats.pdf"),height=4,width=5)
ggplot(sampleCombies_extended,aes(x=V1,y=V2,fill=overlap))+geom_tile()+xlab("")+ylab("")+theme(axis.text.x  = element_text(angle=90, vjust=0.5))
dev.off()
#------------------------------------------------------------------
#differential methylation analysis using limmaP from RnBeads
#------------------------------------------------------------------

#calculate basic statistics for both groups
diffTab=merge2_meth[,c("meth.ID","meta","seqnames","start.x","end.x","start.y","end.y"),with=FALSE]
diffTab[,mean.meth.g1:=rowMeans(merge2_meth[,comp_matrix[[1]]$meth.idx,with=FALSE],na.rm=TRUE),]
diffTab[,mean.meth.g2:=rowMeans(merge2_meth[,comp_matrix[[2]]$meth.idx,with=FALSE],na.rm=TRUE),]
diffTab[,sd.meth.g1:=apply(merge2_meth[,comp_matrix[[1]]$meth.idx,with=FALSE],1,function(x){sd(x,na.rm=TRUE)}),]
diffTab[,sd.meth.g2:=apply(merge2_meth[,comp_matrix[[2]]$meth.idx,with=FALSE],1,function(x){sd(x,na.rm=TRUE)}),]
diffTab[,mean.cov.g1:=rowMeans(merge2_meth[,comp_matrix[[1]]$cov.idx,with=FALSE],na.rm=TRUE),]
diffTab[,mean.cov.g2:=rowMeans(merge2_meth[,comp_matrix[[2]]$cov.idx,with=FALSE],na.rm=TRUE),]
diffTab[,mean.mean.cov:=rowMeans(cbind(mean.cov.g1,mean.cov.g2),na.rm=TRUE),]
diffTab[,mean.mean.cov:=ifelse(is.na(mean.mean.cov),0,mean.mean.cov),]
diffTab[,diff.mean.meth:=mean.meth.g1-mean.meth.g2,]
eps=0.01 #to prevent deviding by 0
diffTab[,quot.mean.meth:=(mean.meth.g1+eps)/(mean.meth.g2+eps),]
diffTab[,quot.log2.mean.meth:=log2(quot.mean.meth),]

#calculate p-values for differential methylation
diffMeth_matrix=as.matrix(merge2_meth[,c(comp_matrix[[1]]$meth.idx,comp_matrix[[2]]$meth.idx),with=FALSE])
diffMeth_matrix=diffMeth_matrix/100
diffTab[,p.val:=limmaP(diffMeth_matrix,inds.g1=c(1:nrow(comp_matrix[[1]])),inds.g2=c((nrow(comp_matrix[[1]])+1):(nrow(comp_matrix[[1]])+nrow(comp_matrix[[2]]))),adjustment.table=NULL,paired=FALSE),]
diffTab[,p.val.adjust:=p.adjust(p.val,method="fdr"),]

#calulate combined rank for differential methylation for each individual motif
diffTab[,rank.diff:=rank(-abs(diff.mean.meth),na.last="keep",ties.method="min"),]
diffTab[,rank.quot:=rank(-abs(quot.log2.mean.meth),na.last="keep",ties.method="min"),]
diffTab[,rank.pval:=rank(p.val.adjust,na.last="keep",ties.method="min"),]
diffTab[,max.rank:=pmax(rank.diff,rank.quot,rank.pval),]
diffTab=diffTab[order(max.rank)]

write.table(diffTab,paste0(out_dir,"/",species,"_diff_meth_motif.tsv"),quote=FALSE,sep="\t",row.names=FALSE)
#diffTab=fread(paste0(out_dir,"/",species,"_diff_meth_motif.tsv"))

#combine pValues per deduced genome fragment
#diffTab_ref=diffTab[,list(sites=.N,p.val.adjust=combineTestPvalsMeth(p.val.adjust,testWeights=mean.mean.cov,correlated=TRUE),meth_hi=as.character(ifelse(mean(quot.mean.meth,na.rm=TRUE)>1,names(comp_matrix)[1],ifelse(mean(quot.mean.meth,na.rm=TRUE)<1,names(comp_matrix)[2],"tie"))),mean.meth.g1=mean(mean.meth.g1, na.rm=TRUE),mean.meth.g2=mean(mean.meth.g2, na.rm=TRUE),mean.cov.g1=mean(mean.cov.g1, na.rm=TRUE),mean.cov.g2=mean(mean.cov.g2, na.rm=TRUE),mean.mean.cov=mean(mean.mean.cov,na.rm=TRUE),diff.mean.meth=mean(diff.mean.meth,na.rm=TRUE),quot.log2.mean.meth=mean(quot.log2.mean.meth,na.rm=TRUE)),by=c("meta","seqnames","start.x","end.x")]

diffTab_ref=diffTab[,list(sites=.N,p.val.adjust=combineTestPvalsMeth(p.val.adjust,testWeights=mean.mean.cov,correlated=TRUE),meth_hi=as.character(ifelse(mean(diff.mean.meth,na.rm=TRUE)>0,names(comp_matrix)[1],ifelse(mean(diff.mean.meth,na.rm=TRUE)<0,names(comp_matrix)[2],"tie"))),mean.meth.g1=mean(mean.meth.g1, na.rm=TRUE),mean.meth.g2=mean(mean.meth.g2, na.rm=TRUE),mean.cov.g1=mean(mean.cov.g1, na.rm=TRUE),mean.cov.g2=mean(mean.cov.g2, na.rm=TRUE),mean.mean.cov=mean(mean.mean.cov,na.rm=TRUE),diff.mean.meth=mean(diff.mean.meth,na.rm=TRUE),quot.log2.mean.meth=mean(quot.log2.mean.meth,na.rm=TRUE)),by=c("meta","seqnames","start.x","end.x")]



#calulate combined rank for differential methylation for each deduced genome fragment
diffTab_ref[,rank.diff:=rank(-abs(diff.mean.meth),na.last="keep",ties.method="min"),]
diffTab_ref[,rank.quot:=rank(-abs(quot.log2.mean.meth),na.last="keep",ties.method="min"),]
diffTab_ref[,rank.pval:=rank(p.val.adjust,na.last="keep",ties.method="min"),]
diffTab_ref[,max.rank:=pmax(rank.diff,rank.quot,rank.pval),]

setnames(diffTab_ref,names(diffTab_ref),c("dedRef_ID","dedRef_concat","start","end","meth.sites","meth.pval_adj","meth.meth_hi","meth.meth_mean_g1","meth.meth_mean_g2","meth.cov_mean_g1","meth.cov_mean_g2","meth.cov_mean","meth.meth_diff","meth.meth_log2ratio","meth.rank_diff","meth.rank_ratio","meth.rank_pval","meth.rank_max"))

diffTab_ref.ordered=diffTab_ref[order(meth.rank_max)]
write.table(diffTab_ref.ordered,paste0(out_dir,"/",species,"_diff_meth.tsv"),quote=FALSE,sep="\t",row.names=FALSE)


#--------------------------------------------------------------------------------------------------------------------------------
#Visualize mean methylation values for the compared groups for individual motifs and for deduced genome fragments(scatter plots)
#--------------------------------------------------------------------------------------------------------------------------------

cov_trsh_l=c(0,8)
cov_trsh_u=c(Inf,200)

comp=c(names(comp_matrix)[1:2])
for (i in 1:length(cov_trsh_l)){
  sub_motif=na.omit(diffTab[mean.cov.g1>=cov_trsh_l[i]&mean.cov.g1>=cov_trsh_l[i]&mean.cov.g2<=cov_trsh_u[i]&mean.cov.g2<=cov_trsh_u[i],])
  sub_frag=na.omit(diffTab_ref.ordered[meth.cov_mean_g1>=cov_trsh_l[i]&meth.cov_mean_g1>=cov_trsh_l[i]&meth.cov_mean_g2<=cov_trsh_u[i]&meth.cov_mean_g2<=cov_trsh_u[i],])


  if (nrow(sub_motif)>0){  
    size=4
    top_diff_data=sub_motif[p.val.adjust<0.05&mean.cov.g1>8&mean.cov.g2>8][1:min((nrow(sub_frag)/10),1000)]#[1:1000]
    
    pdf(paste0(out_dir,"/diffMeth_",comp[1],"-",comp[2],"_",cov_trsh_l[i],"-",cov_trsh_u[i],"_motif.pdf"),width=6,height=5)
    ggp_motif=ggplot(sub_motif,aes(x=mean.meth.g1,y=mean.meth.g2))+ stat_binhex(bins=30,col="white")+geom_point(data=top_diff_data, aes(x=mean.meth.g1,y=mean.meth.g2),col="green",size=2,shape=21) + annotate("text",size=size,x=90,y=10,label=paste0("r = ",round(cor(sub_motif[,mean.meth.g1],sub_motif[,mean.meth.g2]),3),"\n cov>= ",cov_trsh_l[i], "\n cov <=",cov_trsh_u[i] ,"\nN = ",dim(sub_motif)[1])[1])+scale_fill_gradient(limits=c(0,nrow(sub_motif)/200),high="blue",low="white",na.value="blue")+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+xlab(paste0("% methylation ",comp[1]))+ylab(paste0("% methylation ",comp[2]))+xlim(0,100)+ylim(0,100)
    
    if (nrow(na.omit(top_diff_data)) < 10) {
      ggp_motif=ggp_motif+geom_point(data=top_diff_data, aes(x=mean.meth.g1,y=mean.meth.g2),col="green",size=2,shape=21)}
    print(ggp_motif)
    dev.off() 
  }
  if (nrow(sub_frag)>0){
    size=4
    top_diff_data=sub_frag[meth.pval_adj<0.05&meth.cov_mean_g1>8&meth.cov_mean_g2>8][1:min((nrow(sub_frag)/10),500)]#[1:500]
    
    pdf(paste0(out_dir,"/diffMeth_",comp[1],"-",comp[2],"_",cov_trsh_l[i],"-",cov_trsh_u[i],"_frag.pdf"),width=6,height=5)
    ggp_frag=ggplot(sub_frag,aes(x=meth.meth_mean_g1,y=meth.meth_mean_g2))+ stat_binhex(bins=30,col="white") +geom_point(data=top_diff_data, aes(x=meth.meth_mean_g1,y=meth.meth_mean_g2),col="green",size=2,shape=21) + annotate("text",size=size,x=90,y=10,label=paste0("r = ",round(cor(sub_frag[,meth.meth_mean_g1],sub_frag[,meth.meth_mean_g2]),3),"\n cov>= ",cov_trsh_l[i], "\n cov <=",cov_trsh_u[i] ,"\nN = ",dim(sub_frag)[1])[1])+scale_fill_gradient(limits=c(0,nrow(sub_frag)/200),high="blue",low="white",na.value="blue")+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+xlab(paste0("% methylation ",comp[1]))+ylab(paste0("% methylation ",comp[2]))+xlim(0,100)+ylim(0,100)
    if (nrow(na.omit(top_diff_data)) < 10) {
      ggp_frag=ggp_frag+geom_point(data=top_diff_data, aes(x=meth.meth_mean_g1,y=meth.meth_mean_g2),col="green",size=2,shape=21)}
    print(ggp_frag)
    dev.off()
   
  }
}

#---------------------------------------------------------------------------------------------------------
#Save sequences of to differentially methylated fragments for down-stream analysis (motif enrichment etc)
#---------------------------------------------------------------------------------------------------------
#For now only activat for CpG methylation, because in our samples there is no non-cpg differential methylation which leads to errors here
#TODO: activate this also for nonCpG methylation 

if (motif == "cpg" ){

top_threshold=as.numeric(nTopDiffMeth)
cov_thres=2

fasta=readDNAStringSet(paste0("reduced/consensus/",sub("_concat",".fa",RRBSdir)))
groups=unique(diffTab_ref.ordered$meth.meth_hi)


top1_all=diffTab_ref.ordered[meth.meth_hi==groups[1]&meth.cov_mean_g1>=cov_thres& meth.cov_mean_g2>=cov_thres&meth.pval_adj<=0.05]
top1=top1_all[1:min(top_threshold,nrow(top1_all))]$dedRef_ID

top2_all=diffTab_ref.ordered[meth.meth_hi==groups[2]&meth.cov_mean_g1>=cov_thres& meth.cov_mean_g2>=cov_thres&meth.pval_adj<=0.05]
top2=top2_all[1:min(top_threshold,nrow(top2_all))]$dedRef_ID

bottom_all=diffTab_ref.ordered[meth.cov_mean_g1>=cov_thres& meth.cov_mean_g2>=cov_thres&!is.na(meth.rank_max)]
bottom1=bottom_all[(nrow(bottom_all)-min(top_threshold,nrow(top1_all))):nrow(bottom_all)]$dedRef_ID
bottom2=bottom_all[(nrow(bottom_all)-min(top_threshold,nrow(top2_all))):nrow(bottom_all)]$dedRef_ID

system("mkdir -p motifAnalysis")

if (length(na.omit(top1))>0){
  top_fa1=fasta[top1]
  writeXStringSet(top_fa1,paste0("motifAnalysis/","ded_top",top_threshold,"_cov",cov_thres,"_",groups[1],".fa"))
  bottom_fa1=fasta[bottom1]
  writeXStringSet(bottom_fa1,paste0("motifAnalysis/","ded_bottom",top_threshold,"_cov",cov_thres,"_",groups[1],".fa"))}else{message("No significantly hypermethylated fragments for ",comp[1])}
if (length(na.omit(top2))>0){
  top_fa2=fasta[top2]
  writeXStringSet(top_fa2,paste0("motifAnalysis/","ded_top",top_threshold,"_cov",cov_thres,"_",groups[2],".fa"))
  bottom_fa2=fasta[bottom2]
  writeXStringSet(bottom_fa2,paste0("motifAnalysis/","ded_bottom",top_threshold,"_cov",cov_thres,"_",groups[2],".fa"))}else{message("No significantly hypermethylated fragments for ",comp[2])}



stats=rbind(stats,data.frame(assay="diffMeth",condition=c(paste0(groups[1],"_methHi"),paste0(groups[2],"_methHi")),percent=c(nrow(top1_all)/nrow(diffTab_ref.ordered),nrow(top2_all)/nrow(diffTab_ref.ordered)),number=c(nrow(top1_all),nrow(top2_all)),correlation=NA))

writeXStringSet(fasta,paste0("motifAnalysis/","allDeducedGenomeFragments.fa"))
write.table(stats,paste0(out_dir,"/",species,"_diffMeth_stats.tsv"),quote=FALSE,sep="\t",row.names=FALSE)
}

system(paste0("echo '' >",wd,"/",motif,"_diffMeth.done"))

