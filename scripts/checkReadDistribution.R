library(data.table)
library(ggplot2)
theme_set(theme_bw())
args=commandArgs(trailingOnly = TRUE)

wd=args[1]
sample=args[2]
#wd="/scratch/lab_bock/jklughammer/projects/Zebrafinch/toSelf_filtered_0.08mm_final_concat/BSF_0185_H3GLTBBXX_4__G1AV207_2A/biseqMethcalling/"
setwd(wd)

cpgReads=fread(paste0("RRBS_cpgReads_",sample,".bed"),drop=c(4,7,8,9,11,12,13),sep="\t")
setnames(cpgReads,names(cpgReads),c("chr","start","end","meth","orientation","Nmotifs"))
cpgReads[,coverage:=.N,by=c("chr","start")]
cpgReads[,coverage_topped_1000:=ifelse(coverage<1000,coverage,1000),]

cpgReads_bins=cpgReads[,list(readsCoveredByMin1read=length(chr[coverage>0]),readsCoveredByMin5reads=length(chr[coverage>4]),readsCoveredByMin10reads=length(chr[coverage>9]),readsCoveredByMin50reads=length(chr[coverage>49]),readsCoveredByMin100reads=length(chr[coverage>99]),readsCoveredByMin1000reads=length(chr[coverage>999]),readsCoveredByMin10000reads=length(chr[coverage>9990])),]
cpgReads_coverage=cpgReads[,list(Nreads=.N),by=coverage][order(coverage)]
cpgReads_coverage[,Nregions:=Nreads/coverage,]

write.table(cpgReads_bins,"cpg_readDistribution.tsv",sep="\t",quote=FALSE,row.names=FALSE)

pdf("readCoverage.pdf",height=5,width=5)
ggplot(cpgReads,aes(x=coverage))+geom_density()
ggplot(cpgReads,aes(x=coverage_topped_1000))+geom_density()
ggplot(cpgReads,aes(x=coverage_topped_1000))+geom_histogram()
dev.off()
