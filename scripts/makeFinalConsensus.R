#!/usr/bin/env Rscript

#Author: Johanna Klughammer
#Date: 26.07.2015


options(stringsAsFactors=FALSE)
library(data.table)
library(Biostrings)


args=commandArgs(trailingOnly = TRUE)

in_file= args[1]
out_dir=args[2]
wd=args[3]
working_dir=args[4]
consensus_dist=args[5]
cLimit=args[6]

setwd(wd)


#functions
adjustPos=function(mapped_seq,target_pos){
  Ns=lapply(target_pos,function(x){paste0(rep("N",x-1),collapse="")})
  newSeq=paste0(Ns,mapped_seq)
  return(newSeq)
}

makeConsensus=function(forConsensus){
  cm=as.data.frame(t(consensusMatrix(DNAStringSet(forConsensus))))
  raw_cons=colnames(cm)[1:4][apply(cm[,1:4], 1, which.max)]
  Ts=which(raw_cons == "T")
  raw_cons[Ts[cm$C[Ts]/(cm$C[Ts]+cm$T[Ts])>=cLimit]]="C"
  cons=paste0(raw_cons,collapse="")
  return(cons)
}

print(gc())

system.time({
  print("prepare file")
  all=fread(in_file)
  setnames(all,c(names(all)),c("mapped_name","mapped_seq","target_name","target_pos","CIGAR","mismatches","mapped_length","aln_type","flag"))
  print(gc())
  })



#reduce file (only keep first occurrence of a read)
system.time({
  print("reduce file")
#sort by 1. target_n, 2. mapped_n, 3. target_name
  print("sort")
  all=all[,mapped_n:=-.N,by=mapped_name]
  all=all[,target_n:=-.N,by=target_name]
  setkeyv(all,c("target_n","target_name","mapped_n"))
  print(gc())
#Build helper_u table do be able to find the first occurrence of a read in target as well as mapped   
  print("build helper_u")
  all[,rowID:=1:nrow(all),]
  helper_u=rbindlist(list(cbind(all[,list(target_name,rowID)],rep("target",nrow(all))),cbind(all[,list(mapped_name,rowID)],rep("mapped",nrow(all)))))
  setnames(helper_u,c("target_name","V2"),c("name","type"))
  helper_u=helper_u[order(rowID)]
  helper_u[,dupl:=duplicated(name),]
  print(gc())
#Prepare the search table
  print("prepare search table")
  in_dt=all
  print(gc())
  rm(all)
  print(gc())
  in_dt[,targetKeep:=FALSE,]
  in_dt[,mappedKeep:=FALSE,] 
  
#Build helper listing the start and end position for each target read (rows in in_dt occupied by a read)
  print("build helper")
  helper=in_dt[,.N,by=target_name]
  helper[,end:=cumsum(N),]
  helper$start=c(0,helper$end[1:(nrow(helper)-1)])+1
  helper2=copy(helper)
  setkey(helper2,target_name)
  print(gc())
  
#Find indices of targets to keep and mark them in in_dt
  print("tag to keep")
  temp_h=helper2[helper_u[type=="target"&dupl==FALSE,,]$name]
  targetKeep_ind=temp_h[,seq(start,end),by=target_name]$V1
  in_dt[targetKeep_ind,targetKeep:=TRUE,]

#Find indices of mapped to keep and mark them in in_dt
  mappedKeep_ind=helper_u[type=="mapped"&dupl==FALSE,,]$rowID
  in_dt[mappedKeep_ind,mappedKeep:=TRUE,]
  in_dt[mapped_name==target_name&targetKeep==TRUE,mappedKeep:=TRUE,]
  print(gc()) 
#output reduced table
  print("make out table")
  out_dt=in_dt[targetKeep&mappedKeep,,]  #try if this is more memory efficient and still works
  out_dt[,target_n:=.N,by=target_name]
  print(gc())
})  
  
write.table(out_dt,paste0(out_dir,in_file,"_red"),sep="\t",quote=FALSE,row.names=FALSE)
rm(in_dt)
rm(helper_u)
rm(helper)
rm(helper2)
print(gc())

system.time({
  print("adjust for consensus")
  out_dt[bitwAnd(flag,16)==16,mapped_seq:=unlist(lapply(mapped_seq,function(x){as.character(reverseComplement(DNAString(x)))})),]
  out_dt[target_n==1,forConsensus:=mapped_seq,]
  out_dt[target_n>1,forConsensus:=adjustPos(mapped_seq,target_pos),]
  print(gc())
})

write.table(out_dt,paste0(out_dir,in_file,"_forCons"),sep="\t",quote=FALSE,row.names=FALSE)


system.time({
  print("calculate consensus")
  out_dt[target_n==1,consensus:=mapped_seq,]
  out_dt[target_n>1,consensus:=makeConsensus(forConsensus),by=target_name] #C maintaining consensus
  out_dt[,edit_dist:=adist(consensus,forConsensus),by=rowID]
  out_dt[,diffNts:=edit_dist-abs(nchar(consensus)-nchar(mapped_seq)),]
  out_dt[,edit_dist_conv:=adist(gsub("C","T",consensus),gsub("C","T",forConsensus)),by=rowID]
  out_dt[,diffNts_conv:=edit_dist_conv-abs(nchar(consensus)-nchar(mapped_seq)),]
  print(gc())
  })

write.table(out_dt,paste0(out_dir,in_file,"_cons"),sep="\t",quote=FALSE,row.names=FALSE)

system.time({
  print("relax consensus")
  mismatch=consensus_dist
  count=0
  prevScore=out_dt[,.N,by=diffNts_conv][diffNts_conv>0,sum(N),]+1
  while (prevScore>out_dt[,.N,by=diffNts_conv][diffNts_conv>0,sum(N),]){
    count=count+1
    prevScore=out_dt[,.N,by=diffNts_conv][diffNts_conv>0,sum(N),]
    print(paste0(count,": ",prevScore," sequences > 0 mismatches"))
    out_dt[diffNts_conv/nchar(mapped_seq)>mismatch,consensus:=makeConsensus(forConsensus),by=consensus] #C maintaining consensus
    out_dt[diffNts_conv/nchar(mapped_seq)>mismatch,edit_dist:=adist(consensus,forConsensus),by=rowID]
    out_dt[diffNts_conv/nchar(mapped_seq)>mismatch,diffNts:=edit_dist-abs(nchar(consensus)-nchar(mapped_seq)),]
    out_dt[diffNts_conv/nchar(mapped_seq)>mismatch,edit_dist_conv:=adist(gsub("C","T",consensus),gsub("C","T",forConsensus)),by=rowID]
    out_dt[diffNts_conv/nchar(mapped_seq)>mismatch,diffNts_conv:=edit_dist_conv-abs(nchar(consensus)-nchar(mapped_seq)),]
  }

if (nrow(out_dt[diffNts_conv/nchar(mapped_seq)>mismatch]) > 2)  {
  count=0
  prevScore=out_dt[,.N,by=diffNts_conv][diffNts_conv>0,sum(N),]+1
  while (prevScore>out_dt[,.N,by=diffNts_conv][diffNts_conv>0,sum(N),]){
    count=count+1
    prevScore=out_dt[,.N,by=diffNts_conv][diffNts_conv>0,sum(N),]
    print(paste0(count,": ",prevScore," sequences > 0 mismatches"))
    #subdivide consensus groups into highest similarity groups
    out_dt[diffNts_conv/nchar(mapped_seq)>mismatch,simil:=rowID[apply(adist(gsub("C","T",forConsensus),gsub("C","T",forConsensus)),1,function(x){x[which.min(x)]=Inf;which.min(x)})],by="consensus"]  
    out_dt[rowID%in%simil,simil:=rowID,]
    
    out_dt[diffNts_conv/nchar(mapped_seq)>mismatch,consensus:=makeConsensus(forConsensus),by=c("consensus","simil")] #C maintaining consensus
    out_dt[diffNts_conv/nchar(mapped_seq)>mismatch,edit_dist:=adist(consensus,forConsensus),by=rowID]
    out_dt[diffNts_conv/nchar(mapped_seq)>mismatch,diffNts:=edit_dist-abs(nchar(consensus)-nchar(mapped_seq)),]
    out_dt[diffNts_conv/nchar(mapped_seq)>mismatch,edit_dist_conv:=adist(gsub("C","T",consensus),gsub("C","T",forConsensus)),by=rowID]
    out_dt[diffNts_conv/nchar(mapped_seq)>mismatch,diffNts_conv:=edit_dist_conv-abs(nchar(consensus)-nchar(mapped_seq)),]
  } 
}  
  rest=nrow(out_dt[diffNts_conv/nchar(mapped_seq)>mismatch])
  print(paste0(rest," sequences with > ",mismatch," % mismatches to the consensus"))
  if (out_dt[,max(diffNts_conv/nchar(mapped_seq))]>mismatch){
    print("Consensus did not converge. Replacing consensus of non fitting sequenes with mapped sequence.")
    out_dt[diffNts_conv/nchar(mapped_seq)>mismatch,consensus:=mapped_seq,] 
    out_dt[diffNts_conv/nchar(mapped_seq)>mismatch,edit_dist:=adist(consensus,forConsensus),by=rowID]
    out_dt[diffNts_conv/nchar(mapped_seq)>mismatch,diffNts:=edit_dist-abs(nchar(forConsensus)-nchar(mapped_seq)),]
    out_dt[diffNts_conv/nchar(mapped_seq)>mismatch,edit_dist_conv:=adist(gsub("C","T",consensus),gsub("C","T",forConsensus)),by=rowID]
    out_dt[diffNts_conv/nchar(mapped_seq)>mismatch,diffNts_conv:=edit_dist_conv-abs(nchar(forConsensus)-nchar(mapped_seq)),]
  }
  print(gc())    
})
write.table(out_dt,paste0(out_dir,in_file,"_cons_relaxed"),sep="\t",quote=FALSE,row.names=FALSE)


system.time({
  print("reduce output")
  final=out_dt[,list(.N,round(mean(diffNts_conv),1),max(diffNts_conv),min(diffNts_conv)),by=list(target_name,consensus)]
  setnames(final,c("N","V2","V3","V4"),c("consN","mean_diffNts","max_diffNts","min_diffNts"))
  final[,target_nameN:=.N,by="target_name"]
  final[target_nameN>1,target_name:=outer(target_name,c(1:.N),paste,sep="-"),by="target_name"]
  final[,target_nameN:=NULL,]
  print(gc())
})


write.table(final,paste0(out_dir,in_file,"_final_all"),sep="\t",quote=FALSE,col.names=FALSE,row.names=FALSE)

#remove all deducedReferences, that don't contain a CpG
final=final[grepl("CG|GC",consensus)]
write.table(final,paste0(out_dir,in_file,"_final"),sep="\t",quote=FALSE,col.names=FALSE,row.names=FALSE)


system(paste0("echo '' > ",working_dir,"/reduce_redundency.done"))






