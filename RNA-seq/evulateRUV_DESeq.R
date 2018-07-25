#!/usr/bin/Rscript

library(RUVSeq)
library(DESeq2)


evulateERCC = function(erccf,pre){
    #ERCC 
    ercc = read.table(erccf,sep="\t",header=T,row.names=1,check.names=F)
    #conditions 
    condition = gsub("[0-9]","",colnames(ercc))
    design = data.frame(row.names=colnames(ercc),condition=condition) 
    ###DESeq normalization using ercc estimated size factors
    #prepare data
    ercccds = DESeqDataSetFromMatrix(countData=ercc,colData=design,design=~ 1)
    #normalization
    ercccds = estimateSizeFactors(ercccds)
    sf = sizeFactors(ercccds)
    ercc_normalized = as.data.frame(counts(ercccds,normalized=TRUE))
    write.table(ercc_normalized,paste(pre,"deseq.txt",sep="_"),sep="\t",fileEncoding="utf-8",quote=F,eol="\n")
    write.table(sf,paste(pre,"sizefactor",sep="_"),sep="\t",fileEncoding="utf-8",quote=F,col="\n")
    ##RUV
    #ercc_set = newSeqExpressionSet(as.matrix(ercc),phenoData=design)
    ercc_set = newSeqExpressionSet(as.matrix(ercc))
    ercc_ruv = normCounts(RUVg(ercc_set,rownames(ercc),k=1))
    write.table(ercc_ruv,paste(pre,"ruvseq.txt",sep="_"),sep="\t",fileEncoding="utf-8",quote=F,eol="\n")
}

erccf = "../../../5.ParsedCountsMatrix/batches/TL_ERCC_counts_batch1.txt"
evulateERCC(erccf,"batch1")
erccf = "../../../5.ParsedCountsMatrix/batches/TL_ERCC_counts_batch2.txt"
evulateERCC(erccf,"batch2")
