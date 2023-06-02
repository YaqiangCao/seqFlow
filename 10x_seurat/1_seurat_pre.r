library(Seurat)
library(DoubletFinder)
suppressMessages(require(DoubletFinder))


pre10x = function(data, pn,minGenes=500, minCells=10,topFeature=3000,dr=0.1){
    #pre-process replicates of 10x RNA-seq results based on seurat; reference to following
    #1. https://satijalab.org/seurat/articles/pbmc3k_tutorial.html 
    #2. https://nbisweden.github.io/workshop-scRNAseq/labs/compiled/seurat/seurat_01_qc.html#Predict_doublets

    # @param data, a seurat object, may merged from multiple experiments
    # @param pn, project name for each seurat object 
    # @param dr, doublelets expected ratio

    start_time = Sys.time()

    #step 1
    #read in cell ranger processed data and create seurate objects
    #ds = c()
    #cids = c()
    #for (i in seq_along(bcms)){
    #    n = paste(pn,'_rep',i,sep="")
    #    d = CreateSeuratObject(counts=Read10X(bcms[i]), project=n)
    #    d$type = pn
    #    cids = append(cids,  n)
    #    ds = append(ds,d)
    #}

    #write the raw data for following check
    write.table(data@assays$RNA@counts,paste(pn, "raw_counts.txt",sep="_"),col.names=TRUE, row.names=TRUE, quote=FALSE,sep="\t")

    #step 3 quality control of mito and ribosome genes, mice 
    print(paste(pn,"raw dim"))
    print(dim(data))
    data = PercentageFeatureSet(data, "^mt-", col.name = "percent_mito")
    data = PercentageFeatureSet(data, "^Rp[sl]", col.name = "percent_ribo")
    pdf(paste(pn, "1_raw_qc.pdf",sep="_"))
    feats = c("Xist","nFeature_RNA", "nCount_RNA", "percent_mito", "percent_ribo")
    VlnPlot(data, group.by = "orig.ident", features = feats, pt.size = 0, ncol = 5) + NoLegend()
    dev.off()
    
    #start filtering
    selected_c = WhichCells(data, expression = nFeature_RNA > minGenes)
    selected_f = rownames(data)[Matrix::rowSums(data) > minCells]
    data.filt = subset(data, features = selected_f, cells = selected_c)
    print(paste(pn,"filter minimal genes and cells dim"))
    print(dim(data.filt))
    selected_mito = WhichCells(data.filt, expression = percent_mito < 5) #already converted to percentage
    selected_ribo = WhichCells(data.filt, expression = percent_ribo > 5)
    data.filt = subset(data.filt, cells = selected_mito)
    data.filt = subset(data.filt, cells = selected_ribo)
    print(paste(pn,"filter mito and ribosome genes dim", dim(data.filt)))
    #filter most expressed genes are Malat1, Tmsb4x, Actb
    data.filt = data.filt[!grepl("Malat1",rownames(data.filt)),]
    data.filt = data.filt[!grepl("Tmsb4x",rownames(data.filt)),]
    data.filt = data.filt[!grepl("Actb",rownames(data.filt)),]
    
    #step 3 quality control of mito and ribosome genes, mice 
    #normalize data 
    data.filt = NormalizeData( data.filt )
    #feature selection
    data.filt = FindVariableFeatures(data.filt, selection.method = "vst", nfeatures = topFeature)
    dim(data.filt)
    top20 = head(VariableFeatures(data.filt), 20)

    #show top variable features 
    pdf(paste(pn,"2_top20_var.pdf",sep="_"))
    plot1 = VariableFeaturePlot(data.filt)
    plot2 = LabelPoints(plot = plot1, points = top20, repel = TRUE)
    dev.off()

    #scaling 
    data.filt = ScaleData(data.filt, features=rownames(data.filt))
    #PCA
    data.filt = RunPCA(data.filt, verbose = F, npcs = 20, features = VariableFeatures(object = data.filt))
    #shwo PCA
    pdf(paste(pn,"3_pca.pdf",sep="_"))
    VizDimLoadings(data.filt, dims = 1:2, reduction = "pca")
    DimPlot(data.filt, reduction = "pca")
    DimHeatmap(data.filt, dims = 1:20, cells = 1000, balanced = TRUE)
    dev.off()
    data.filt = JackStraw(data.filt, num.replicate = 100)
    data.filt = ScoreJackStraw(data.filt, dims = 1:20)
    pdf(paste(pn,"4_selPc.pdf",sep="_"))
    JackStrawPlot(data.filt, dims = 1:20)
    ElbowPlot(data.filt)
    dev.off()

    pdf(paste(pn,"5_pca.pdf",sep="_"))
    data.filt = RunUMAP(data.filt, dims = 1:20, verbose = F)
    DimPlot(data.filt, reduction="umap", group.by="orig.ident") + NoAxes() 
    dev.off()

    #remove potential doublets 
    nExp = round(ncol(data.filt) * dr)  # expect 8% doublets
    data.filt = doubletFinder_v3(data.filt, pN = 0.25, pK = 0.09, nExp = nExp, PCs = 1:20)
    DF.name = colnames(data.filt@meta.data)[grepl("DF.classification", colnames(data.filt@meta.data))]
    pdf(paste(pn,"6_doublets.pdf",sep="_"))
    cowplot::plot_grid(ncol = 2, DimPlot(data.filt, group.by = "orig.ident") + NoAxes(), DimPlot(data.filt, group.by = DF.name) + NoAxes())
    VlnPlot(data.filt, features = "nFeature_RNA", group.by = DF.name, pt.size = 0)
    dev.off()
    data.filt = data.filt[, data.filt@meta.data[, DF.name] == "Singlet"]
    
    #save 
    write.table(data.filt@assays$RNA@data,paste(pn, "normalizedFiltered_counts.txt",sep="_") ,col.names=TRUE, row.names=TRUE, quote=FALSE,sep="\t")
    umap_coordinates = Embeddings(data.filt, reduction = "umap")
    write.table(umap_coordinates,paste(pn,"normalizedFiltered_umap.txt",sep="_"),col.names=TRUE, row.names=TRUE, quote=FALSE,sep="\t")
    saveRDS(data.filt, file = paste(pn,".rds",sep=""))
    
    end_time = Sys.time()
    print(paste(pn, "Running time:", end_time-start_time))
}


#only for test
#R1.data = Read10X("../1.cellRanger/R1/filtered_feature_bc_matrix")
#R6.data = Read10X("../1.cellRanger/R6/filtered_feature_bc_matrix")
#R8.data = Read10X("../1.cellRanger/R8/filtered_feature_bc_matrix")
#R1 = CreateSeuratObject(counts=R1.data, project="TregIgGRep1")
#R6 = CreateSeuratObject(counts=R6.data, project="TregCtla4Rep2")
#R8 = CreateSeuratObject(counts=R8.data, project="TregCtla4Rep4")
#R1$type = "IgG"
#R6$type = "Ctla4"
#R8$type = "Ctla4"
#data = merge( R1, y=c(R6,R8), add.cell.ids=c("IgG_rep1","Ctla4_rep2","Ctla4_rep4") )
#pre10x(data,"init")


#read in
R1.data = Read10X("../1.cellRanger/R1/filtered_feature_bc_matrix")
R2.data = Read10X("../1.cellRanger/R2/filtered_feature_bc_matrix")
R3.data = Read10X("../1.cellRanger/R3/filtered_feature_bc_matrix")
R4.data = Read10X("../1.cellRanger/R4/filtered_feature_bc_matrix")
R5.data = Read10X("../1.cellRanger/R5/filtered_feature_bc_matrix")
R6.data = Read10X("../1.cellRanger/R6/filtered_feature_bc_matrix")
R7.data = Read10X("../1.cellRanger/R7/filtered_feature_bc_matrix")
R8.data = Read10X("../1.cellRanger/R8/filtered_feature_bc_matrix")

#seurat objects
R1 = CreateSeuratObject(counts=R1.data, project="IgGRep1")
R2 = CreateSeuratObject(counts=R2.data, project="IgGRep2")
R3 = CreateSeuratObject(counts=R3.data, project="IgGRep3")
R4 = CreateSeuratObject(counts=R4.data, project="IgGRep4")
R5 = CreateSeuratObject(counts=R5.data, project="Ctla4Rep1")
R6 = CreateSeuratObject(counts=R6.data, project="Ctla4Rep2")
R7 = CreateSeuratObject(counts=R7.data, project="Ctla4Rep3")
R8 = CreateSeuratObject(counts=R8.data, project="Ctla4Rep4")

#add metadata
R1$type = "IgG"
R2$type = "IgG"
R3$type = "IgG"
R4$type = "IgG"
R5$type = "Ctla4"
R6$type = "Ctla4"
R7$type = "Ctla4"
R8$type = "Ctla4"

#combine all data
IgG = merge( R1, y=c(R2,R3,R4), add.cell.ids=c("IgG_rep1","IgG_rep2","IgG_rep3","IgG_rep4") )
Ctla4 = merge( R5, y=c(R6,R7,R8), add.cell.ids=c("Ctla4_rep1","Ctla4_rep2","Ctla4_rep3","Ctla4_rep4") )
pre10x(IgG,"IgG" )
pre10x(Ctla4,"Ctla4")
