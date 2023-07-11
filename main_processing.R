library(DoubletFinder)
library(Seurat) #dont use v5, DoubletFinder does not support Assay5

RDoublet <- function(tmp){
    sweep.res.list <- paramSweep_v3(tmp, PCs = 1:20, sct = FALSE)
    sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
    bcmvn <- find.pK(sweep.stats)
    pKopt <- as.numeric(as.character(bcmvn$pK[bcmvn$BCmetric == max(bcmvn$BCmetric)]))
    pKopt <- pKopt[order(pKopt, decreasing = TRUE) ]
    pKopt <- pKopt[1]
    homotypic.prop <- modelHomotypic(tmp$seurat_clusters)
    nExp_poi <- round(0.1*length(colnames(tmp)))  ## Assuming 10% doublet formation rate
    nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
    tmp <- doubletFinder_v3(tmp, PCs = 1:20, pN = 0.25, pK = pKopt, nExp = nExp_poi, reuse.pANN = FALSE)
    tmp <- doubletFinder_v3(tmp, PCs = 1:20, pN = 0.25, pK = pKopt, nExp = nExp_poi.adj, reuse.pANN = paste("pANN_0.25",pKopt,nExp_poi, sep="_"))
    return (tmp)
}

ProcessSeu <- function(Seurat){
    Seurat <- NormalizeData(Seurat)
    Seurat <- FindVariableFeatures(Seurat, selection.method = "vst", nfeatures = 2000)
    Seurat <- ScaleData(Seurat)
    Seurat <- RunPCA(Seurat)
    Seurat <- FindNeighbors(Seurat, dims = 1:20)
    Seurat <- FindClusters(Seurat, resolution = 0.5)
    Seurat <- RunUMAP(Seurat, dims = 1:20)
    Seurat <- RunTSNE(Seurat,  dims.use = 1:20 )
    DimPlot(object = Seurat, reduction = "umap")
    return (Seurat)
}

s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes
library(Seurat)

processSeuratObjects <- function(startIndex, endIndex) {
    seuratObjects <- list()
    
    for (i in startIndex:endIndex) {
        # Define the filename and project name
        filename <- paste0("E://scRNAseq_PBMC/cellranger_run/D23-", i, "/outs/filtered_feature_bc_matrix")
        projectName <- as.character(i)
        
        # Read the data
        data <- Read10X(data.dir = filename)
        
        # Create Seurat object
        seuratObj <- CreateSeuratObject(counts = data, project = projectName, min.cells = 3, min.features = 200)
        
        # Calculate percentage of ribosomal genes
        seuratObj[["percent.rb"]] <- PercentageFeatureSet(seuratObj, pattern = "^RPS|^RPL|^MRPS|^MRPL", assay = 'RNA')
        
        # Perform cell cycle scoring
        seuratObj <- CellCycleScoring(seuratObj, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)
        
        # Calculate percentage of mitochondrial genes
        seuratObj[["percent.mt"]] <- PercentageFeatureSet(seuratObj, pattern = "^MT-")
        
        # Generate violin plots
        VlnPlot(seuratObj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", 'percent.rb'), ncol = 4)
        
        # Subset cells based on specified criteria
        seuratObj <- subset(seuratObj, subset = nCount_RNA > 300 & nCount_RNA < 10000 & nFeature_RNA > 600 & nFeature_RNA < 4000 & percent.mt < 15 & percent.rb < 40)
        
        # Scale the data
        seuratObj <- ScaleData(seuratObj, verbose = TRUE, vars.to.regress = c('nCount_RNA', 'nFeature_RNA', 'percent.mt', "percent.rb","S.Score","G2M.Score"))
        seuratObj <- ProcessSeu(seuratObj)
        # Perform doublet detection
        seuratObj <- RDoublet(seuratObj)
        seuratObj <- subset(seuratObj , cells = colnames(seuratObj )[which(seuratObj [[]][13] == 'Singlet')])
        seuratObj <- subset(seuratObj , cells = colnames(seuratObj )[which(seuratObj [[]][12] == 'Singlet')])
        seuratObj <- ProcessSeu(seuratObj)
        # Store the Seurat object in the list
        seuratObjects[[paste0("s", i)]] <- seuratObj
    }
    
    return(seuratObjects)
}

# Usage example:
startIndex <- 3805
endIndex <- 3812
seuratObjs <- processSeuratObjects(startIndex, endIndex)

# Access the individual Seurat objects
s3805 <- seuratObjs[["s3805"]]
s3806 <- seuratObjs[["s3806"]]
s3807 <- seuratObjs[["s3807"]]
s3808 <- seuratObjs[["s3808"]]
s3809 <- seuratObjs[["s3809"]]
s3810 <- seuratObjs[["s3810"]]
s3811 <- seuratObjs[["s3811"]]
s3812 <- seuratObjs[["s3812"]]

s3805$condition <- 'Control'
s3806$condition <- 'Control'
s3807$condition <- 'Control'
s3808$condition <- 'Glaucoma'
s3809$condition <- 'Glaucoma'
s3810$condition <- 'Glaucoma'
s3811$condition <- 'Glaucoma'
s3812$condition <- 'Glaucoma'

ProcessInt <- function(data.integrated){
    data.integrated <- ScaleData(data.integrated, verbose = FALSE)
    data.integrated <- RunPCA(data.integrated, npcs = 30, verbose = FALSE)
    data.integrated <- FindNeighbors(data.integrated, dims = 1:20)
    data.integrated <- FindClusters(data.integrated, resolution = 0.5)
    data.integrated <- RunUMAP(data.integrated, reduction = "pca", dims = 1:20)
    data.integrated <- RunTSNE(data.integrated,  dims.use = 1:20 )
}

integration_list <- list(s3805,s3806,s3807,s3808,s3809,s3810,s3811,s3812)

features <- SelectIntegrationFeatures(object.list = integration_list)
data.anchors <- FindIntegrationAnchors(object.list = integration_list, anchor.features = features)

data.combined <- IntegrateData(anchorset = data.anchors)

shuhong.combined <- ProcessInt(data.combined)

shuhong.combined.markers <- FindAllMarkers(shuhong.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

shuhong.combined.markers %>%
    group_by(cluster) %>%
    slice_max(n = 2, order_by = avg_log2FC)

shuhong.combined.markers %>%
    group_by(cluster) %>%
    top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(shuhong.combined, features = top10$gene) + NoLegend()

SaveH5Seurat(shuhong.combined, 'E://scRNAseq_PBMC/combined.h5Seurat')
