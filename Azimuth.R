remotes::install_github('satijalab/azimuth', ref = 'master')
library(Azimuth)
shuhong.combined <- RunAzimuth(shuhong.combined, reference = "pbmcref")
DimPlot(shuhong.combined, group.by = "predicted.celltype.l2", label = TRUE, label.size = 3, label.box = TRUE, repel = TRUE) + NoLegend()
DoHeatmap(shuhong.combined, features = c('FCGR3A','NCAM1','CD19','CD14','CD4','CD8A','CD3E','CD3G','FCGR1A','CD45','TRG','TRD','TRA','TRB','PTPRC','B3GAT1','CD28','HLA-DRB1',
                                         'CD27','CD45','IGHD','CR2','CD38','CD24','IGHM','THBD','LIN28A','IL3RA','ITGAX','CLEC4C','CD16','IFNA6'),
          group.by = 'seurat_clusters') + NoLegend()

genes <- c("CD4", "IL2RA", "FOXP3", "CTLA4", "TNFRSF18", "PTPRC",
           "ITGA2B", "GP1BA", "ITGB3",
           "CD19", "MS4A1", "CD38", "CD27", "SDC1",
           "IL3RA", "CLEC4C", "CLEC4D", "ITGAX", "HLA-DRA",
           "CD3D", "NCAM1", "FCGR3A", "KLRD1", "NCR1", "MKI67",
           "TRAV1-2", "CLEC2B", "TRGV3",
           "IL7R", "ITGAE", "KLRB1",
           "CD34", "LIN28A",
           "TRDC",  "CD8A",
           "LY6G6D", "TFR1", "GYPA",
           "CD1C", "CD14",  "NRP1",
           "CCR7", 
            "IGHD")


DotPlot(shuhong.combined, features = genes,
          group.by = 'predicted.celltype.l2') + NoLegend()
#PTPRCRA = PTPRC, CD138 = SDC1, CD71 = TFR1
