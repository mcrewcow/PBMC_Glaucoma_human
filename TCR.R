devtools::install_github("ncborcherding/scRepertoire")
suppressMessages(library(scRepertoire))
suppressMessages(library(Seurat))

library(scRepertoire)

S1 <- read.csv("E://scRNAseq_PBMC/cellranger_run/D23-3695/outs/filtered_contig_annotations.csv")
S1$condition <- 'Control'
S2 <- read.csv("E://scRNAseq_PBMC/cellranger_run/D23-3696/outs/filtered_contig_annotations.csv")
S2$condition <- 'Control'
S3 <- read.csv("E://scRNAseq_PBMC/cellranger_run/D23-3697/outs/filtered_contig_annotations.csv")
S3$condition <- 'Control'

S4 <- read.csv("E://scRNAseq_PBMC/cellranger_run/D23-3698/outs/filtered_contig_annotations.csv")
S4$condition <- 'Glaucoma'
S5 <- read.csv("E://scRNAseq_PBMC/cellranger_run/D23-3699/outs/filtered_contig_annotations.csv")
S5$condition <- 'Glaucoma'
S6 <- read.csv("E://scRNAseq_PBMC/cellranger_run/D23-3700/outs/filtered_contig_annotations.csv")
S6$condition <- 'Glaucoma'
S7 <- read.csv("E://scRNAseq_PBMC/cellranger_run/D23-3701/outs/filtered_contig_annotations.csv")
S7$condition <- 'Glaucoma'
S8 <- read.csv("E://scRNAseq_PBMC/cellranger_run/D23-3702/outs/filtered_contig_annotations.csv")
S8$condition <- 'Glaucoma'

contig_list <- list(S1, S2, S3, S4,S5,S6,S7,S8)
combined <- combineTCR(contig_list, 
                       samples = c("Control",'Control',"Control",'Glaucoma','Glaucoma','Glaucoma','Glaucoma','Glaucoma'), 
                       ID = c("1",'2',"3",'4','5','6','7','8')) #or repeat

quantContig(combined, cloneCall="gene+nt", chain = "TRB", group.by = 'ID', scale = F)
abundanceContig(combined, cloneCall = "gene", scale = F)
lengthContig(combined, cloneCall="aa", chain = "both", group.by = 'condition') #or clone nt

lengthContig(combined, cloneCall="nt", chain = "TRA")

compareClonotypes(combined, 
                  numbers = 10, 
                  samples = c("Control_Control", "Glaucoma_Glaucoma"), 
                  cloneCall="aa", 
                  graph = "alluvial")

vizGenes(combined, gene = "V", 
         chain = "TRA", 
         plot = "bar", 
         order = "variance", 
         scale = TRUE)

vizGenes(combined[c(1:3)], 
         gene = "V", 
         chain = "TRB", 
         y.axis = "J", 
         plot = "heatmap", 
         scale = TRUE, 
         order = "gene")

vizGenes(combined[c(4:8)], 
         gene = "V", 
         chain = "TRB", 
         y.axis = "J", 
         plot = "heatmap", 
         scale = TRUE, 
         order = "gene")

clonalHomeostasis(combined, cloneCall = "gene", 
                  cloneTypes = c(Rare = 1e-04, 
                                 Small = 0.001, 
                                 Medium = 0.01, 
                                 Large = 0.1, 
                                 Hyperexpanded = 1))
clonalHomeostasis(combined, cloneCall = "aa")
clonalProportion(combined, cloneCall = "gene",
                 split = c(10, 100, 1000, 10000, 30000, 1e+05)) 
clonalProportion(combined, cloneCall = "nt") 

scatterClonotype(combined, 
                 cloneCall ="gene", 
                 x.axis = "PY_P", 
                 y.axis = "PY_T",
                 dot.size = "total",
                 graph = "proportion")
