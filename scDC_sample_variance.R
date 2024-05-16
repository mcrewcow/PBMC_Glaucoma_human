library(scDC)
library(broom.mixed)
exprsMat <- GetAssayData(object = PBMC, assay = "RNA", slot = "data")
cellTypes <- PBMC$predicted.celltype.l2
subject <- PBMC$sample
cond <- PBMC$condition 
dim(exprsMat)
table(subject, cellTypes)
table(cond, cellTypes)
res_scDC_noClust <- scDC_noClustering(cellTypes, subject, calCI = TRUE,
calCI_method = c("percentile", "BCa", "multinom"),
nboot = 50)
barplotCI(res_scDC_noClust, c("Control","Control","Control",
                              "Control","Control","Control",
                              "Control","Control","Control",
                              "Control","Control","Control",
                              "Control","Control","Control",
                              "Control","Control","Control",
                              "Control","Control","Control",
                              "Control","Control","Control",
                              "Control","Control","Control",
                              "Glaucoma","Glaucoma","Glaucoma","Glaucoma","Glaucoma",
                              "Glaucoma","Glaucoma","Glaucoma","Glaucoma","Glaucoma",
                              "Glaucoma","Glaucoma","Glaucoma","Glaucoma","Glaucoma",
                              "Glaucoma","Glaucoma","Glaucoma","Glaucoma","Glaucoma",
                              "Glaucoma","Glaucoma","Glaucoma","Glaucoma","Glaucoma",
                              "Glaucoma","Glaucoma","Glaucoma","Glaucoma","Glaucoma",
                              "Glaucoma","Glaucoma","Glaucoma","Glaucoma","Glaucoma",
                              "Glaucoma","Glaucoma","Glaucoma","Glaucoma","Glaucoma",
                              "Glaucoma","Glaucoma","Glaucoma","Glaucoma","Glaucoma"))
res_GLM <- fitGLM(res_scDC_noClust, c("Control","Control","Control",
                                      "Control","Control","Control",
                                      "Control","Control","Control",
                                      "Control","Control","Control",
                                      "Control","Control","Control",
                                      "Control","Control","Control",
                                      "Control","Control","Control",
                                      "Control","Control","Control",
                                      "Control","Control","Control",
                                      "Glaucoma","Glaucoma","Glaucoma","Glaucoma","Glaucoma",
                                      "Glaucoma","Glaucoma","Glaucoma","Glaucoma","Glaucoma",
                                      "Glaucoma","Glaucoma","Glaucoma","Glaucoma","Glaucoma",
                                      "Glaucoma","Glaucoma","Glaucoma","Glaucoma","Glaucoma",
                                      "Glaucoma","Glaucoma","Glaucoma","Glaucoma","Glaucoma",
                                      "Glaucoma","Glaucoma","Glaucoma","Glaucoma","Glaucoma",
                                      "Glaucoma","Glaucoma","Glaucoma","Glaucoma","Glaucoma",
                                      "Glaucoma","Glaucoma","Glaucoma","Glaucoma","Glaucoma",
                                      "Glaucoma","Glaucoma","Glaucoma","Glaucoma","Glaucoma"),
pairwise = F)
summary(res_GLM$pool_res_fixed)
summary(res_GLM$pool_res_random)


data <- data.frame(
  cond = rep(c("Control", "Glaucoma"), each = 9),
  cellType = rep(c("CD4 CTL", "CD4 Naive", "CD4 Proliferating", "CD4 TCM", "CD4 TEM", "dnT", "gdT", "MAIT", "Treg"), 2),
  count = c(17, 82, 7, 1354, 76, 14, 29, 45, 29, 74, 222, 10, 2912, 167, 54, 30, 38, 24)
)

fixed_pvalues <- data.frame(
  cellType = c('CD4 CTL',"CD4 Naive", "CD4 Proliferating", "CD4 TCM", "CD4 TEM", "dnT", "gdT", "MAIT", "Treg"),
  p.value = c(4.545449e-05, 1.261986e-04, 2.728365e-01, 8.709151e-08, 4.462173e-03, 1.415184e-01, 3.385058e-01, 5.852680e-01, 1.904133e-01)
)


get_significance <- function(p) {
  if (p < 0.001) {
    return("***")
  } else if (p < 0.01) {
    return("**")
  } else if (p < 0.05) {
    return("*")
  } else {
    return("NS")
  }
}

p <- ggplot(data, aes(x = cellType, y = count, fill = cond)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  theme_minimal() +
  labs(title = "Cell Counts for Control and Glaucoma Conditions",
       x = "Cell Type", y = "Count") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# Calculate significance labels
fixed_significance <- fixed_pvalues %>%
  mutate(label = sapply(p.value, get_significance))

# Calculate the y positions for the annotations
annotation_data <- data %>%
  group_by(cellType) %>%
  summarise(max_count = max(count)) %>%
  inner_join(fixed_significance, by = "cellType") %>%
  mutate(y_position = max_count + 100)

# Add fixed effect p-values with increased size and bold font
p <- p + geom_text(data = annotation_data, aes(x = cellType, y = y_position, label = label), 
                   position = position_dodge(width = 0.9), size = 4, fontface = "bold", inherit.aes = FALSE)

# Show plot
print(p)

data <- data.frame(
  cellType = c('CD4 CTL', 'CD4 Naive', 'CD4 Proliferating', 'CD4 TCM', 'CD4 TEM', 'dnT', 'gdT', 'MAIT', 'Treg',
               'CD4 CTL', 'CD4 Naive', 'CD4 Proliferating', 'CD4 TCM', 'CD4 TEM', 'dnT', 'gdT', 'MAIT', 'Treg',
               'CD4 CTL', 'CD4 Naive', 'CD4 Proliferating', 'CD4 TCM', 'CD4 TEM', 'dnT', 'gdT', 'MAIT', 'Treg',
               'CD4 CTL', 'CD4 Naive', 'CD4 Proliferating', 'CD4 TCM', 'CD4 TEM', 'dnT', 'gdT', 'MAIT', 'Treg',
               'CD4 CTL', 'CD4 Naive', 'CD4 Proliferating', 'CD4 TCM', 'CD4 TEM', 'dnT', 'gdT', 'MAIT', 'Treg',
               'CD4 CTL', 'CD4 Naive', 'CD4 Proliferating', 'CD4 TCM', 'CD4 TEM', 'dnT', 'gdT', 'MAIT', 'Treg',
               'CD4 CTL', 'CD4 Naive', 'CD4 Proliferating', 'CD4 TCM', 'CD4 TEM', 'dnT', 'gdT', 'MAIT', 'Treg',
               'CD4 CTL', 'CD4 Naive', 'CD4 Proliferating', 'CD4 TCM', 'CD4 TEM', 'dnT', 'gdT', 'MAIT', 'Treg'),
  subject = rep(1:8, each = 9),
  value = c(0, 17, 0, 225, 5, 1, 2, 0, 24,
            2, 23, 0, 302, 12, 2, 1, 3, 0,
            15, 42, 7, 827, 59, 11, 26, 42, 5,
            27, 8, 3, 394, 57, 15, 6, 5, 4,
            4, 44, 0, 575, 45, 7, 11, 9, 4,
            1, 59, 0, 566, 5, 12, 1, 23, 9,
            39, 1, 2, 385, 41, 2, 1, 0, 1,
            3, 110, 5, 992, 19, 18, 11, 1, 6)
)

# Add condition based on subject
data <- data %>%
  mutate(cond = ifelse(subject %in% 1:3, "Control", "Glaucoma"))

# P-values for cell types
fixed_pvalues <- data.frame(
  cellType = c('CD4 CTL', 'CD4 Naive', 'CD4 Proliferating', 'CD4 TCM', 'CD4 TEM', 'dnT', 'gdT', 'MAIT', 'Treg'),
  p.value = c(3.143458e-03, 8.739365e-04, 2.728365e-01, 1.203813e-12, 4.462173e-03, 8.530561e-01, 3.385058e-01, 2.063392e-02, 1.904133e-01)
)

# Define function to get significance label
get_significance <- function(p) {
  if (p < 0.001) {
    return("***")
  } else if (p < 0.01) {
    return("**")
  } else if (p < 0.05) {
    return("*")
  } else {
    return("NS")
  }
}

# Create ggplot object
p <- ggplot(data, aes(x = cellType, y = value, fill = cond)) +
  geom_boxplot(position = position_dodge()) +
  geom_jitter(position = position_dodge(0.75), size = 1) +
  theme_minimal() +
  labs(title = "Cell Counts for Control and Glaucoma Conditions",
       x = "Cell Type", y = "Count") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# Calculate significance labels
fixed_significance <- fixed_pvalues %>%
  mutate(label = sapply(p.value, get_significance))

# Calculate the y positions for the annotations
annotation_data <- data %>%
  group_by(cellType) %>%
  summarise(max_value = max(value)) %>%
  inner_join(fixed_significance, by = "cellType") %>%
  mutate(y_position = max_value + 50)  # Adjust the offset as needed

# Add fixed effect p-values with increased size and bold font
p <- p + geom_text(data = annotation_data, aes(x = cellType, y = y_position, label = label), 
                   position = position_dodge(width = 0.9), size = 4, fontface = "bold", inherit.aes = FALSE)

# Show plot
print(p)

PBMC <- SetIdent(PBMC, value = 'condition')
markers <- FindMarkers(PBMC, ident.1 = "Control", ident.2 = "Glaucoma",
                       min.cells.feature = 1,
                       min.cells.group = 1, assay = 'RNA', 
                       logfc.threshold = 0, slot = 'data')

# List of genes to check
genes_to_check <- c('CD3E', 'LAT', 'LCK', 'CD247', 'CDC42','LCP2')
FeaturePlot(PBMC, 'LCP2')

# Filter results for specified genes
significant_genes <- markers[rownames(markers) %in% genes_to_check, ]

# Check significance
significant_genes <- significant_genes %>%
  mutate(significance = ifelse(p_val_adj < 0.05, "Significant", "Not Significant"))

# Print results
print(significant_genes)
