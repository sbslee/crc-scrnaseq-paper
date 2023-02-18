library(dplyr)
library(tidyr)
library(Seurat)
library(patchwork)
library(plyr)
library(celldex)
library(SingleR)
library(ggplot2)
library(ggpubr)
library(stringr)
library(EnhancedVolcano)
library(ggraph)
library(scRepertoire)
library(infercnv)
# library(devtools)
# install_github("ncborcherding/scRepertoire@a7b3c79f5fa574e705e16bece47007788bbba6fb")
set.seed(1)

##################
# Pre-processing #
##################

data.dir <- '/home/sbslee/scRNAseq/multi/aggr/outs/count/filtered_feature_bc_matrix'
data <- Read10X(data.dir = data.dir)
asan1 <- CreateSeuratObject(counts = data, project = 'asan', min.cells = 3, min.features = 200)
asan1[['percent.mt']] <- PercentageFeatureSet(asan1, pattern = '^MT-')
saveRDS(asan1, './dat/asan1.rds') # 31,658 features and 550,102 cells

# asan1 <- readRDS('./dat/asan1.rds')

pdf('./out/qc-metrics.pdf', width = 10, height = 7)
VlnPlot(asan1, features = c('nFeature_RNA', 'nCount_RNA', 'percent.mt'), ncol = 3, pt.size = 0)
dev.off()

pdf('./out/qc-filter.pdf', width = 10, height = 7)
p1 <- FeatureScatter(asan1, feature1 = 'nCount_RNA', feature2 = 'percent.mt') + geom_hline(yintercept = 20, linetype = 'dashed') + theme(legend.position = 'none')
p2 <- FeatureScatter(asan1, feature1 = 'nCount_RNA', feature2 = 'nFeature_RNA') + geom_hline(yintercept = 5000, linetype = 'dashed') + theme(legend.position = 'none')
p1 + p2
dev.off()

#############
# Filtering #
#############

# asan1 <- readRDS('./dat/asan1.rds')

asan2 <- subset(asan1, subset = nFeature_RNA < 5000 & percent.mt < 20)

saveRDS(asan2, './dat/asan2.rds') # 31,658 features and 456,779 cells

##############
# Processing #
##############

# asan2 <- readRDS('./dat/asan2.rds')

asan3 <- NormalizeData(asan2)
asan3 <- FindVariableFeatures(asan3, selection.method = 'vst', nfeatures = 2000)

pdf('./out/qc-variable-genes.pdf', width = 10, height = 7)
VariableFeaturePlot(asan3)
dev.off()

asan3 <- ScaleData(asan3)
asan3 <- RunPCA(asan3, npcs = 50)
asan3 <- RunUMAP(asan3, dims = 1:50)
asan3 <- FindNeighbors(asan3, dims = 1:50)
asan3 <- FindClusters(asan3, resolution = 0.05)

batch <- sub('.*-', '', names(asan3$orig.ident))
aggr <- read.csv('/home/sbslee/scRNAseq/multi/aggr-info.csv')

## Add batch information

sample <- as.character(aggr$sample_id)
names(sample) <- as.character(rownames(aggr))
sample = revalue(x = batch, sample)

patient <- as.character(aggr$donor)
names(patient) <- as.character(rownames(aggr))
patient = revalue(x = batch, patient)

tissue <- as.character(aggr$origin)
tissue <- str_sub(tissue, start = -1, end = -1)
names(tissue) <- as.character(rownames(aggr))
tissue = revalue(x = batch, tissue)

names(batch) <- rownames(x = asan3@meta.data)
names(sample) <- rownames(x = asan3@meta.data)
names(patient) <- rownames(x = asan3@meta.data)
names(tissue) <- rownames(x = asan3@meta.data)

asan3[['batch']] <- batch
asan3[['sample']] <- sample
asan3[['patient']] <- patient
asan3[['tissue']] <- tissue

pdf('./out/5gex-umap-batch.pdf', width = 10, height = 7)
DimPlot(asan3, reduction = 'umap', group.by = 'sample')
dev.off()

pdf('./out/5gex-umap-leiden.pdf', width = 10, height = 7)
DimPlot(asan3, reduction = 'umap', label = T)
dev.off()

## Cell type identification

# Epithelial cells: EPCAM, KRT8, KRT18
# T-cells: CD3D
# B-cells: CD79A
# Myeloid: LYZ
# Fibroblast: COL1A1
# Endothelial: CLDN5

pdf('./out/5gex-canonical-markers.pdf', width = 10, height = 7)
VlnPlot(asan3, group.by = 'seurat_clusters', features = c('EPCAM', 'KRT18', 'KRT8', 'CD3D', 'CD79A', 'LYZ', 'COL1A1', 'CLDN5'), pt.size = 0, ncol = 4)
dev.off()

asan3[['celltype']] <- mapvalues(asan3@meta.data$seurat_clusters, from=c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13), to=c('T-cells', 'B-cells', 'Epithelial', 'B-cells', 'Myeloid', 'Fibroblast', 'B-cells', 'Endothelial', 'Fibroblast', 'Epithelial', 'Other', 'Other', 'Other', 'Other'))

pdf('./out/5gex-umap-celltype.pdf', width = 10, height = 7)
DimPlot(asan3, reduction = 'umap', label = T, group.by = 'celltype')
dev.off()

## Add clinical data

df1 <- read.csv('/home/sbslee/scRNAseq/analysis/patient-table-20221203.csv')
df1$PatientID <- paste(df1$PatientID, 'N', sep='')
df2 <- df1
df2$PatientID <- gsub('N', 'T', df2$PatientID)
df3 <- rbind(df1, df2)
df4 <- merge(asan3[[]], df3, by.x='sample', by.y='PatientID')

cols <- c('Sex', 'Age', 'DifferentiationLevel', 'Location', 'LymphovascularInvasion', 'VenousInvasion', 'PerineuralInvasion', 'LymphNodeMetastasis', 'MSI', 'TNM', 'TILsIntra', 'TILsPeri', 'TumorDensity')

for (col in cols){
  asan3 <- AddMetaData(asan3, unlist(df4[col]), col)
}

saveRDS(asan3, './dat/asan3.rds')

# Tumor-infiltrating lymphocyte

asan3 <- AddMetaData(asan3, ifelse(asan3@meta.data$TumorDensity > 22.345, 'High', 'Low'), 'BinTumorDensity')
asan3 <- AddMetaData(asan3, ifelse(asan3@meta.data$TILsPeri > 9.1565, 'High', 'Low'), 'BinTILsPeri')
asan3 <- AddMetaData(asan3, ifelse(asan3@meta.data$TILsIntra > 6.532, 'High', 'Low'), 'BinTILsIntra')

####################
# Epithelial cells #
####################

# asan3 <- readRDS('./dat/asan3.rds')
# asan3.epi <- readRDS('./dat/asan3-epi.rds')

asan3.epi <- subset(x = asan3, subset = celltype == 'Epithelial')
asan3.epi <- RunPCA(asan3.epi, npcs = 50)
asan3.epi <- RunUMAP(asan3.epi, dims = 1:50)
asan3.epi <- FindNeighbors(asan3.epi, dims = 1:50)
asan3.epi <- FindClusters(asan3.epi, resolution = 0.05)

saveRDS(asan3.epi, './dat/asan3-epi.rds') # 46,086 cells

pdf('./out/epi-umap-tissue.pdf', width = 10, height = 7)
DimPlot(asan3.epi, reduction = 'umap', group.by = 'tissue')
dev.off()

pdf('./out/epi-umap-patient.pdf', width = 10, height = 7)
DimPlot(asan3.epi, reduction = 'umap', group.by = 'patient')
dev.off()

pdf('./out/epi-umap-msi.pdf', width = 10, height = 7)
DimPlot(asan3.epi, reduction = 'umap', group.by = 'MSI', split = 'tissue')
dev.off()

# For InferCNV
counts_matrix <- asan3.epi@assays$RNA@counts[,colnames(asan3.epi)]
cells <- rownames(asan3.epi@meta.data)
annot <- asan3.epi@meta.data$tissue
df <- cbind(cells, annot)
write.table(df, './annot-file.tsv', sep = '\t', row.names = F, quote = F, col.names = F)

# DEG analysis
markers <- FindMarkers(asan3.epi, ident.1 = 'T', ident.2 = 'N', group.by = 'tissue')
write.csv(markers, './out/epi-deg-table.csv')

pdf('./out/epi-umap-deg.pdf')
FeaturePlot(asan3.epi, features = c('HES4', 'PADI2', 'ID3', 'RPL11', 'IFI6', 'TMEM54'))
dev.off()

# Tumor-infiltrating lymphocyte
pdf('./out/epi-til-intra.pdf', width = 10, height = 7)
FeaturePlot(asan3.epi, features = 'TILsIntra', order = T, split = 'tissue')
dev.off()

pdf('./out/epi-til-peri.pdf', width = 10, height = 7)
FeaturePlot(asan3.epi, features = 'TILsPeri', order = T, split = 'tissue')
dev.off()

pdf('./out/epi-til-density.pdf', width = 10, height = 7)
FeaturePlot(asan3.epi, features = 'TumorDensity', order = T, split = 'tissue')
dev.off()

# Immuno-oncology targets
pdf('./out/epi-immune-ici.pdf')
FeaturePlot(asan3.epi, features = c('CD274', 'PDCD1', 'CTLA4'))
dev.off()

pdf('./out/epi-immune-checkpoint.pdf')
FeaturePlot(asan3.epi, features = c('LAG3', 'HAVCR2', 'TIGIT', 'VSIR', 'CD276', 'BTLA'))
dev.off()

pdf('./out/epi-immune-mage.pdf')
FeaturePlot(asan3.epi, features = c('MAGEA10', 'MAGEA11', 'MAGEA12', 'MAGEA2B', 'MAGEA3', 'MAGEA4', 'MAGEA4-AS1', 'MAGEA6', 'MAGEA8', 'MAGEA8-AS1', 'MAGEB17', 'MAGEB2', 'MAGEB5', 'MAGEB6'))
dev.off()

##############
# Fibroblast #
##############

# asan3 <- readRDS('./dat/asan3.rds')
# asan3.caf <- readRDS('./dat/asan3-caf.rds')

asan3.caf <- subset(x = asan3, subset = celltype == 'Fibroblast')
asan3.caf <- RunPCA(asan3.caf, npcs = 50)
asan3.caf <- RunUMAP(asan3.caf, dims = 1:50)
asan3.caf <- FindNeighbors(asan3.caf, dims = 1:50)
asan3.caf <- FindClusters(asan3.caf, resolution = 0.05)

saveRDS(asan3.caf, './dat/asan3-caf.rds') # 22,758 cells

pdf('./out/caf-umap-tissue.pdf', width = 10, height = 7)
DimPlot(asan3.caf, reduction = 'umap', group.by = 'tissue')
dev.off()

pdf('./out/caf-umap-patient.pdf', width = 10, height = 7)
DimPlot(asan3.caf, reduction = 'umap', group.by = 'patient')
dev.off()

# DEG analysis
markers <- FindMarkers(asan3.caf, ident.1 = 'T', ident.2 = 'N', group.by = 'tissue')
write.csv(markers, './out/caf-deg-table.csv')

pdf('./out/caf-umap-deg.pdf')
FeaturePlot(asan3.caf, features = c('HES4', 'FBLIM1', 'TINAGL1', 'TXNIP', 'DCN', 'MGP'))
dev.off()

#################
# Clinical data #
#################

asan3 <- readRDS('./dat/asan3.rds')

asan3 <- AddMetaData(asan3, ifelse(asan3@meta.data$Age < 65, '<65', '>=65'), 'Age65')

pdf('./out/5gex-umap-age.pdf', width = 10, height = 7)
FeaturePlot(asan3, features = 'Age', order = T, split = 'tissue')
dev.off()

asan3.tumor <- subset(x = asan3, subset = tissue == 'T')
asan3.normal <- subset(x = asan3, subset = tissue == 'N')

saveRDS(asan3.tumor, './dat/asan3-tumor.rds') # 221,704 cells
saveRDS(asan3.normal, './dat/asan3-normal.rds') # 235,075 cells

df.all <- as.data.frame(table(asan3$celltype))
df.tumor <- as.data.frame(table(asan3.tumor$celltype))
df.normal <- as.data.frame(table(asan3.normal$celltype))
df.all['Tissue'] = 'All'
df.tumor['Tissue'] = 'T'
df.normal['Tissue'] = 'N'
df = rbind(df.all, df.tumor, df.normal)
colnames(df) = c('Type', 'Fraction', 'Tissue')

pdf('./out/5gex-cell-prop.pdf', width = 10, height = 7)
ggplot(df, aes(fill = Type, y = Fraction, x = Tissue)) + geom_bar(position = 'fill', stat = 'identity')
dev.off()

plot_one <- function(object, split.by, legend=F) {
  obj.list <- SplitObject(object, split.by = split.by)
  for (i in 1:length(obj.list)) {
    obj.list[[i]] <- as.data.frame(table(obj.list[[i]]$celltype))
    obj.list[[i]][split.by] <- names(obj.list[i])
  }
  df <- bind_rows(obj.list)
  colnames(df) <- c('Type', 'Fraction', split.by)
  p <- ggplot(df, aes_string(fill = 'Type', y = 'Fraction', x = split.by)) + geom_bar(position = 'fill', stat = 'identity')
  if (!legend) {
    p <- p + theme(legend.position='none')
  }
  return(p)
}

t.msi <- plot_one(asan3.tumor, 'MSI')
t.tnm <- plot_one(asan3.tumor, 'TNM')
t.sex <- plot_one(asan3.tumor, 'Sex')
t.dif <- plot_one(asan3.tumor, 'DifferentiationLevel')
t.loc <- plot_one(asan3.tumor, 'Location')
t.lvi <- plot_one(asan3.tumor, 'LymphovascularInvasion')
t.vni <- plot_one(asan3.tumor, 'VenousInvasion')
t.pni <- plot_one(asan3.tumor, 'PerineuralInvasion')
t.lnm <- plot_one(asan3.tumor, 'LymphNodeMetastasis')
t.tia <- plot_one(asan3.tumor, 'BinTILsIntra')
t.tpi <- plot_one(asan3.tumor, 'BinTILsPeri')
t.tey <- plot_one(asan3.tumor, 'BinTumorDensity')
t.age <- plot_one(asan3.tumor, 'Age65', legend = T)

n.msi <- plot_one(asan3.normal, 'MSI')
n.tnm <- plot_one(asan3.normal, 'TNM')
n.sex <- plot_one(asan3.normal, 'Sex')
n.dif <- plot_one(asan3.normal, 'DifferentiationLevel')
n.loc <- plot_one(asan3.normal, 'Location')
n.lvi <- plot_one(asan3.normal, 'LymphovascularInvasion')
n.vni <- plot_one(asan3.normal, 'VenousInvasion')
n.pni <- plot_one(asan3.normal, 'PerineuralInvasion')
n.lnm <- plot_one(asan3.normal, 'LymphNodeMetastasis')
n.tia <- plot_one(asan3.normal, 'BinTILsIntra')
n.tpi <- plot_one(asan3.normal, 'BinTILsPeri')
n.tey <- plot_one(asan3.normal, 'BinTumorDensity')
n.age <- plot_one(asan3.normal, 'Age65', legend = T)

pdf('./out/5gex-clin-tumor.pdf', width = 9, height = 12)
ggarrange(t.msi, t.tnm, t.sex, t.dif, t.loc, t.lvi, t.vni, t.pni, t.lnm, t.tia, t.tpi, t.tey, t.age, ncol = 4, nrow = 4, common.legend = T, legend = 'right')
dev.off()

pdf('./out/5gex-clin-normal.pdf', width = 9, height = 12)
ggarrange(n.msi, n.tnm, n.sex, n.dif, n.loc, n.lvi, n.vni, n.pni, n.lnm, n.tia, n.tpi, n.tey, n.age, ncol = 4, nrow = 4, common.legend = T, legend = 'right')
dev.off()

################
# TCR analysis #
################

aggr <- read.csv('/home/sbslee/scRNAseq/multi/aggr-info.csv')
aggr$vdj_t <- paste(aggr$sample_outs, '/vdj_t/filtered_contig_annotations.csv', sep = '')
contig_list <- list()

for (file in aggr$vdj_t) {
  df <- read.csv(file)
  contig_list <- append(contig_list, list(df))
}

combined <- combineTCR(
  contig_list,
  samples = aggr$donor,
  ID = substr(aggr$sample_id, 3, 3),
  cells = 'T-AB'
)

for (i in seq_along(combined)) {
  combined[[i]]$barcode <- sapply(strsplit(combined[[i]]$barcode, '_'), function(f) f[[3]])
  combined[[i]]$barcode <- paste(sapply(strsplit(combined[[i]]$barcode, '-'), function(f) f[[1]]), i, sep = '-')
  combined[[i]]$sample_id <- aggr$sample_id[[i]]
}

asan3.imm <- readRDS('./dat/asan3-imm.rds')

asan3.tcr <- combineExpression(
  combined,
  asan3.imm,
  cloneCall = 'aa',
  proportion = F,
  filterNA = T,
  cloneTypes = c(Single = 1, Small = 5, Medium = 20, Large = 100, Hyperexpanded = 2000)
)

targets <- c('Treg', 'NK_CD56bright', 'NK Proliferating', 'NK', 'MAIT', 'gdT', 'dnT', 'CD8 TEM', 'CD8 TCM', 'CD8 Proliferating', 'CD8 Naive', 'CD4 TEM', 'CD4 TCM', 'CD4 Proliferating', 'CD4 Native', 'CD4 CTL')

asan3.tcr <- subset(asan3.tcr, subset = final.celltype %in% targets)

saveRDS(asan3.tcr, './dat/asan3-tcr.rds') # 138,963 cells

asan3.tcr <- readRDS('./dat/asan3-tcr.rds')

# Filter out cells using gene expression

for (i in seq_along(combined)) {
  keep <- combined[[i]]$barcode %in% colnames(asan3.tcr)
  combined[[i]] <- combined[[i]][keep, ]
}

n <- 0

for (i in seq_along(combined)) {
  n <- n + dim(combined[[i]])[1]
}

saveRDS(combined, './dat/combined.rds') #138,963 cells

# Compare the proportion of subtypes between tumor and normal 

asan3.tcr.t <- subset(x = asan3.tcr, subset = tissue == 'T')
asan3.tcr.n <- subset(x = asan3.tcr, subset = tissue == 'N')

df.t <- as.data.frame(asan3.tcr.t@meta.data %>% group_by(as.factor(final.celltype), as.factor(cloneType), .drop = F) %>% dplyr::count())
df.n <- as.data.frame(asan3.tcr.n@meta.data %>% group_by(as.factor(final.celltype), as.factor(cloneType), .drop = F) %>% dplyr::count())

colnames(df.t) = c('final.celltype', 'cloneType', 'count')
colnames(df.n) = c('final.celltype', 'cloneType', 'count')

write.csv(df.t, './out/tcr-prop-tumor.csv', row.names = F)
write.csv(df.n, './out/tcr-prop-normal.csv', row.names = F)

legend.order <- c('Hyperexpanded (100 < X <= 2000)', 'Large (20 < X <= 100)', 'Medium (5 < X <= 20)', 'Small (1 < X <= 5)', 'Single (0 < X <= 1)')

df.t$cloneType <- factor(df.t$cloneType, levels = legend.order)
df.n$cloneType <- factor(df.n$cloneType, levels = legend.order)

pdf('./out/tcr-prop-subtypes.pdf', width = 15, height = 6)
p1 <- ggplot(df.n, aes(fill = cloneType, y = count, x = final.celltype)) + geom_bar(position = 'fill', stat = 'identity') + coord_flip() + ggtitle('Normal') + theme(plot.title = element_text(size = 20, hjust = 0.5), axis.title.x = element_blank(), axis.title.y = element_blank())
p2 <- ggplot(df.t, aes(fill = cloneType, y = count, x = final.celltype)) + geom_bar(position = 'fill', stat = 'identity') + coord_flip() + ggtitle('Tumor') + theme(plot.title = element_text(size = 20, hjust = 0.5), axis.title.x = element_blank(), axis.title.y = element_blank())
ggarrange(p1, p2, ncol = 2, nrow = 1, common.legend = T, legend = 'right')
dev.off()

# Make scRepertoire plots

pdf('./out/tcr-quantContig.pdf', width = 10, height = 10)
quantContig(combined, cloneCall = 'aa', scale = F)
dev.off()

pdf('./out/tcr-abundanceContig.pdf', width = 10, height = 10)
abundanceContig(combined, cloneCall = 'aa', scale = F)
dev.off()

pdf('./out/tcr-lengthContig.pdf', width = 10, height = 10)
lengthContig(combined, cloneCall = 'aa')
dev.off()

pdf('./out/tcr-compareClonotypes-combined.pdf', width = 20, height = 10)
p.15 <- compareClonotypes(combined, numbers = 10, samples = c('15_N', '15_T'), cloneCall = 'aa', graph = 'alluvial') + theme(legend.position='none')
p.16 <- compareClonotypes(combined, numbers = 10, samples = c('16_N', '16_T'), cloneCall = 'aa', graph = 'alluvial') + theme(legend.position='none')
p.35 <- compareClonotypes(combined, numbers = 10, samples = c('35_N', '35_T'), cloneCall = 'aa', graph = 'alluvial') + theme(legend.position='none')
p.41 <- compareClonotypes(combined, numbers = 10, samples = c('41_N', '41_T'), cloneCall = 'aa', graph = 'alluvial') + theme(legend.position='none')
p.47 <- compareClonotypes(combined, numbers = 10, samples = c('47_N', '47_T'), cloneCall = 'aa', graph = 'alluvial') + theme(legend.position='none')
p.48 <- compareClonotypes(combined, numbers = 10, samples = c('48_N', '48_T'), cloneCall = 'aa', graph = 'alluvial') + theme(legend.position='none')
ggarrange(p.15, p.16, p.35, p.41, p.47, p.48, ncol = 3, nrow = 2)
dev.off()

# Make individual plots

for (donor in unique(aggr$donor)) {
  normal <- paste(donor, '_N', sep = '')
  tumor <- paste(donor, '_T', sep = '')
  pdf(paste('./out/tcr-compareClonotypes-', donor,'.pdf', sep = ''), width = 5, height = 5)
  print(compareClonotypes(combined, numbers = 10, samples = c(normal, tumor), cloneCall = 'aa', graph = 'alluvial') + theme(legend.position='none'))
  dev.off()
}

pdf('./out/tcr-clonalHomeostasis.pdf', width = 15, height = 5)
clonalHomeostasis(combined, cloneCall = 'aa') + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

pdf('./out/tcr-clonesizeDistribution.pdf', width = 10, height = 10)
clonesizeDistribution(combined, cloneCall = 'aa', method = 'ward.D2')
dev.off()

pdf('./out/tcr-clonalOverlap.pdf', width = 10, height = 10)
clonalOverlap(combined, cloneCall = 'aa', method = 'morisita')
dev.off()

# pdf('./out/tcr-clonalDiversity.pdf', width = 10, height = 10)
# clonalDiversity(combined, cloneCall = 'aa', n.boots = 100, group.by = 'sample',  x.axis = 'ID')
# dev.off()

# pdf('./out/tcr-scatterClonotype.pdf', width = 10, height = 10)
# scatterClonotype(combined, cloneCall = 'aa', x.axis = 'N', y.axis = 'T', graph = 'proportion')
# dev.off()

################
# Immune cells #
################

library(Seurat)
library(Azimuth)
library(SeuratData)
library(patchwork)

# conda create -n azimuth2 -c conda-forge r-seuratdisk r-hdf5r bioconductor-glmgampoi r-seurat r-devtools
# asan3 <- readRDS('./dat/asan3.rds')
# asan3.imm <- subset(x = asan3, subset = celltype %in% c('T-cells', 'B-cells', 'Myeloid'))
# saveRDS(asan3.imm, './dat/asan3-imm.rds') # 381,488 cells

asan3.imm <- readRDS('./dat/asan3-imm.rds')

# asan3.imm <- RunAzimuth(asan3.imm, reference = 'pbmcref')

asan3.imm <- AddMetaData(asan3.imm, ifelse(asan3.imm@meta.data$TumorDensity > 22.345, 'High', 'Low'), 'BinTumorDensity')
asan3.imm <- AddMetaData(asan3.imm, ifelse(asan3.imm@meta.data$TILsPeri > 9.1565, 'High', 'Low'), 'BinTILsPeri')
asan3.imm <- AddMetaData(asan3.imm, ifelse(asan3.imm@meta.data$TILsIntra > 6.532, 'High', 'Low'), 'BinTILsIntra')
asan3.imm <- AddMetaData(asan3.imm, ifelse(asan3.imm@meta.data$Age < 65, '<65', '>=65'), 'Age65')

pdf('./out/imm-umap-together-l2.pdf', width = 10, height = 7)
p1 <- DimPlot(asan3.imm, reduction = 'umap', group.by = 'predicted.celltype.l2')
LabelClusters(plot = p1, id = 'predicted.celltype.l2')
dev.off()

pdf('./out/imm-umap-split-l2.pdf', width = 10, height = 7)
DimPlot(asan3.imm, reduction = 'umap', split.by = 'predicted.celltype.l2', ncol = 6, group.by = 'predicted.celltype.l2') + NoLegend()
dev.off()

double.positive <- GetAssayData(object = asan3.imm)['SELL', ] > 0 & GetAssayData(object = asan3.imm)['CCR7', ] > 0

asan3.imm <- AddMetaData(asan3.imm, double.positive, 'double.positive')

# Manually assign TCM vs. TEM
one_row <- function(r) {
  if (r['predicted.celltype.l2'] == 'CD4 TCM' || r['predicted.celltype.l2'] == 'CD4 TEM') {
    if (r['double.positive']) {
      'CD4 TCM'
    } else {
      'CD4 TEM'
    }
  } else if (r['predicted.celltype.l2'] == 'CD8 TCM' || r['predicted.celltype.l2'] == 'CD8 TEM') {
    if (r['double.positive']) {
      'CD8 TCM'
    } else {
      'CD8 TEM'
    }
  } else {
    r['predicted.celltype.l2']
  }
}

asan3.imm <- AddMetaData(asan3.imm, apply(asan3.imm@meta.data, 1, one_row), 'final.celltype')

pdf('./out/imm-double-positive.pdf', width = 10, height = 7)
VlnPlot(asan3.imm, group.by = 'final.celltype', features = c('CCR7', 'SELL'), pt.size = 0, ncol = 1)
dev.off()

pdf('./out/imm-umap-tissue.pdf', width = 10, height = 7)
DimPlot(asan3.imm, reduction = 'umap', split.by = 'tissue', ncol = 6, group.by = 'final.celltype')
dev.off()

pdf('./out/imm-umap-split-final.pdf', width = 10, height = 7)
DimPlot(asan3.imm, reduction = 'umap', split.by = 'final.celltype', ncol = 6, group.by = 'final.celltype') + NoLegend()
dev.off()

pdf('./out/imm-umap-together-final.pdf', width = 10, height = 7)
p1 <- DimPlot(asan3.imm, reduction = 'umap', group.by = 'final.celltype')
LabelClusters(plot = p1, id = 'final.celltype')
dev.off()

asan3.imm.tumor <- subset(x = asan3.imm, subset = tissue == 'T')
asan3.imm.normal <- subset(x = asan3.imm, subset = tissue == 'N')

df.all <- as.data.frame(table(asan3.imm$final.celltype))
df.tumor <- as.data.frame(table(asan3.imm.tumor$final.celltype))
df.normal <- as.data.frame(table(asan3.imm.normal$final.celltype))
df.all['Tissue'] = 'All'
df.tumor['Tissue'] = 'T'
df.normal['Tissue'] = 'N'
df = rbind(df.all, df.tumor, df.normal)
colnames(df) = c('Type', 'Fraction', 'Tissue')

write.csv(df.all, './out/imm-prop-all.csv', row.names=F)
write.csv(df.tumor, './out/imm-prop-tumor.csv', row.names=F)
write.csv(df.normal, './out/imm-prop-normal.csv', row.names=F)

pdf('./out/imm-cell-prop.pdf', width = 10, height = 7)
ggplot(df, aes(fill = Type, y = Fraction, x = Tissue)) + geom_bar(position = 'fill', stat = 'identity')
dev.off()

plot_one <- function(object, split.by, legend=F) {
  obj.list <- SplitObject(object, split.by = split.by)
  for (i in 1:length(obj.list)) {
    obj.list[[i]] <- as.data.frame(table(obj.list[[i]]$final.celltype))
    obj.list[[i]][split.by] <- names(obj.list[i])
  }
  df <- bind_rows(obj.list)
  colnames(df) <- c('Type', 'Fraction', split.by)
  p <- ggplot(df, aes_string(fill = 'Type', y = 'Fraction', x = split.by)) + geom_bar(position = 'fill', stat = 'identity')
  if (!legend) {
    p <- p + theme(legend.position='none')
  }
  return(p)
}

t.msi <- plot_one(asan3.imm.tumor, 'MSI')
t.tnm <- plot_one(asan3.imm.tumor, 'TNM')
t.sex <- plot_one(asan3.imm.tumor, 'Sex')
t.dif <- plot_one(asan3.imm.tumor, 'DifferentiationLevel')
t.loc <- plot_one(asan3.imm.tumor, 'Location')
t.lvi <- plot_one(asan3.imm.tumor, 'LymphovascularInvasion')
t.vni <- plot_one(asan3.imm.tumor, 'VenousInvasion')
t.pni <- plot_one(asan3.imm.tumor, 'PerineuralInvasion')
t.lnm <- plot_one(asan3.imm.tumor, 'LymphNodeMetastasis')
t.tia <- plot_one(asan3.imm.tumor, 'BinTILsIntra')
t.tpi <- plot_one(asan3.imm.tumor, 'BinTILsPeri')
t.tey <- plot_one(asan3.imm.tumor, 'BinTumorDensity')
t.age <- plot_one(asan3.imm.tumor, 'Age65', legend = T)

n.msi <- plot_one(asan3.imm.normal, 'MSI')
n.tnm <- plot_one(asan3.imm.normal, 'TNM')
n.sex <- plot_one(asan3.imm.normal, 'Sex')
n.dif <- plot_one(asan3.imm.normal, 'DifferentiationLevel')
n.loc <- plot_one(asan3.imm.normal, 'Location')
n.lvi <- plot_one(asan3.imm.normal, 'LymphovascularInvasion')
n.vni <- plot_one(asan3.imm.normal, 'VenousInvasion')
n.pni <- plot_one(asan3.imm.normal, 'PerineuralInvasion')
n.lnm <- plot_one(asan3.imm.normal, 'LymphNodeMetastasis')
n.tia <- plot_one(asan3.imm.normal, 'BinTILsIntra')
n.tpi <- plot_one(asan3.imm.normal, 'BinTILsPeri')
n.tey <- plot_one(asan3.imm.normal, 'BinTumorDensity')
n.age <- plot_one(asan3.imm.normal, 'Age65', legend = T)

pdf('./out/imm-clin-tumor.pdf', width = 9, height = 12)
ggarrange(t.msi, t.tnm, t.sex, t.dif, t.loc, t.lvi, t.vni, t.pni, t.lnm, t.tia, t.tpi, t.tey, t.age, ncol = 4, nrow = 4, common.legend = T, legend = 'right')
dev.off()

pdf('./out/imm-clin-normal.pdf', width = 9, height = 12)
ggarrange(n.msi, n.tnm, n.sex, n.dif, n.loc, n.lvi, n.vni, n.pni, n.lnm, n.tia, n.tpi, n.tey, n.age, ncol = 4, nrow = 4, common.legend = T, legend = 'right')
dev.off()
