library(dplyr)
library(Seurat)
library(patchwork)
library(stringr)


# Load the dataset
brain.data <- read.table("Brain_Neurons-counts.csv", header = TRUE, 
                         sep=',', 
                         row.names = 1)
# Load metadata 
metadata_file <- read.csv("metadata_FACS.csv")
metadata_file <- na.omit(metadata_file)


# Initialize the Seurat object with the raw (non-normalized data)
brain.counts <- CreateSeuratObject(counts = brain.data,
                                   project = "brain",
                                   names.field = 1,
                                   min.cells = 3,
                                   min.features = 200)

##### QC

# seems there's no any mt genes in the initial counts dataset
brain.counts[["percent.mt"]] <- PercentageFeatureSet(brain.counts, pattern = "^mt-")

# look at RNA and Genes quantities distributions
VlnPlot(brain.counts, features = c("nFeature_RNA", "nCount_RNA"))
FeatureScatter(brain.counts, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

# here we see approximate droplets thresholds
hist(brain.counts@meta.data[,'nFeature_RNA'], breaks=40)
# here we see approximate doublets thresholds
hist(brain.counts@meta.data[,'nCount_RNA'], breaks=40)

# filtering "by eye"
brain.counts.filtered <- subset(brain.counts, 
               subset = nFeature_RNA > 700 & nFeature_RNA < 7500 & nCount_RNA < 6.0e+06)


# here we observe there are no outliers anymore
FeatureScatter(brain.counts.filtered, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
hist(brain.counts.filtered@meta.data[,'nCount_RNA'], breaks=40)
hist(brain.counts.filtered@meta.data[,'nFeature_RNA'], breaks=40)


# extracting plate.barcode from column name
brain.counts.filtered@meta.data$plate.barcode <- str_extract(colnames(brain.counts.filtered), "(?<=\\.)[A-Z0-9]{6,9}")
# storing cells.names as merge() removes indexes (cell.names)
cells.names <- rownames(brain.counts.filtered@meta.data) 
# adding subtissue and sex to each cell in Seurat object metadata
brain.counts.filtered@meta.data <- left_join(brain.counts.filtered@meta.data, 
                                metadata_file[, c('plate.barcode','subtissue','mouse.sex')], 
                                by = 'plate.barcode')
# adding back cells.names after merge
rownames(brain.counts.filtered@meta.data) <- cells.names


#### Normalization (log10k)
brain.counts.norm <- NormalizeData(brain.counts.filtered)

#### Feature selection
brain.counts.fs <- FindVariableFeatures(brain.counts.norm, selection.method = "vst", nfeatures = 5000)

top10 <- head(VariableFeatures(brain.counts.fs), 10)
plot1 <- VariableFeaturePlot(brain.counts.fs)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

#### Scaling + PCA
all.genes <- rownames(brain.counts.fs)
brain.counts.scaled <- ScaleData(brain.counts.fs, features = all.genes)
brain.counts.pca <- RunPCA(brain.counts.scaled, 
                           features = VariableFeatures(object = brain.counts.scaled))

DimPlot(brain.counts.pca, reduction = "pca")
# see how many PCs really contribute to variance
ElbowPlot(brain.counts.pca)

##### Clustering

brain.counts.clustered <- FindNeighbors(brain.counts.pca, dims = 1:11)
brain.counts.clustered <- FindClusters(brain.counts.clustered, resolution = 0.5)
Idents(brain.counts.clustered) # see the number of clusters

#### UMAP

brain.counts.umap <- RunUMAP(brain.counts.clustered, dims = 1:11)
plot1 <- DimPlot(brain.counts.umap, reduction = "umap", group.by = "mouse.sex", split.by = )
plot2 <- DimPlot(brain.counts.umap, reduction = "umap", group.by = "subtissue")
plot3 <- DimPlot(brain.counts.umap, reduction = "umap", group.by = "seurat_clusters")
plot1 + plot3

#### Barplot for M/F distribution among clusters

counts <- table(brain.counts.umap@meta.data$mouse.sex,
                brain.counts.umap@meta.data$seurat_clusters)
prop.counts <- prop.table(counts, margin = 2)

barplot(prop.counts,
        xlim=c(-0.5, ncol(prop.counts)*3)+2,
        col=c('red','blue'),
        beside=T,
        args.legend = list(
          x=ncol(prop.counts)*3+4,
          y=1.05),
        legend=c('F','M'))

# We see that 10th clusters have abundant female cells while other clusters consist of primarily male cells

#### Differential Expression

Idents(brain.counts.umap) <- brain.counts.umap$subtissue # reset identity clusters to find DE between subtissues
hippo_cortex.markers <- FindMarkers(brain.counts.umap, ident.1 = "Hippocampus", ident.2 = "Cortex", min.pct = 0.25)
hippo_cortex.markers %>% slice_max(n = 30, order_by = avg_log2FC) %>% rownames -> top30
DoHeatmap(brain.counts.umap, features = top30) + NoLegend()

# On the Heatmap we see 3 more or less well distinguishable Hippocampus cell groups that significantly differ from Cortex cells by the selected gene set expression (upregulation)
# Possibly, these cells are subclusters of hippocampus cells wielding such an upregulated genes expression for performing any function.
# For instance:
#   Ly6c1 - lymphocyte antigen 6 complex, surface markers of monocytes/macrophages
#   Cldn5 - Claudin-5, component of tight junction strands (reported to be expressed in endothelial cells to form blood-brain-barrier)
#   Jun - regulation of cell proliferation
#   Ly6a - marker of murine (HSCs), also plays an important role in regulating T-cell proliferation in mice 

# Initially it seems that some inflammatory process takes place in hippocampus, involving angiogenesis, cell proliferation and leukocytes presence