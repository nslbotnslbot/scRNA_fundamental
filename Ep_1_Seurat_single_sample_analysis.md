# scRNA_fundamental

## Step 1: Install the necessary packages

Lists of packages I used in this episode.

```r
rm(list = ls())

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("multtest")
BiocManager::install("Seurat")
BiocManager::install("SeuratData")
BiocManager::install("dplyr")
BiocManager::install("mindr")
BiocManager::install("tidyverse")
BiocManager::install("remotes")
BiocManager::install("R.utils")
remotes::install_github('satijalab/seurat-data')

library(multtest)
library(Seurat)
library(dplyr)
library(mindr)
library(tidyverse)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(R.utils)
```

10x genomics have many samples to practice, we chose Peripheral Blood Mononuclear Cells (PBMCs).

```r
download.file('https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz','pbmc3k.gz')
untar(gunzip("pbmc3k.gz"))
```

## Step 2 Download sample data or local data#

### Peripheral Blood Mononuclear Cells (PBMCs) —— 外周血单个核细胞
```r
### We try the custom way, Local 10x data##

pbmc3k.data <- Read10X(data.dir = "./filtered_gene_bc_matrices/hg19/")
pbmc3k <- CreateSeuratObject(counts = pbmc3k.data, project = "pbmc3k", min.cells = 3, min.features = 200)

### min.cells = how many cell types a gene is expressed in at least, min.features = how many genes a cell expresses at least. Only when the conditions are met will the gene be retained

### Have a glimpse
pbmc3k

###An object of class Seurat 
###13714 features across 2700 samples within 1 assay 
###Active assay: RNA (13714 features, 0 variable features)
###1 layer present: counts

ncol(pbmc3k)
### 2700

ncol(pbmc3k.data)
### 2700

pbmc3k_express_matrix <- as.data.frame(GetAssayData(pbmc3k[["RNA"]], slot = "counts"))
### saving the matrix

write.table(pbmc3k_express_matrix,'testcount.txt', sep = '\t')
### write into the 'txt' file, but wasting space

pbmc3k[["percent.mt"]] <- PercentageFeatureSet(pbmc3k, pattern = "^MT-")
### Lowercase "mt" for mouse, Uppercase "MT" for human

head(pbmc3k@meta.data,5)

###                 orig.ident nCount_RNA nFeature_RNA
###AAACATACAACCAC-1     pbmc3k       2419          779
###AAACATTGAGCTAC-1     pbmc3k       4903         1352
###AAACATTGATCAGC-1     pbmc3k       3147         1129
###AAACCGTGCTTCCG-1     pbmc3k       2639          960
###AAACCGTGTATGCG-1     pbmc3k        980          521
```

```r
VlnPlot(pbmc3k, features = c('nFeature_RNA',"nCount_RNA","percent.mt"), ncol = 3)
### Normally before VlnPlot, you should run "NormalzeData()" firstly...
```
![Fig_1](figs/Fig_1.png "violin plot of pbmc3k")

```r
plot1 <- FeatureScatter(pbmc3k, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc3k, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots
```
![Fig_2](figs/Fig_2.png)

```r
pbmc3k <- subset (pbmc3k, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
### 单细胞基因太少可能质量较差，测序深度或者细胞活性差；细胞基因数量过多可能是双细胞多细胞混合体，并没有形成单细胞悬液

ncol(pbmc3k)
### 2638
### Checking the results after filtering
```


## Step 3 Normalization
```r
pbmc3k <- NormalizeData(pbmc3k, normalization.method = "LogNormalize", scale.factor = 2000)
### Like TPM for normal RNA-seq, Eliminate the influence of sequencing depth (library size) on the 
### expression level between cells to make the expression levels of different cells comparable.
```

## Step 4 Find variable genes
```r
pbmc3k <- FindVariableFeatures(pbmc3k, selection.method = "vst", nfeatures = 2000)
top10_pbmc3k <- head(VariableFeatures(pbmc3k), 10)
```

```r
### Visualization variable genes
plot1 <- VariableFeaturePlot(pbmc3k)
plot2 <- LabelPoints(plot = plot1, points = top10_pbmc3k, repel = TRUE)
plot1 + plot2
```
![Fig_3](figs/Fig_3.png)

### Optional Step: Scaledata

```
### The Scaledata() is the standardization of each gene in all cells prevents errors caused by overly outlier expression of a certain gene during the subsequent cell clustering process.
pbmc3k <- ScaleData(pbmc3k, features = rownames(pbmc3k)) ### all
pbmc3k <- ScaleData(pbmc3k) ### 高变基因
```

## Step 5 PCA dimension reduction

```r
pbmc3k <- RunPCA(pbmc3k, features = VariableFeatures(object = pbmc3k))
print(pbmc3k[["pca"]],dims = 1:5, nfeatures = 5)
VizDimLoadings(pbmc3k, dims = 1:2, reduction = "pca")
DimPlot(pbmc3k, reduction = "pca")
DimHeatmap(pbmc3k, dims = 1 , cells = 500, balanced = TRUE)
DimHeatmap(pbmc3k, dims = 1:15, cells = 500, balanced = TRUE)
```

## Step 6 Perform "hypothesis tests" on each principal component of PCA to generate P-values
```r
pbmc3k <- JackStraw(pbmc3k, num.replicate = 100)
pbmc3k <- ScoreJackStraw(pbmc3k, dims = 1:20)
JackStrawPlot(pbmc3k, dims = 1:20)
ElbowPlot(pbmc3k)
### 发现维度10左右就稳定了
```

## Step 7 Find neighbors
```r
pbmc3k <- FindNeighbors(pbmc3k, dims = 1:10) 
pbmc3k <- FindClusters(pbmc3k, resolution = 0.5)
```r
### 寻找近邻和集群，resolution 0.4-1.2，越大细胞cluster越多。

pbmc3k <- RunUMAP(pbmc3k, dims = 1:10)  
DimPlot(pbmc3k, reduction = "umap")
### UMAP降维


pbmc3k <- RunTSNE(pbmc3k, dims = 1:10)
DimPlot(pbmc3k, reduction = "tsne")
### T-SNE降维,更复杂，分类更好，适合小数据

cluster5.markers <- FindMarkers(pbmc3k, ident.1 = 5, ident.2 = c(0,3), min.pct = 0.25)
### 计算所有细胞

pbmc3k.markers <- FindAllMarkers(pbmc3k, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
### 寻找所有markergene
pbmc3k.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
### 按照计算出的高表达的markergene传递给管道符号，选出表达量top n的基因
FeaturePlot(pbmc3k, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ"))

top10 <-pbmc3k.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(pbmc3k, features = top10$gene) + NoLegend()
### 用来参考细胞高表达基因的网络

new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono",
                     "NK", "DC","Platelet")
### id 重命名

names(new.cluster.ids) <- levels(pbmc3k)
pbmc3k <- RenameIdents(pbmc3k, new.cluster.ids)
DimPlot(pbmc3k, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

sessionInfo()
### 环境确认，可以用

### 用RDS二进制文件保存，节省空间
```
