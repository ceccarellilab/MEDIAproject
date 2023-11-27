---
title: "SingleCell data pre-processing and data analyses"
author: "Maria Claudia Costa"
date: "2023"
output:
  html_document:
    toc: true
    toc_depth: 2
    toc_float: true
    fig.width: 15
    fig.height: 7
---
```{r setup, include=FALSE}
knitr::opts_chunk$set(dpi=300)
```
## Loading libraries and set directories

```{r, message=FALSE}
library(dplyr)
library(apeglm)
library(org.Hs.eg.db)
library(Seurat)
library(sctransform)
library(SeuratObject)
library(DESeq2)
library(EnsDb.Hsapiens.v79)
library(enrichplot)
library(clusterProfiler)
library(AnnotationDbi)
library(SingleCellExperiment)
library(muscat)
library(ggplot2)
library(gprofiler2)
library(ggrepel)
library(msigdbr)


MEDIApath="/home/mclaudia/MEDIA Project/"
MEDIAgraph="/home/mclaudia/MEDIA Project/Codici_finali/Paper_figures/"
MEDIAsave="/home/mclaudia/MEDIA Project/Codici_finali/"
wd="/home/mclaudia/MEDIA Project/GSE152805_RAW_scell/"
```
# SingleCell Dataset from Chou et al., 2020 (doi:10.1038/s41598-020-67730-y): Dataset 1

## Build merged Seurat Object of the three patients: s113, s116, s118

### MT: medial region of the tibia (affected)
### OLT: outer lateral region of the tibia (healthy)

```{r, message=FALSE, warning=FALSE, eval=FALSE}

MT1.data <-Read10X(data.dir =paste0(wd,"MT_113"))
MT1<- CreateSeuratObject(counts = MT1.data, 
                         project = "MT1", min.cells = 300, min.features = 500)

MT2.data <-Read10X(data.dir = paste0(wd,"MT_116"))
MT2<- CreateSeuratObject(counts = MT2.data, 
                         project = "MT2", min.cells = 300, min.features = 500)

MT3.data <-Read10X(data.dir = paste0(wd,"MT_118"))
MT3<- CreateSeuratObject(counts = MT3.data, 
                         project = "MT3", min.cells = 300, min.features = 500)

olt1.data <-Read10X(data.dir = paste0(wd,"olt_113"))
olt1<- CreateSeuratObject(counts = olt1.data, 
                          project = "olt1", min.cells = 300, min.features = 500)

olt2.data <-Read10X(data.dir = paste0(wd,"olt_116"))
olt2<- CreateSeuratObject(counts = olt2.data, 
                          project = "olt2", min.cells = 300, min.features = 500)

olt3.data <-Read10X(data.dir = paste0(wd,"olt_118"))
olt3<- CreateSeuratObject(counts = olt3.data, 
                          project = "olt3", min.cells = 300, min.features = 500)

tot113 <- merge(MT1, y = olt1, add.cell.ids = c("MT_s113", "Olt_s113"), project = "tot113")
length(Cells(tot113))

tot116 <- merge(MT2, y = olt2, add.cell.ids = c("MT_s116", "Olt_s116"), project = "tot116")
length(Cells(tot116))

tot118 <- merge(MT3, y = olt3, add.cell.ids = c("MT_s118", "Olt_s118"), project = "tot118")
length(Cells(tot118))

pt_list <-list(tot113,tot116,tot118)
names(pt_list) <-c("tot113","tot116","tot118")
```

## Pre-processing

```{r, message=FALSE, eval=FALSE}
for( i in 1:length(pt_list)){
  
  counts <- GetAssayData(object = pt_list[[i]], slot = "counts")
  
  nonzero <- counts > 0
  keep_genes <- Matrix::rowSums(nonzero) >= 100
  filtered_counts <- counts[keep_genes, ]
  filtered_seurat <- CreateSeuratObject(filtered_counts, meta.data = pt_list[[i]]@meta.data)
  
  tmp <-filtered_seurat
  
  tmp[["percent.mt"]] <- PercentageFeatureSet(tmp , pattern = "^MT-")

  VlnPlot(tmp, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  
  plot1 <- FeatureScatter(tmp, feature1 = "nCount_RNA", feature2 = "percent.mt")
  plot2 <- FeatureScatter(tmp, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
 # plot1 + plot2, visualize to set nFeature_RNA, nCount_RNA and percent.mt thresholds
  
  if(i == 1) tmp <- subset(tmp, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 4 & nCount_RNA < 50000)
  if(i == 2) tmp <- subset(tmp, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 5 & nCount_RNA < 50000)
  if(i == 3) tmp <- subset(tmp, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 3 & nCount_RNA < 40000)
  
  tmp<- NormalizeData(tmp)
  tmp1<- FindVariableFeatures(tmp, selection.method = "vst", nfeatures = 8500)
  pt_list[[i]] <-tmp1
}
```

## Integrate patient data

```{r, message=FALSE,eval=FALSE}
features <- SelectIntegrationFeatures(pt_list, nfeatures = 8500)

anchors <- FindIntegrationAnchors(pt_list, anchor.features=features)

OAtot <-IntegrateData(anchors, new.assay.name = "integrate")
```

```{r, message=FALSE,echo=FALSE}
load(paste0(MEDIAsave,"OA_batch_anchors_nor.RData"))

```

## Examine and visualize PCA

```{r, message=FALSE}
all.genes <-  rownames(OAtot)

OAtot <- ScaleData(OAtot, features = all.genes)
OAtot <- RunPCA(OAtot, features = VariableFeatures(object = OAtot))
print(OAtot[["pca"]], dims = 1:5, nfeatures = 5)
```

## Run UMAP: 10 PC chosen

```{r, message=FALSE,fig.width=17,fig.height=10}

VizDimLoadings(OAtot, dims = 1:2, reduction = "pca")
DimPlot(OAtot, reduction = "pca")

DimHeatmap(OAtot, dims = 1:15, cells = 500, balanced = TRUE)

ElbowPlot(OAtot) ###9-12 PC

OAtot <- FindNeighbors(OAtot, dims = 1:10)
OAtot <- FindClusters(OAtot, resolution = 0.5) 

OAtot <- RunUMAP(OAtot, dims = 1:10)

```

### Save UMAP coordinates for reproducibility

```{r, message=FALSE, eval=FALSE}
umapCoord <- Embeddings(object = OAtot[["umap"]])
save(umapCoord,file=paste0(MEDIAsave,"CoordUMAP.RData"))
```

## Visualize UMAP for: clusterization, class, single patient
```{r,eval=FALSE}
load(paste0(MEDIAsave,"CoordUMAP.RData"))
```

```{r, message=FALSE,echo=FALSE}
load("~/MEDIA Project/CoordUMAP.RData")
```

```{r, message=FALSE}
OAtot@meta.data$type <-sapply(strsplit(rownames(OAtot@meta.data),"_"), function(x) x[1])
OAtot@meta.data$patient <-sapply(strsplit(rownames(OAtot@meta.data),"_"), function(x) x[2])

OAtot@reductions[["umap"]]@cell.embeddings  <- umapCoord

DimPlot(OAtot, reduction = "umap", label = TRUE, label.size = 7)+
   ggtitle("Clustering")+
   theme(plot.title = element_text(hjust = 0.5))

DimPlot(OAtot, reduction = "umap", group.by = "type", label = TRUE, 
        label.size = 7,cols= c("turquoise","orangered"))+
        ggtitle("Class")

DimPlot(OAtot, reduction = "umap", group.by = "patient")+
        ggtitle("Patient")
```

```{r,eval=FALSE,echo=FALSE}

knitr::opts_chunk$set(dpi=150)

pdf(file = paste0(MEDIAgraph,"UMAP_clusters.pdf"),  width = 15,height = 11)
DimPlot(OAtot, reduction = "umap", label = TRUE, label.size = 7,group.by = "seurat_clusters")
dev.off()

pdf(file = paste0(MEDIAgraph,"UMAP_type.pdf"),  width = 15,height = 11)
DimPlot(OAtot, reduction = "umap", 
        group.by = "type", label = TRUE, 
        label.size = 7,cols= c("turquoise","orangered"))
dev.off()

pdf(file = paste0(MEDIAgraph,"UMAP_patient.pdf"),  width = 15,height = 11)
DimPlot(OAtot, reduction = "umap", group.by = "patient")
dev.off()
```
## Pseudo-bulk singleCell data trasformation

```{r, message=FALSE}
ei <-OAtot@meta.data[,c(1,6,7,8)]

head(ei) 
ei$clus <-0  #merge in unique

OA_sce <-as.SingleCellExperiment(OAtot, assay = "RNA") 
OA_sce$clus <-0

(OA_sce2 <- prepSCE(OA_sce, 
                    kid = "seurat_clusters",
                    gid = "type",  
                    sid = "orig.ident",   
                    drop = FALSE))  

OA_sce2$cond <-paste0(OA_sce2$sample_id,"_",OA_sce2$patient)

experiment_info <- data.frame(sample_id=as.factor(c("MT1","MT2","MT3","olt1","olt2","olt3")),
                              cond=as.factor(c("MT1_s113 ","MT2_s116",
                                               "MT3_s118","olt1_s113","olt2_s116","olt3_s118")),
                              patient=as.factor(rep(c("s113","s116","s118"),2)),
                              group_id=as.factor(c("MT","MT","MT","Olt","Olt","Olt")),
                              n_cells=as.numeric(table(ei$orig.ident)))


OA_sce2@metadata <- list(experiment_info)
names(OA_sce2@metadata) <-"experiment_info"


table(OA_sce2$cluster_id,OA_sce2$sample_id)
```

```{r, message=FALSE}

pb <- aggregateData(OA_sce2, assay = "counts", fun = "sum")
colData(pb)
```

```{r, message=FALSE}

plotPCA(OA_sce2, colour_by = "group_id")

pbMDS(pb)
```

```{r, message=FALSE}
colData(pb)$patient <-as.factor(colData(pb)$patient)
colData(pb)$cond <-as.factor(colData(pb)$cond)

ei <- metadata(pb)$experiment_info

pb$group_id <-relevel(pb$group_id,ref="Olt")
ei$group_id <-relevel(ei$group_id,ref="Olt")


ass <-matrix(0, nrow = nrow(pb@int_elementMetadata), ncol = 6)
colnames(ass) <-c("MT1","MT2","MT3","olt1","olt2","olt3")

for(i in c(1:9)){
  as <-as.data.frame(as.matrix(pb@assays@data@listData[[i]]))
  ass <-ass+as
}

ass <- ass+1  #add dummy value = 1

cd <- colData(pb)
y <- DESeqDataSetFromMatrix(ass, cd, design=~patient+group_id)

dds <- DESeq(y)
resultsNames(dds)
res <-results(dds,alpha=.9,name="group_id_MT_vs_Olt",pAdjustMethod = "fdr")

res3 <-lfcShrink(dds, coef="group_id_MT_vs_Olt",res =res,type = "apeglm")
res3 <- as.data.frame(res3)
res3 <-res3[order(res3$padj, decreasing = FALSE),]

res3[1:10,]
```

#### Save files

```{r, message=FALSE,eval=FALSE}

save(res3,file=paste0(MEDIAsave,"DE_DESeq2_SCell.RData"))

```

## Volcano Plot of DE genes

```{r, message=FALSE, warning=FALSE}
 
annot <-res3[,c(2,5)]
annot$gene_names <-rownames(res3)

annot$col <-"gray50"
annot[annot$log2FoldChange <= -9 & annot$padj <=0.05, "col"] <- "springgreen3"
annot[annot$log2FoldChange  >= 9 & annot$padj <=0.05, "col"] <- "red"

annot$label <- NA
annot[annot$log2FoldChange  >= 10 & annot$padj <=0.05, "label"] <- annot[annot$log2FoldChange  >= 10 & annot$padj <=0.05, "gene_names"]
annot[annot$log2FoldChange  <= -10 & annot$padj <=0.05, "label"] <- annot[annot$log2FoldChange   <= -10 & annot$padj <=0.05, "gene_names"]

ggplot(data=annot, aes(x=log2FoldChange, y=-log10(padj), color=col,label=label)) +
  geom_point(size=1) +
  geom_label_repel(size=4,
                   min.segment.length = 0,direction = "both",
                   max.overlaps =40,
                   color="black",
                   segment.linetype=2) +
  scale_shape_manual(values=c(8, 19))+
  scale_color_manual(values = c("red","gray50","springgreen3"),
                     breaks = c("red","gray50","springgreen3"),
                     labels = c("UP","notDE","DOWN")) +
  labs(color="Legend")+
  xlab("log2(FC)")+
  ylab("-log10(FDR)")+
  theme_bw()+
  ggtitle("Volcano plot of DE genes - Dataset 1")+
  geom_vline(xintercept=c(-9, 9), col="gray70", linetype=2)+
  geom_hline(yintercept=-log10(0.05), col="gray70", linetype=2)
```

## Enrichment with GSEA: BP and REACTOME

```{r, message=FALSE,warning=FALSE}
hs=get("org.Hs.eg.db")

genes1 <- annot[, c("log2FoldChange","gene_names")]
genes <- genes1$log2FoldChange
names(genes) <-genes1$gene_names

genes <-sort(genes, decreasing = TRUE)

gsebp <- gseGO(geneList=genes, 
               ont ="BP", 
               keyType = "SYMBOL", 
               pvalueCutoff = 0.01,
               eps = 0,
               verbose = TRUE, 
               OrgDb = hs, 
               pAdjustMethod = "fdr")

dotplot(gsebp, showCategory = 20, x="NES",split=".sign") + facet_grid(.~.sign)

rc <-msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME")
sig <-rc[,c(3,4)]
sig$gs_name <-gsub("^REACTOME_","",sig$gs_name)

gsreac <- clusterProfiler::GSEA(genes,TERM2GENE = sig,
                             eps = 0,
                             verbose = TRUE, 
                             minGSSize =5, pAdjustMethod = "fdr",
                             pvalueCutoff =0.1)
dotplot(gsreac, showCategory = 10, x="NES",split=".sign") + facet_grid(.~.sign)+
  theme(axis.text.y = element_text(size = 8))
```
