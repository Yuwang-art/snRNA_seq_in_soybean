library(dplyr)
library(Seurat)
library(patchwork)
library(Matrix)
library(data.table)
library(ggplot2)
library(tidyverse)
library(parallel)
library(pheatmap)
#Reads the information of homologous genes into environment####
#Orthogroups.csv is produced by "orthofinder” (version 2.2.7, under the parameter: -S diamond)
orthoRaw <- as.data.frame(fread('Orthogroups.csv',header = T))
head(orthoRaw)
dim(orthoRaw)
colnames(orthoRaw) <- c("group","M",'S')
orthoOverlapMtSoybean <- orthoRaw[orthoRaw$group %in% orthoCount[rowSums(orthoCount[,2:3]>0)==2,]$group,]
#Sum the UMI counts for all homologous genes within each group to obtain a total UMI count####
orthoOverlapMtSoybean$S <- gsub('_','-',orthoOverlapMtSoybean$S)
orthoOverlapMtSoybean$MtRoot <- orthoOverlapMtSoybean$M
orthoOverlapMtSoybean$MtNodule <- orthoOverlapMtSoybean$M
#listUMI is the list including the matrices of soybean and Medicago
funUMIGeneGroup <- function(x){
  dat1 <- listUMI[[x]]
  dat2 <- orthoOverlapMtSoybean[,c("group",x)]
  funGroupSumUMI <- function(i){
    gene <- unlist(strsplit(as.character(dat2[i,2]),split = ", "))
    gene <- gene[gene %in% rownames(dat1)]
    if (length(gene)>0){
      df <- data.frame(stringsAsFactors = F)
      if (length(gene) ==1){
        df[1,1:ncol(dat1)] <- dat1[gene,]
      } else {
        df[1,1:ncol(dat1)] <- colSums(dat1[gene,])
      }
      colnames(df) <- colnames(dat1)
      rownames(df) <- dat2$group[i]
      return(df)
    }
  }
  out <- do.call('rbind',mclapply(1:nrow(dat2),function(i){funGroupSumUMI(i)},mc.cores = 70))
  return(out)
}
listUMIGeneGroup <-  mclapply(names(listUMI),function(x){funUMIGeneGroup(x)},mc.cores = length(names(listUMI)))
names(listUMIGeneGroup) <- names(listUMI)
#Integration####
all_integ <- merge(listUMIGeneGroup[[1]],listUMIGeneGroup[[2]],by="row.names")
all_integ <- merge(all_integ, listUMIGeneGroup[[3]],by.x="Row.names",by.y='row.names')
all_integ[1:5,1:5]
rownames(all_integ) <- all_integ$Row.names
all_integ <- all_integ[,-1]
dim(all_integ)
sam <- separate(data = data.frame(cell=colnames(all_integ)), col = "cell",
                into = c("barcode", "sample"), sep = "-")$sample
table(sam)
mydata <- CreateSeuratObject(counts = all_integ, project = "mydata_scRNAseq")
mydata@meta.data$sample <- sam
table(mydata@meta.data$sam)
plant.list <- SplitObject(mydata, split.by = "sample")
listSCTplant <- list()
for (x in names(plant.list)){
  ls <- plant.list[[x]]
  listSCTplant[[x]] <- SCTransform(ls, verbose = T, variable.features.n = 3000)
}
plant.features <- SelectIntegrationFeatures(object.list = listSCTplant,nfeatures = 3000)
listSCTplant <- PrepSCTIntegration(object.list = listSCTplant, anchor.features = plant.features,
                                   verbose = T)
reference_dataset <- which(names(listSCTplant) == "N.group")
noduleSpecific.anchors <- FindIntegrationAnchors(object.list = listSCTplant, normalization.method = "SCT",
                                       anchor.features = plant.features, reference = reference_dataset)
threeSample.integrated <- IntegrateData(anchorset = noduleSpecific.anchors, normalization.method = "SCT")
threeSample.integrated <- RunPCA(object = threeSample.integrated, verbose = FALSE, npcs = 30)
threeSample.integrated <- RunUMAP(object = threeSample.integrated, dims = 1:30)
DimPlot(threeSample.integrated, reduction = "umap",group.by = 'sample')
DefaultAssay(threeSample.integrated) <- "RNA"
threeSample.integrated <- NormalizeData(threeSample.integrated, normalization.method = "LogNormalize",
                                        scale.factor = 1e6)
#Integration of nodule-specific cells####
all_integ <- merge(listUMIGeneGroup[[2]],
                   listUMIGeneGroup[[3]][,cellID],#cellID is a vector including the cell barcodes of three nodule-specific cell clusters in N group
                   by="row.names")
all_integ[1:5,1:5]
rownames(all_integ) <- all_integ$Row.names
all_integ <- all_integ[,-1]
dim(all_integ)
sam <- separate(data = data.frame(cell=colnames(all_integ)), col = "cell",
                into = c("barcode", "sample"), sep = "-")$sample
mydata <- CreateSeuratObject(counts = all_integ, project = "mydata_scRNAseq")
mydata@meta.data$sample <- sam
table(mydata@meta.data$sam)
plant.list <- SplitObject(mydata, split.by = "sample")
listSCTplant <- list()
for (x in names(plant.list)){
  ls <- plant.list[[x]]
  listSCTplant[[x]] <- SCTransform(ls, verbose = T, variable.features.n = 3000)
}
plant.features <- SelectIntegrationFeatures(object.list = listSCTplant,nfeatures = 3000)
listSCTplant <- PrepSCTIntegration(object.list = listSCTplant, anchor.features = plant.features,
                                   verbose = T)
noduleSpecific.anchors <- FindIntegrationAnchors(object.list = listSCTplant, normalization.method = "SCT",
                                       anchor.features = plant.features)#, reference = reference_dataset
noduleSpecific.integrated <- IntegrateData(anchorset = noduleSpecific.anchors, normalization.method = "SCT")
noduleSpecific.integrated <- RunPCA(object = noduleSpecific.integrated, verbose = FALSE, npcs = 30)
noduleSpecific.integrated <- RunUMAP(object = noduleSpecific.integrated, dims = 1:30)
DimPlot(noduleSpecific.integrated, reduction = "umap")
noduleSpecific.integrated <- FindNeighbors(noduleSpecific.integrated, reduction = "pca", dims = 1:30)
noduleSpecific.integrated <- FindClusters(noduleSpecific.integrated, resolution = 0.2)
DimPlot(noduleSpecific.integrated, reduction = "umap",label = T)
DefaultAssay(noduleSpecific.integrated) <- "RNA"
noduleSpecific.integrated <- NormalizeData(noduleSpecific.integrated, normalization.method = "LogNormalize",
                                 scale.factor = 1e6)
DefaultAssay(noduleSpecific.integrated) <- "integrated"
