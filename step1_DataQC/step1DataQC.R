## Data QC
# (0) 加载R包
library(dplyr)
library(Seurat)
library(patchwork)
#library(clustree)
library(ggsci)
library(ggplot2) 
options(future.globals.maxSize = 20000 * 1024^2)


# （1）读取数据
Ctr1 <- Read10X(data.dir = "/szrmyy/wangjgLab/scRNA/chenjh/P0004_mouseLiver/00.data/Ctr1/outs/filtered_feature_bc_matrix")
Ctr2 <-Read10X(data.dir="/szrmyy/wangjgLab/scRNA/chenjh/P0004_mouseLiver/00.data/Ctr2/outs/filtered_feature_bc_matrix")
Ctr3 <-Read10X(data.dir="/szrmyy/wangjgLab/scRNA/chenjh/P0004_mouseLiver/00.data/Ctr3/outs/filtered_feature_bc_matrix")

PFOA252 <- Read10X(data.dir = "/szrmyy/wangjgLab/scRNA/chenjh/P0004_mouseLiver/00.data/PFOA252/outs/filtered_feature_bc_matrix")
PFOA253 <- Read10X(data.dir = "/szrmyy/wangjgLab/scRNA/chenjh/P0004_mouseLiver/00.data/PFOA253/outs/filtered_feature_bc_matrix")
PFOA254 <- Read10X(data.dir = "/szrmyy/wangjgLab/scRNA/chenjh/P0004_mouseLiver/00.data/PFOA254/outs/filtered_feature_bc_matrix")

Ctr1 <- CreateSeuratObject(counts = Ctr1, project = "Ctr1", min.cells = 3, min.features = 200)
Ctr2 <- CreateSeuratObject(counts = Ctr2, project = "Ctr2", min.cells = 3, min.features = 200)
Ctr3 <- CreateSeuratObject(counts = Ctr3, project = "Ctr3", min.cells = 3, min.features = 200)

PFOA252 <- CreateSeuratObject(counts = PFOA252, project = "PFOA252", min.cells = 3, min.features = 200)
PFOA253 <- CreateSeuratObject(counts = PFOA253, project = "PFOA253", min.cells = 3, min.features = 200)
PFOA254 <- CreateSeuratObject(counts = PFOA254, project = "PFOA254", min.cells = 3, min.features = 200)

Ctr1 # 20198 features across 11708 samples
Ctr1$sample <- "Ctr1"
Ctr1$type <- "Control"

#Ctr1T # 19991 features across 11046 samples
Ctr2$sample <- "Ctr2"
Ctr2$type <- "Control"

Ctr3$sample <- "Ctr3"
Ctr3$type <- "Control"

PFOA252 # 20144 features across 11527 samples
PFOA252$sample <- "PFOA252"
PFOA252$type <- "Treatment"

PFOA253 # 19037 features across 11121 samples
PFOA253$sample <- "PFOA253"
PFOA253$type <- "Treatment"

PFOA254 # 18914 features across 9400 samples
PFOA254$sample <- "PFOA254"
PFOA254$type <- "Treatment"


# （2）统计线粒体相关基因，
Ctr1[["percent.mt"]] <- PercentageFeatureSet(Ctr1, pattern = "mt-") 
Ctr2[["percent.mt"]] <- PercentageFeatureSet(Ctr2, pattern = "mt-") 
Ctr3[["percent.mt"]] <- PercentageFeatureSet(Ctr3, pattern = "mt-")
PFOA252[["percent.mt"]] <- PercentageFeatureSet(PFOA252, pattern = "mt-") 
PFOA253[["percent.mt"]] <- PercentageFeatureSet(PFOA253, pattern = "mt-") 
PFOA254[["percent.mt"]] <- PercentageFeatureSet(PFOA254, pattern = "mt-") 
pdf("beforeFilter.pdf")
VlnPlot(Ctr1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) # 10068 samples
VlnPlot(Ctr2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) # 9326 samples
VlnPlot(Ctr3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) # 9326 samples
VlnPlot(PFOA252, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) # 9561 samples
VlnPlot(PFOA253, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) # 8281 samples
VlnPlot(PFOA254, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) # 6272 samples
dev.off()
# （3）提取子集：对检测基因数,reads数以及线粒体含量并进行过滤
Ctr1 <- subset(Ctr1, subset = nFeature_RNA > 200 & nFeature_RNA < 4500 & nCount_RNA > 500 & nCount_RNA < 18000 & percent.mt < 15)
Ctr2 <- subset(Ctr2, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & nCount_RNA > 500 & nCount_RNA < 18000 & percent.mt < 15)
Ctr3 <- subset(Ctr3, subset = nFeature_RNA > 200 & nFeature_RNA < 4500 & nCount_RNA > 500 & nCount_RNA < 18000 & percent.mt < 15)
PFOA252 <- subset(PFOA252, subset = nFeature_RNA > 200 & nFeature_RNA < 4800 & nCount_RNA > 500 & nCount_RNA < 20000 & percent.mt < 15)
PFOA253 <- subset(PFOA253, subset = nFeature_RNA > 200 & nFeature_RNA < 4500 & nCount_RNA > 500 & nCount_RNA < 20000 & percent.mt < 15)
PFOA254 <- subset(PFOA254, subset = nFeature_RNA > 200 & nFeature_RNA < 4500 & nCount_RNA > 500 & nCount_RNA < 20000 & percent.mt < 15)

# （3）绘制小提琴图：检测基因数、reads数以及线粒体含量
pdf("afterFilter.pdf")
VlnPlot(Ctr1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) # 10068 samples
VlnPlot(Ctr2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) # 9326 samples
VlnPlot(Ctr3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) # 9326 samples
VlnPlot(PFOA252, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) # 9561 samples
VlnPlot(PFOA253, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) # 8281 samples
VlnPlot(PFOA254, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) # 6272 samples
dev.off()

# （4） 运行 SCT 和 CCA 对多个样本进行数据整合
Ctr1 <- SCTransform(Ctr1, return.only.var.genes = FALSE,assay = "RNA", verbose = FALSE)
Ctr2 <- SCTransform(Ctr2, return.only.var.genes = FALSE,assay = "RNA", verbose = FALSE)
Ctr3 <- SCTransform(Ctr3, return.only.var.genes = FALSE,assay = "RNA", verbose = FALSE)
PFOA252 <- SCTransform(PFOA252, return.only.var.genes = FALSE,assay = "RNA", verbose = FALSE)
PFOA253 <- SCTransform(PFOA253, return.only.var.genes = FALSE,assay = "RNA", verbose = FALSE)
PFOA254 <- SCTransform(PFOA254, return.only.var.genes = FALSE,assay = "RNA", verbose = FALSE)

object_list = list(Ctr1,Ctr2,Ctr3,PFOA252,PFOA253,PFOA254)
selfeatures <- SelectIntegrationFeatures(object.list = object_list, nfeatures = 3000)
scc.list <- PrepSCTIntegration(object.list = object_list, anchor.features = selfeatures, verbose = FALSE)
scc.anchors <- FindIntegrationAnchors(object.list = scc.list, normalization.method = "SCT",anchor.features = selfeatures, verbose = FALSE)
scc_integrated <- IntegrateData(anchorset = scc.anchors, normalization.method = "SCT",verbose = FALSE)

remove(object_list,scc.anchors,scc.list)
scc_integrated # 45152 features across 52211 samples within 3 assays
head(scc_integrated@meta.data) #查看meta信息
scc_integrated@assays$RNA@data #查看归一化后数据

pdf("SCCVlnPlot.pdf")
VlnPlot(scc_integrated, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0)
dev.off()

saveRDS(scc_integrated,file="SCT_combine.sample.rds")
saveRDS(c(Ctr1,Ctr2,Ctr3,PFOA252,PFOA253,PFOA254),file="SCT_each.sample.rds")

##  cellType annotation
data<-readRDS("SCT_combine.sample40_0.8.rds")
pdf("3-clusters.pdf")
DimPlot(data,label=TRUE)+ theme(plot.title = element_text(vjust = 0.5,hjust = 0.5,size = 15)) + labs(title =  paste0("Cluster id"))
DimPlot(data,label = F,group.by = "type",cols = c('Control' = '#FFCC33', 'Treatment' = '#006633')) + theme(plot.title = element_text(vjust = 0.5,hjust = 0.5,size = 15)) + labs(title =  paste0("Dataset Type"))
DimPlot(data,label = F,group.by = "orig.ident")+ theme(plot.title = element_text(vjust = 0.5,hjust = 0.5,size = 15)) + labs(title =  paste0("Sample Id"))    #+ scale_color_nejm()
dev.off()

tt = as.character(data@meta.data$seurat_clusters)
tt[tt== "0" | tt =="25" | tt =="27" | tt == "29" | tt == "34"] = "B cell"
tt[tt== "30"] = "cholangiocytes"
tt[tt== "26"] = "DC"
tt[tt== "23"] = "LSEC"
tt[tt== "19"] = "Hepatocytes"
tt[tt== "28" | tt== "32"] = "HSC"
tt[tt== "8" | tt== "9" | tt == "12" | tt == "35"] = "Kc"
tt[tt== "1" | tt== "2" | tt== "3" | tt== "4" | tt== "5" | tt=="6" | tt =="7" | tt=="11"] = "LSEC"
tt[tt== "18" | tt=="24"] = "Macrophages"
tt[tt== "15" | tt== "16"] = "MEC"
tt[tt== "21"] = "Monocytes"
tt[tt== "10" | tt=="36"] = "Neutrophils"
tt[tt== "20" | tt== "22"] = "NK"
tt[tt== "17"] = "pDC"
tt[tt== "13" | tt== "31"] = "T cell"
tt[tt== "14"] = "T cell"
tt[tt== "33" | tt== "37"] = "unknown"

data@meta.data$celltype <-tt
saveRDS(data,file="SCT_combine.sample_pca36_0.8Anno.rds")
Idents(data)<-"celltype"
pdf("2-celltypeDimplot.pdf")
p<-DimPlot(data,label=TRUE)+ theme(plot.title = element_text(vjust = 0.5,hjust = 0.5,size = 15)) + labs(title =  paste0("Cluster id"))
DimPlot(data,label = F,group.by = "type",cols = c('Control' = "#84BD00FF", 'Treatment' = "#CC0C00FF")) + theme(plot.title = element_text(vjust = 0.5,hjust = 0.5,size = 15)) + labs(title =  paste0("Dataset Type"))
DimPlot(data,label = F,group.by = "orig.ident")+ theme(plot.title = element_text(vjust = 0.5,hjust = 0.5,size = 15)) + labs(title =  paste0("Sample Id"))    #+ scale_color_nejm()
dev.off()