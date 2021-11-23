data<-readRDS("../HepchoHSCRemoveCellsAnno.RDS")

rm(list = ls())
library(SCENIC)
packageVersion("SCENIC")
library(SCopeLoomR)
scenicLoomPath='Samplef_loom_path_scenic_output.loom'
loom <- open_loom(scenicLoomPath)
regulons_incidMat <- get_regulons(loom, column.attr.name="Regulons")
regulons <- regulonsToGeneLists(regulons_incidMat)
regulonAUC <- get_regulons_AUC(loom, column.attr.name="RegulonsAUC")
pyScenicDir <- "pySCENIC_example/output"
library(SCENIC)
library(SCopeLoomR)
loom <- open_loom("/szrmyy/wangjgLab/scRNA/chenjh/P0004_mouseLiver/05.subgroupAnalysis/06.CelltypeSet/03.Hephsc/02.TF/Samplef_loom_path_scenic_output.loom", mode="r")

regulons_incidMat <- get_regulons(loom, column.attr.name='Regulons')
regulons <- regulonsToGeneLists(regulons_incidMat)
regulonsAUC <- get_regulons_AUC(loom, column.attr.name='RegulonsAUC')
regulonsAucThresholds <- get_regulon_thresholds(loom)
embeddings <- get_embeddings(loom)

exprMat <- get_dgem(loom)
cellInfo <- get_cell_annotation(loom)              #注释信息？获取不到
clusterings <- get_clusterings_with_name(loom)     # cluster 信息
close_loom(loom)

cellInfo<-read.csv(file="/szrmyy/wangjgLab/scRNA/chenjh/P0004_mouseLiver/05.subgroupAnalysis/06.CelltypeSet/03.Hephsc/02.TF/hepmeta.data.csv",check.name=TRUE,row.names=1)
regulonAUC <- regulonsAUC[onlyNonDuplicatedExtended(rownames(regulonsAUC)),]
regulonActivity_byCellType <- sapply(split(rownames(cellInfo), cellInfo$subType),
                                     function(cells) rowMeans(getAUC(regulonAUC)[,cells]))
regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))
pdf("6-Tf_AUCregulated.pdf")  # 热图展现AUC
pheatmap::pheatmap(regulonActivity_byCellType_Scaled,color=colorRampPalette(c("blue","white","red"))(100),
                   angle_col=45,cellwidth = 30,cellheight = 20,breaks=seq(-3, 3, length.out = 100),treeheight_row=10, treeheight_col=10, border_color=NA)
dev.off()