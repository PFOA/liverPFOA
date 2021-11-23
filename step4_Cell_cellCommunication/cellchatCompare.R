library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)
library(Seurat)
library(reticulate)
#use_python("/szrmyy/wangjgLab/scRNA/chenjh/platform/software/project/cellchat/bin/python")

NCcellchat<-readRDS("Control_cellchat.RDS")
BPDcellchat <-readRDS("Treatment_cellchat.RDS")


pdf("1-compareBPD_NC.pdf")
object.list <- list(Control = NCcellchat,Treatment = BPDcellchat)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2
dev.off()
pdf("1-compareBPD_NCweight.pdf")
par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")
dev.off()
pdf("1-compareBPD_NCheatmap.pdf")
gg1 <- netVisual_heatmap(cellchat)
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
gg1 + gg2
dev.off()
future::plan("multiprocess", workers = 60)

pdf("1-compareBPD_NCcirclesplit.pdf")
weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}
dev.off()
#====  main sources and targets
pdf("2-majorSourcesTargets.pdf",height=6,width=8)
num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
}
patchwork::wrap_plots(plots = gg)
dev.off()

# 鉴定cDC和fibro的差异的通路
#gg1 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Inflam. DC", signaling.exclude = "MIF")
#gg2 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "cDC1", signaling.exclude = c("MIF"))
#patchwork::wrap_plots(plots = list(gg1,gg2))

# Identify the conserved and context-specific signaling pathways
# 信号通路的相似性 功能
pdf("3-functionalsimilarity.pdf")
cellchat <- computeNetSimilarityPairwise(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat,umap.method = 'uwot', type = "functional")
cellchat <- netClustering(cellchat, type = "functional")
netVisual_embeddingPairwise(cellchat, type = "functional", label.size = 3.5)
dev.off()

pdf("3-structuresimilarity.pdf")
cellchat <- computeNetSimilarityPairwise(cellchat, type = "structural")
cellchat <- netEmbedding(cellchat, umap.method = 'uwot',type = "structural")
cellchat <- netClustering(cellchat, type = "structural")
netVisual_embeddingPairwise(cellchat, type = "structural", label.size = 3.5)
dev.off()
#netVisual_embeddingPairwiseZoomIn(cellchat, type = "structural", nCol = 2)
#  signaling networks with larger (or less) difference  Larger distance implies larger difference of the communication networks between two datasets
# in terms of either functional or structure

pdf("4-rankSimilarity.pdf")
rankSimilarity(cellchat, type = "functional")   #  problem? umap
dev.off()

pdf("4-netVisual_Weightheatmap.pdf",width=12,height=8)
gg1 <- netVisual_heatmap(cellchat,cluster.rows=T,cluster.cols=F)
gg2 <- netVisual_heatmap(cellchat, cluster.rows=T,cluster.cols=F,width=10,height=6,measure = "weight")
gg2
dev.off()


# conserved and context-specific signaling pathways. The top signaling pathways colored red are enriched in NL skin, and these colored green were enriched in the LS skin.
pdf("4-specificsignalingpathways.pdf",height=12,width=12)
gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, font.size =8, return.data =T,do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F,  font.size =8,do.stat = TRUE)
gg1 + gg2
dev.off()
# Compare outgoing (or incoming) signaling associated with each cell population
library(ComplexHeatmap)
i = 1
# combining all the identified signaling pathways from different datasets 
pdf("4-specificsignalingpathwaysoutgoing.pdf",width=12,height=13)
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, font.size =6,title = names(object.list)[i], width = 9, height = 16)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, font.size =6,title = names(object.list)[i+1], width = 9, height = 16)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.off()
pdf("4-specificsignalingpathwaysincoming.pdf",width=12,height=13)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "incoming", signaling = pathway.union, font.size=6,title = names(object.list)[i], width = 9, height = 16, color.heatmap = "GnBu")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "incoming", signaling = pathway.union, font.size=6,title = names(object.list)[i+1], width = 9, height = 16, color.heatmap = "GnBu")
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.off()
pdf("4-specificsignalingpathwaysAll.pdf",width=12,height=13)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "all", signaling = pathway.union, font.size=6,title = names(object.list)[i], width = 9, height = 16, color.heatmap = "OrRd")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "all", signaling = pathway.union, font.size=6,title = names(object.list)[i+1], width = 9, height = 16, color.heatmap = "OrRd")
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.off()

# Part III   Identify the upgulated and down-regulated signaling ligand-receptor pairs
#3.1 通信概率识别通信信号 signaling.LSIncreased = gg1$data.

pdf("5-signalingligandreceptorpairs_HSCtoOtherCelltype.pdf",width=8,height=12)
netVisual_bubble(cellchat, sources.use = c("HSC"), targets.use = c("Kc","cholangiocytes","Neutrophils","pDC","MEC","LSEC"),  comparison = c(1, 2), pairLR.use = pairLR.use,angle.x = 45)
dev.off()
#3.2 上调和下调
pdf("5-signalingligandreceptorpairs2_upDown_1.pdf",width=6,height=10)
gg1 <- netVisual_bubble(cellchat, sources.use = c("HSC"),targets.use =c("Kc","cholangiocytes","Neutrophils","pDC","MEC","LSEC"),comparison = c(1, 2), max.dataset = 2, font.size=6,title.name = "Increased signaling in Treatment", angle.x = 45, remove.isolate = T)
print(gg1)
dev.off()
pdf("5-signalingligandreceptorpairs2_upDown_2.pdf",width=6,height=10)
gg2 <- netVisual_bubble(cellchat, sources.use = c("HSC"),targets.use =c("Kc","cholangiocytes","Neutrophils","pDC","MEC","LSEC"),comparison = c(1, 2), max.dataset = 1, font.size=6,title.name = "Decreased signaling in Treatment", angle.x = 45, remove.isolate = T)
print(gg2)
dev.off()



#dysfunctional signaling by using differential expression analysis
pos.dataset = "Treatment"
features.name = pos.dataset
cellchat <- identifyOverExpressedGenes(cellchat, group.dataset = "datasets", pos.dataset = pos.dataset, features.name = features.name, only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.1, thresh.p = 1)
net <- netMappingDEG(cellchat, features.name = features.name)
net.up <- subsetCommunication(cellchat, net = net, datasets = "Treatment",ligand.logFC = 0.2, receptor.logFC = NULL)
net.down <- subsetCommunication(cellchat, net = net, datasets = "Control",ligand.logFC = -0.1, receptor.logFC = -0.1)
gene.up <- extractGeneSubsetFromPair(net.up, cellchat)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat)

# ==================  rankNet ===============================
Treatment<-c("TGFb","MK","SN","CD80","IL16","CD34","GALECTIN","VEGF","PROS","THBS","APP","COLLAGEN" 
,"OSM","SEMA4","BST2","CCL","APRIL","CD200","CD45")
pairLR.use <- extractEnrichedLR(cellchat, signaling = Treatment)
pdf("5-rankNet.pdf",width=8,height=12)
gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE)
gg1
dev.off()


pairLR.use.up = net.up[, "interaction_name", drop = F]
pdf("5-singnalingligandrecptorexpression.pdf")
gg1 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.up, sources.use = c(2,7,12), targets.use = c(1,3:6,8,9:11,13), comparison = c(1, 2),  angle.x = 45, remove.isolate = T,title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
pairLR.use.down = net.down[, "interaction_name", drop = F]
gg2 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.down, sources.use = c(2,7,12),targets.use = c(1,3:6,8,9:11,13), comparison = c(1, 2),  angle.x = 45, remove.isolate = T,title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
gg1 + gg2
dev.off()
pdf("5-singnalingligandChord.pdf")
par(mfrow = c(1,2), xpd=TRUE)
netVisual_chord_gene(object.list[[2]], sources.use = c(1:4), targets.use = c(5:11), slot.name = 'net', net = net.up, lab.cex = 1.2, small.gap = 3.8, title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
netVisual_chord_gene(object.list[[1]], sources.use = c(1:4), targets.use = c(5:11), slot.name = 'net', net = net.down, lab.cex = 1.2, small.gap = 3.8, title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
dev.off()

# Part IV Visually compare cell-cell communication using Hierarchy plot, Circle plot or Chord diagram
pdf("6-CXClaggregate.pdf")
pathways.show <- c("CXCL") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}
dev.off()
pdf("6-CXClheatmap.pdf")
pathways.show <- c("CXCL") 
par(mfrow = c(1,2), xpd=TRUE)
ht <- list()
for (i in 1:length(object.list)) {
  ht[[i]] <- netVisual_heatmap(object.list[[i]], signaling = pathways.show, color.heatmap = "Reds",title.name = paste(pathways.show, "signaling ",names(object.list)[i]))
}
ComplexHeatmap::draw(ht[[1]] + ht[[2]], ht_gap = unit(0.5, "cm"))
dev.off()
pdf("6-CXClaggregatechord.pdf")
pathways.show <- c("CXCL") 
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "chord", signaling.name = paste(pathways.show, names(object.list)[i]))
}
dev.off()
#Part V: Compare the signaling gene expression 
cellchat@meta$datasets = factor(cellchat@meta$datasets, levels = c("BPD", "NC")) # set factor level
pdf("7-CXCBPD_NC_expression.pdfL")
plotGeneExpression(cellchat, signaling = "CXCL", split.by = "datasets", colors.ggplot = T)
dev.off()
saveRDS(cellchat, file = "cellchat_comparisonAnalysis_Ratlung_BPD_vs_NC.rds")


# pathway analysis  
pathways.show <- c("BMP")
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
pdf(paste("6-",pathways.show,"pathwayCircle.pdf",sep=""),height=8, width=12)
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}
dev.off()

pdf(paste("6-",pathways.show,"pathwayHeatmap.pdf",sep=""),height=8, width=12)
par(mfrow = c(1,2), xpd=TRUE)
ht <- list()
for (i in 1:length(object.list)) {
  ht[[i]] <- netVisual_heatmap(object.list[[i]], signaling = pathways.show, color.heatmap = "Reds",title.name = paste(pathways.show, "signaling ",names(object.list)[i]))
}
ComplexHeatmap::draw(ht[[1]] + ht[[2]], ht_gap = unit(0.5, "cm"))
dev.off()
cellchat<-readRDS("cellchat_comparisonAnalysis_Ratlung_BPD_vs_NC.rds")
pdf(paste("6-",pathways.show,"pathwayLRDEG.pdf",sep=""),height=8, width=12)
cellchat@meta$datasets = factor(cellchat@meta$datasets, levels = c("BPD", "NC")) # set factor level
plotGeneExpression(cellchat, signaling = pathways.show, split.by = "datasets", colors.ggplot = T)
dev.off()

#调节role
pdf(paste("6-",pathways.show,"pathwayCelltypeRole.pdf",sep=""),height=8, width=12)
cellchat <- netAnalysis_computeCentrality(object.list[[1]], slot.name = "netP")
p1<-netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 9, height = 2.5, font.size = 6)
cellchat <- netAnalysis_computeCentrality(object.list[[2]], slot.name = "netP")
p2<-netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 9, height = 2.5, font.size = 6)
p1+p2
dev.off()
