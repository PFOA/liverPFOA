## dimplot of seurat object
Idents(data)<-"celltype"
pdf("2-celltypeDimplot.pdf")
p<-DimPlot(data,label=TRUE)+ theme(plot.title = element_text(vjust = 0.5,hjust = 0.5,size = 15)) + labs(title =  paste0("Cluster id"))
DimPlot(data,label = F,group.by = "type",cols = c('Control' = "#84BD00FF", 'Treatment' = "#CC0C00FF")) + theme(plot.title = element_text(vjust = 0.5,hjust = 0.5,size = 15)) + labs(title =  paste0("Dataset Type"))
DimPlot(data,label = F,group.by = "orig.ident")+ theme(plot.title = element_text(vjust = 0.5,hjust = 0.5,size = 15)) + labs(title =  paste0("Sample Id"))    #+ scale_color_nejm()
dev.off()

## Changes in cell composition

# logFC
count<-table(data$type,data$celltype)
percentage<-apply(count,2,function(x){x/rowSums(count)})
p<-percentage[c("Treatment","Control"),]
logfc<-percentage["Treatment",]/percentage["Control",]
mege<-rbind(p,abs(log2(logfc)))
rownames(mege)<- c("Treatment","Control","logfc(abs)")

# t-test
sampleRep<-table(data$celltype,data$sample)
celltype<-as.character(unique(data$celltype))
dat<-data.frame(celltype="",pvalue=0)
for(i  in celltype){
   pvalue<-t.test(sampleRep[i,c("PFOA252","PFOA253","PFOA254")],sampleRep[i,c("Ctr1","Ctr2","Ctr3")])$p.value
   pvalue<-round(pvalue,4)
   aveExr<-data.frame(celltype=i,pvalue=pvalue)
   dat=rbind(dat,aveExr)
  
}
dat<-dat[-1,]
write.csv(file="celltypePvalue.csv",dat)

# Visualization of cellular components 
data$orig.ident <-factor(data$sample,levels=c("Ctr1","Ctr2","Ctr3","PFOA252","PFOA253","PFOA254"))
dtype<-as.data.frame(table(data@meta.data[,c('celltype','type',"orig.ident")]))
pdf("2-Frac_celltype.pdf")
ggplot(dtype,aes(x=celltype,y=Freq,fill=type))+geom_bar(stat = "identity",position = 'fill')+theme_classic()+labs(x="",y="Fraction of different clusters")+coord_flip()+scale_fill_manual(values=c("Treatment" = pal_startrek()(4)[1], "Control" = pal_startrek()(4)[3]))
dev.off()
pdf("2-Frac_2_celltype.pdf")
ggplot(dtype,aes(x=celltype,y=Freq,fill=orig.ident))+geom_bar(stat = "identity",position = 'fill')+theme_classic()+labs(x="",y="Fraction of different clusters")+coord_flip()+scale_fill_manual(values=c("Ctr1" = pal_locuszoom()(7)[1], "Ctr2" = pal_locuszoom()(7)[2],"Ctr3" = pal_locuszoom()(7)[3],"PFOA252" = pal_locuszoom()(7)[4],"PFOA253" = pal_locuszoom()(7)[5], "PFOA254" = pal_locuszoom()(8)[6]))
dev.off()

# Visualization of Maker genes
data@meta.data$celltype=factor(data@meta.data$celltype,levels = c("Kc","cholangiocytes","HSC","DC","Monocytes","Macrophages","pDC","MEC","TNK","Neutrophils","Hepatocytes","LSEC","B cell","unknown" ))
gene<-c("C1qa","Epcam","Dcn","Xcr1","S100a4","Bcl11a","Pecam1","Trbc2","Csf3r","Pck1","Kdr","Cd79b")
gene<-c("Cd79b","Kdr","Pecam1","S100a4","C1qa","Nkg7","Pck1","Csf3r","Bcl11a","Xcr1","Cx3cr1","Dcn","Epcam")  # order gene
library(MySeuratWrappers)
pdf("2-MySeuratWrappersVlnPlot.pdf")
DefaultAssay(data)<-"RNA"
VlnPlot(data, features = gene,group.by="celltype",stacked=T,pt.size=0,features.face = "bold.italic")+theme(axis.ticks.y = element_blank(),axis.text.y = element_blank())+labs(x="",y="")
dev.off()

pdf("4-markers.pdf")
DoHeatmap(data, features = top10$gene) + NoLegend()+theme(text=element_text(size=5))
DotPlot(data,features=unique(top10$gene))+NoLegend()+theme(axis.text.x=element_text(angle=45,hjust=0.8,size=5))
dev.off()


## Differential gene analysis in certain celltype (Treatment VS Control)
data <- readRDS("/szrmyy/wangjgLab/scRNA/chenjh/P0004_mouseLiver/02.celltypeAnno/SCT_combine.sample40_0.8AnnoTNK.rds")
Celltype <- unique(data$celltype)

Idents(data)<-"celltype"
celltype = as.character(unique(data$celltype))
dat = data.frame(gene="",logFC=0,adj.P=0,Celltype="",State="")
for (i in celltype) {
 subdata=subset(data,idents=i)
 Idents(subdata)<-"type"
 marker<-FindMarkers(subdata,ident.1="Treatment",min.pct=0.25,logfc.threshold=0)
 submarker = data.frame(gene=row.names(marker),logFC=marker$avg_logFC,adj.P=marker$p_val_adj,Celltype=i,State=ifelse(marker$p_val_adj<0.05,ifelse(marker$avg_logFC>0.25,"Up",ifelse(marker$avg_logFC< -0.25,"Down","No")),"No"))
 dat=rbind(dat,submarker)
}
dat<-dat[-1,]
write.csv(dat,"Celltype_Sampletype_marker.csv",quote=F,row.names=F)

## Functional enrichment analysis of DEGs
u <-subset(deg, State =="Up")
d <-subset(deg, State =="Down")
x<- u$gene
if (type=="SYMBOL") {
  eg <- bitr(x, fromType="SYMBOL", toType=c("ENTREZID","ENSEMBL"), OrgDb="org.Mm.eg.db")
}else if(type=="ENSEMBL"){
  eg <- bitr(x, fromType="ENSEMBL", toType=c("ENTREZID","SYMBOL"), OrgDb="org.Mm.eg.db")
}
ego_MF <- enrichGO(OrgDb="org.Mm.eg.db", gene = eg[,2], pvalueCutoff = 0.05, ont = "MF", readable=TRUE)
ego_CC <- enrichGO(OrgDb="org.Mm.eg.db", gene = eg[,2], pvalueCutoff = 0.05, ont = "CC", readable=TRUE)
ego_BP <- enrichGO(OrgDb="org.Mm.eg.db", gene = eg[,2], pvalueCutoff = 0.05, ont = "BP", readable=TRUE)
library(clusterProfiler)
ids=bitr(markers$gene,'SYMBOL','ENTREZID','org.Mm.eg.db')
markers=merge(markers,ids,by.x='gene',by.y='SYMBOL')
gcSample=split(markers$ENTREZID, markers$cluster)
kegg <- compareCluster(gcSample, fun="enrichKEGG",organism="mmu",pvalueCutoff=0.05)
mf <- compareCluster(gcSample,fun="enrichGO", OrgDb="org.Mm.eg.db", ont= "BP")

# GO 不同亚群的功能，用bar图展现，注意不同亚群的GO term
library(RColorBrewer)
color <- brewer.pal(14,"Dark2")
color<-c(color,color)
colorl <- rep(color,each=12)
ggplot(top10) +
        aes(x = ID, y = -log10(p.adjust), fill = celltype) +
        geom_bar(stat = "identity",colour="black") +
        #scale_fill_hue() +
        scale_fill_manual(values =color)+
        theme(
                axis.title=element_text(size=15,face="plain",color="black"),
                axis.text = element_text(size=12,face="plain",color="black"),
                axis.text.x = element_text(angle = 90,colour = colorl,hjust=1,vjust=0.6),
                axis.title.x = element_blank(),
                legend.title = element_blank(),
                legend.text = element_text(size = 15, face = "bold"),
                legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
                legend.direction = "horizontal",
                legend.position = c(0.8,0.9),
                legend.background = element_blank(),
                panel.background = element_rect(fill = "transparent",colour = "black"),

                plot.background = element_blank()
        )


