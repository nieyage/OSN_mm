.libPaths("/data/R02/nieyg/ori/biosoft/conda/envs/r4-base/lib/R/library")

library(BSgenome)
library(Matrix)
library(Signac)
library(Seurat)
library(BSgenome.Mmusculus.UCSC.mm10)
library(tidyverse)
library(cowplot)
library(patchwork)
library(sctransform)
library(scater)
library(celda)
library(hdf5r)

OSN_mm<-readRDS("./02_All_celltype/WNN_OSN_integrated_all_celltype.rds")

Neutrophils<- subset(OSN_mm,idents="Neutrophils")

# ATAC analysis
# We exclude the first dimension as this is typically correlated with sequencing depth
DefaultAssay(Neutrophils) <- "ATAC"
Neutrophils <- RunUMAP(Neutrophils, reduction = 'integrated_LSI', dims = 2:40, reduction.name = "umap.atac", reduction.key = "atacUMAP_")

# tune the resolution
Neutrophils <- FindNeighbors(object = Neutrophils, reduction = 'integrated_LSI', dims = 2:40)

Neutrophils <- FindClusters(Neutrophils, resolution =1, algorithm = 3, verbose = FALSE)

my47colors <- c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
         '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
         '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398',
         '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
         '#712820', '#DCC1DD', '#CCE0F5', '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
         '#968175', '#161853', '#FF9999', '#344CB7', '#FFCC1D', '#116530', '#678983',
         '#678983', '#A19882', '#FFBCBC', '#24A19C', '#FF9A76')
pdf("./03_Neutrophils/Neutrophils_cluster_ATAC.pdf",width=15,height=5)
###cluster
p1 <- DimPlot(Neutrophils, cols=my47colors, reduction = "umap.atac", group.by = "seurat_clusters", label = TRUE, label.size = 3.5, repel = TRUE) + ggtitle("seurat_clusters")
p2 <- DimPlot(Neutrophils, cols=my47colors, reduction = "umap.atac", group.by = "orig.ident", shuffle=TRUE,label = F, label.size = 5, repel = TRUE) + ggtitle("sample")
p3 <- FeaturePlot(Neutrophils, reduction = 'umap.atac',features = "atac_S_M.pulmonis_ReadsPercentage", ncol = 1)
p1 + p2 +p3 
dev.off()


Idents(Neutrophils)<-Neutrophils$seurat_clusters
#####further annotation########
Neutrophils <- RenameIdents(
  object = Neutrophils,
  '0' = 'C1',
  '1' = 'C2',
  '2' = 'C2',
  '3' = 'C1',
  '4' = 'C1',
  '5' = 'C1',
  '6' = 'C1',
  '7' = 'C3',
  '8' = 'C3'
  )
Neutrophils@meta.data$subtype<-Idents(Neutrophils)
table(Neutrophils$subtype,Neutrophils$orig.ident)

Neutrophils$subtype<-factor(Neutrophils$subtype,
  levels=c("C1","C2","C3"))
Idents(Neutrophils)<-Neutrophils$subtype;

# Dimplot with annotation
pdf("./03_Neutrophils/Neutrophils_cluster_ATAC_subtype.pdf",width=10,height=8)
###cluster
p1 <- DimPlot(Neutrophils, cols=my47colors[-1:-2], reduction = "umap.atac", group.by = "seurat_clusters", label = TRUE, label.size = 3.5, repel = TRUE) + ggtitle("seurat_clusters")
p2 <- DimPlot(Neutrophils, cols=my47colors[-1:-2], reduction = "umap.atac", group.by = "subtype", label = TRUE, label.size = 3.5, repel = TRUE) + ggtitle("subtype")
pp1<- (p1 + p2 )
p3 <- DimPlot(Neutrophils, cols=my47colors, reduction = "umap.atac", group.by = "orig.ident", shuffle=TRUE,label = F, label.size = 5, repel = TRUE) + ggtitle("sample")
p4 <- FeaturePlot(Neutrophils, reduction = 'umap.atac',features = "atac_S_M.pulmonis_ReadsPercentage", ncol = 1)
pp2<- (p3 + p4 )
pp1/pp2
dev.off()


# cell proportion distribution 
pdf("./03_Neutrophils/Neutrophils_subtype_proportion.pdf",width=3,height=5)
df <- as.data.frame(Neutrophils@meta.data)
df_ct <- df %>%
    group_by(orig.ident, subtype) %>%
    summarise(counts = n()) %>%
    mutate(cell_proportion = counts / sum(counts))

p <- ggplot(df_ct, aes(orig.ident, cell_proportion, fill=subtype)) + 
    geom_bar(stat="identity", position = position_stack(reverse = TRUE)) +
    scale_fill_manual(values = my47colors[-1:-2]) +
    theme_cowplot() +
    xlab("") + ylab("") +
    theme(axis.text.x = element_text(angle = 60, hjust = 1))
p
dev.off()


pdf("./03_Neutrophils/Neutrophils_subtype_M.pulmonis_violinplot.pdf",width=4,height=8)

p1 <- VlnPlot(subset(Neutrophils,atac_S_M.pulmonis_ReadsPercentage>0), "atac_S_M.pulmonis_ReadsPercentage",add.noise=FALSE, pt.size = 0,col=my47colors[-1:-2],group.by = "subtype") +
geom_boxplot(width=0.1,outlier.size=0.01, position = position_dodge(0.9))

p2 <- VlnPlot(Neutrophils, "atac_S_M.pulmonis_ReadsPercentage", pt.size = 0.01,col=my47colors[-1:-2],group.by = "subtype") +
geom_boxplot(width=0.1,outlier.size=0.01, position = position_dodge(0.9))
p1/p2
dev.off()

# Find DEP for C1,2,3
# recall peak 
library(BSgenome.Mmusculus.UCSC.mm10)
# ORN recall peak for subcluster 
DefaultAssay(Neutrophils)<-"ATAC"
peak<-CallPeaks(
       Neutrophils,
       group.by = "subtype",
       macs2.path = "/public/home/nieyg/biosoft/conda/bin/macs2",
       broad = FALSE,
       format = "BED",
       fragment.tempdir = tempdir(),
       effective.genome.size = 1.87e+09, #mm10 size
       outdir="./03_Neutrophils/",
       combine.peaks=TRUE
)
macs2_counts <- FeatureMatrix(
     fragments = Fragments(Neutrophils),
     features = peak,
     cells = colnames(Neutrophils)
     )     
Neutrophils[["peaks_subtype"]] <- CreateChromatinAssay(
  counts = macs2_counts,
  #fragments = fragpath[[i]],
  fragments = Fragments(Neutrophils),
  annotation = Annotation(Neutrophils)
)

# DEG for each subtype 
DefaultAssay(Neutrophils)<-"ACTIVITY"
#Neutrophils<- SCTransform(Neutrophils)
Idents(Neutrophils)<- Neutrophils$subtype
Neutrophils <- ScaleData(Neutrophils,features=rownames(Neutrophils))

#DefaultAssay(Neutrophils)<-"SCT"
markers <- FindAllMarkers(Neutrophils, only.pos = TRUE, min.pct = 0.05, logfc.threshold = 0.1)
table(markers$cluster)
write.csv(markers,"./03_Neutrophils/Neutrophils_FindAllMarkers_DEG_ACTIVITY.csv")

# verification
myUmapcolors <- c(  '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
         '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
         '#9FA3A8', '#E0D4CA', '#5F3D69', '#58A4C3', '#AA9A59', '#E63863', '#E39A35', 
         '#C1E6F3', '#6778AE', '#B53E2B', '#712820', '#DCC1DD', '#CCE0F5', '#625D9E', 
         '#68A180', '#3A6963', '#968175', '#161853', '#FF9999', '#344CB7', '#FFCC1D', 
         '#116530', '#678983', '#A19882', '#FFBCBC', '#24A19C', '#FF9A76', "#8DD3C7",
         "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", 
         "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F", "#E41A1C", "#377EB8", "#4DAF4A", 
         "#FF7F00", "#FFFF33", "#A65628", "#F781BF")
signif_markers <- markers[markers$p_val_adj<0.05,] 
top10 <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC);
blueYellow = c("1"="#352A86","2"="#343DAE","3"="#0262E0","4"="#1389D2","5"="#2DB7A3","6"="#A5BE6A","7"="#F8BA43","8"="#F6DA23","9"="#F8FA0D")
solarExtra<- c("#3361A5", "#248AF3", "#14B3FF" ,"#88CEEF" ,"#C1D5DC", "#EAD397" ,"#FDB31A" ,"#E42A2A","#A31D1D")

DefaultAssay(Neutrophils)<- "ACTIVITY"
pdf("./03_Neutrophils/Neutrophils_DEG_heatmap_ACTIVITY.pdf",width=7,height=7)
DoHeatmap(object = Neutrophils,disp.min = -0.7,disp.max = 0.7,features=top10$gene,group.colors =myUmapcolors,size = 2,group.by = "subtype") +scale_fill_gradientn(colors = solarExtra[3:8])
DoHeatmap(object = Neutrophils,disp.min = -0.7,disp.max = 0.7,features=markers$gene, group.colors =myUmapcolors,size = 2,group.by = "subtype") +scale_fill_gradientn(colors = solarExtra[3:8])
dev.off();

saveRDS(Neutrophils,"./03_Neutrophils/Neutrophils_subtype.rds")

.libPaths("/data/R02/nieyg/ori/biosoft/conda/envs/R/lib/R/library")
library(BSgenome)
library(Matrix)
library(Signac)
library(Seurat)
library(BSgenome.Mmusculus.UCSC.mm10)
library(tidyverse)
library(cowplot)
library(patchwork)
library(sctransform)
library(scater)

Neutrophils<- readRDS("./03_Neutrophils/Neutrophils_subtype.rds")


markers<- read.csv("./03_Neutrophils/Neutrophils_FindAllMarkers_DEG_ACTIVITY.csv")
# GO and KEGG for CM subtype
C1<-markers[markers$cluster=="C1",8]
C2<-markers[markers$cluster=="C2",8]
C3<-markers[markers$cluster=="C3",8]

library(pheatmap)
library(clusterProfiler)
library(org.Mm.eg.db)

CAV_list<- list(C1,C2,C3)
subtype<- c("C1","C2","C3")
all_ego<- data.frame()
pdf("./03_Neutrophils/Neutrophils_subtype_DEG_GO_BP_ACTIVITY.pdf")
for(i in 1:length(subtype)){
  gene<- CAV_list[[i]];
  gene.df <- bitr(gene, fromType = "SYMBOL",
                toType = c("ENSEMBL", "ENTREZID"),
                OrgDb = org.Mm.eg.db)
    ego <- enrichGO(gene.df$ENTREZID,
                keyType = 'ENTREZID',
                OrgDb = org.Mm.eg.db,
                ont = "BP", ###BP,MF,CC
                pAdjustMethod  = "BH",
                pvalueCutoff = 0.5,
                qvalueCutoff = 0.5,
                readable = TRUE)
    p<- barplot(ego, showCategory=20,,label_format=50)
    print(p)
    write.csv(ego,paste0("./03_Neutrophils/",subtype[i],"-GO-BP.csv",sep=""))
    ego<- as.data.frame(ego)
    ego$celltype<- subtype[i];
    all_ego<- rbind(all_ego,ego);
}
dev.off()

write.csv(all_ego,"./03_Neutrophils/Neutrophils_Allsubtype_GO.csv")


.libPaths("/data/R02/nieyg/ori/biosoft/conda/envs/r43/lib/R/library")
DefaultAssay(Neutrophils)<-"peaks_subtype"
da_peaks <- FindAllMarkers(
  object = Neutrophils,
  test.use = 'LR',
  logfc.threshold = 0.1,
  min.pct = 0.1,
  latent.vars = 'nCount_peaks'
)
write.csv(da_peaks,"./03_Neutrophils/Neutrophils_FindAllMarkers_DEP.csv")
table(da_peaks$cluster)

da_peaks2<- da_peaks[which(abs(da_peaks$avg_log2FC)>0.5),]
table(da_peaks2$cluster)

#verification

C1_top.da.peak <- rownames(da_peaks[da_peaks$cluster=="C1" & da_peaks$p_val < 0.005 & da_peaks$avg_log2FC > 0.5, ])
C2_top.da.peak <- rownames(da_peaks[da_peaks$cluster=="C2" & da_peaks$p_val < 0.005 & da_peaks$avg_log2FC > 0.5, ])
C3_top.da.peak <- rownames(da_peaks[da_peaks$cluster=="C3" & da_peaks$p_val < 0.005 & da_peaks$avg_log2FC > 0.5, ])


peak2show<- c(C1_top.da.peak,C2_top.da.peak,C3_top.da.peak)
Neutrophils<- ScaleData(Neutrophils,features=peak2show)
pdf("./03_Neutrophils/Neutrophils_DEP_heatmap.pdf",width=10,height=10)
DoHeatmap(object = Neutrophils,disp.min = -0.5,disp.max = 0.4,features=peak2show, 
  group.colors =myUmapcolors,size = 4,group.by = "subtype") +scale_fill_gradientn(colors = blueYellow[3:8])
#DoHeatmap(object = Neutrophils,features=peak2show, 
#  group.colors =myUmapcolors,size = 4,group.by = "subtype") #+scale_fill_gradientn(colors = blueYellow[3:8])

dev.off();

# get top differentially accessible peaks
write.csv(da_peaks,"./03_Neutrophils/Neutrophils_FindAllMarkers_DEP.csv")

da_peaks<- read.csv("./03_Neutrophils/Neutrophils_FindAllMarkers_DEP.csv")
C1_top.da.peak <- da_peaks[da_peaks$cluster=="C1" & da_peaks$p_val < 0.005 & da_peaks$avg_log2FC > 0.5, 1]
C2_top.da.peak <- da_peaks[da_peaks$cluster=="C2" & da_peaks$p_val < 0.005 & da_peaks$avg_log2FC > 0.5, 1]
C3_top.da.peak <- da_peaks[da_peaks$cluster=="C3" & da_peaks$p_val < 0.005 & da_peaks$avg_log2FC > 0.5, 1]


library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Mmusculus.UCSC.mm10)
library(patchwork)
# Get a list of motif position frequency matrices from the JASPAR database
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)
DefaultAssay(Neutrophils)<- "peaks_subtype"
# add motif information
Neutrophils <- AddMotifs(
  object = Neutrophils,
  genome = BSgenome.Mmusculus.UCSC.mm10,
  pfm = pfm
)
C1_enriched.motifs <- FindMotifs(object = Neutrophils,features = C1_top.da.peak)
C2_enriched.motifs <- FindMotifs(object = Neutrophils,features = C2_top.da.peak)
C3_enriched.motifs <- FindMotifs(object = Neutrophils,features = C3_top.da.peak)

C1_enriched.motifs$rank <- 1:nrow(C1_enriched.motifs)
C2_enriched.motifs$rank <- 1:nrow(C2_enriched.motifs)
C3_enriched.motifs$rank <- 1:nrow(C3_enriched.motifs)
library(ggrepel)
convert2_mouse_symbol <- function(symbol){
  first_char <- toupper(substr(symbol,1,1))
  other_chars <- tolower(substr(symbol,2,nchar(symbol)))
  coverted_symbol <- str_c(first_char,other_chars)
  return(coverted_symbol)
}
pdf("./03_Neutrophils/Neutrophils_DEP_C123_Top_MotifPlot.pdf",width=6,height=3)
MotifPlot(object = Neutrophils,motifs = head(rownames(C1_enriched.motifs)))
MotifPlot(object = Neutrophils,motifs = head(rownames(C2_enriched.motifs)))
MotifPlot(object = Neutrophils,motifs = head(rownames(C3_enriched.motifs)))
dev.off()
out_dir<- "/data/R02/nieyg/project/OSN_mm/joint/03_Neutrophils/"
## Figure ranking plot of enriched motifs
C1_enriched_motifs_df <- C1_enriched.motifs
C1_enriched_motifs_df$motif.name <- gsub("\\(var\\.\\d\\)","",C1_enriched_motifs_df$motif.name)
C1_enriched_motifs_df$motif.name <- gsub("\\_02","",C1_enriched_motifs_df$motif.name)
C1_enriched_motifs_df <- C1_enriched_motifs_df[-grep("\\:\\:",C1_enriched_motifs_df$motif.name),]

C1_enriched_motifs_df$mouse_symbol <- sapply(1:nrow(C1_enriched_motifs_df),function(i){
  mouse_symbol <- convert2_mouse_symbol(C1_enriched_motifs_df$motif.name[i])
  mouse_symbol
  })
C1_enriched_motifs_df$label <- ifelse(C1_enriched_motifs_df$mouse_symbol %in% c(C1_enriched_motifs_df$mouse_symbol[1:10],"Lhx2"),C1_enriched_motifs_df$mouse_symbol,"")
C1_enriched_motifs_df$rank <- 1:nrow(C1_enriched_motifs_df)
C1_enriched_motifs_df$mlog10pvalue <- -log10(C1_enriched_motifs_df$pvalue)
pdf(str_c(out_dir,"C1_peaks_enriched_motifs_ranking_plot.pdf"),height=7.5)
ggplot(data=C1_enriched_motifs_df[1:50,],aes(x=rank,y=mlog10pvalue))+
  geom_point(aes(color=fold.enrichment,size=percent.observed))+
  labs(x="Motif rank",y="-log10(p_value)",color="Enrichment",size="Percentage of peaks with corresponding motif",title="Enriched motifs of peaks in C1")+
  geom_text_repel(aes(label=label),size=5,point.padding = 0.2,
    nudge_x = .15,
    nudge_y = .5,
    segment.curvature = -1e-20,
    arrow = arrow(length = unit(0.015, "npc")))+
  scale_colour_gradient(low = "yellow", high = "red")+
  theme_classic()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title = element_text(hjust = 0.5,size=rel(1.6)),axis.title=element_text(size=rel(1.6)),axis.text.y = element_text(color="black",size=rel(1.5)), axis.text.x = element_text(color="black",size=rel(1.5)),axis.line = element_line(colour="black",size = 1),legend.position="bottom")
dev.off()
write.table(C1_enriched_motifs_df,str_c(out_dir,"C1_peaks_enriched_motifs.txt"),sep="\t",quote=FALSE,col.names=TRUE,row.names=FALSE)

C2_enriched_motifs_df <- C2_enriched.motifs
C2_enriched_motifs_df$motif.name <- gsub("\\(var\\.\\d\\)","",C2_enriched_motifs_df$motif.name)
C2_enriched_motifs_df$motif.name <- gsub("\\_02","",C2_enriched_motifs_df$motif.name)
C2_enriched_motifs_df <- C2_enriched_motifs_df[-grep("\\:\\:",C2_enriched_motifs_df$motif.name),]

C2_enriched_motifs_df$mouse_symbol <- sapply(1:nrow(C2_enriched_motifs_df),function(i){
  mouse_symbol <- convert2_mouse_symbol(C2_enriched_motifs_df$motif.name[i])
  mouse_symbol
  })
C2_enriched_motifs_df$label <- ifelse(C2_enriched_motifs_df$mouse_symbol %in% c(C2_enriched_motifs_df$mouse_symbol[1:10],"Lhx2"),C2_enriched_motifs_df$mouse_symbol,"")
C2_enriched_motifs_df$rank <- 1:nrow(C2_enriched_motifs_df)
C2_enriched_motifs_df$mlog10pvalue <- -log10(C2_enriched_motifs_df$pvalue)
pdf(str_c(out_dir,"C2_peaks_enriched_motifs_ranking_plot.pdf"),height=7.5)
ggplot(data=C2_enriched_motifs_df[1:50,],aes(x=rank,y=mlog10pvalue))+
  geom_point(aes(color=fold.enrichment,size=percent.observed))+
  labs(x="Motif rank",y="-log10(p_value)",color="Enrichment",size="Percentage of peaks with corresponding motif",title="Enriched motifs of peaks in C2")+
  geom_text_repel(aes(label=label),size=5,point.padding = 0.2,
    nudge_x = .15,
    nudge_y = .5,
    segment.curvature = -1e-20,
    arrow = arrow(length = unit(0.015, "npc")))+
  scale_colour_gradient(low = "yellow", high = "red")+
  theme_classic()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title = element_text(hjust = 0.5,size=rel(1.6)),axis.title=element_text(size=rel(1.6)),axis.text.y = element_text(color="black",size=rel(1.5)), axis.text.x = element_text(color="black",size=rel(1.5)),axis.line = element_line(colour="black",size = 1),legend.position="bottom")
dev.off()
write.table(C2_enriched_motifs_df,str_c(out_dir,"C2_peaks_enriched_motifs.txt"),sep="\t",quote=FALSE,col.names=TRUE,row.names=FALSE)


C3_enriched_motifs_df <- C3_enriched.motifs
C3_enriched_motifs_df$motif.name <- gsub("\\(var\\.\\d\\)","",C3_enriched_motifs_df$motif.name)
C3_enriched_motifs_df$motif.name <- gsub("\\_02","",C3_enriched_motifs_df$motif.name)
C3_enriched_motifs_df <- C3_enriched_motifs_df[-grep("\\:\\:",C3_enriched_motifs_df$motif.name),]

C3_enriched_motifs_df$mouse_symbol <- sapply(1:nrow(C3_enriched_motifs_df),function(i){
  mouse_symbol <- convert2_mouse_symbol(C3_enriched_motifs_df$motif.name[i])
  mouse_symbol
  })
C3_enriched_motifs_df$label <- ifelse(C3_enriched_motifs_df$mouse_symbol %in% c(C3_enriched_motifs_df$mouse_symbol[1:10],"Lhx2"),C3_enriched_motifs_df$mouse_symbol,"")
C3_enriched_motifs_df$rank <- 1:nrow(C3_enriched_motifs_df)
C3_enriched_motifs_df$mlog10pvalue <- -log10(C3_enriched_motifs_df$pvalue)
pdf(str_c(out_dir,"C3_peaks_enriched_motifs_ranking_plot.pdf"),height=7.5)
ggplot(data=C3_enriched_motifs_df[1:50,],aes(x=rank,y=mlog10pvalue))+
  geom_point(aes(color=fold.enrichment,size=percent.observed))+
  labs(x="Motif rank",y="-log10(p_value)",color="Enrichment",size="Percentage of peaks with corresponding motif",title="Enriched motifs of peaks in C3")+
  geom_text_repel(aes(label=label),size=5,point.padding = 0.2,
    nudge_x = .15,
    nudge_y = .5,
    segment.curvature = -1e-20,
    arrow = arrow(length = unit(0.015, "npc")))+
  scale_colour_gradient(low = "yellow", high = "red")+
  theme_classic()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title = element_text(hjust = 0.5,size=rel(1.6)),axis.title=element_text(size=rel(1.6)),axis.text.y = element_text(color="black",size=rel(1.5)), axis.text.x = element_text(color="black",size=rel(1.5)),axis.line = element_line(colour="black",size = 1),legend.position="bottom")
dev.off()
write.table(C3_enriched_motifs_df,str_c(out_dir,"C3_peaks_enriched_motifs.txt"),sep="\t",quote=FALSE,col.names=TRUE,row.names=FALSE)


## chromVAR
# Compute a per-cell motif activity score by running chromVAR
DefaultAssay(Neutrophils) <- "ATAC"
Neutrophils <- RunChromVAR(
  object = Neutrophils,
  genome = BSgenome.Mmusculus.UCSC.mm10
)

# gene expression of Top TF 
motifs <- head(C2_enriched_motifs_df$motif)
motif_names <- head(C2_enriched_motifs_df$mouse_symbol)


library(RColorBrewer)
gene_exp_palettes <- c(brewer.pal(8,"Pastel2")[8],brewer.pal(9,"Greens")[3:8])
DefaultAssay(Neutrophils) <- "ACTIVITY"
gene_exp_UMAP_plot.ls <- lapply(1:length(motif_names),function(i){
  gene_exp_df <- FetchData(Neutrophils,vars = c(motif_names[i],"atacUMAP_1","atacUMAP_2"),slot = "data")
  colnames(gene_exp_df)[1] <- "gene"
  p <- ggplot(data=gene_exp_df,aes(x=atacUMAP_1,y=atacUMAP_2,color=gene))+
    geom_point(size=0.1)+
    scale_colour_gradientn(colours = gene_exp_palettes,limits=c(0,quantile(gene_exp_df$gene,seq(0,1,0.01))[96]),oob = scales::squish,breaks=c(0,quantile(gene_exp_df$gene,seq(0,1,0.01))[96]),labels = c("0", round(quantile(gene_exp_df$gene,seq(0,1,0.01))[96],1)),guide = guide_colorbar(frame.colour = "black",frame.linewidth= 1,title.position="top", title.hjust = 0.5)) + 
    scale_x_continuous(limits=c(-12.5,7.5),expand = c(0,0))+
    scale_y_continuous(limits=c(-7.5,8),expand = c(0,0))+
    labs(title=motif_names[i],color="Gene ACTIVITY",x="UMAP 1",y="UMAP 2") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background=element_rect(fill='transparent', color='black',linetype="solid",size=3),plot.title = element_text(hjust = 0.5,size=rel(3),face="bold"),axis.text=element_blank(),axis.ticks=element_blank(),axis.title.x=element_text(size=rel(1.5),vjust=0.5),axis.title.y=element_text(size=rel(1.5),vjust=0.5),plot.margin=unit(c(0.3, 0.2, 0.7, 0.5), "cm"),legend.title = element_text(size=15),legend.text=element_text(size=10),legend.position=c(0.2,0.15),legend.direction="horizontal",legend.key.size = unit(0.8, 'cm'))
    p
  })

violin_palettes <- c(brewer.pal(9,"Greens")[7],brewer.pal(8,"Set2")[8])
#violin_palettes <- colorspace::lighten(violin_palettes,amount = 0.3)
fmt_dcimals <- function(decimals=0){
    function(x) format(x,nsmall = decimals,scientific = FALSE)
  }
library(ggpubr)

Neutrophils$S_M.pulmonis <- ifelse(Neutrophils$subtype=="C2","S_M.pulmonis","other")
Neutrophils$S_M.pulmonis <- factor(Neutrophils$S_M.pulmonis,levels=c("S_M.pulmonis","other"))
Idents(Neutrophils)<- Neutrophils$S_M.pulmonis

markers<- motif_names
Cd36_vs_other_exp_df <- FindMarkers(Neutrophils,features=markers,ident.1="S_M.pulmonis",logfc.threshold=-Inf,min.pct=-Inf)
gene_exp_violin_plot.ls <- lapply(1:length(markers),function(i){
  gene_exp_df <- FetchData(Neutrophils,vars = c(markers[i],"atacUMAP_1","atacUMAP_2","S_M.pulmonis"),slot = "data")
  colnames(gene_exp_df)[1] <- "gene" 
  gene_exp_df$S_M.pulmonis <- factor(gene_exp_df$S_M.pulmonis,levels=c("S_M.pulmonis","other"),labels=c("S_M.pulmonis","other"))
  Cd36_vs_other_p_val_adj <- signif(Cd36_vs_other_exp_df$p_val_adj[which(rownames(Cd36_vs_other_exp_df)==markers[i])],3)
  stat_df <- tibble::as_tibble(data.frame(group1="S_M.pulmonis",group2="other",p.adj=Cd36_vs_other_p_val_adj,y.position=round(max(gene_exp_df$gene),1)+0.7))
  stat_df$p.signif <- ifelse(stat_df$p.adj<0.05,stat_df$p.adj,"ns")
  p <- ggplot() + 
    geom_violin(data=gene_exp_df,aes(x=S_M.pulmonis,y=gene,fill=S_M.pulmonis),scale="width",trim=FALSE)+
    geom_boxplot(data=gene_exp_df,aes(x=S_M.pulmonis,y=gene,fill=S_M.pulmonis),width=0.05, fill="white",outlier.size=0.8)+
    scale_fill_manual(values=violin_palettes)+
    labs(y="Gene expression",title=markers[i]) + 
    scale_y_continuous(labels = fmt_dcimals(0),limits=c(round(min(gene_exp_df$gene),1),round(max(gene_exp_df$gene),1)+1))+
    stat_pvalue_manual(stat_df, label = "p = {p.adj}",label.size=5) + 
    theme_classic()+
    theme(
          plot.title = element_text(hjust = 0.5,size=rel(3),face="bold"),
          axis.title.x = element_blank(),
          axis.title.y = element_text(color="black",size=rel(1.8)),
          axis.text.y = element_text(color="black",size=rel(1.5)), 
          axis.text.x = element_text(size=rel(1.5),color="black",vjust=1),
          axis.line = element_line(colour="black",size=1.5),
          axis.ticks = element_line(),
          plot.margin=unit(c(0, 0.2, 0.7, 0.2), "cm"),
          legend.position = " none")
  p
  })

pdf(str_c(out_dir,"S_M.pulmonis_Top_Motif_gene_expression_UMAP_VlnPlot_1.pdf"),width=12,height=21)
(gene_exp_UMAP_plot.ls[[1]] + gene_exp_violin_plot.ls[[1]] + plot_layout(widths=c(5.5,5))) / 
(gene_exp_UMAP_plot.ls[[2]] + gene_exp_violin_plot.ls[[2]]+plot_layout(widths=c(5.5,5))) / 
(gene_exp_UMAP_plot.ls[[3]] + gene_exp_violin_plot.ls[[3]] + plot_layout(widths=c(5.5,5))) + plot_annotation(title="Gene ACTIVITY of TFs",theme = theme(plot.title = element_text(hjust = 0.5,size = rel(2.8))))
dev.off()
pdf(str_c(out_dir,"S_M.pulmonis_Top_Motif_gene_expression_UMAP_VlnPlot_2.pdf"),width=12,height=21)
(gene_exp_UMAP_plot.ls[[4]] + gene_exp_violin_plot.ls[[4]] + plot_layout(widths=c(5.5,5))) / 
(gene_exp_UMAP_plot.ls[[5]] + gene_exp_violin_plot.ls[[5]]+plot_layout(widths=c(5.5,5))) / 
(gene_exp_UMAP_plot.ls[[6]] + gene_exp_violin_plot.ls[[6]] + plot_layout(widths=c(5.5,5))) + plot_annotation(title="Gene ACTIVITY of TFs",theme = theme(plot.title = element_text(hjust = 0.5,size = rel(2.8))))
dev.off()

## Figure 3F,H,J - ChromVar Z-score for Top TF (UMAP + Violin plot)
library(ArchR)
chromvar_UMAP_plot.ls <- lapply(1:length(motifs),function(i){
  chromvar_df <- as.data.frame(Neutrophils@assays$chromvar@data) %>%
    tibble::rownames_to_column("motif")  %>%
    dplyr::filter(motif %in% motifs[i])
  chromvar_df <- reshape2::melt(chromvar_df,id="motif",variable.name="cell",value.name ="chromvar_zscore")
  meta_df <- FetchData(Neutrophils,vars=c("atacUMAP_1","atacUMAP_2"))
  chromvar_df <- merge(chromvar_df,meta_df,by.x="cell",by.y="row.names")
  p <- ggplot(data=chromvar_df,aes(x=atacUMAP_1,y=atacUMAP_2,color=chromvar_zscore))+
    geom_point(size=0.4)+
    scale_colour_gradientn(colours = as.vector(ArchRPalettes[["solarExtra"]]),limits=c(-1.5,1.5),oob = scales::squish,breaks=c(-1.5,0,1.5),guide = guide_colorbar(frame.colour = "black",frame.linewidth= 1,title.position="top", title.hjust = 0.5)) + 
    labs(title=motif_names[i],x="UMAP 1",y="UMAP 2",color="ChromVAR z-score")+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background=element_rect(fill='transparent', color='black',linetype="solid",size=2.5),plot.title = element_text(hjust = 0.3,size=rel(3),face="bold"),axis.text=element_blank(),axis.ticks=element_blank(),axis.title.x=element_text(size=rel(1.5),vjust=0.5),axis.title.y=element_text(size=rel(1.5),vjust=0.5),legend.direction="horizontal",legend.position=c(0.8,0.1),legend.key.size = unit(0.65, 'cm'),legend.title=element_text(size=rel(1.1))) 
  p
  })
Cd36_vs_other_differential.activity<-  FindMarkers(
  object = Neutrophils,
  ident.1 = 'S_M.pulmonis',
  ident.2 = 'other',
  only.pos = TRUE,
  mean.fxn = rowMeans,
  fc.name = "avg_diff"
)

DefaultAssay(Neutrophils) <- "chromvar"
newpalette <- c(brewer.pal(9,"Greens")[5],brewer.pal(9,"Pastel1")[9])
Cd36_vs_other_differential.activity$Row.names<- rownames(Cd36_vs_other_differential.activity)
chromvar_violin_plot.ls <- lapply(1:length(motifs),function(i){
  chromvar_df <- FetchData(Neutrophils,vars=c(motifs[i],"S_M.pulmonis"))
  colnames(chromvar_df)[1] <- "motif"
  chromvar_df$S_M.pulmonis <- factor(chromvar_df$S_M.pulmonis,levels=c("S_M.pulmonis","other"),labels=c("S_M.pulmonis","other"))
  Cd36_vs_other_p_val_adj <- signif(Cd36_vs_other_differential.activity$p_val_adj[which(Cd36_vs_other_differential.activity$Row.names==motifs[i])],3)
  stat_df <- tibble::as_tibble(data.frame(group1="S_M.pulmonis",group2="other",p.adj=Cd36_vs_other_p_val_adj,y.position=round(max(chromvar_df$motif),1)+1.5))
  p <- ggplot() + 
    geom_violin(data=chromvar_df,aes(x=S_M.pulmonis,y=motif,fill=S_M.pulmonis),scale="width",trim=FALSE)+
    geom_boxplot(data=chromvar_df,aes(x=S_M.pulmonis,y=motif,fill=S_M.pulmonis),width=0.05, fill="white",outlier.size=0.8)+
    scale_fill_manual(values=newpalette)+
    scale_y_continuous(labels = fmt_dcimals(1),limits=c(round(min(chromvar_df$motif),1),round(max(chromvar_df$motif),1)+2)) + 
    labs(y="ChromVAR z-score",title=motif_names[i]) +
    stat_pvalue_manual(stat_df, label = "p = {p.adj}",label.size=4.5) +
    theme_classic()+
    theme(
          plot.title = element_text(hjust = 0.5,size=rel(3),face="bold"),
          axis.title.x = element_blank(),
          axis.title.y = element_text(color="black",size=rel(1.8)),
          axis.text.y = element_text(color="black",size=rel(1.5)), 
          axis.text.x = element_text(size=rel(1.5),color="black",vjust=1),
          axis.line = element_line(colour="black",size=1.5),
          axis.ticks = element_line(),
          plot.margin=unit(c(0, 0.2, 0.7, 0.2), "cm"),
          legend.position = " none") 
  p 
  }) 


pdf(str_c(out_dir,"S_M.pulmonis_Top_Motif_chromVAR_zscore_UMAP_VlnPlot_1.pdf"),width=12,,height=21)
(chromvar_UMAP_plot.ls[[1]] + chromvar_violin_plot.ls[[1]]+plot_layout(widths=c(5.5,5))) / 
(chromvar_UMAP_plot.ls[[2]] + chromvar_violin_plot.ls[[2]] + plot_layout(widths=c(5.5,5))) / 
(chromvar_UMAP_plot.ls[[3]] + chromvar_violin_plot.ls[[3]]+plot_layout(widths=c(5.5,5))) + plot_annotation(title="Accessibility of cis-elements",theme = theme(plot.title = element_text(hjust = 0.5,size = rel(2.8))))
dev.off()


pdf(str_c(out_dir,"S_M.pulmonis_Top_Motif_chromVAR_zscore_UMAP_VlnPlot_2.pdf"),width=12,,height=21)
(chromvar_UMAP_plot.ls[[4]] + chromvar_violin_plot.ls[[4]]+plot_layout(widths=c(5.5,5))) / 
(chromvar_UMAP_plot.ls[[5]] + chromvar_violin_plot.ls[[5]] + plot_layout(widths=c(5.5,5))) / 
(chromvar_UMAP_plot.ls[[6]] + chromvar_violin_plot.ls[[6]]+plot_layout(widths=c(5.5,5))) + plot_annotation(title="Accessibility of cis-elements",theme = theme(plot.title = element_text(hjust = 0.5,size = rel(2.8))))
dev.off()

# Top 20 Motif 
TopMotif<- C2_enriched_motifs_df$motif[1:20]
TopMotif_name<- C2_enriched_motifs_df$mouse_symbol[1:20]

DefaultAssay(Neutrophils)<- "ACTIVITY"
Cd36_vs_other_exp_df <- FindMarkers(Neutrophils,features=TopMotif_name,ident.1="S_M.pulmonis",logfc.threshold=-Inf,min.pct=-Inf)


# differential gene expression vs differential TF activity
DefaultAssay(Neutrophils) <- 'chromvar'
diff_chromvar_df <- FindMarkers(
  object = Neutrophils,
  ident.1 = 'S_M.pulmonis',
  ident.2 = 'other',
  mean.fxn = rowMeans,
  fc.name = "avg_diff",
  logfc.threshold=-Inf,
  min.pct=-Inf
)
pfm_df <- c()
for (i in 1:length(pfm)){
  pfm_df <- data.frame(ID=pfm[[i]]@ID,name=pfm[[i]]@name,species=pfm[[i]]@tags$species) %>% rbind(pfm_df,.)
}
diff_chromvar_df <- merge(diff_chromvar_df,pfm_df,by.x="row.names",by.y="ID",all.x=TRUE)
diff_chromvar_df$name <- gsub("\\(var\\.\\d\\)","",diff_chromvar_df$name)
diff_chromvar_df <- diff_chromvar_df[-grep("\\:\\:",diff_chromvar_df$name),]
convert2_mouse_symbol <- function(symbol){
  first_char <- toupper(substr(symbol,1,1))
  other_chars <- tolower(substr(symbol,2,nchar(symbol)))
  coverted_symbol <- str_c(first_char,other_chars)
  return(coverted_symbol)
}
diff_chromvar_df$mouse_symbol <- sapply(1:nrow(diff_chromvar_df),function(i){
  mouse_symbol <- convert2_mouse_symbol(diff_chromvar_df$name[i])
  mouse_symbol
  })

DefaultAssay(Neutrophils) <- "ACTIVITY"
genes <- unique(diff_chromvar_df$mouse_symbol[which(diff_chromvar_df$mouse_symbol %in% rownames(Neutrophils))])
diff_exp_df <- FindMarkers(Neutrophils,ident.1="S_M.pulmonis",features=genes,logfc.threshold=-Inf,min.pct=-Inf)

merged_df <- merge(diff_exp_df,diff_chromvar_df,by.x="row.names",by.y="mouse_symbol",all.x=TRUE)
merged_df <- merged_df %>% 
  arrange(abs(avg_diff))

high_label<- merged_df[which(abs(merged_df$avg_log2FC)>0.15&abs(merged_df$avg_diff)>0.7),"Row.names"]

show_label<- c(high_label,C2_enriched_motifs_df$mouse_symbol[1:6])
merged_df$label <- ifelse(merged_df$Row.names %in% show_label,merged_df$Row.names,"")
merged_df$chromvar_mlog10pvalue <- -log10(merged_df$p_val.y)
merged_df_01 <- merged_df %>%
  filter(label=="")
merged_df_02 <- merged_df %>%
  filter(label!="")

df_01 <- merged_df %>%
  filter(avg_log2FC<=-log2(1.5),p_val.x<0.05,avg_diff<0,p_val.y<0.05)
# Tshz2
df_02 <- merged_df %>%
  filter(avg_log2FC<0,p_val.x<0.05,avg_diff>0,p_val.y<0.05)
df_03 <- merged_df %>%
  filter(avg_log2FC>0,avg_diff<0,p_val.y<0.05)

# Figure 3C
pdf(str_c(out_dir,"S_M.pulmonis_vs_other_TF_expression_vs_activity.pdf"),width=10)
ggplot()+
  geom_rect(data=data.frame(xmin=0,xmax=Inf,ymin=0,ymax=Inf),aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),fill = brewer.pal(9,"Reds")[3] , alpha = 0.2)+
  geom_rect(data=data.frame(xmin=-Inf,xmax=0,ymin=-Inf,ymax=0),aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),fill = brewer.pal(9,"Reds")[3] , alpha = 0.2)+
  geom_rect(data=data.frame(xmin=0,xmax=Inf,ymin=-Inf,ymax=0),aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),fill = brewer.pal(9,"Blues")[3] , alpha = 0.2)+
  geom_rect(data=data.frame(xmin=-Inf,xmax=0,ymin=0,ymax=Inf),aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),fill = brewer.pal(9,"Blues")[3] , alpha = 0.2)+
  geom_point(data=merged_df_01,aes(x=avg_diff,y=avg_log2FC,size=pmax(pmin(chromvar_mlog10pvalue, 60), 0)),color="lightgrey")+
  geom_point(data=merged_df_02,aes(x=avg_diff,y=avg_log2FC,size=pmax(pmin(chromvar_mlog10pvalue, 60), 0),fill=avg_log2FC),pch=21,color="black")+
  scale_size_continuous(range=c(1,8),breaks = seq(0,60,20),limits=c(0,60))+
  scale_fill_gradientn(colours =c(brewer.pal(9,"Blues")[3],"#FFFFFF",brewer.pal(9,"Reds")[c(4,5,6,8)]) ,limits=c(-0.2,0.8),oob = scales::squish,guide = guide_colorbar(frame.colour = "black",frame.linewidth= 0.5,title.position="top",title.hjust = 0.5))+
  labs(y="Differential Gene ACTIVITY log2FoldChange",x="Differential TF chromVAR activity",title="S_M.pulmonis vs others",size="-log10(Differential TF chromVAR activity p value)",fill="log2(expression FoldChange)")+
  geom_text_repel(data=merged_df_02,aes(x=avg_diff,y=avg_log2FC,label=label),size=5,point.padding=unit(1.6, "lines"),arrow = arrow(length=unit(0.01, "npc")),max.overlaps=500,force=3,segment.color = "#cccccc")+
  ylim(c(-0.3,0.3))+
  xlim(c(-1.5,1.5))+
  geom_hline(yintercept=0,linetype="longdash")+
  geom_vline(xintercept=0,linetype="longdash")+
  theme_classic()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title = element_text(hjust = 0.5,size=rel(1.6)),axis.title=element_text(size=rel(1.6)),axis.text.y = element_text(color="black",size=rel(1.5)), axis.text.x = element_text(color="black",size=rel(1.5)),axis.line = element_line(colour="black",size = 1))
dev.off()



markers<- show_label
DefaultAssay(Neutrophils)<- "ACTIVITY"
Cd36_vs_other_exp_df <- FindMarkers(Neutrophils,features=markers,ident.1="S_M.pulmonis",logfc.threshold=-Inf,min.pct=-Inf)
gene_exp_violin_plot.ls <- lapply(1:length(markers),function(i){
  gene_exp_df <- FetchData(Neutrophils,vars = c(markers[i],"atacUMAP_1","atacUMAP_2","S_M.pulmonis"),slot = "data")
  colnames(gene_exp_df)[1] <- "gene" 
  gene_exp_df$S_M.pulmonis <- factor(gene_exp_df$S_M.pulmonis,levels=c("S_M.pulmonis","other"),labels=c("S_M.pulmonis","other"))
  Cd36_vs_other_p_val_adj <- signif(Cd36_vs_other_exp_df$p_val_adj[which(rownames(Cd36_vs_other_exp_df)==markers[i])],3)
  stat_df <- tibble::as_tibble(data.frame(group1="S_M.pulmonis",group2="other",p.adj=Cd36_vs_other_p_val_adj,y.position=round(max(gene_exp_df$gene),1)+0.7))
  stat_df$p.signif <- ifelse(stat_df$p.adj<0.05,stat_df$p.adj,"ns")
  p <- ggplot() + 
    geom_violin(data=gene_exp_df,aes(x=S_M.pulmonis,y=gene,fill=S_M.pulmonis),scale="width",trim=FALSE)+
    geom_boxplot(data=gene_exp_df,aes(x=S_M.pulmonis,y=gene,fill=S_M.pulmonis),width=0.05, fill="white",outlier.size=0.8)+
    scale_fill_manual(values=violin_palettes)+
    labs(y="Gene expression",title=markers[i]) + 
    scale_y_continuous(labels = fmt_dcimals(0),limits=c(round(min(gene_exp_df$gene),1),round(max(gene_exp_df$gene),1)+1))+
    stat_pvalue_manual(stat_df, label = "p = {p.adj}",label.size=5) + 
    theme_classic()+
    theme(
          plot.title = element_text(hjust = 0.5,size=rel(3),face="bold"),
          axis.title.x = element_blank(),
          axis.title.y = element_text(color="black",size=rel(1.8)),
          axis.text.y = element_text(color="black",size=rel(1.5)), 
          axis.text.x = element_text(size=rel(1.5),color="black",vjust=1),
          axis.line = element_line(colour="black",size=1.5),
          axis.ticks = element_line(),
          plot.margin=unit(c(0, 0.2, 0.7, 0.2), "cm"),
          legend.position = " none")
  p
  })

library(RColorBrewer)
gene_exp_palettes <- c(brewer.pal(8,"Pastel2")[8],brewer.pal(9,"Greens")[3:8])
DefaultAssay(Neutrophils) <- "ACTIVITY"
gene_exp_UMAP_plot.ls <- lapply(1:length(markers),function(i){
  gene_exp_df <- FetchData(Neutrophils,vars = c(markers[i],"atacUMAP_1","atacUMAP_2"),slot = "data")
  colnames(gene_exp_df)[1] <- "gene"
  p <- ggplot(data=gene_exp_df,aes(x=atacUMAP_1,y=atacUMAP_2,color=gene))+
    geom_point(size=0.1)+
    scale_colour_gradientn(colours = gene_exp_palettes,limits=c(0,quantile(gene_exp_df$gene,seq(0,1,0.01))[96]),oob = scales::squish,breaks=c(0,quantile(gene_exp_df$gene,seq(0,1,0.01))[96]),labels = c("0", round(quantile(gene_exp_df$gene,seq(0,1,0.01))[96],1)),guide = guide_colorbar(frame.colour = "black",frame.linewidth= 1,title.position="top", title.hjust = 0.5)) + 
    scale_x_continuous(limits=c(-12.5,7.5),expand = c(0,0))+
    scale_y_continuous(limits=c(-7.5,8),expand = c(0,0))+
    labs(title=markers[i],color="Gene ACTIVITY",x="UMAP 1",y="UMAP 2") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background=element_rect(fill='transparent', color='black',linetype="solid",size=3),plot.title = element_text(hjust = 0.5,size=rel(3),face="bold"),axis.text=element_blank(),axis.ticks=element_blank(),axis.title.x=element_text(size=rel(1.5),vjust=0.5),axis.title.y=element_text(size=rel(1.5),vjust=0.5),plot.margin=unit(c(0.3, 0.2, 0.7, 0.5), "cm"),legend.title = element_text(size=15),legend.text=element_text(size=10),legend.position=c(0.2,0.15),legend.direction="horizontal",legend.key.size = unit(0.8, 'cm'))
    p
  })

pdf(str_c(out_dir,"S_M.pulmonis_Etv6_Batf_gene_expression_UMAP_VlnPlot.pdf"),width=12,height=14)
(gene_exp_UMAP_plot.ls[[1]] + gene_exp_violin_plot.ls[[1]] + plot_layout(widths=c(5.5,5))) / 
(gene_exp_UMAP_plot.ls[[2]] + gene_exp_violin_plot.ls[[2]] + plot_layout(widths=c(5.5,5))) + plot_annotation(title="Gene ACTIVITY of TFs",theme = theme(plot.title = element_text(hjust = 0.5,size = rel(2.8))))
dev.off()

## Figure 3F,H,J - ChromVar Z-score for Top TF (UMAP + Violin plot)
library(ArchR)
chromvar_UMAP_plot.ls <- lapply(1:length(motifs),function(i){
  chromvar_df <- as.data.frame(Neutrophils@assays$chromvar@data) %>%
    tibble::rownames_to_column("motif")  %>%
    dplyr::filter(motif %in% motifs[i])
  chromvar_df <- reshape2::melt(chromvar_df,id="motif",variable.name="cell",value.name ="chromvar_zscore")
  meta_df <- FetchData(Neutrophils,vars=c("atacUMAP_1","atacUMAP_2"))
  chromvar_df <- merge(chromvar_df,meta_df,by.x="cell",by.y="row.names")
  p <- ggplot(data=chromvar_df,aes(x=atacUMAP_1,y=atacUMAP_2,color=chromvar_zscore))+
    geom_point(size=0.4)+
    scale_colour_gradientn(colours = as.vector(ArchRPalettes[["solarExtra"]]),limits=c(-1.5,1.5),oob = scales::squish,breaks=c(-1.5,0,1.5),guide = guide_colorbar(frame.colour = "black",frame.linewidth= 1,title.position="top", title.hjust = 0.5)) + 
    labs(title=show_label[i],x="UMAP 1",y="UMAP 2",color="ChromVAR z-score")+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background=element_rect(fill='transparent', color='black',linetype="solid",size=2.5),plot.title = element_text(hjust = 0.3,size=rel(3),face="bold"),axis.text=element_blank(),axis.ticks=element_blank(),axis.title.x=element_text(size=rel(1.5),vjust=0.5),axis.title.y=element_text(size=rel(1.5),vjust=0.5),legend.direction="horizontal",legend.position=c(0.8,0.1),legend.key.size = unit(0.65, 'cm'),legend.title=element_text(size=rel(1.1))) 
  p
  })
DefaultAssay(Neutrophils)<- "chromvar"
# Cd36_vs_other_differential.activity<-  FindMarkers(
#   object = Neutrophils,
#   ident.1 = 'S_M.pulmonis',
#   ident.2 = 'other',
#   only.pos = TRUE,
#   mean.fxn = rowMeans,
#   fc.name = "avg_diff"
# )
Cd36_vs_other_differential.activity<- diff_chromvar_df
show_label_motif<- merged_df[match(show_label,merged_df$Row.names),]$Row.names.y
motifs<- show_label_motif[1:2]
DefaultAssay(Neutrophils) <- "chromvar"
newpalette <- c(brewer.pal(9,"Greens")[5],brewer.pal(9,"Pastel1")[9])
#Cd36_vs_other_differential.activity$Row.names<- rownames(Cd36_vs_other_differential.activity)
chromvar_violin_plot.ls <- lapply(1:length(motifs),function(i){
  chromvar_df <- FetchData(Neutrophils,vars=c(motifs[i],"S_M.pulmonis"))
  colnames(chromvar_df)[1] <- "motif"
  chromvar_df$S_M.pulmonis <- factor(chromvar_df$S_M.pulmonis,levels=c("S_M.pulmonis","other"),labels=c("S_M.pulmonis","other"))
  Cd36_vs_other_p_val_adj <- signif(Cd36_vs_other_differential.activity$p_val_adj[which(Cd36_vs_other_differential.activity$Row.names==show_label_motif[i])],3)
  stat_df <- tibble::as_tibble(data.frame(group1="S_M.pulmonis",group2="other",p.adj=Cd36_vs_other_p_val_adj,y.position=round(max(chromvar_df$motif),1)+1.5))
  p <- ggplot() + 
    geom_violin(data=chromvar_df,aes(x=S_M.pulmonis,y=motif,fill=S_M.pulmonis),scale="width",trim=FALSE)+
    geom_boxplot(data=chromvar_df,aes(x=S_M.pulmonis,y=motif,fill=S_M.pulmonis),width=0.05, fill="white",outlier.size=0.8)+
    scale_fill_manual(values=newpalette)+
    scale_y_continuous(labels = fmt_dcimals(1),limits=c(round(min(chromvar_df$motif),1),round(max(chromvar_df$motif),1)+2)) + 
    labs(y="ChromVAR z-score",title=show_label[i]) +
    stat_pvalue_manual(stat_df, label = "p = {p.adj}",label.size=4.5) +
    theme_classic()+
    theme(
          plot.title = element_text(hjust = 0.5,size=rel(3),face="bold"),
          axis.title.x = element_blank(),
          axis.title.y = element_text(color="black",size=rel(1.8)),
          axis.text.y = element_text(color="black",size=rel(1.5)), 
          axis.text.x = element_text(size=rel(1.5),color="black",vjust=1),
          axis.line = element_line(colour="black",size=1.5),
          axis.ticks = element_line(),
          plot.margin=unit(c(0, 0.2, 0.7, 0.2), "cm"),
          legend.position = " none") 
  p 
  }) 


pdf(str_c(out_dir,"S_M.pulmonis_Etv6_Batf_chromVAR_zscore_UMAP_VlnPlot.pdf"),width=12,,height=14)
(chromvar_UMAP_plot.ls[[1]] + chromvar_violin_plot.ls[[1]]+plot_layout(widths=c(5.5,5))) / 
(chromvar_UMAP_plot.ls[[2]] + chromvar_violin_plot.ls[[2]]+plot_layout(widths=c(5.5,5))) + plot_annotation(title="Accessibility of cis-elements",theme = theme(plot.title = element_text(hjust = 0.5,size = rel(2.8))))
dev.off()


