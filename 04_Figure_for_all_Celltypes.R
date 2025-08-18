
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


library(RColorBrewer)
celltype_palette <- c(brewer.pal(9,"YlGn")[c(9,5,2)],
                      brewer.pal(12,"Set3")[8],
                      brewer.pal(9,"RdPu")[c(2,4,5,7)],
                      brewer.pal(9,"Blues")[c(2,4,6,8,9)],
                      brewer.pal(9,"BuPu")[5],brewer.pal(9,"BuPu")[7])


# WNN UMAP with annotation and cellnumber
cell_number_df <- OSN_mm@meta.data[,"cell_type",drop=FALSE] %>%
  group_by(cell_type) %>%
  summarise(counts=n())
cell_number_df$log10_counts <- log10(cell_number_df$counts)
cell_number_df$cell_type <- factor(cell_number_df$cell_type,levels=rev(levels(OSN_mm)))
cell_type_df <- cell_number_df[,"cell_type",drop=FALSE]
cell_type_df$x <- 1
pdf("./02_All_celltype/OSN_WNNUMAP_with_cell_types.pdf",width=9.5)
p1 <- DimPlot(OSN_mm, reduction = "wnn.umap",label=TRUE,cols=celltype_palette,repel=TRUE,label.size=5)+
  labs(title="WNNUMAP")+
  #geom_segment(aes(x = -14.641, y = -16.31245, xend = -12.2, yend = -16.31245),arrow = arrow(length = unit(0.3, "cm")))+
  #geom_segment(aes(x = -14.641, y = -16.31245, xend = -14.641, yend = -13.31245),arrow = arrow(length = unit(0.3, "cm")))+
  #annotate(geom = "text", x = -13.5, y = -17, label = "UMAP_1", color = "black") +
  #annotate(geom = "text", x = -15.2, y = -15, label = "UMAP_2", color = "black",angle = 90) +
  theme(plot.title = element_text(hjust = 0.5,size=rel(1.5),face="bold"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background=element_rect(fill='transparent', color='black',linetype="solid",size=1.5),axis.ticks=element_blank(),axis.text=element_blank(),axis.title=element_text(size=rel(1)),legend.position="none",plot.margin=unit(c(0.2, 0, 0.2, 0.3), "cm"))
p2 <- ggplot(data=cell_number_df,aes(x=log10_counts,y=cell_type,fill=cell_type))+
  geom_bar(stat = "identity",width=0.8)+
  labs(x="Log10 (# cells)")+
  #scale_x_continuous(limits=c(0,11100),breaks = c(1000,11000))+
  #scale_y_discrete(position = "right")+
  scale_fill_manual(values=rev(celltype_palette))+
  theme(plot.title = element_text(hjust = 0.5,size=rel(1.5),face="bold"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background=element_rect(fill='transparent', color='black',linetype="solid",size=1.5),axis.ticks.y=element_blank(),axis.text.y=element_blank(),axis.title.y=element_blank(),legend.position="none",plot.margin=unit(c(0.2, 0, 0.2, 0), "cm"),axis.text.x=element_text(size=rel(1),color="black"),axis.title.x=element_text(size=rel(1)))
p3 <- ggplot(data=cell_type_df,aes(x=x,y=cell_type))+
  geom_tile(aes(fill=cell_type),color="black")+
  scale_y_discrete(position = "right",expand=c(0,0))+
  scale_fill_manual(values=rev(celltype_palette))+
  theme(plot.title = element_text(hjust = 0.5,size=rel(1.5),face="bold"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background=element_rect(fill='transparent', color='black',linetype="solid",size=1.5),axis.ticks.x=element_blank(),axis.text.x=element_blank(),axis.title.x=element_blank(),legend.position="none",plot.margin=unit(c(0.2, 0, 0.2, 0.2), "cm"),axis.text.y=element_text(size=rel(1.5),color="black"),axis.title.y=element_blank())
p1 + p2 + p3 + plot_layout(widths=c(8,1.5,0.2))
dev.off()



library(RColorBrewer)
celltype_palette <- c(brewer.pal(9,"YlGn")[c(9,5,2)],
                      brewer.pal(12,"Set3")[8],
                      brewer.pal(9,"RdPu")[c(2,4,5,7)],
                      brewer.pal(9,"Blues")[c(2,4,6,8,9)],
                      brewer.pal(9,"BuPu")[5],brewer.pal(9,"BuPu")[7])


# ATAC UMAP with annotation and cellnumber
cell_number_df <- OSN_mm@meta.data[,"cell_type",drop=FALSE] %>%
  dplyr::group_by(cell_type) %>%
  dplyr::summarise(counts=n())

cell_number_df$log10_counts <- log10(cell_number_df$counts)
cell_number_df$cell_type <- factor(cell_number_df$cell_type,levels=rev(levels(OSN_mm)))
cell_type_df <- cell_number_df[,"cell_type",drop=FALSE]
cell_type_df$x <- 1
pdf("./02_All_celltype/OSN_ATAC_UMAP_with_cell_types.pdf",width=9.5)
p1 <- DimPlot(OSN_mm, reduction = "umap.atac",label=TRUE,cols=celltype_palette,repel=TRUE,label.size=5)+
  labs(title="ATACUMAP")+
  #geom_segment(aes(x = -14.641, y = -16.31245, xend = -12.2, yend = -16.31245),arrow = arrow(length = unit(0.3, "cm")))+
  #geom_segment(aes(x = -14.641, y = -16.31245, xend = -14.641, yend = -13.31245),arrow = arrow(length = unit(0.3, "cm")))+
  #annotate(geom = "text", x = -13.5, y = -17, label = "UMAP_1", color = "black") +
  #annotate(geom = "text", x = -15.2, y = -15, label = "UMAP_2", color = "black",angle = 90) +
  theme(plot.title = element_text(hjust = 0.5,size=rel(1.5),face="bold"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background=element_rect(fill='transparent', color='black',linetype="solid",size=1.5),axis.ticks=element_blank(),axis.text=element_blank(),axis.title=element_text(size=rel(1)),legend.position="none",plot.margin=unit(c(0.2, 0, 0.2, 0.3), "cm"))
p2 <- ggplot(data=cell_number_df,aes(x=log10_counts,y=cell_type,fill=cell_type))+
  geom_bar(stat = "identity",width=0.8)+
  labs(x="Log10 (# cells)")+
  #scale_x_continuous(limits=c(0,11100),breaks = c(1000,11000))+
  #scale_y_discrete(position = "right")+
  scale_fill_manual(values=rev(celltype_palette))+
  theme(plot.title = element_text(hjust = 0.5,size=rel(1.5),face="bold"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background=element_rect(fill='transparent', color='black',linetype="solid",size=1.5),axis.ticks.y=element_blank(),axis.text.y=element_blank(),axis.title.y=element_blank(),legend.position="none",plot.margin=unit(c(0.2, 0, 0.2, 0), "cm"),axis.text.x=element_text(size=rel(1),color="black"),axis.title.x=element_text(size=rel(1)))
p3 <- ggplot(data=cell_type_df,aes(x=x,y=cell_type))+
  geom_tile(aes(fill=cell_type),color="black")+
  scale_y_discrete(position = "right",expand=c(0,0))+
  scale_fill_manual(values=rev(celltype_palette))+
  theme(plot.title = element_text(hjust = 0.5,size=rel(1.5),face="bold"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background=element_rect(fill='transparent', color='black',linetype="solid",size=1.5),axis.ticks.x=element_blank(),axis.text.x=element_blank(),axis.title.x=element_blank(),legend.position="none",plot.margin=unit(c(0.2, 0, 0.2, 0.2), "cm"),axis.text.y=element_text(size=rel(1.5),color="black"),axis.title.y=element_blank())
p1 + p2 + p3 + plot_layout(widths=c(8,1.5,0.2))
dev.off()


pdf('./02_All_celltype/atac_S_M.pulmonis_ReadsPercentage_FeaturePlot_split_by_sample.pdf', width=12, height=5)
FeaturePlot(OSN_mm, reduction = 'umap.atac',split.by="orig.ident",keep.scale = "all",features = "atac_S_M.pulmonis_ReadsPercentage")& theme(legend.position = "right")
dev.off()
markers <- c("Syt1", #神经元
              "Omp",# Mature ORNs
              "Nqo1","Ncam2","Cd36",
              "Gap43","Gng8",#Immature ORNs
              "Trp63","Krt14",#HBCs
              "Sox2","Ermn","Cyp2g1","Cyp1a2",#支持细胞  Sustentacular cells
              "Sox9","Sox10",#Bowman's gland
              "Pebp1","Calb2", #球周细胞Periglomerular cells
              "Ascl3","Cftr",#Microvillar cells
              "Krt18","Trpm5", #Brush cells
              "Cd37","Cd79a",#B cells
              "S100a9","S100a8",#中性粒细胞 Neutrophils
              "Lyz2","S100a4", #Monocytes
              "Hbb-bs","Hbb-bt",#红细胞,Erythrocytes
              "Mcpt8","Ccl4","Ptprc","Itga2", #Basophils
              "C1qa","Ms4a7"#Macrophages
)
DefaultAssay(OSN_mm)<- "ACTIVITY"
pdf("./02_All_celltype/All_Marker_gene_Activity_Umap.pdf",height=6,width=6)
for(i in markers){
  print(i)
  p1 <- FeaturePlot(OSN_mm,order=TRUE, reduction = 'umap.atac',features =i)
print(p1)}
dev.off()



# cell proportion distribution 
pdf("./02_All_celltype/OSN_celltypetype_proportion.pdf",width=4,height=5)
df <- as.data.frame(OSN_mm@meta.data)
df_ct <- df %>%
    group_by(orig.ident, cell_type) %>%
    summarise(counts = n()) %>%
    mutate(cell_proportion = counts / sum(counts))

p <- ggplot(df_ct, aes(orig.ident, cell_proportion, fill=cell_type)) + 
    geom_bar(stat="identity", position = position_stack(reverse = TRUE)) +
    scale_fill_manual(values = celltype_palette) +
    theme_cowplot() +
    xlab("") + ylab("") +
    theme(axis.text.x = element_text(angle = 60, hjust = 1))
p
dev.off()

# atac_S_M.pulmonis_ReadsPercentage violin plot 

pdf("./02_All_celltype/OSN_celltypetype_M.pulmonis_violinplot.pdf",width=10,height=10)

p1 <- VlnPlot(subset(OSN_mm,rna_S_M.pulmonis_ReadsPercentage>0), "rna_S_M.pulmonis_ReadsPercentage",add.noise=FALSE, pt.size = 0,col=celltype_palette,group.by = "cell_type") +
geom_boxplot(width=0.2,outlier.size=0.01, position = position_dodge(0.9))
p2 <- VlnPlot(subset(OSN_mm,atac_S_M.pulmonis_ReadsPercentage>0), "atac_S_M.pulmonis_ReadsPercentage",add.noise=FALSE, pt.size = 0,col=celltype_palette,group.by = "cell_type") +
geom_boxplot(width=0.2,outlier.size=0.01, position = position_dodge(0.9))
p1/p2
p1 <- VlnPlot(OSN_mm, "rna_S_M.pulmonis_ReadsPercentage",add.noise=FALSE, pt.size = 0.01,col=celltype_palette,group.by = "cell_type") +
geom_boxplot(width=0.2,outlier.size=0.01, position = position_dodge(0.9))
p2 <- VlnPlot(OSN_mm, "atac_S_M.pulmonis_ReadsPercentage",add.noise=FALSE, pt.size = 0.01,col=celltype_palette,group.by = "cell_type") +
geom_boxplot(width=0.2,outlier.size=0.01, position = position_dodge(0.9))
p1/p2
dev.off()

## marker gene dotplot for show 
DefaultAssay(OSN_mm)<-"RNA"
Idents(OSN_mm)<- OSN_mm$cell_type
OSN_mm <- ScaleData(OSN_mm,features=rownames(OSN_mm))
markers <- FindAllMarkers(OSN_mm, only.pos = TRUE, min.pct = 0.05, logfc.threshold = 0.5)
table(markers$cluster)


# figure S1c - Dot plot showing the expression of marker genes for all cell types
markers <- c(              
             "Trp63","Krt14",#HBCs
             "Ncam2","Gng8",#Immature ORNs
             "Syt1", #神经元
             "Omp",# Mature ORNs
             "Cyp2g1","Ermn",#支持细胞
             "Ascl3","Cftr",#Microvillar cells
              "Atp1a2","Sox10",#Bowman's gland Malat1
              "Pebp1","Calb2", #球周细胞Periglomerular cells
              "Atp2a3","Trpm5", #Brush cells
              "Ptprc","Itga2", #Basophils
              "Cd37","Cd79a",#B cells
              "Hbb-bs","Hbb-bt",#红细胞,Erythrocytes
              "Ms4a7","Cdkn1c", #Monocytes
              "C1qa","Lyz2",#。
              "S100a9","S100a8"#中性粒细胞 Neutrophils
)

Idents(OSN_mm) <- factor(Idents(OSN_mm),levels=rev(levels(OSN_mm)))
DefaultAssay(OSN_mm) <- "RNA"
p <- DotPlot(OSN_mm,features=markers,cols =c("lightblue","red"))+RotatedAxis()
DotPlot_df <- p$data
pdf("./02_All_celltype/OSN_scRNA_cell_type_markers_dotplot.pdf",width=10,height=6)
ggplot(data=DotPlot_df,aes(x=features.plot,y=id))+
  geom_point(aes(size=pct.exp,color=avg.exp.scaled),)+
  scale_size_continuous(range=c(1,6),breaks = c(0,25,50,75),labels=c("0","25","50","75"),limits=c(0,100))+
  scale_color_gradientn(colours =c(brewer.pal(9,"Blues")[2],brewer.pal(9,"Reds")[c(2,4,6)]) ,limits=c(-1,2.5),oob = scales::squish)+
  labs(x="",y="",size="Percent Expressed",color="Scaled average expression")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5,size=rel(1.5),face="bold"),panel.background=element_rect(fill='transparent', color='black',linetype="solid",size=1.5),axis.text.x = element_text(angle =45,hjust=1,color='black',size=11),axis.text.y=element_text(color="black",size=11),axis.title=element_blank(),legend.position="bottom",legend.direction="horizontal")
dev.off()

## marker gene dotplot for show 
DefaultAssay(OSN_mm)<-"ACTIVITY"
Idents(OSN_mm)<- OSN_mm$cell_type
OSN_mm <- ScaleData(OSN_mm,features=rownames(OSN_mm))
markers <- FindAllMarkers(OSN_mm, only.pos = TRUE, min.pct = 0.05, logfc.threshold = 0.1)
table(markers$cluster)

markers <- c(              
             "Trp63","Krt14",#HBCs
             "Ncam2","Gng8",#Immature ORNs
             "Omp",# Mature ORNs
             "Cyp2g1","Ermn",#支持细胞
             "Dgki","Cftr",#Microvillar cells
              "Grik3","Sox10",#Bowman's gland
              "Pebp1","Calb2", #球周细胞Periglomerular cells
              "Krt18","Trpm5", #Brush cells
              "Ptprc","Itga2", #Basophils
              "Cd37","Cd79a",#B cells
              "Hbb-bs","Hbb-bt",#红细胞,Erythrocytes
              "Cx3cr1","Csf1", #Monocytes
              "C1qa","Lyz2",#Macrophages
              "S100a9","S100a8"#中性粒细胞 Neutrophils
)

DefaultAssay(OSN_mm) <- "ACTIVITY"
p <- DotPlot(OSN_mm,features=markers)+RotatedAxis()
DotPlot_df <- p$data
pdf("./02_All_celltype/OSN_scATAC_cell_type_markers_dotplot.pdf",width=10,height=6)
ggplot(data=DotPlot_df,aes(x=features.plot,y=id))+
  geom_point(aes(size=pct.exp,color=avg.exp.scaled),)+
  scale_size_continuous(range=c(1,6),breaks = c(0,25,50,75),labels=c("0","25","50","75"),limits=c(0,100))+
  scale_color_gradientn(colours =c(brewer.pal(9,"Blues")) ,limits=c(-1,2.5),oob = scales::squish)+
  labs(x="",y="",size="Percent Expressed",color="Scaled average expression")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5,size=rel(1.5),face="bold"),panel.background=element_rect(fill='transparent', color='black',linetype="solid",size=1.5),axis.text.x = element_text(angle =45,hjust=1,color='black',size=11),axis.text.y=element_text(color="black",size=11),axis.title=element_blank(),legend.position="bottom",legend.direction="horizontal")
dev.off()


## marker peaks of all cell types
markers <- c("Trp63","Krt14",#HBCs
             "Ncam2","Gng8",#Immature ORNs
             "Omp",# Mature ORNs
             "Cyp2g1","Ermn",#支持细胞
             "Ascl3","Cftr",#Microvillar cells
              "Sox9","Sox10",#Bowman's gland
              "Pebp1","Calb2", #球周细胞Periglomerular cells
              "Krt18","Trpm5", #Brush cells
              "Ptprc","Itga2", #Basophils
              "Cd37","Cd79a",#B cells
              "Hbb-bs","Hbb-bt",#红细胞,Erythrocytes
              "Cx3cr1","Ms4a7", #Monocytes
              "Mtss1","S100a12","Serpin1","Siglec10","Aplp2","Mpeg1","Fcgr3a","Cdkn1c","Cx3cr1","Cd11b","Iba1","Cd68","Cx3cr1",
              "C1qa","Lyz2",#Macrophages
              "S100a9","S100a8"#中性粒细胞 Neutrophils
)
pdf("./02_All_celltype/OSN_scATAC_cell_type_markers_CoveragePlot.pdf")
for (i in 1:length(markers)){
  p <- CoveragePlot(
    object = OSN_mm,
    region = markers[i],
    extend.upstream = 1000,
    extend.downstream = 1000,
    peaks=FALSE
  )
  print(p)
}
dev.off()

DefaultAssay(OSN_mm)<- "peaks_All_cluster"
annotation <- Annotation(OSN_mm)
genes_in_annotation <- unique(annotation$gene_name)
gene<- c("Cx3cr1","Csf1")

gene<-gene[which(gene%in%genes_in_annotation)]

pdf("./02_All_celltype/Monocyte_OSN_scATAC_cell_type_markers_CoveragePlot.pdf")
for (i in 1:length(gene)){
  print(i)
  p <- CoveragePlot(
    object = OSN_mm,
    region = gene[i],
    extend.upstream = 1000,
    extend.downstream = 1000,
    peaks=FALSE
  )
  print(p)
}
dev.off()

regions <- c(
  "chr11-100207000-100209000",#Krt14
  "chr7-16894000-16895000",# Gng8
  "chr7-98144500-98146500",# Omp
  "chr2-58052364-58053864",# Ermn
  "chr6-18290000-18298000",# Cftr
  "chr15-79163000-79166000",# Sox10
  "chr8-110167500-110169000",# Calb2
  "chr7-143093500-143095000",# Trpm5
  "chr13-114844000-114848000", # Itga2
  "chr7-24897000-24897700",# Cd79a
  "chr7-103827500-103828300",# Hbb−bs
  "chr3-107754000-107756000",# Csf1
  "chr4-136898000-136899999",# C1qa
  "chr3-90695400-90696200"# S100a9
  )
names <- c("Krt14","Gng8","Omp","Ermn","Cftr","Sox10","Calb2","Trpm5","Itga2","Cd79a","Hbb-bs","Csf1","C1qa","S100a9")
cov_plot.ls <- lapply(1:length(regions),function(i){
  cov_plot <- CoveragePlot(object = OSN_mm,region=regions[i],peaks=FALSE,annotation = FALSE)
  cov_df <- cov_plot$data 
  cov_plot <- ggplot(data = cov_df,mapping = aes(x = position, y = coverage, fill = group)) +
    geom_area(stat = "identity") +
    geom_hline(yintercept = 0, size = 0.1) +
    facet_wrap(facets = ~group, strip.position = "left", ncol = 1) +
    ylim(c(0, signif(max(cov_df$coverage, na.rm = TRUE), digits = 2))) +
    theme_bw() + 
    labs(title=names[i]) +
    theme(plot.title=element_text(face = "bold",size = rel(1.5),hjust = 0.5),panel.background=element_rect(fill='transparent', color='black',linetype="solid"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.position = "none",strip.text.y.left = element_blank(), strip.background = element_blank(),axis.text.y = element_blank(),panel.spacing.y = unit(x = 0, units = "line"),plot.margin=unit(c(0, 0,0,0), "cm"),axis.title.y = element_blank(),axis.ticks.y = element_blank(),axis.text.x = element_blank(),axis.title.x = element_blank(),axis.ticks=element_blank()) +
    scale_fill_manual(values = celltype_palette)
  cov_plot
  })

gene_plot.ls <- lapply(1:length(regions),function(i){
  gene_plot <- AnnotationPlot(object = OSN_mm,region = regions[i])
  gene_plot <- gene_plot + theme(panel.background=element_rect(fill='transparent', color='black',linetype="solid"),axis.text.x = element_blank(),axis.title.x = element_blank(),axis.ticks.x = element_blank(),axis.title.y=element_blank(),plot.margin=unit(c(0.2,0.03,0,0), "cm"))
  gene_plot
  })

pdf("./02_All_celltype/OSN_all_cell_type_markers_combined_CoveragePlot.pdf",width=18,height=6.5)
((cov_plot.ls[[1]] + gene_plot.ls[[1]] + plot_layout(heights=c(7,0.6))) |
(cov_plot.ls[[2]] + gene_plot.ls[[2]] + plot_layout(heights=c(7,0.6))) |
(cov_plot.ls[[3]] + gene_plot.ls[[3]] + plot_layout(heights=c(7,0.6))) |
(cov_plot.ls[[4]] + gene_plot.ls[[4]] + plot_layout(heights=c(7,0.6))) |
(cov_plot.ls[[5]] + gene_plot.ls[[5]] + plot_layout(heights=c(7,0.6))) |
(cov_plot.ls[[6]] + gene_plot.ls[[6]] + plot_layout(heights=c(7,0.6))) |
(cov_plot.ls[[7]] + gene_plot.ls[[7]] + plot_layout(heights=c(7,0.6))) |
(cov_plot.ls[[8]] + gene_plot.ls[[8]] + plot_layout(heights=c(7,0.6))) |
(cov_plot.ls[[9]] + gene_plot.ls[[9]] + plot_layout(heights=c(7,0.6))) |
(cov_plot.ls[[10]] + gene_plot.ls[[10]] + plot_layout(heights=c(7,0.6))) |
(cov_plot.ls[[11]] + gene_plot.ls[[11]] + plot_layout(heights=c(7,0.6))) |
(cov_plot.ls[[12]] + gene_plot.ls[[12]] + plot_layout(heights=c(7,0.6))) |
(cov_plot.ls[[13]] + gene_plot.ls[[13]] + plot_layout(heights=c(7,0.6))) | 
(cov_plot.ls[[14]] + gene_plot.ls[[14]] + plot_layout(heights=c(7,0.6)))) +  plot_layout(widths=c(1.1,1,1,1,1,1,1,1,1,1,1,1,1,1,1))
dev.off()

## Figure S4B - UMAP split by sample
newpalette <- c('#E5D2DD', '#53A85F')
pdf("./02_All_celltype/OSN_split_by_sample_UMAP.pdf",width=7,height=12)
DimPlot(OSN_mm, reduction = "wnn.umap",label=FALSE,cols=newpalette,group.by="orig.ident",split.by="orig.ident",ncol=1)+labs(title="WNN")+theme(plot.title = element_text(hjust = 0.5,size=rel(1.5),face="bold"))
dev.off()


## Figure S4B - UMAP split by sample
newpalette <- c('#E5D2DD', '#53A85F')
pdf("./02_All_celltype/OSN_split_by_sample_ATACUMAP.pdf",width=12,height=6)
DimPlot(OSN_mm, reduction = "umap.atac",label=FALSE,cols=newpalette,group.by="orig.ident",split.by="orig.ident",ncol=2)+labs(title="ATAC")+theme(plot.title = element_text(hjust = 0.5,size=rel(1.5),face="bold"))
dev.off()



