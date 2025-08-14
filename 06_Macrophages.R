
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

Macrophages<- subset(OSN_mm,idents="Macrophages")

# ATAC analysis
# We exclude the first dimension as this is typically correlated with sequencing depth
DefaultAssay(Macrophages) <- "ATAC"
Macrophages <- RunUMAP(Macrophages, reduction = 'integrated_LSI', dims = 2:30, reduction.name = "umap.atac", reduction.key = "atacUMAP_")
Macrophages <- FindNeighbors(object = Macrophages, reduction = 'integrated_LSI', dims = 2:30)
Macrophages <- FindClusters(Macrophages, resolution =0.5, algorithm = 3, verbose = FALSE)

my47colors <- c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
         '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
         '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398',
         '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
         '#712820', '#DCC1DD', '#CCE0F5', '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
         '#968175', '#161853', '#FF9999', '#344CB7', '#FFCC1D', '#116530', '#678983',
         '#678983', '#A19882', '#FFBCBC', '#24A19C', '#FF9A76')
pdf("./04_Macrophages/Macrophages_cluster_ATAC.pdf",width=15,height=5)
###cluster
p1 <- DimPlot(Macrophages, cols=my47colors, reduction = "umap.atac", group.by = "seurat_clusters", label = TRUE, label.size = 3.5, repel = TRUE) + ggtitle("seurat_clusters")
p2 <- DimPlot(Macrophages, cols=my47colors, reduction = "umap.atac", group.by = "orig.ident", shuffle=TRUE,label = F, label.size = 5, repel = TRUE) + ggtitle("sample")
p3 <- FeaturePlot(Macrophages, reduction = 'umap.atac',features = "atac_S_M.pulmonis_ReadsPercentage", ncol = 1)
p1 + p2 +p3 
dev.off()


Idents(Macrophages)<-Macrophages$seurat_clusters

Macrophages@meta.data$subtype<-Macrophages$seurat_clusters
table(Macrophages$subtype,Macrophages$orig.ident)
Idents(Macrophages)<-Macrophages$subtype;


# C3 vs other 


Idents(Macrophages)<-Macrophages$seurat_clusters
#####further annotation########
Macrophages <- RenameIdents(
  object = Macrophages,
  '0' = 'other',
  '1' = 'other',
  '2' = 'other',
  '3' = 'C3',
  '4' = 'other',
  '5' = 'other',
  '6' = 'other',
  '7' = 'other',
  '8' = 'other'
  )
Macrophages@meta.data$subtype<-Idents(Macrophages)
table(Macrophages$subtype,Macrophages$orig.ident)

Macrophages$subtype<-factor(Macrophages$subtype,
  levels=c("other","C3"))
Idents(Macrophages)<-Macrophages$subtype;


# Dimplot with annotation
pdf("./04_Macrophages/Macrophages_cluster_ATAC_subtype.pdf",width=10,height=8)
###cluster
p1 <- DimPlot(Macrophages, cols=my47colors[-1:-6], reduction = "umap.atac", group.by = "seurat_clusters", label = TRUE, label.size = 3.5, repel = TRUE) + ggtitle("seurat_clusters")
p2 <- DimPlot(Macrophages, cols=my47colors[-1:-6], reduction = "umap.atac", group.by = "subtype", label = TRUE, label.size = 3.5, repel = TRUE) + ggtitle("subtype")
pp1<- (p1 + p2 )
p3 <- DimPlot(Macrophages, cols=my47colors, reduction = "umap.atac", group.by = "orig.ident", shuffle=TRUE,label = F, label.size = 5, repel = TRUE) + ggtitle("sample")
p4 <- FeaturePlot(Macrophages, reduction = 'umap.atac',features = "atac_S_M.pulmonis_ReadsPercentage", ncol = 1)
pp2<- (p3 + p4 )
pp1/pp2
dev.off()


# cell proportion distribution 
pdf("./04_Macrophages/Macrophages_subtype_proportion.pdf",width=3,height=5)
df <- as.data.frame(Macrophages@meta.data)
df_ct <- df %>%
    group_by(orig.ident, subtype) %>%
    summarise(counts = n()) %>%
    mutate(cell_proportion = counts / sum(counts))

p <- ggplot(df_ct, aes(orig.ident, cell_proportion, fill=subtype)) + 
    geom_bar(stat="identity", position = position_stack(reverse = TRUE)) +
    scale_fill_manual(values = my47colors[-1:-6]) +
    theme_cowplot() +
    xlab("") + ylab("") +
    theme(axis.text.x = element_text(angle = 60, hjust = 1))
p
dev.off()


pdf("./04_Macrophages/Macrophages_subtype_M.pulmonis_violinplot.pdf",width=6,height=8)

p1 <- VlnPlot(subset(Macrophages,atac_S_M.pulmonis_ReadsPercentage>0), "atac_S_M.pulmonis_ReadsPercentage",add.noise=FALSE, pt.size = 0,col=my47colors[-1:-2],group.by = "subtype") +
geom_boxplot(width=0.1,outlier.size=0.01, position = position_dodge(0.9))

p2 <- VlnPlot(Macrophages, "atac_S_M.pulmonis_ReadsPercentage", pt.size = 0.01,col=my47colors[-1:-2],group.by = "subtype") +
geom_boxplot(width=0.1,outlier.size=0.01, position = position_dodge(0.9))
p1/p2
dev.off()


# 1. The DEG between C3 and other ;
log2FC = 0.3
padj = 0.05 
DefaultAssay(Macrophages)<- "ACTIVITY";
Macrophages <- ScaleData(Macrophages, features =rownames(Macrophages),verbose = FALSE)
IZ_CAvsCTRL_A_snRNA_markers <- FindMarkers(Macrophages,ident.1="C3",ident.2="other", min.pct = 0, logfc.threshold = 0)
markers<- IZ_CAvsCTRL_A_snRNA_markers
markers$threshold="ns";
markers[which(markers$avg_log2FC  > log2FC & markers$p_val_adj <padj),]$threshold="up_in_C3";
markers[which(markers$avg_log2FC  < (-log2FC) & markers$p_val_adj < padj),]$threshold="down_in_C3";
markers$threshold=factor(markers$threshold, levels=c('down_in_C3','up_in_C3','ns'))

blueYellow = c("1"="#352A86","2"="#343DAE","3"="#0262E0","4"="#1389D2","5"="#2DB7A3","6"="#A5BE6A","7"="#F8BA43","8"="#F6DA23","9"="#F8FA0D")
solarExtra<- c("#3361A5", "#248AF3", "#14B3FF" ,"#88CEEF" ,"#C1D5DC", "#EAD397" ,"#FDB31A" ,"#E42A2A","#A31D1D")

tmp<- markers[which(markers$threshold%in%c('down_in_C3','up_in_C3')),]
tmp<- tmp[order(tmp$avg_log2FC,decreasing=FALSE),]
write.csv(tmp,"./04_Macrophages/Macrophages_subtype_heatmapshow_DEG.csv")
gene<- rownames(tmp)
Macrophages<- ScaleData(Macrophages)
pdf("./04_Macrophages/Macrophages_DEG_heatmap_ACTIVITY.pdf",width=10,height=10)
DoHeatmap(object = Macrophages,disp.min = -1.5,disp.max = 1.5,features=gene, group.colors =my47colors,size = 4,group.by = "subtype") +scale_fill_gradientn(colors = solarExtra[3:8])
DoHeatmap(object = Macrophages,features=gene, group.colors =my47colors,size = 4,group.by = "subtype") +scale_fill_gradientn(colors = solarExtra)
dev.off();



# Find DEP for C3 vs other 
# recall peak 
library(BSgenome.Mmusculus.UCSC.mm10)
# ORN recall peak for subcluster 
DefaultAssay(Macrophages)<-"ATAC"
peak<-CallPeaks(
       Macrophages,
       group.by = "subtype",
       macs2.path = "/public/home/nieyg/biosoft/conda/bin/macs2",
       broad = FALSE,
       format = "BED",
       fragment.tempdir = tempdir(),
       effective.genome.size = 1.87e+09, #mm10 size
       outdir="./04_Macrophages/",
       combine.peaks=TRUE
)
macs2_counts <- FeatureMatrix(
     fragments = Fragments(Macrophages),
     features = peak,
     cells = colnames(Macrophages)
     )     
Macrophages[["peaks_subtype"]] <- CreateChromatinAssay(
  counts = macs2_counts,
  #fragments = fragpath[[i]],
  fragments = Fragments(Macrophages),
  annotation = Annotation(Macrophages)
)

saveRDS(Macrophages,"./04_Macrophages/Macrophages_subtype.rds")

Macrophages<- readRDS("./04_Macrophages/Macrophages_subtype.rds")


## Figure 1B - Heatmap showing DEGs between Cd36+ ORNs and non-Cd36+ ORNs
filtered_Ifi27_cells_mature_ORNs<- filtered_mature_ORNs_scATAC
Idents(filtered_Ifi27_cells_mature_ORNs) <- filtered_Ifi27_cells_mature_ORNs$Cd36_or_not
DefaultAssay(filtered_Ifi27_cells_mature_ORNs) <- "ACTIVITY"

Cd36_ORNs_markers <- FindMarkers(filtered_Ifi27_cells_mature_ORNs,logfc.threshold = 0.12,ident.1="S_M.pulmonis",assay="ACTIVITY")
Cd36_ORNs_markers <- Cd36_ORNs_markers %>%
  dplyr::filter(p_val<0.05,avg_log2FC >= log2(1)) %>%
  dplyr::arrange(desc(avg_log2FC))

other_ORNs_markers <- FindMarkers(filtered_Ifi27_cells_mature_ORNs,only.pos =T,logfc.threshold = 0.12,ident.1="other",assay="ACTIVITY")
other_ORNs_markers <- other_ORNs_markers %>%
  dplyr::filter(p_val<0.05,avg_log2FC >= log2(1)) %>%
  dplyr::arrange(desc(avg_log2FC))

set.seed(100)
scRNA_sampled_df <- c()
for (i in cell_types){
  cells <- rownames(filtered_Ifi27_cells_mature_ORNs@meta.data)[which(filtered_Ifi27_cells_mature_ORNs$Cd36_or_not==i)]
  if (length(cells)>=200){
    sampled_cells <- sample(cells,200,replace=FALSE)
  } else{
    sampled_cells <- sample(cells,200,replace=TRUE)
  }
  scRNA_sampled_df <- data.frame(idx=1:200,cell_type=i,sampled_cell=sampled_cells) %>% rbind(scRNA_sampled_df ,.)
}
DefaultAssay(filtered_Ifi27_cells_mature_ORNs) <- "ACTIVITY"
scRNA_df <- FetchData(object = filtered_Ifi27_cells_mature_ORNs, vars = c("Cd36_or_not",rownames(Cd36_ORNs_markers),rownames(other_ORNs_markers)),slot = "data",cells=scRNA_sampled_df$sampled_cell)
scRNA_df$Cd36_or_not <- factor(scRNA_df$Cd36_or_not,levels=c("S_M.pulmonis","other"),labels=c("S_M.pulmonis","other"))
scaled_scRNA_df <- t(apply(scRNA_df[,-1],2,scale))
colnames(scaled_scRNA_df) <- rownames(scRNA_df)
rescaled_scRNA_df <- scaled_scRNA_df
rescaled_scRNA_df[rescaled_scRNA_df>=1.5] <- 1.5
rescaled_scRNA_df[rescaled_scRNA_df<=-1.5] <- -1.5

library(ComplexHeatmap)
cell_type_cols <- c("#00A087FF",brewer.pal(8,"Set2")[8])
names(cell_type_cols) <- cell_types
highlight_genes <-c(rownames(Cd36_ORNs_markers[order(Cd36_ORNs_markers$avg_log2FC,decreasing=T),])[1:5],
  rownames(other_ORNs_markers[order(other_ORNs_markers$avg_log2FC,decreasing=T),])[1:5])


#highlight_genes <- c("Cd36","Cyp4v3","Abcd2","Mboat1","Snca","Rab3gap1","Macrod1","Hdac9","Tshz1","Zfp407","Cd9","Robo2","Cd55","Klkb1","Sgms2")
pdf(str_c(out_dir,"S_M.pulmonis_vs_other_MPs_DEGs_heatmap.pdf"),width=3.8,height=6)
ha1 <- HeatmapAnnotation(
  df=data.frame(row.names=rownames(scRNA_df),cell_type=scRNA_df$Cd36_or_not),
  col=list(cell_type=cell_type_cols),
  annotation_legend_param = list(cell_type=list(title="Cell type")),
  #gp = gpar(col = "black"),
  show_legend=FALSE,
  show_annotation_name=FALSE
  )
ha2 <- rowAnnotation(gene = anno_mark(at = which(rownames(rescaled_scRNA_df) %in% highlight_genes), labels = rownames(rescaled_scRNA_df)[which(rownames(rescaled_scRNA_df) %in% highlight_genes)],labels_gp=gpar(fontsize = 10),side = "left"))
h1 <- Heatmap(rescaled_scRNA_df,
  name="RNA z-score",
  col=as.vector(ArchRPalettes[["solarExtra"]]),
  cluster_columns=FALSE,
  column_split=factor(as.vector(scRNA_df$Cd36_or_not),levels=cell_types),
  column_title_gp = gpar(fontsize = 11),
  show_column_names=FALSE,
  show_row_names=FALSE,
  row_title_rot = 0,
  row_title_gp = gpar(fontsize = 11),
  column_gap = unit(0, "mm"),
  cluster_rows=FALSE,
  cluster_column_slices=FALSE,
  top_annotation=ha1,
  left_annotation=ha2,
  heatmap_legend_param = list(direction = "horizontal",legend_width = unit(4, "cm"), title_position = "topcenter",border="black")
)
draw(h1, heatmap_legend_side = "bottom")
dev.off()

write.csv(Cd36_ORNs_markers,"./04_Macrophages/S_M.pulmonis_MP_up_DEG_ACTIVITY.csv")
write.csv(other_ORNs_markers,"./04_Macrophages/S_M.pulmonis_MP_down_DEG_ACTIVITY.csv")

## GO enrichment of Cd36+ ORNs markers
.libPaths("/data/R02/nieyg/ori/biosoft/conda/envs/R/lib/R/library")

library(clusterProfiler)
library(org.Mm.eg.db)
Cd36_ORNs_markers<- read.csv("./04_Macrophages/S_M.pulmonis_MP_up_DEG_ACTIVITY.csv")
other_ORNs_markers<- read.csv("./04_Macrophages/S_M.pulmonis_MP_down_DEG_ACTIVITY.csv")
GO_ORA <- function(gene_sets){
  ENTREZID_id <- bitr(gene_sets,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Mm.eg.db",drop=TRUE)
  BP_GO_ORA <- enrichGO(gene=ENTREZID_id$ENTREZID,
    OrgDb=org.Mm.eg.db,
    ont= "BP",
    readable= TRUE)
  BP_GO_ORA <- simplify(BP_GO_ORA)
  MF_GO_ORA <- enrichGO(gene=ENTREZID_id$ENTREZID,
    OrgDb=org.Mm.eg.db,
    ont= "MF",
    readable= TRUE)
  MF_GO_ORA <- simplify(MF_GO_ORA)
  CC_GO_ORA <- enrichGO(gene=ENTREZID_id$ENTREZID,
    OrgDb=org.Mm.eg.db,
    ont= "CC",
    readable= TRUE)
  CC_GO_ORA <- simplify(CC_GO_ORA)
  BP_GO_ORA_result <- as.data.frame(BP_GO_ORA)
  MF_GO_ORA_result <- as.data.frame(MF_GO_ORA)
  CC_GO_ORA_result <- as.data.frame(CC_GO_ORA)
  BP_GO_ORA_result <- BP_GO_ORA_result[order(BP_GO_ORA_result$pvalue),]
  MF_GO_ORA_result <- MF_GO_ORA_result[order(MF_GO_ORA_result$pvalue),]
  CC_GO_ORA_result <- CC_GO_ORA_result[order(CC_GO_ORA_result$pvalue),]
  results <- list(BP_GO_ORA=BP_GO_ORA_result,
    MF_GO_ORA=MF_GO_ORA_result,
    CC_GO_ORA=CC_GO_ORA_result)
  return(results)
}

Cd36_ORNs_markers_GO <- GO_ORA(Cd36_ORNs_markers$X)
write.table(Cd36_ORNs_markers_GO[[1]],str_c(out_dir,"S_M.pulmonis_up_DEGs_BP_GO_ORA.txt"),row.names=FALSE,col.names=TRUE,sep="\t",quote=FALSE)
write.table(Cd36_ORNs_markers_GO[[2]],str_c(out_dir,"S_M.pulmonis_up_DEGs_MF_GO_ORA.txt"),row.names=FALSE,col.names=TRUE,sep="\t",quote=FALSE)
write.table(Cd36_ORNs_markers_GO[[3]],str_c(out_dir,"S_M.pulmonis_up_DEGs_CC_GO_ORA.txt"),row.names=FALSE,col.names=TRUE,sep="\t",quote=FALSE)

other_ORNs_markers_GO <- GO_ORA(other_ORNs_markers$X)
write.table(other_ORNs_markers_GO[[1]],str_c(out_dir,"S_M.pulmonis_down_DEGs_BP_GO_ORA.txt"),row.names=FALSE,col.names=TRUE,sep="\t",quote=FALSE)
write.table(other_ORNs_markers_GO[[2]],str_c(out_dir,"S_M.pulmonis_down_DEGs_MF_GO_ORA.txt"),row.names=FALSE,col.names=TRUE,sep="\t",quote=FALSE)
write.table(other_ORNs_markers_GO[[3]],str_c(out_dir,"S_M.pulmonis_down_DEGs_CC_GO_ORA.txt"),row.names=FALSE,col.names=TRUE,sep="\t",quote=FALSE)

library(pheatmap)
library(clusterProfiler)
library(org.Mm.eg.db)

CAV_list<- list(Cd36_ORNs_markers$X,other_ORNs_markers$X)
subtype<- c("S_M.pulmonis","others")
all_ego<- data.frame()
pdf("./04_Macrophages/Macrophages_subtype_DEG_GO_BP_ACTIVITY.pdf")
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
    write.csv(ego,paste0("./04_Macrophages/",subtype[i],"-GO-BP.csv",sep=""))
    ego<- as.data.frame(ego)
    ego$celltype<- subtype[i];
    all_ego<- rbind(all_ego,ego);
}
dev.off()

write.csv(all_ego,"./04_Macrophages/Macrophages_Allsubtype_GO.csv")


## Find Cd36+ vs Cd36- OSNs differentially accessible peaks 
DefaultAssay(Macrophages) <- "ATAC"
filtered_mature_ORNs_scATAC<- Macrophages
filtered_mature_ORNs_scATAC$S_M.pulmonis <- ifelse(filtered_mature_ORNs_scATAC$subtype=="C3","S_M.pulmonis","other")
filtered_mature_ORNs_scATAC$S_M.pulmonis <- factor(filtered_mature_ORNs_scATAC$S_M.pulmonis,levels=c("S_M.pulmonis","other"))
Idents(filtered_mature_ORNs_scATAC) <- filtered_mature_ORNs_scATAC$S_M.pulmonis
DefaultAssay(filtered_mature_ORNs_scATAC) <- "peaks_subtype"

Cd36_vs_other_da_peaks <- FindMarkers(
  object = filtered_mature_ORNs_scATAC,
  ident.1 = "S_M.pulmonis",
  ident.2 = "other",
  test.use = 'LR',
  latent.vars = 'nCount_ATAC'
)

Cd36_up_peaks <- Cd36_vs_other_da_peaks %>%
  dplyr::filter(p_val < 0.05,avg_log2FC>=log2(1.5)) %>%
  dplyr::arrange(desc(avg_log2FC)) 
Cd36_up_peaks_ann <-  ClosestFeature(object=filtered_mature_ORNs_scATAC,
  regions=rownames(Cd36_up_peaks))
strig2GRanges <- function(string) {
  chr <- sapply(1:length(string),function(i){
    unlist(strsplit(string[i],"-"))[1]
    })
  start <- sapply(1:length(string),function(i){
    as.numeric(unlist(strsplit(string[i],"-"))[2])
    })
  end <- sapply(1:length(string),function(i){
    as.numeric(unlist(strsplit(string[i],"-"))[3])
    })
  gr <- GRanges(chr,IRanges(start = start, end = end), strand = "*")
  gr
}
Cd36_up_peaks.gr <- strig2GRanges(rownames(Cd36_up_peaks))
Cd36_up_peaks_ann_df <- merge(Cd36_up_peaks,Cd36_up_peaks_ann[,c("gene_name","query_region","distance")],by.x="row.names",by.y="query_region",all.x=TRUE,sort=FALSE)
out_dir<- "/data/R02/nieyg/project/OSN_mm/joint/04_Macrophages/"
write.table(Cd36_up_peaks_ann_df,str_c(out_dir,"S_M.pulmonis_enrich_up_peaks_with_annotation.txt"),row.names=FALSE,col.names=TRUE,sep="\t",quote=FALSE)


Cd36_down_peaks <- Cd36_vs_other_da_peaks %>%
  dplyr::filter(p_val < 0.05,avg_log2FC<=-log2(1.5)) %>%
  dplyr::arrange(avg_log2FC) 
Cd36_down_peaks_ann <- ClosestFeature(object=filtered_mature_ORNs_scATAC,
  regions=rownames(Cd36_down_peaks))
Cd36_down_peaks_ann_df <-  merge(Cd36_down_peaks,Cd36_down_peaks_ann[,c("gene_name","query_region","distance")],by.x="row.names",by.y="query_region",all.x=TRUE,sort=FALSE)
write.table(Cd36_down_peaks_ann_df,str_c(out_dir,"S_M.pulmonis_enrich_down_peaks_with_annotation.txt"),row.names=FALSE,col.names=TRUE,sep="\t",quote=FALSE)


## Figure 3B - Heatmap showing DARs

filtered_mature_ORNs_scATAC$Cd36_or_not<- filtered_mature_ORNs_scATAC$S_M.pulmonis
filtered_mature_ORNs_scATAC$Cd36_or_not<- factor(filtered_mature_ORNs_scATAC$Cd36_or_not,levels=c("S_M.pulmonis","other"))
DefaultAssay(filtered_mature_ORNs_scATAC) <- "peaks_subtype"
set.seed(100)
cell_types <- c("S_M.pulmonis","other")
sampled_df <- c()
for (i in cell_types){
  cells <- rownames(filtered_mature_ORNs_scATAC@meta.data)[which(filtered_mature_ORNs_scATAC$Cd36_or_not==i)]
  if (length(cells)>=200){
    sampled_cells <- sample(cells,200,replace=FALSE)
  } else{
    sampled_cells <- sample(cells,200,replace=TRUE)
  }
  sampled_df <- data.frame(idx=1:200,cell_type=i,sampled_cell=sampled_cells) %>% rbind(sampled_df ,.)
}
scATAC_df <- FetchData(object = filtered_mature_ORNs_scATAC, vars = unique(c("Cd36_or_not",rownames(Cd36_up_peaks),rownames(Cd36_down_peaks))),slot = "data",cells=sampled_df$sampled_cell)
scaled_scATAC_df <- t(apply(scATAC_df[,-1],2,scale))
colnames(scaled_scATAC_df) <- rownames(scATAC_df)
rescaled_scATAC_df <- scaled_scATAC_df
 rescaled_scATAC_df[which(rescaled_scATAC_df>=1)] <- 1
 rescaled_scATAC_df[which(rescaled_scATAC_df<= -1)] <- -1

cell_type_cols <- c("#00A087FF",brewer.pal(8,"Set2")[8])
names(cell_type_cols) <- cell_types
library(ComplexHeatmap)
library(ArchR)
pdf(str_c(out_dir,"S_M.pulmonis_vs_other_MPs_DARs_heatmap.pdf"),width=3,height=6)
ha2 <- HeatmapAnnotation(
  df=data.frame(row.names=rownames(scATAC_df),subtype=scATAC_df$Cd36_or_not),
  col=list(subtype=cell_type_cols),
  annotation_legend_param = list(subtype=list(title="Cell type")),
  show_legend=FALSE,
  show_annotation_name=FALSE
  )
ha3<- HeatmapAnnotation(subtype = c(rep("#00A087FF",200),rep(brewer.pal(8,"Set2")[8],200)))
h2 <- Heatmap(rescaled_scATAC_df,
  name="ATAC z-score",
  col=as.vector(ArchRPalettes[["blueYellow"]]),
  cluster_columns=FALSE,
  column_split=factor(as.vector(scATAC_df$Cd36_or_not),levels=cell_types),
  column_title_gp = gpar(fontsize = 11),
  show_column_names=FALSE,
  show_row_names=FALSE,
  row_title_rot = 0,
  row_title_gp = gpar(fontsize = 11),
  column_gap = unit(0, "mm"),
  cluster_rows=FALSE,
  cluster_column_slices=FALSE,
  top_annotation=ha2,
  heatmap_legend_param = list(direction = "horizontal",legend_width = unit(4, "cm"), title_position = "topcenter",border="black")
)
draw(h2, heatmap_legend_side = "bottom")

dev.off()


## motif enrichment of DARs
library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Mmusculus.UCSC.mm10)
library(patchwork)
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)

pfm_df <- c()
for (i in 1:length(pfm)){
  pfm_df <- data.frame(ID=pfm[[i]]@ID,name=pfm[[i]]@name,species=pfm[[i]]@tags$species) %>% rbind(pfm_df,.)
}
DefaultAssay(filtered_mature_ORNs_scATAC)<-"peaks_subtype"
filtered_mature_ORNs_scATAC <- AddMotifs(
  object = filtered_mature_ORNs_scATAC,
  genome = BSgenome.Mmusculus.UCSC.mm10,
  pfm = pfm
)

Cd36_up_enriched.motifs <- FindMotifs(
  object = filtered_mature_ORNs_scATAC,
  features = rownames(Cd36_up_peaks)
) 
Cd36_up_enriched.motifs$rank <- 1:nrow(Cd36_up_enriched.motifs)

Cd36_down_enriched.motifs <- FindMotifs(
  object = filtered_mature_ORNs_scATAC,
  features = rownames(Cd36_down_peaks)
)
Cd36_down_enriched.motifs$rank <- 1:nrow(Cd36_down_enriched.motifs)
out_dir<- "/data/R02/nieyg/project/OSN_mm/joint/04_Macrophages/"

## Figure S5D - ranking plot of Cd36+ OSNs down peaks enriched motifs
library(ggrepel)
Cd36_down_enriched_motifs_df <- Cd36_down_enriched.motifs
Cd36_down_enriched_motifs_df$motif.name <- gsub("\\(var\\.\\d\\)","",Cd36_down_enriched_motifs_df$motif.name)
Cd36_down_enriched_motifs_df$motif.name <- gsub("\\_02","",Cd36_down_enriched_motifs_df$motif.name)
Cd36_down_enriched_motifs_df <- Cd36_down_enriched_motifs_df[-grep("\\:\\:",Cd36_down_enriched_motifs_df$motif.name),]
convert2_mouse_symbol <- function(symbol){
  first_char <- toupper(substr(symbol,1,1))
  other_chars <- tolower(substr(symbol,2,nchar(symbol)))
  coverted_symbol <- str_c(first_char,other_chars)
  return(coverted_symbol)
}
Cd36_down_enriched_motifs_df$mouse_symbol <- sapply(1:nrow(Cd36_down_enriched_motifs_df),function(i){
  mouse_symbol <- convert2_mouse_symbol(Cd36_down_enriched_motifs_df$motif.name[i])
  mouse_symbol
  })
Cd36_down_enriched_motifs_df$label <- ifelse(Cd36_down_enriched_motifs_df$mouse_symbol %in% c(Cd36_down_enriched_motifs_df$mouse_symbol[1:10],"Lhx2"),Cd36_down_enriched_motifs_df$mouse_symbol,"")
Cd36_down_enriched_motifs_df$rank <- 1:nrow(Cd36_down_enriched_motifs_df)
Cd36_down_enriched_motifs_df$mlog10pvalue <- -log10(Cd36_down_enriched_motifs_df$pvalue)
pdf(str_c(out_dir,"S_M.pulmonis_down_peaks_enriched_motifs_ranking_plot.pdf"),height=7.5)
ggplot(data=Cd36_down_enriched_motifs_df[1:50,],aes(x=rank,y=mlog10pvalue))+
  geom_point(aes(color=fold.enrichment,size=percent.observed))+
  labs(x="Motif rank",y="-log10(p_value)",color="Enrichment",size="Percentage of peaks with corresponding motif",title="Enriched motifs of down-regulated peaks in S_M.pulmonis MPs")+
  geom_text_repel(aes(label=label),size=5,point.padding = 0.2,
    nudge_x = .15,
    nudge_y = .5,
    segment.curvature = -1e-20,
    arrow = arrow(length = unit(0.015, "npc")))+
  scale_colour_gradient(low = "yellow", high = "red")+
  theme_classic()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title = element_text(hjust = 0.5,size=rel(1.6)),axis.title=element_text(size=rel(1.6)),axis.text.y = element_text(color="black",size=rel(1.5)), axis.text.x = element_text(color="black",size=rel(1.5)),axis.line = element_line(colour="black",size = 1),legend.position="bottom")
dev.off()
write.table(Cd36_down_enriched_motifs_df,str_c(out_dir,"S_M.pulmonis_down_peaks_enriched_motifs.txt"),sep="\t",quote=FALSE,col.names=TRUE,row.names=FALSE)

# Figure S5C - ranking plot of Cd36+ OSNs up peaks enriched motifs
Cd36_up_enriched_motifs_df <- Cd36_up_enriched.motifs
Cd36_up_enriched_motifs_df$motif.name <- gsub("\\(var\\.\\d\\)","",Cd36_up_enriched_motifs_df$motif.name)
Cd36_up_enriched_motifs_df$motif.name <- gsub("\\_02","",Cd36_up_enriched_motifs_df$motif.name)
Cd36_up_enriched_motifs_df <- Cd36_up_enriched_motifs_df[-grep("\\:\\:",Cd36_up_enriched_motifs_df$motif.name),]
convert2_mouse_symbol <- function(symbol){
  first_char <- toupper(substr(symbol,1,1))
  other_chars <- tolower(substr(symbol,2,nchar(symbol)))
  coverted_symbol <- str_c(first_char,other_chars)
  return(coverted_symbol)
}
Cd36_up_enriched_motifs_df$mouse_symbol <- sapply(1:nrow(Cd36_up_enriched_motifs_df),function(i){
  mouse_symbol <- convert2_mouse_symbol(Cd36_up_enriched_motifs_df$motif.name[i])
  mouse_symbol
  })
Cd36_up_enriched_motifs_df$label <- ifelse(Cd36_up_enriched_motifs_df$mouse_symbol %in% c(Cd36_up_enriched_motifs_df$mouse_symbol[1:10],"Tshz1","Atf4"),Cd36_up_enriched_motifs_df$mouse_symbol,"")
Cd36_up_enriched_motifs_df$rank <- 1:nrow(Cd36_up_enriched_motifs_df)
Cd36_up_enriched_motifs_df$mlog10pvalue <- -log10(Cd36_up_enriched_motifs_df$pvalue)
pdf(str_c(out_dir,"S_M.pulmonis_up_peaks_enriched_motifs_ranking_plot.pdf"),height=7.5)
ggplot(data=Cd36_up_enriched_motifs_df[1:50,],aes(x=rank,y=mlog10pvalue))+
  geom_point(aes(color=fold.enrichment,size=percent.observed))+
  labs(x="Motif rank",y="-log10(p_value)",color="Enrichment",size="Percentage of peaks with corresponding motif",title="Enriched motifs of up-regulated peaks in S_M.pulmonis MPs")+
  geom_text_repel(aes(label=label),size=5,point.padding = 0.2,
    nudge_x = .15,
    nudge_y = .5,
    segment.curvature = -1e-20,
    arrow = arrow(length = unit(0.015, "npc")),max.overlaps=500)+
  scale_colour_gradient(low = "yellow", high = "red")+
  theme_classic()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title = element_text(hjust = 0.5,size=rel(1.6)),axis.title=element_text(size=rel(1.6)),axis.text.y = element_text(color="black",size=rel(1.5)), axis.text.x = element_text(color="black",size=rel(1.5)),axis.line = element_line(colour="black",size = 1),legend.position="bottom")
dev.off()
write.table(Cd36_up_enriched_motifs_df,str_c(out_dir,"S_M.pulmonis_up_peaks_enriched_motifs.txt"),sep="\t",quote=FALSE,col.names=TRUE,row.names=FALSE)

pdf("./04_Macrophages/Macrophages_DEP_S_M.pulmonis_upordown_Top_MotifPlot.pdf",width=6,height=3)
MotifPlot(object = filtered_mature_ORNs_scATAC,motifs = head(rownames(Cd36_up_enriched_motifs_df)))
MotifPlot(object = filtered_mature_ORNs_scATAC,motifs = head(rownames(Cd36_down_enriched_motifs_df)))
dev.off()


## chromVAR
# Compute a per-cell motif activity score by running chromVAR
DefaultAssay(filtered_mature_ORNs_scATAC) <- "ATAC"
filtered_mature_ORNs_scATAC <- AddMotifs(
  object = filtered_mature_ORNs_scATAC,
  genome = BSgenome.Mmusculus.UCSC.mm10,
  pfm = pfm
)
filtered_mature_ORNs_scATAC <- RunChromVAR(
  object = filtered_mature_ORNs_scATAC,
  genome = BSgenome.Mmusculus.UCSC.mm10
)

# differential gene expression vs differential TF activity
DefaultAssay(filtered_mature_ORNs_scATAC) <- 'chromvar'
diff_chromvar_df <- FindMarkers(
  object = filtered_mature_ORNs_scATAC,
  ident.1 = 'S_M.pulmonis',
  ident.2 = 'other',
  mean.fxn = rowMeans,
  fc.name = "avg_diff",
  logfc.threshold=-Inf,
  min.pct=-Inf
)
diff_chromvar_df <- merge(diff_chromvar_df,pfm_df,by.x="row.names",by.y="ID",all.x=TRUE)
#diff_chromvar_df <- diff_chromvar_df[-grep("Tshz1-01",diff_chromvar_df$Row.names),]
#diff_chromvar_df$name[1:5] <- c("Trps1","Tshz1","Tshz2","Tshz2","Tshz2")
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

DefaultAssay(filtered_mature_ORNs_scATAC) <- "ACTIVITY"
genes <- unique(diff_chromvar_df$mouse_symbol[which(diff_chromvar_df$mouse_symbol %in% rownames(filtered_mature_ORNs_scATAC))])
diff_exp_df <- FindMarkers(filtered_mature_ORNs_scATAC,ident.1="S_M.pulmonis",features=genes,logfc.threshold=-Inf,min.pct=-Inf)

merged_df <- merge(diff_exp_df,diff_chromvar_df,by.x="row.names",by.y="mouse_symbol",all.x=TRUE)
merged_df <- merged_df %>% 
  arrange(abs(avg_diff))

high_label<- merged_df[which(abs(merged_df$avg_log2FC)>0.15&abs(merged_df$avg_diff)>0.7),"Row.names"]
show_label<- c(high_label,Cd36_up_enriched_motifs_df$mouse_symbol[1:6],Cd36_down_enriched_motifs_df$mouse_symbol[1:6])

merged_df$label <- ifelse(merged_df$Row.names %in% show_label,merged_df$Row.names,"")

merged_df$chromvar_mlog10pvalue <- -log10(merged_df$p_val.y)

merged_df<- merged_df[-441,]

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
pdf(str_c(out_dir,"S_M.pulmonis_vs_other_MPs_TF_expression_vs_activity.pdf"),width=10)
ggplot()+
  geom_rect(data=data.frame(xmin=0,xmax=Inf,ymin=0,ymax=Inf),aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),fill = brewer.pal(9,"Reds")[3] , alpha = 0.2)+
  geom_rect(data=data.frame(xmin=-Inf,xmax=0,ymin=-Inf,ymax=0),aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),fill = brewer.pal(9,"Reds")[3] , alpha = 0.2)+
  geom_rect(data=data.frame(xmin=0,xmax=Inf,ymin=-Inf,ymax=0),aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),fill = brewer.pal(9,"Blues")[3] , alpha = 0.2)+
  geom_rect(data=data.frame(xmin=-Inf,xmax=0,ymin=0,ymax=Inf),aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),fill = brewer.pal(9,"Blues")[3] , alpha = 0.2)+
  geom_point(data=merged_df_01,aes(x=avg_diff,y=avg_log2FC,size=pmax(pmin(chromvar_mlog10pvalue, 60), 0)),color="lightgrey")+
  geom_point(data=merged_df_02,aes(x=avg_diff,y=avg_log2FC,size=pmax(pmin(chromvar_mlog10pvalue, 60), 0),fill=avg_log2FC),pch=21,color="black")+
  scale_size_continuous(range=c(1,8),breaks = seq(0,60,20),limits=c(0,60))+
  scale_fill_gradientn(colours =c(brewer.pal(9,"Blues")[3],"#FFFFFF",brewer.pal(9,"Reds")[c(4,5,6,8)]) ,limits=c(-0.2,0.8),oob = scales::squish,guide = guide_colorbar(frame.colour = "black",frame.linewidth= 0.5,title.position="top",title.hjust = 0.5))+
  labs(y="Differential ACTIVITY log2FoldChange",x="Differential TF chromVAR activity",title="S_M.pulmonis vs others",size="-log10(Differential TF chromVAR activity p value)",fill="log2(expression FoldChange)")+
  geom_text_repel(data=merged_df_02,aes(x=avg_diff,y=avg_log2FC,label=label),size=5,point.padding=unit(1.6, "lines"),arrow = arrow(length=unit(0.01, "npc")),max.overlaps=500,force=3,segment.color = "#cccccc")+
  ylim(c(-0.35,0.35))+
  xlim(c(-4,4))+
  geom_hline(yintercept=0,linetype="longdash")+
  geom_vline(xintercept=0,linetype="longdash")+
  theme_classic()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title = element_text(hjust = 0.5,size=rel(1.6)),axis.title=element_text(size=rel(1.6)),axis.text.y = element_text(color="black",size=rel(1.5)), axis.text.x = element_text(color="black",size=rel(1.5)),axis.line = element_line(colour="black",size = 1))
dev.off()


## Figure 3E,G,I - gene expression of Mef2a,Lhx2 and Tshz1 (UMAP + Violin plot)


markers<- show_label[-12]
DefaultAssay(filtered_mature_ORNs_scATAC)<- "ACTIVITY"
Cd36_vs_other_exp_df <- FindMarkers(filtered_mature_ORNs_scATAC,features=markers,ident.1="S_M.pulmonis",logfc.threshold=-Inf,min.pct=-Inf)
gene_exp_violin_plot.ls <- lapply(1:length(markers),function(i){
  gene_exp_df <- FetchData(filtered_mature_ORNs_scATAC,vars = c(markers[i],"atacUMAP_1","atacUMAP_2","S_M.pulmonis"),slot = "data")
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
DefaultAssay(filtered_mature_ORNs_scATAC) <- "ACTIVITY"
gene_exp_UMAP_plot.ls <- lapply(1:length(markers),function(i){
  gene_exp_df <- FetchData(filtered_mature_ORNs_scATAC,vars = c(markers[i],"atacUMAP_1","atacUMAP_2"),slot = "data")
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


pdf(str_c(out_dir,"S_M.pulmonis_Tcf4_Bach2_gene_expression_UMAP_VlnPlot.pdf"),width=12,height=21)
(gene_exp_UMAP_plot.ls[[1]] + gene_exp_violin_plot.ls[[1]] + plot_layout(widths=c(5.5,5))) / 
(gene_exp_UMAP_plot.ls[[2]] + gene_exp_violin_plot.ls[[2]] + plot_layout(widths=c(5.5,5))) / 
(gene_exp_UMAP_plot.ls[[3]] + gene_exp_violin_plot.ls[[3]] + plot_layout(widths=c(5.5,5))) + plot_annotation(title="Gene ACTIVITY of TFs",theme = theme(plot.title = element_text(hjust = 0.5,size = rel(2.8))))
dev.off()

## Figure 3F,H,J - ChromVar Z-score for Top TF (UMAP + Violin plot)
library(ArchR)
motifs<- show_label_motif[1:3]
chromvar_UMAP_plot.ls <- lapply(1:length(motifs),function(i){
  chromvar_df <- as.data.frame(filtered_mature_ORNs_scATAC@assays$chromvar@data) %>%
    tibble::rownames_to_column("motif")  %>%
    dplyr::filter(motif %in% motifs[i])
  chromvar_df <- reshape2::melt(chromvar_df,id="motif",variable.name="cell",value.name ="chromvar_zscore")
  meta_df <- FetchData(filtered_mature_ORNs_scATAC,vars=c("atacUMAP_1","atacUMAP_2"))
  chromvar_df <- merge(chromvar_df,meta_df,by.x="cell",by.y="row.names")
  p <- ggplot(data=chromvar_df,aes(x=atacUMAP_1,y=atacUMAP_2,color=chromvar_zscore))+
    geom_point(size=0.4)+
    scale_colour_gradientn(colours = as.vector(ArchRPalettes[["solarExtra"]]),limits=c(-1.5,1.5),oob = scales::squish,breaks=c(-1.5,0,1.5),guide = guide_colorbar(frame.colour = "black",frame.linewidth= 1,title.position="top", title.hjust = 0.5)) + 
    labs(title=show_label[i],x="UMAP 1",y="UMAP 2",color="ChromVAR z-score")+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background=element_rect(fill='transparent', color='black',linetype="solid",size=2.5),plot.title = element_text(hjust = 0.3,size=rel(3),face="bold"),axis.text=element_blank(),axis.ticks=element_blank(),axis.title.x=element_text(size=rel(1.5),vjust=0.5),axis.title.y=element_text(size=rel(1.5),vjust=0.5),legend.direction="horizontal",legend.position=c(0.8,0.1),legend.key.size = unit(0.65, 'cm'),legend.title=element_text(size=rel(1.1))) 
  p
  })
DefaultAssay(filtered_mature_ORNs_scATAC)<- "chromvar"
# Cd36_vs_other_differential.activity<-  FindMarkers(
#   object = filtered_mature_ORNs_scATAC,
#   ident.1 = 'S_M.pulmonis',
#   ident.2 = 'other',
#   only.pos = TRUE,
#   mean.fxn = rowMeans,
#   fc.name = "avg_diff"
# )
Cd36_vs_other_differential.activity<- diff_chromvar_df
show_label_motif<- merged_df[match(show_label,merged_df$Row.names),]$Row.names.y

DefaultAssay(filtered_mature_ORNs_scATAC) <- "chromvar"
newpalette <- c(brewer.pal(9,"Greens")[5],brewer.pal(9,"Pastel1")[9])
#Cd36_vs_other_differential.activity$Row.names<- rownames(Cd36_vs_other_differential.activity)
chromvar_violin_plot.ls <- lapply(1:length(motifs),function(i){
  chromvar_df <- FetchData(filtered_mature_ORNs_scATAC,vars=c(motifs[i],"S_M.pulmonis"))
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


pdf(str_c(out_dir,"S_M.pulmonis_Tcf4_Bach2_chromVAR_zscore_UMAP_VlnPlot.pdf"),width=12,,height=21)
(chromvar_UMAP_plot.ls[[1]] + chromvar_violin_plot.ls[[1]]+plot_layout(widths=c(5.5,5))) / 
(chromvar_UMAP_plot.ls[[2]] + chromvar_violin_plot.ls[[2]]+plot_layout(widths=c(5.5,5))) / 
(chromvar_UMAP_plot.ls[[3]] + chromvar_violin_plot.ls[[3]]+plot_layout(widths=c(5.5,5))) + plot_annotation(title="Accessibility of cis-elements",theme = theme(plot.title = element_text(hjust = 0.5,size = rel(2.8))))
dev.off()



