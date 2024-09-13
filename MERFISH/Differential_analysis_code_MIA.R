##This script was used to run differential gene expression using MAST for comparing MIAM vs PBSM and MIAF vs PBSF

#load libraries
library(dplyr)
library(Seurat)
library(patchwork)
library(methods)
library(Matrix)
library(MAST)
library(DT)

#read MIA seurat dataset
mia <- readRDS("mia_dataset.rds")
DefaultAssay(mia) <- "Vizgen"
mia <- JoinLayers(mia)

cell_types<- unique(mia@meta.data$celltype)
 
mia$group_celltype <-paste(mia$Condition, mia$celltype, sep="_")
des<-list()

Idents(mia)<-mia$group_celltype

#Run MAST on each celltype for ech comparison
for(celltype in cell_types){
    des[[celltype]]<-list()

    MIAM<-paste0("MIAM_", celltype)

    PBSM<-paste0("PBSM_", celltype)
    
    comparisons<-list(MIAM_PBSM=c(MIAM, PBSM))

    for(comparison in comparisons){

      name<-paste(comparison, collapse="_vs_")

      mast<-FindMarkers(mia, ident.1 = comparison[1], ident.2 = comparison[2],

                        verbose = T, test.use = "MAST", pseudocount.use = 1, logfc.threshold = 0, min.pct = 0.01,

                        slot = "counts")

      des[[celltype]][[name]]<-mast

    }
}	


# get some annatations and re-calculate the p adjust, using <0.05 for significane,

des_collapsed<-list()

for(celltype in names(des)){

  dat<-des[[celltype]]

  annotated_dat<-list()

  for(comparison in names(dat)){

    comp<-dat[[comparison]]

    comp$comparison<-gsub(paste0("_", celltype), "", comparison)

    comp$gene<-rownames(comp)

    comp$FDR<-p.adjust(comp$p_val, method = "fdr")

    comp$sig<-ifelse(comp$FDR<0.05 & abs(comp$avg_log2FC) >log(0.1, base = 2), T, F)

    comp$celltype<-celltype

    subclasses<-unlist(strsplit(x = comparison, split="_vs_"))

    cells<-subset(mia@meta.data, idents = subclasses)

    expression<-as.data.frame(log2(AverageExpression(cells, verbose = F)$Vizgen))

   expression$gene<-rownames(as.data.frame(expression))

   colnames(expression)[1:2]<-c("exp.1", "exp.2")

   expression[mapply(is.infinite, expression)]<-0

   comp<-left_join(comp, expression, by="gene")

    annotated_dat[[comparison]]<-comp

  }

  annotated_dat<-do.call("rbind", annotated_dat)

  rownames(annotated_dat)<-c(1:nrow(annotated_dat))

  des_collapsed[[celltype]]<-annotated_dat

}

des_collapsed<-do.call("rbind", des_collapsed)

rownames(des_collapsed)<-c(1:nrow(des_collapsed))
write.csv(des_collapsed, "MIA_males_DEG.csv")




######Following code was used to make volcano plots#####
#load libraries
library(gridExtra)
library(grid)
library(ggrepel)

#read DEG files of males and female comparison as csv
df1 <- read.csv("MIA_males_DEG.csv")
df2<- read.csv("MIA_females_DEG.csv")

#Code for volcano plot of Ventral_RP
df4 <- subset(df2, celltype == "Ventral_RP")
df3 <- subset(df1, celltype == "Ventral_RP")
colnames(df4)[colnames(df4) == "avg_log2FC"] <- "avg_log2FC_female"                                                                            
colnames(df3)[colnames(df3) == "avg_log2FC"] <- "avg_log2FC_male"
merged_df <- merge(df3,df4, by = 'gene', all=TRUE)

df3$gene_label <- rownames(df3)
  df3$group <- with(df3, ifelse(avg_log2FC_male > 0 & sig == TRUE, "up_sig",
                              ifelse(avg_log2FC_male > 0 & sig == FALSE, "up_non_sig",
                                     ifelse(avg_log2FC_male <= 0 & sig == TRUE, "down_sig", "down_non_sig"))))

 df3$neg_log_padj <- -log10(df3$FDR)
replacement_number <- 500
df3$neg_log_padj <- replace(df3$neg_log_padj, is.infinite(df3$neg_log_padj), replacement_number)
volcano1 <- ggplot(df3, aes(x = avg_log2FC_male, y = neg_log_padj, color = group)) +
  geom_point(size = 1.5) +  
  theme(legend.position = "none") +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black")
  ) +
  ggtitle("MIAM vs PBSM") +
  ylab("-log10(FDR)") +
  xlab("log2 Fold Change") +
  scale_color_manual(values = c("up_sig" = "#cc6864", "down_sig" = "#4B8CFF", "up_non_sig" = "gray", "down_non_sig" = "gray")) +
  guides(color = FALSE) +
  geom_hline(yintercept = -log10(0.05), alpha = 0.5, linetype = "dashed", color = "gray") +
  geom_vline(xintercept = c(-0.1, 0.1), alpha = 0.5, linetype = "dashed", color = "gray") +
  coord_cartesian(ylim = c(0, 520)) +  
  geom_text_repel(
    data = subset(df3, gene %in% c("Bcan","Cmtm7","Fgf17","Il2","Kitl","Wnt5a","Wnt7b","Ly6d","Myc","Bmp7","Il18bp","Ifnar1","Gli2","Gli3","Fgf11","Fgf13","Il6st","Kif7","Vcan","Notch2","Vegfa","Notch1","Notch3","Mki67","Sox8","Tmeff1","Lrp6","Igfbp2","Ncan","Fzd2","Ptch1","Apc","Wnt7a","Ifngr1")),  # Highlight specific genes
    aes(label = gene),
    color = "black",
    box.padding = unit(.7, "lines"), max.overlaps=Inf,
    hjust= 0.30
  )   

df4$gene_label <- rownames(df4)
  df4$group <- with(df4, ifelse(avg_log2FC_female > 0 & sig == TRUE, "up_sig",
                              ifelse(avg_log2FC_female > 0 & sig == FALSE, "up_non_sig",
                                     ifelse(avg_log2FC_female <= 0 & sig == TRUE, "down_sig", "down_non_sig"))))

 df4$neg_log_padj <- -log10(df4$FDR)
replacement_number <- 500
df4$neg_log_padj <- replace(df4$neg_log_padj, is.infinite(df4$neg_log_padj), replacement_number)
 
volcano2 <- ggplot(df4, aes(x = avg_log2FC_female, y = neg_log_padj, color = group)) +
  geom_point(size = 1.5) +  
  theme(legend.position = "none") +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black")
  ) +
  ggtitle("MIAF vs PBSF") +
  ylab("-log10(FDR)") +
  xlab("log2 Fold Change") +
  scale_color_manual(values = c("up_sig" = "#ccb8e8", "down_sig" = "#1e968a", "up_non_sig" = "gray", "down_non_sig" = "gray")) +
  guides(color = FALSE) +
  geom_hline(yintercept = -log10(0.05), alpha = 0.5, linetype = "dashed", color = "gray") +
  geom_vline(xintercept = c(-0.1, 0.1), alpha = 0.5, linetype = "dashed", color = "gray") +
  coord_cartesian(ylim = c(0, 520)) +  
  geom_text_repel(
    data = subset(df4, gene %in% c("Bcan","Fgf17","Il2","Kitl","Wnt5a","Wnt7b","Ly6d","Myc","Bmp7","Il1f10","Il18bp","Ifnar1","Gli2","Fgf11","Fgf13","Il6st","Kif7","Vcan","Notch2","Vegfa","Notch1","Notch3","Mki67","Tmeff1","Lrp6","Igfbp2","Fgfr1","Gli1","Ncan","Fzd2","Ptch1","Fgf15", "Apc","Wnt7a","Ifngr1")),  # Highlight specific genes
    aes(label = gene),
    color = "black",
    box.padding = unit(.7, "lines"),max.overlaps=Inf,
    hjust= 0.30
  )  
pdf("mia_VRP_volcano_plot.pdf", width = 10)
grid.arrange(volcano1, volcano2, ncol = 2) 
dev.off()

########################################
#Code for volcano plot of Dorsal RP

df4 <- subset(df2, celltype == "Dorsal_RP")
df3 <- subset(df1, celltype == "Dorsal_RP")
colnames(df4)[colnames(df4) == "avg_log2FC"] <- "avg_log2FC_female"                                                                            
colnames(df3)[colnames(df3) == "avg_log2FC"] <- "avg_log2FC_male"
merged_df <- merge(df3,df4, by = 'gene', all=TRUE)

df3$gene_label <- rownames(df3)
  df3$group <- with(df3, ifelse(avg_log2FC_male > 0 & sig == TRUE, "up_sig",
                              ifelse(avg_log2FC_male > 0 & sig == FALSE, "up_non_sig",
                                     ifelse(avg_log2FC_male <= 0 & sig == TRUE, "down_sig", "down_non_sig"))))

 df3$neg_log_padj <- -log10(df3$FDR)
replacement_number <- 500
df3$neg_log_padj <- replace(df3$neg_log_padj, is.infinite(df3$neg_log_padj), replacement_number)
volcano1 <- ggplot(df3, aes(x = avg_log2FC_male, y = neg_log_padj, color = group)) +
  geom_point(size = 1.5) +  
  theme(legend.position = "none") +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black")
  ) +
  ggtitle("MIAM vs PBSM") +
  ylab("-log10(FDR)") +
  xlab("log2 Fold Change") +
  scale_color_manual(values = c("up_sig" = "#cc6864", "down_sig" = "#4B8CFF", "up_non_sig" = "gray", "down_non_sig" = "gray")) +
  guides(color = FALSE) +
  geom_hline(yintercept = -log10(0.05), alpha = 0.5, linetype = "dashed", color = "gray") +
  geom_vline(xintercept = c(-0.1, 0.1), alpha = 0.5, linetype = "dashed", color = "gray") +
  coord_cartesian(ylim = c(0, 220)) +  
  geom_text_repel(
    data = subset(df3, gene %in% c("Bcan","Cmtm6","Fgf11","Il20rb","Il12b","Ifngr1","Fgf13","Bmp7","Notch3","Il18bp","Ifnar1","Gli2","Gli3","Wnt5a","Vcan","Notch2","Arx","Ctnnb1","Mki67","Nes","Txlna","Vegfa","Sirpa","Csnk1a1","Neurog2","Fzd2","Il1f10","Igfbp2","Il2","Tnfrsf4","Axin1","Ackr3","Ifnar2","Fgfr2")),  # Highlight specific genes
    aes(label = gene),
    color = "black",
    box.padding = unit(.7, "lines"), max.overlaps=Inf,
    hjust= 0.30
  )   

df4$gene_label <- rownames(df4)
  df4$group <- with(df4, ifelse(avg_log2FC_female > 0 & sig == TRUE, "up_sig",
                              ifelse(avg_log2FC_female > 0 & sig == FALSE, "up_non_sig",
                                     ifelse(avg_log2FC_female <= 0 & sig == TRUE, "down_sig", "down_non_sig"))))

 df4$neg_log_padj <- -log10(df4$FDR)
replacement_number <- 500
df4$neg_log_padj <- replace(df4$neg_log_padj, is.infinite(df4$neg_log_padj), replacement_number)
 
volcano2 <- ggplot(df4, aes(x = avg_log2FC_female, y = neg_log_padj, color = group)) +
  geom_point(size = 1.5) +  
  theme(legend.position = "none") +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black")
  ) +
  ggtitle("MIAF vs PBSF") +
  ylab("-log10(FDR)") +
  xlab("log2 Fold Change") +
  scale_color_manual(values = c("up_sig" = "#ccb8e8", "down_sig" = "#1e968a", "up_non_sig" = "gray", "down_non_sig" = "gray")) +
  guides(color = FALSE) +
  geom_hline(yintercept = -log10(0.05), alpha = 0.5, linetype = "dashed", color = "gray") +
  geom_vline(xintercept = c(-0.1, 0.1), alpha = 0.5, linetype = "dashed", color = "gray") +
  coord_cartesian(ylim = c(0, 200)) +  
  geom_text_repel(
    data = subset(df4, gene %in% c("Bcan","Cmtm6","Il20rb","Il12b","Ifngr1","Notch3","Wnt5a","Il6st","Vcan","Notch2", "Vegfa","Csnk1a1","Fgfr1","Fzd2","Il1f10","Igfbp2","Il2","Tnfrsf4","Axin1","Ackr3","Ifnar2")),  # Highlight specific genes
    aes(label = gene),
    color = "black",
    box.padding = unit(.7, "lines"),max.overlaps=Inf,
    hjust= 0.30
  )  
pdf("mia_dRP_volcano_plot.pdf", width = 10)
grid.arrange(volcano1, volcano2, ncol = 2) 
dev.off()

######################################################
#Code for volcano plot of SCPN

df4 <- subset(df2, celltype == "SCPN")
df3 <- subset(df1, celltype == "SCPN")
colnames(df4)[colnames(df4) == "avg_log2FC"] <- "avg_log2FC_female"                                                                            
colnames(df3)[colnames(df3) == "avg_log2FC"] <- "avg_log2FC_male"
merged_df <- merge(df3,df4, by = 'gene', all=TRUE)

df3$gene_label <- rownames(df3)
  df3$group <- with(df3, ifelse(avg_log2FC_male > 0 & sig == TRUE, "up_sig",
                              ifelse(avg_log2FC_male > 0 & sig == FALSE, "up_non_sig",
                                     ifelse(avg_log2FC_male <= 0 & sig == TRUE, "down_sig", "down_non_sig"))))

 df3$neg_log_padj <- -log10(df3$FDR)
replacement_number <- 500
df3$neg_log_padj <- replace(df3$neg_log_padj, is.infinite(df3$neg_log_padj), replacement_number)
volcano1 <- ggplot(df3, aes(x = avg_log2FC_male, y = neg_log_padj, color = group)) +
  geom_point(size = 1.5) +  
  theme(legend.position = "none") +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black")
  ) +
  ggtitle("MIAM vs PBSM") +
  ylab("-log10(FDR)") +
  xlab("log2 Fold Change") +
  scale_color_manual(values = c("up_sig" = "#cc6864", "down_sig" = "#4B8CFF", "up_non_sig" = "gray", "down_non_sig" = "gray")) +
  guides(color = FALSE) +
  geom_hline(yintercept = -log10(0.05), alpha = 0.5, linetype = "dashed", color = "gray") +
  geom_vline(xintercept = c(-0.1, 0.1), alpha = 0.5, linetype = "dashed", color = "gray") +
  coord_cartesian(ylim = c(0, 250)) +  
  geom_text_repel(
    data = subset(df3, gene %in% c("Ly6d","Ctnnb1","Fgfr2","Mdk","Il12b","Wnt5a","Hgf","Cntn6","Vcan","Angpt1","Il18bp","Flt3","Dll1","Cxcr4","Ccnd1","Fgf12","Fgf13","Wnt7a","Tubb3","Tmeff1","Hmgcs1","Il1b","Txlna","Vegfa","Wnt7b","Ncan","Il12b","Dlg4","Kitl","Sirpa","Igfbp2","Sparc","Fzd2","Il1f10","Kit","Ctnnb1","Cnr1","Nrtn","Il16","Fgf21","Trp53")),  # Highlight specific genes
    aes(label = gene),
    color = "black",
    box.padding = unit(.7, "lines"), max.overlaps=Inf,
    hjust= 0.30
  )   


df4$gene_label <- rownames(df4)
  df4$group <- with(df4, ifelse(avg_log2FC_female > 0 & sig == TRUE, "up_sig",
                              ifelse(avg_log2FC_female > 0 & sig == FALSE, "up_non_sig",
                                     ifelse(avg_log2FC_female <= 0 & sig == TRUE, "down_sig", "down_non_sig"))))

 df4$neg_log_padj <- -log10(df4$FDR)
replacement_number <- 500
df4$neg_log_padj <- replace(df4$neg_log_padj, is.infinite(df4$neg_log_padj), replacement_number)
 
volcano2 <- ggplot(df4, aes(x = avg_log2FC_female, y = neg_log_padj, color = group)) +
  geom_point(size = 1.5) +  
  theme(legend.position = "none") +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black")
  ) +
  ggtitle("MIAF vs PBSF") +
  ylab("-log10(FDR)") +
  xlab("log2 Fold Change") +
  scale_color_manual(values = c("up_sig" = "#ccb8e8", "down_sig" = "#1e968a", "up_non_sig" = "gray", "down_non_sig" = "gray")) +
  guides(color = FALSE) +
  geom_hline(yintercept = -log10(0.05), alpha = 0.5, linetype = "dashed", color = "gray") +
  geom_vline(xintercept = c(-0.1, 0.1), alpha = 0.5, linetype = "dashed", color = "gray") +
  coord_cartesian(ylim = c(0, 100)) +  
  geom_text_repel(
    data = subset(df4, gene %in% c("Ly6d","Mdk","Notch2","Wnt5a","Myc","Cntn6","Vcan","Angpt1","Il18bp","Ifnar1","Csnk1a1","Cxcr4","Kit","Fgf12","Ptch1","Wnt7a","Tubb3","Tmeff1","Hmgcs1","Il1b","Txlna","Vegfa","Ncan","Il12b","Dlg4","Kitl","Sirpa","Igfbp2","Sparc","Fzd2","Il1f10")),  # Highlight specific genes
    aes(label = gene),
    color = "black",
    box.padding = unit(.7, "lines"),max.overlaps=Inf,
    hjust= 0.30
  )  
pdf("mia_SCPN_diff.pdf", width = 10) 
grid.arrange(volcano1, volcano2, ncol = 2)
dev.off() 

#####################################################

#Code for volcano plot of Striatal_Inhibitory neurons
df4 <- subset(df2, celltype == "Striatal_Inh")
df3 <- subset(df1, celltype == "Striatal_Inh")
colnames(df4)[colnames(df4) == "avg_log2FC"] <- "avg_log2FC_female"                                                                            
colnames(df3)[colnames(df3) == "avg_log2FC"] <- "avg_log2FC_male"
merged_df <- merge(df3,df4, by = 'gene', all=TRUE)

df3$gene_label <- rownames(df3)
  df3$group <- with(df3, ifelse(avg_log2FC_male > 0 & sig == TRUE, "up_sig",
                              ifelse(avg_log2FC_male > 0 & sig == FALSE, "up_non_sig",
                                     ifelse(avg_log2FC_male <= 0 & sig == TRUE, "down_sig", "down_non_sig"))))

 df3$neg_log_padj <- -log10(df3$FDR)
replacement_number <- 500
df3$neg_log_padj <- replace(df3$neg_log_padj, is.infinite(df3$neg_log_padj), replacement_number)
volcano1 <- ggplot(df3, aes(x = avg_log2FC_male, y = neg_log_padj, color = group)) +
  geom_point(size = 1.5) +  
  theme(legend.position = "none") +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black")
  ) +
  ggtitle("MIAM vs PBSM") +
  ylab("-log10(FDR)") +
  xlab("log2 Fold Change") +
  scale_color_manual(values = c("up_sig" = "#cc6864", "down_sig" = "#4B8CFF", "up_non_sig" = "gray", "down_non_sig" = "gray")) +
  guides(color = FALSE) +
  geom_hline(yintercept = -log10(0.05), alpha = 0.5, linetype = "dashed", color = "gray") +
  geom_vline(xintercept = c(-0.1, 0.1), alpha = 0.5, linetype = "dashed", color = "gray") +
  coord_cartesian(ylim = c(0, 200)) +  
  geom_text_repel(
    data = subset(df3, gene %in% c("Ctnnb1","Tubb3","Arx","Slc32a1","Tmeff1","Hmgcs1","Il1b","Txlna","Vegfa","Foxj1","Hgf","Lrp6","Ncan","Hdgf","Flt3","Wnt7a","Tubb3","Il18bp","Vcan","Dlg4","Il12b","Fgf13","Ptch1","Cnr1")),  # Highlight specific genes
    aes(label = gene),
    color = "black",
    box.padding = unit(.7, "lines"), max.overlaps=Inf,
    hjust= 0.30
  )   


df4$gene_label <- rownames(df4)
  df4$group <- with(df4, ifelse(avg_log2FC_female > 0 & sig == TRUE, "up_sig",
                              ifelse(avg_log2FC_female > 0 & sig == FALSE, "up_non_sig",
                                     ifelse(avg_log2FC_female <= 0 & sig == TRUE, "down_sig", "down_non_sig"))))

 df4$neg_log_padj <- -log10(df4$FDR)
replacement_number <- 500
df4$neg_log_padj <- replace(df4$neg_log_padj, is.infinite(df4$neg_log_padj), replacement_number)
 
volcano2 <- ggplot(df4, aes(x = avg_log2FC_female, y = neg_log_padj, color = group)) +
  geom_point(size = 1.5) +  
  theme(legend.position = "none") +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black")
  ) +
  ggtitle("MIAF vs PBSF") +
  ylab("-log10(FDR)") +
  xlab("log2 Fold Change") +
  scale_color_manual(values = c("up_sig" = "#ccb8e8", "down_sig" = "#1e968a", "up_non_sig" = "gray", "down_non_sig" = "gray")) +
  guides(color = FALSE) +
  geom_hline(yintercept = -log10(0.05), alpha = 0.5, linetype = "dashed", color = "gray") +
  geom_vline(xintercept = c(-0.1, 0.1), alpha = 0.5, linetype = "dashed", color = "gray") +
  coord_cartesian(ylim = c(0, 170)) +  
  geom_text_repel(
    data = subset(df4, gene %in% c("Sparc","Il1b","Dcx","Ctnnb1","Ackr3","Ncan","Gad2","Cxcr4","Fstl1","Lrp6","Tubb3","Axin1","Il12b","Arx","Kitl","Il16","Fgfr1","Gsk3b","Apc","Fgfr2","Tgfbr1","Scg2","Ccnd1","Rgs8")),  # Highlight specific genes
    aes(label = gene),
    color = "black",
    box.padding = unit(.7, "lines"),max.overlaps=Inf,
    hjust= 0.30
  )  
pdf("mia_Straital_inh_diff.pdf", width = 10) 
grid.arrange(volcano1, volcano2, ncol = 2)
dev.off() 

##################################################

#Code for volcano plot of IPC
df4 <- subset(df2, celltype == "IPC")
df3 <- subset(df1, celltype == "IPC")
colnames(df4)[colnames(df4) == "avg_log2FC"] <- "avg_log2FC_female"                                                                            
colnames(df3)[colnames(df3) == "avg_log2FC"] <- "avg_log2FC_male"
merged_df <- merge(df3,df4, by = 'gene', all=TRUE)

df3$gene_label <- rownames(df3)
  df3$group <- with(df3, ifelse(avg_log2FC_male > 0 & sig == TRUE, "up_sig",
                              ifelse(avg_log2FC_male > 0 & sig == FALSE, "up_non_sig",
                                     ifelse(avg_log2FC_male <= 0 & sig == TRUE, "down_sig", "down_non_sig"))))

 df3$neg_log_padj <- -log10(df3$FDR)
replacement_number <- 500
df3$neg_log_padj <- replace(df3$neg_log_padj, is.infinite(df3$neg_log_padj), replacement_number)
volcano1 <- ggplot(df3, aes(x = avg_log2FC_male, y = neg_log_padj, color = group)) +
  geom_point(size = 1.5) +  
  theme(legend.position = "none") +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black")
  ) +
  ggtitle("MIAM vs PBSM") +
  ylab("-log10(FDR)") +
  xlab("log2 Fold Change") +
  scale_color_manual(values = c("up_sig" = "#cc6864", "down_sig" = "#4B8CFF", "up_non_sig" = "gray", "down_non_sig" = "gray")) +
  guides(color = FALSE) +
  geom_hline(yintercept = -log10(0.05), alpha = 0.5, linetype = "dashed", color = "gray") +
  geom_vline(xintercept = c(-0.1, 0.1), alpha = 0.5, linetype = "dashed", color = "gray") +
  coord_cartesian(ylim = c(0, 520)) +  
  geom_text_repel(
    data = subset(df3, gene %in% c("Bcan","Aldoc","Ctnnb1","Mki67","Vegfa","Fgfr2","Txlna","Trp53","Tubb3","Sirpa","Il18bp","Lrig1","Cntfr","Tgfbr1","Vcan","Fgf13","Gli3","Il12b","Fgf11","Neurod1","Gli2","Axin1","Ifngr1")),  # Highlight specific genes
    aes(label = gene),
    color = "black",
    box.padding = unit(.7, "lines"), max.overlaps=Inf,
    hjust= 0.30
  )   

df4$gene_label <- rownames(df4)
  df4$group <- with(df4, ifelse(avg_log2FC_female > 0 & sig == TRUE, "up_sig",
                              ifelse(avg_log2FC_female > 0 & sig == FALSE, "up_non_sig",
                                     ifelse(avg_log2FC_female <= 0 & sig == TRUE, "down_sig", "down_non_sig"))))

 df4$neg_log_padj <- -log10(df4$FDR)
replacement_number <- 500
df4$neg_log_padj <- replace(df4$neg_log_padj, is.infinite(df4$neg_log_padj), replacement_number)
 
volcano2 <- ggplot(df4, aes(x = avg_log2FC_female, y = neg_log_padj, color = group)) +
  geom_point(size = 1.5) +  
  theme(legend.position = "none") +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black")
  ) +
  ggtitle("MIAF vs PBSF") +
  ylab("-log10(FDR)") +
  xlab("log2 Fold Change") +
  scale_color_manual(values = c("up_sig" = "#ccb8e8", "down_sig" = "#1e968a", "up_non_sig" = "gray", "down_non_sig" = "gray")) +
  guides(color = FALSE) +
  geom_hline(yintercept = -log10(0.05), alpha = 0.5, linetype = "dashed", color = "gray") +
  geom_vline(xintercept = c(-0.1, 0.1), alpha = 0.5, linetype = "dashed", color = "gray") +
  coord_cartesian(ylim = c(0, 520)) +  
  geom_text_repel(
    data = subset(df4, gene %in% c("Bcan","Cnr1","Kit","Neurod2","Csnk1a1","Nusap1","Ctnnb1","Kif7","Wnt7a","Fzd8","Fgfr2","Vegfa","Igfbp2","Axin1","Ifnar2","Il6st","Notch2","Ifngr1","Fzd2","Fgfrl1","Wnt5a","Il2","Il13","Il1f10")),  # Highlight specific genes
    aes(label = gene),
    color = "black",
    box.padding = unit(.7, "lines"),max.overlaps=Inf,
    hjust= 0.30
  )  
pdf("mia_IPC_volcano_plot.pdf", width = 10)
grid.arrange(volcano1, volcano2, ncol = 2) 
dev.off()

#########################################
#Code for volcano plot of TAPs

df4 <- subset(df2, celltype == "TAPs")
df3 <- subset(df1, celltype == "TAPs")
colnames(df4)[colnames(df4) == "avg_log2FC"] <- "avg_log2FC_female"                                                                            
colnames(df3)[colnames(df3) == "avg_log2FC"] <- "avg_log2FC_male"
merged_df <- merge(df3,df4, by = 'gene', all=TRUE)


df3$gene_label <- rownames(df3)
  df3$group <- with(df3, ifelse(avg_log2FC_male > 0 & sig == TRUE, "up_sig",
                              ifelse(avg_log2FC_male > 0 & sig == FALSE, "up_non_sig",
                                     ifelse(avg_log2FC_male <= 0 & sig == TRUE, "down_sig", "down_non_sig"))))

 df3$neg_log_padj <- -log10(df3$FDR)
replacement_number <- 500
df3$neg_log_padj <- replace(df3$neg_log_padj, is.infinite(df3$neg_log_padj), replacement_number)
volcano1 <- ggplot(df3, aes(x = avg_log2FC_male, y = neg_log_padj, color = group)) +
  geom_point(size = 1.5) +  
  theme(legend.position = "none") +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black")
  ) +
  ggtitle("MIAM vs PBSM") +
  ylab("-log10(FDR)") +
  xlab("log2 Fold Change") +
  scale_color_manual(values = c("up_sig" = "#cc6864", "down_sig" = "#4B8CFF", "up_non_sig" = "gray", "down_non_sig" = "gray")) +
  guides(color = FALSE) +
  geom_hline(yintercept = -log10(0.05), alpha = 0.5, linetype = "dashed", color = "gray") +
  geom_vline(xintercept = c(-0.1, 0.1), alpha = 0.5, linetype = "dashed", color = "gray") +
  coord_cartesian(ylim = c(0, 520)) +  
  geom_text_repel(
    data = subset(df3, gene %in% c("Flt3","Il12b","Ly6d","Ifnlr1","Fgf13","Trp53","Ifngr1","Nusap1","Fgf11","Ackr3","Il13","Il6st","Gli3","Gli1","Fgfrl1","Il18","Fzd1","Gli2","Il1f10","Fzd2","Axin2","Fgfr1","Tgfbr1","Lef1","Txlna","Dlg4","Fgfr2","Il18bp","Vegfa")),  # Highlight specific genes
    aes(label = gene),
    color = "black",
    box.padding = unit(.7, "lines"), max.overlaps=Inf,
    hjust= 0.30
  )   

df4$gene_label <- rownames(df4)
  df4$group <- with(df4, ifelse(avg_log2FC_female > 0 & sig == TRUE, "up_sig",
                              ifelse(avg_log2FC_female > 0 & sig == FALSE, "up_non_sig",
                                     ifelse(avg_log2FC_female <= 0 & sig == TRUE, "down_sig", "down_non_sig"))))

 df4$neg_log_padj <- -log10(df4$FDR)
replacement_number <- 500
df4$neg_log_padj <- replace(df4$neg_log_padj, is.infinite(df4$neg_log_padj), replacement_number)
 
volcano2 <- ggplot(df4, aes(x = avg_log2FC_female, y = neg_log_padj, color = group)) +
  geom_point(size = 1.5) +  
  theme(legend.position = "none") +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black")
  ) +
  ggtitle("MIAF vs PBSF") +
  ylab("-log10(FDR)") +
  xlab("log2 Fold Change") +
  scale_color_manual(values = c("up_sig" = "#ccb8e8", "down_sig" = "#1e968a", "up_non_sig" = "gray", "down_non_sig" = "gray")) +
  guides(color = FALSE) +
  geom_hline(yintercept = -log10(0.05), alpha = 0.5, linetype = "dashed", color = "gray") +
  geom_vline(xintercept = c(-0.1, 0.1), alpha = 0.5, linetype = "dashed", color = "gray") +
  coord_cartesian(ylim = c(0, 520)) +  
  geom_text_repel(
    data = subset(df4, gene %in% c("Gli1","Fzd8","Fgf15","Gli2","Fzd5","Il20rb","Fgf13","Wnt7a","Gli3","Apc","Kif7","Ctnnb1","Fzd1","Ifnar1","Fgfr2","Igfbp2","Angptl2","Rgs8","Tnfrsf4","Kitl","Il2","Ifnar2","Il1f10","Cdkn1a","Wnt5a","Il18","Il13","Wnt7b","Bcan")),  # Highlight specific genes
    aes(label = gene),
    color = "black",
    box.padding = unit(.7, "lines"),max.overlaps=Inf,
    hjust= 0.30
  )  
pdf("mia_TAPS_volcano_plot.pdf", width = 10)
grid.arrange(volcano1, volcano2, ncol = 2) 
dev.off()

######################################

#Code for volcano plot of Inhibitory_NB

df4 <- subset(df2, celltype == "Inhibitory_NB")
df3 <- subset(df1, celltype == "Inhibitory_NB")
colnames(df4)[colnames(df4) == "avg_log2FC"] <- "avg_log2FC_female"                                                                            
colnames(df3)[colnames(df3) == "avg_log2FC"] <- "avg_log2FC_male"
merged_df <- merge(df3,df4, by = 'gene', all=TRUE)

df3$gene_label <- rownames(df3)
  df3$group <- with(df3, ifelse(avg_log2FC_male > 0 & sig == TRUE, "up_sig",
                              ifelse(avg_log2FC_male > 0 & sig == FALSE, "up_non_sig",
                                     ifelse(avg_log2FC_male <= 0 & sig == TRUE, "down_sig", "down_non_sig"))))

 df3$neg_log_padj <- -log10(df3$FDR)
replacement_number <- 500
df3$neg_log_padj <- replace(df3$neg_log_padj, is.infinite(df3$neg_log_padj), replacement_number)
volcano1 <- ggplot(df3, aes(x = avg_log2FC_male, y = neg_log_padj, color = group)) +
  geom_point(size = 1.5) +  
  theme(legend.position = "none") +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black")
  ) +
  ggtitle("MIAM vs PBSM") +
  ylab("-log10(FDR)") +
  xlab("log2 Fold Change") +
  scale_color_manual(values = c("up_sig" = "#cc6864", "down_sig" = "#4B8CFF", "up_non_sig" = "gray", "down_non_sig" = "gray")) +
  guides(color = FALSE) +
  geom_hline(yintercept = -log10(0.05), alpha = 0.5, linetype = "dashed", color = "gray") +
  geom_vline(xintercept = c(-0.1, 0.1), alpha = 0.5, linetype = "dashed", color = "gray") +
  coord_cartesian(ylim = c(0, 520)) +  
  geom_text_repel(
    data = subset(df3, gene %in% c("Hgf","Il12b","Flt3","Il1b","Fgf13","Trp53","Ctnnb1","Cnr1","Fgf11","Nusap1","Wnt7a","Vegfb","Apc","Angpt1","Lrp6","Il1f10","Axin1","Il6st","Fzd2","Fzd8","Txlna","Lrig1","Il18bp","Fgfr2","Vegfa")),  # Highlight specific genes
    aes(label = gene),
    color = "black",
    box.padding = unit(.7, "lines"), max.overlaps=Inf,
    hjust= 0.30
  )   

df4$gene_label <- rownames(df4)
  df4$group <- with(df4, ifelse(avg_log2FC_female > 0 & sig == TRUE, "up_sig",
                              ifelse(avg_log2FC_female > 0 & sig == FALSE, "up_non_sig",
                                     ifelse(avg_log2FC_female <= 0 & sig == TRUE, "down_sig", "down_non_sig"))))

 df4$neg_log_padj <- -log10(df4$FDR)
replacement_number <- 500
df4$neg_log_padj <- replace(df4$neg_log_padj, is.infinite(df4$neg_log_padj), replacement_number)
 
volcano2 <- ggplot(df4, aes(x = avg_log2FC_female, y = neg_log_padj, color = group)) +
  geom_point(size = 1.5) +  
  theme(legend.position = "none") +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black")
  ) +
  ggtitle("MIAF vs PBSF") +
  ylab("-log10(FDR)") +
  xlab("log2 Fold Change") +
  scale_color_manual(values = c("up_sig" = "#ccb8e8", "down_sig" = "#1e968a", "up_non_sig" = "gray", "down_non_sig" = "gray")) +
  guides(color = FALSE) +
  geom_hline(yintercept = -log10(0.05), alpha = 0.5, linetype = "dashed", color = "gray") +
  geom_vline(xintercept = c(-0.1, 0.1), alpha = 0.5, linetype = "dashed", color = "gray") +
  coord_cartesian(ylim = c(0, 520)) +  
  geom_text_repel(
    data = subset(df4, gene %in% c("Angpt1","Cnr1","Fzd8","Wnt7a","Fgfbp3","Ncan","Tgfbr1","Lrp6","Trp53","Il6st","Fgfr2","Angptl2","Fzd2","Igfbp2","Ifnar2","Kif7","Fgfr1","Axin1","Fgfrl1","Kit","Notch2","Il1b","Il18","Il1f10","Il2","Fzd1")),  # Highlight specific genes
    aes(label = gene),
    color = "black",
    box.padding = unit(.7, "lines"),max.overlaps=Inf,
    hjust= 0.30
  )  
pdf("mia_Inhibitory_NB_volcano_plot.pdf", width = 10)
grid.arrange(volcano1, volcano2, ncol = 2) 
dev.off()


###########################################


#Code for volcano plot of MGE_InhNB
df4 <- subset(df2, celltype == "MGE_InhNB")
df3 <- subset(df1, celltype == "MGE_InhNB")
colnames(df4)[colnames(df4) == "avg_log2FC"] <- "avg_log2FC_female"                                                                            
colnames(df3)[colnames(df3) == "avg_log2FC"] <- "avg_log2FC_male"
merged_df <- merge(df3,df4, by = 'gene', all=TRUE)


df3$gene_label <- rownames(df3)
  df3$group <- with(df3, ifelse(avg_log2FC_male > 0 & sig == TRUE, "up_sig",
                              ifelse(avg_log2FC_male > 0 & sig == FALSE, "up_non_sig",
                                     ifelse(avg_log2FC_male <= 0 & sig == TRUE, "down_sig", "down_non_sig"))))

 df3$neg_log_padj <- -log10(df3$FDR)
replacement_number <- 500
df3$neg_log_padj <- replace(df3$neg_log_padj, is.infinite(df3$neg_log_padj), replacement_number)
volcano1 <- ggplot(df3, aes(x = avg_log2FC_male, y = neg_log_padj, color = group)) +
  geom_point(size = 1.5) +  
  theme(legend.position = "none") +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black")
  ) +
  ggtitle("MIAM vs PBSM") +
  ylab("-log10(FDR)") +
  xlab("log2 Fold Change") +
  scale_color_manual(values = c("up_sig" = "#cc6864", "down_sig" = "#4B8CFF", "up_non_sig" = "gray", "down_non_sig" = "gray")) +
  guides(color = FALSE) +
  geom_hline(yintercept = -log10(0.05), alpha = 0.5, linetype = "dashed", color = "gray") +
  geom_vline(xintercept = c(-0.1, 0.1), alpha = 0.5, linetype = "dashed", color = "gray") +
  coord_cartesian(ylim = c(0, 520)) +  
  geom_text_repel(
    data = subset(df3, gene %in% c("Il12b","Il1b","Hgf","Fgf13","Ctnnb1","Tmeff1","Trp53","Wnt7b","Wnt7a","Fgf11","Apc","Il18bp","Txlna","Vegfa","Igf2")),  # Highlight specific genes
    aes(label = gene),
    color = "black",
    box.padding = unit(.7, "lines"), max.overlaps=Inf,
    hjust= 0.30
  )   

df4$gene_label <- rownames(df4)
  df4$group <- with(df4, ifelse(avg_log2FC_female > 0 & sig == TRUE, "up_sig",
                              ifelse(avg_log2FC_female > 0 & sig == FALSE, "up_non_sig",
                                     ifelse(avg_log2FC_female <= 0 & sig == TRUE, "down_sig", "down_non_sig"))))

 df4$neg_log_padj <- -log10(df4$FDR)
replacement_number <- 500
df4$neg_log_padj <- replace(df4$neg_log_padj, is.infinite(df4$neg_log_padj), replacement_number)
 
volcano2 <- ggplot(df4, aes(x = avg_log2FC_female, y = neg_log_padj, color = group)) +
  geom_point(size = 1.5) +  
  theme(legend.position = "none") +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black")
  ) +
  ggtitle("MIAF vs PBSF") +
  ylab("-log10(FDR)") +
  xlab("log2 Fold Change") +
  scale_color_manual(values = c("up_sig" = "#ccb8e8", "down_sig" = "#1e968a", "up_non_sig" = "gray", "down_non_sig" = "gray")) +
  guides(color = FALSE) +
  geom_hline(yintercept = -log10(0.05), alpha = 0.5, linetype = "dashed", color = "gray") +
  geom_vline(xintercept = c(-0.1, 0.1), alpha = 0.5, linetype = "dashed", color = "gray") +
  coord_cartesian(ylim = c(0, 520)) +  
  geom_text_repel(
    data = subset(df4, gene %in% c("Fgf15","Angpt1","Fzd8","Nrg1","Ccnd1","Jag2","Cnr1","Hmgcs1","Cntfr","Vcan","Ntrk3","Nfkb1","Lmo4","Wnt7b","Vegfa","Ctnnb1","Fgfr2","Fgfr1","Fzd2","Il12b","Fzd1","Ifnar2","Il1b","Fstl1","Lrp5","H2afx","Angptl2","Kif7","Fgfrl1","Kit","Il1f10","Fgf21","Il18","Flt3","Tnfsf9","Tgfb3","Kitl","Bcan")),  # Highlight specific genes
    aes(label = gene),
    color = "black",
    box.padding = unit(.7, "lines"),max.overlaps=Inf,
    hjust= 0.30
  )  
pdf("mia_MGE_NB_volcano_plot.pdf", width = 10)
grid.arrange(volcano1, volcano2, ncol = 2) 
dev.off()


####code to make dot plot fig 7a on filtered DEG cav having expression and FDR values for Ackr3, Cxcl12 and Cxcr4

gene <- c("Ackr3", "Cxcl12", "Cxcr4")
markers <- mia$gene %>% unique()
mia <- mia %>%
  filter(gene %in% markers) %>%
  mutate(`-log10(padj)` = -log10(p_val_adj),
         comparison = ifelse(comparison == "MIAM_vs_MIAF", "MIAM_vs_MIAF",
                            ifelse(comparison == "MIAF_vs_PBSF", "MIAF_vs_PBSF", "MIAM_vs_PBSM")))

ggplot(mia, aes(x = celltype, y = gene, size = `-log10(padj)`)) +
  geom_point(aes(color = avg_log2FC)) +
  scale_color_viridis(option = "magma") +
  scale_size_continuous(range = c(1, 10), guide = guide_legend(title = "-log10(padj)")) +
  facet_wrap(~ comparison, scales = "free_y", ncol = 1) +
  theme_minimal() +
  theme(axis.line = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.ticks = element_blank(),
        plot.title = element_text(hjust = 0.5),
        strip.text = element_text(face = "bold", size = 8),
        strip.background = element_blank(),
        panel.spacing = unit(0, "lines")) +
  ylab('') +
  guides(size = guide_legend(title = "-log10(padj)"))
  
