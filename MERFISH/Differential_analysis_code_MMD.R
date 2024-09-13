##This script was used to run differential gene expression using MAST for comparing MMD vs PBS, MMDM vs PBSM, MMDF vs PBSF

#!/usr/bin/Rscript
library(dplyr)
library(Seurat)
library(patchwork)
library(methods)
library(Matrix)
library(MAST)
library(DT)


MMD <- readRDS("MMD_dataset.rds")
DefaultAssay(MMD) <- "Vizgen"
MMD <- JoinLayers(MMD)

cell_types<- unique(MMD@meta.data$broad_celltype)
 
MMD$group_celltype <-paste(MMD$Condition, MMD$broad_celltype, sep="_")
des<-list()

Idents(MMD)<-MMD$group_celltype

#Run MAST on each celltype for ech comparison
for(celltype in cell_types){
    des[[celltype]]<-list()

    MIAM<-paste0("MMD_", celltype)

    PBSM<-paste0("PBS_", celltype)
    
    comparisons<-list(MMD_PBS=c(MMD, PBS))

    for(comparison in comparisons){

      name<-paste(comparison, collapse="_vs_")

      mast<-FindMarkers(MMD, ident.1 = comparison[1], ident.2 = comparison[2],

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

    cells<-subset(MMD@meta.data, idents = subclasses)

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
write.csv(des_collapsed, "MMD_all_DEG.csv")


######Following code was used to make volcano plots#####
#load libraries
library(gridExtra)
library(grid)
library(ggrepel)
#read DEG files of males and female comparison as csv
df1 <- read.csv("MMD_all_DEG.csv")
#Code for volcano plot of Ventral_RP
df3 <- subset(df1, celltype == "Ventral_RP")
                                                                         
df3$gene_label <- rownames(df3)
  df3$group <- with(df3, ifelse(avg_log2FC > 0 & sig == TRUE, "up_sig",
                              ifelse(avg_log2FC > 0 & sig == FALSE, "up_non_sig",
                                     ifelse(avg_log2FC <= 0 & sig == TRUE, "down_sig", "down_non_sig"))))

 df3$neg_log_padj <- -log10(df3$FDR)
replacement_number <- 500

df3$neg_log_padj <- replace(df3$neg_log_padj, is.infinite(df3$neg_log_padj), replacement_number)
#Plot volcano
volcano1 <- ggplot(df3, aes(x = avg_log2FC, y = neg_log_padj, color = group)) +
  geom_point(size = 1.5) +  
  theme(legend.position = "none") +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black")
  ) +
  ggtitle("MMD vs PBS") +
  ylab("-log10(FDR)") +
  xlab("log2 Fold Change") +
  scale_color_manual(values = c("up_sig" = "#F39B7F", "down_sig" = "#4DBBD5", "up_non_sig" = "gray", "down_non_sig" = "gray")) +
  guides(color = FALSE) +
  geom_hline(yintercept = -log10(0.05), alpha = 0.5, linetype = "dashed", color = "gray") +
  geom_vline(xintercept = c(-0.1, 0.1), alpha = 0.5, linetype = "dashed", color = "gray") +
  coord_cartesian(ylim = c(0, 520)) +  
  geom_text_repel(
    data = subset(df3, gene %in% c("Kitl","Lrp4","Fzd8","Arx","Lrp5","Kif7","Tgfbr1","Csnk1a1","Gmfb","Hmgcs1","Fgf15","Tle4","Fgfr2","Ccnd1","Tmeff1","Ctnnb1","Mki67","Ncan","Tnc","Ackr3","Sparc","Vcan","Fstl1","Ptch1","Lrp6","Lamp1","Dlx2","Hdgf","Notch2","Nusap1","Apc")),  # Highlight specific genes
    aes(label = gene),
    color = "black",
    box.padding = unit(.7, "lines"), max.overlaps=Inf,
    hjust= 0.30
  )
pdf("MMD_PBS_volcanoplot_ventral_RP.pdf") 
volcano1  
dev.off()


#######################################

#Code for volcano plot of Dorsal_RP
df3 <- subset(df1, celltype == "Dorsal_RP")
df3$gene_label <- rownames(df3)
  df3$group <- with(df3, ifelse(avg_log2FC > 0 & sig == TRUE, "up_sig",
                              ifelse(avg_log2FC > 0 & sig == FALSE, "up_non_sig",
                                     ifelse(avg_log2FC <= 0 & sig == TRUE, "down_sig", "down_non_sig"))))

 df3$neg_log_padj <- -log10(df3$FDR)
replacement_number <- 500
df3$neg_log_padj <- replace(df3$neg_log_padj, is.infinite(df3$neg_log_padj), replacement_number)
#plot volcano
volcano1 <- ggplot(df3, aes(x = avg_log2FC, y = neg_log_padj, color = group)) +
  geom_point(size = 1.5) +  
  theme(legend.position = "none") +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black")
  ) +
  ggtitle("MMD vs PBS") +
  ylab("-log10(FDR)") +
  xlab("log2 Fold Change") +
  scale_color_manual(values = c("up_sig" = "#F39B7F", "down_sig" = "#4DBBD5", "up_non_sig" = "gray", "down_non_sig" = "gray")) +
  guides(color = FALSE) +
  geom_hline(yintercept = -log10(0.05), alpha = 0.5, linetype = "dashed", color = "gray") +
  geom_vline(xintercept = c(-0.1, 0.1), alpha = 0.5, linetype = "dashed", color = "gray") +
  coord_cartesian(ylim = c(0, 520)) +  
  geom_text_repel(
    data = subset(df3, gene %in% c("Tgfb2","Wnt7b","Tgfbr1","Cxcr4","Gmfb","Vcan","Fn1","Lef1","Vegfc","Hes1","Tle4","Axin2","Fgfr2","Ccnd1","Tmeff1","Nes","Ncan","Fzd8","Tnc","Lrig1","Neurog2","Notch1","Neurod2","Sstr2","Lamp1","Dvl1","Ackr3","Fgfr1","Notch2","Apc","Il18bp","Foxp1","Fzd1","Fstl1","Angptl2")),  # Highlight specific genes
    aes(label = gene),
    color = "black",
    box.padding = unit(.7, "lines"), max.overlaps=Inf,
    hjust= 0.30
  )
pdf("MMD_PBS_volcanoplot_dorsal_RP.pdf")  
volcano1
dev.off()
############################################

#Code for volcano plot of IPC
pdf("MMD_PBS_volcanoplot_IPC.pdf")
df3 <- subset(df1, celltype == "IPC")
df3$gene_label <- rownames(df3)
  df3$group <- with(df3, ifelse(avg_log2FC > 0 & sig == TRUE, "up_sig",
                              ifelse(avg_log2FC > 0 & sig == FALSE, "up_non_sig",
                                     ifelse(avg_log2FC <= 0 & sig == TRUE, "down_sig", "down_non_sig"))))

 df3$neg_log_padj <- -log10(df3$FDR)
replacement_number <- 500
df3$neg_log_padj <- replace(df3$neg_log_padj, is.infinite(df3$neg_log_padj), replacement_number)
#plot volcano
volcano1 <- ggplot(df3, aes(x = avg_log2FC, y = neg_log_padj, color = group)) +
  geom_point(size = 1.5) +  
  theme(legend.position = "none") +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black")
  ) +
  ggtitle("MMD vs PBS") +
  ylab("-log10(FDR)") +
  xlab("log2 Fold Change") +
  scale_color_manual(values = c("up_sig" = "#F39B7F", "down_sig" = "#4DBBD5", "up_non_sig" = "gray", "down_non_sig" = "gray")) +
  guides(color = FALSE) +
  geom_hline(yintercept = -log10(0.05), alpha = 0.5, linetype = "dashed", color = "gray") +
  geom_vline(xintercept = c(-0.1, 0.1), alpha = 0.5, linetype = "dashed", color = "gray") +
  coord_cartesian(ylim = c(0, 520)) +  
  geom_text_repel(
    data = subset(df3, gene %in% c("Tgfbr1","Fgfr2","Hmgcs1","Ccnd1","Ctnnb1","Ncan","Vcan","Mki67","Notch2","Gli3","Ptch1","Fgf13","Lef1","Fzd8","Lrp6","Ptn","Dvl1","Apc","Nes","Fgfrl1","Nusap1","Fzd2","Fgfbp3","Arx","Fgfr1","Lrp5","Il18bp","Angptl2","Gli2","Fzd7","Ifnar1","Il2","Ltbr","Lifr")),  # Highlight specific genes
    aes(label = gene),
    color = "black",
    box.padding = unit(.7, "lines"), max.overlaps=Inf,
    hjust= 0.30
  )
volcano1
dev.off()


####################################################
#Code for volcano plot of TAPs
pdf("MMD_PBS_volcanoplot_Taps.pdf")
df3 <- subset(df1, celltype == "TAPs")
df3$gene_label <- rownames(df3)
  df3$group <- with(df3, ifelse(avg_log2FC > 0 & sig == TRUE, "up_sig",
                              ifelse(avg_log2FC > 0 & sig == FALSE, "up_non_sig",
                                     ifelse(avg_log2FC <= 0 & sig == TRUE, "down_sig", "down_non_sig"))))

 df3$neg_log_padj <- -log10(df3$FDR)
replacement_number <- 500
df3$neg_log_padj <- replace(df3$neg_log_padj, is.infinite(df3$neg_log_padj), replacement_number)
#plot volcano
volcano1 <- ggplot(df3, aes(x = avg_log2FC, y = neg_log_padj, color = group)) +
  geom_point(size = 1.5) +  
  theme(legend.position = "none") +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black")
  ) +
  ggtitle("MMD vs PBS") +
  ylab("-log10(FDR)") +
  xlab("log2 Fold Change") +
  scale_color_manual(values = c("up_sig" = "#F39B7F", "down_sig" = "#4DBBD5", "up_non_sig" = "gray", "down_non_sig" = "gray")) +
  guides(color = FALSE) +
  geom_hline(yintercept = -log10(0.05), alpha = 0.5, linetype = "dashed", color = "gray") +
  geom_vline(xintercept = c(-0.1, 0.1), alpha = 0.5, linetype = "dashed", color = "gray") +
  coord_cartesian(ylim = c(0, 520)) +  
  geom_text_repel(
    data = subset(df3, gene %in% c("Tgfbr1","Fgfr2","Hmgcs1","Ccnd1","Ctnnb1","Ncan","Vcan","Mki67","Notch2","Gli3","Ptch1","Fgf13","Lef1","Lrp6","Ptn","Dvl1","Apc","Nes","Fgfrl1","Nusap1","Fzd2","Arx","Fgfr1","Lrp5","Il18bp","Angptl2","Gli2","Fzd7","Il2","Ltbr","Lifr")),  # Highlight specific genes
    aes(label = gene),
    color = "black",
    box.padding = unit(.7, "lines"), max.overlaps=Inf,
    hjust= 0.30
  )
volcano1
dev.off()

####code to make dot plot fig 7a on filtered DEG cav having expression and FDR values for Ackr3, Cxcl12 and Cxcr4

gene <- c("Ackr3", "Cxcl12", "Cxcr4")
markers <- mia$gene %>% unique()
MMD <- MMD %>%
  filter(gene %in% markers) %>%
  mutate(`-log10(padj)` = -log10(p_val_adj),
         comparison = ifelse(comparison == "MMDM_vs_MMDF", "MMDM_vs_MMDF",
                            ifelse(comparison == "MMDF_vs_PBSF", "MMDF_vs_PBSF", "MMDM_vs_PBSM")))

ggplot(MMD, aes(x = celltype, y = gene, size = `-log10(padj)`)) +
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
  