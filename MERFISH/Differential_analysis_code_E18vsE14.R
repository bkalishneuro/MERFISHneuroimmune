##This script was used to run differential gene expression using MAST for comparing MIAM vs PBSM and MIAF vs PBSF

#!/usr/bin/Rscript
library(dplyr)
library(Seurat)
library(patchwork)
library(methods)
library(Matrix)
library(MAST)
library(DT)

#read developmental dataset (merged seurat object for E18 and E14 males PBS) 
dev <- readRDS("dev_dataset.rds")
DefaultAssay(dev) <- "Vizgen"
dev <- JoinLayers(dev)

cell_types<- unique(dev@meta.data$broad_celltype)
 
dev$group_celltype <-paste(dev$Condition, dev$broad_celltype, sep="_")
des<-list()

Idents(dev)<-dev$group_celltype

#Run MAST on each celltype for ech comparison
for(celltype in cell_types){
    des[[celltype]]<-list()

    MIAM<-paste0("E18M_", celltype)

    PBSM<-paste0("E14M_", celltype)
    
    comparisons<-list(E18M_E14M=c(E18M, E14M))

    for(comparison in comparisons){

      name<-paste(comparison, collapse="_vs_")

      mast<-FindMarkers(dev, ident.1 = comparison[1], ident.2 = comparison[2],

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

    cells<-subset(dev@meta.data, idents = subclasses)

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
write.csv(des_collapsed, "E14_males_DEG.csv")




######Following code was used to make volcano plots#####
#load libraries
library(gridExtra)
library(grid)
library(ggrepel)

df <- read.csv("E18_E14_cell_expression_spreadsheet/E18_vs_E14_NB.csv")

##
df$gene_label <- rownames(df)
  df$group <- with(df, ifelse(avg_log2FC > 0 & sig == TRUE, "up_sig",
                              ifelse(avg_log2FC > 0 & sig == FALSE, "up_non_sig",
                                     ifelse(avg_log2FC <= 0 & sig == TRUE, "down_sig", "down_non_sig"))))

 df$neg_log_padj <- -log10(df$FDR)
replacement_number <- 500
df$neg_log_padj <- replace(df$neg_log_padj, is.infinite(df$neg_log_padj), replacement_number)
volcano <- ggplot(df, aes(x = avg_log2FC, y = neg_log_padj, color = group)) +
  geom_point(size = 1.5) +  
  theme(legend.position = "none") +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black")
  ) +
  ggtitle("Cell") +
  ylab("-log10(FDR)") +
  xlab("log2 Fold Change") +
  scale_color_manual(values = c("up_sig" = "#1e968a", "down_sig" = "#fccd56", "up_non_sig" = "gray", "down_non_sig" = "gray")) +
  guides(color = FALSE) +
  geom_hline(yintercept = -log10(0.05), alpha = 0.5, linetype = "dashed", color = "gray") +
  geom_vline(xintercept = c(-0.1, 0.1), alpha = 0.5, linetype = "dashed", color = "gray") +
  coord_cartesian(ylim = c(0, 520)) +  
  geom_text_repel(
    data = subset(df, gene %in% c("Igf2","Csf1", "Igfbp2","Dll1", "Il12b","Il13", "Il16", "Il17a", "Il17a","Il17ra", "Il18","DII1","Cxcl12", "Fgf21", "Bmp2", "Gli2")),  # Highlight specific genes
    aes(label = gene),
    color = "black",
    box.padding = unit(.7, "lines"),
    hjust= 0.30
  )  
  
volcano 
###code for volcano plot of IPC

df <- read.csv("E18_E14_cell_expression_spreadsheet/E18_vs_E14_IPC.csv")


df$gene_label <- rownames(df)
  df$group <- with(df, ifelse(avg_log2FC > 0 & sig == TRUE, "up_sig",
                              ifelse(avg_log2FC > 0 & sig == FALSE, "up_non_sig",
                                     ifelse(avg_log2FC <= 0 & sig == TRUE, "down_sig", "down_non_sig"))))

 df$neg_log_padj <- -log10(df$FDR)
replacement_number <- 500
df$neg_log_padj <- replace(df$neg_log_padj, is.infinite(df$neg_log_padj), replacement_number)
volcano <- ggplot(df, aes(x = avg_log2FC, y = neg_log_padj, color = group)) +
  geom_point(size = 1) +  
  theme(legend.position = "none") +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black")
  ) +
  ggtitle("Cell") +
  ylab("-log10(FDR)") +
  xlab("log2 Fold Change") +
  scale_color_manual(values = c("up_sig" = "#D71D93", "down_sig" = "#4B8CFF", "up_non_sig" = "gray", "down_non_sig" = "gray")) +
  guides(color = FALSE) +
  geom_hline(yintercept = -log10(0.05), alpha = 0.5, linetype = "dashed", color = "gray") +
  geom_vline(xintercept = c(-0.1, 0.1), alpha = 0.5, linetype = "dashed", color = "gray") +
  coord_cartesian(ylim = c(0, 520)) +  
  geom_text_repel(
    data = subset(df, gene %in% c("Gli1", "Egfr", "Gria2","Cxcl12", "Cxcr4","Neurod1", "Cnr1", "Fgf11","Gabbr1", "Neurog1","Bmp7","Notch3", "Fzd8", "Wnt7a", "Dll1", "Notch2", "Aldoc", "Ccnd1", "Sox2")),  # Highlight specific genes
    aes(label = gene),
    color = "black",
    box.padding = unit(.7, "lines"),
    hjust= 0.30
  )
volcano  
  
####code to make dot plot fig 7a on filtered DEG cav having expression and FDR values for Ackr3, Cxcl12 and Cxcr4

 
mia1_filtered <- mia1 %>%
  filter(gene %in% c("Ackr3", "Cxcl12","Cxcr4"), comparison == "E18_vs_E14") %>%
  mutate(`-log10(padj)` = ifelse(p_val_adj == 0, 300, -log10(p_val_adj)))  # Set a max value for p-value = 0

# Create the plot
ggplot(mia1_filtered, aes(x = celltype, y = gene, size = `-log10(padj)`, color = avg_log2FC)) +
  geom_point() +
  scale_color_viridis(option = "magma") +
  scale_size_continuous(range = c(1, 10)) +
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