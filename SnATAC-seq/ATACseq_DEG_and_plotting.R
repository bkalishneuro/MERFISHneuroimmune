##Code for differential accessibilty for snATAC-seq data

#read atacseq dataset
combined <- readRDS("E12_ATACseq.rds")
DefaultAssay(combined)<-"ATAC"

# cell type comparisons all_manual_labels

# manual labels
combined$group_manual_labels<-paste(combined$condition, combined$cell_type, sep="_")
Idents(combined)<-combined$group_manual_labels

celltypes<-na.omit(unique(combined$cell_type))

counts<-plyr::count(combined@meta.data[, c("condition", "group_manual_labels")])

des<-list()
##run differential accessibilty on all celltypes

for(celltype in celltypes){
    des[[celltype]]<-list()
    Female<-paste0("Female_", celltype)
    Male<-paste0("Male_", celltype)
    comparisons<-list(Female_Male=c(Female, Male))
    for(comparison in comparisons){
         if(isTRUE(counts$freq[counts$group_manual_labels==comparison[1]] > 50 
                   & counts$freq[counts$group_manual_labels==comparison[2]] > 50)){
             name<-paste(comparison, collapse="_vs_")
             diffs<-FindMarkers(combined, ident.1 = comparison[1], ident.2 = comparison[2], 
                         latent.vars = 'nCount_ATAC', test.use = "LR", logfc.threshold=0.1)
             diffs$FDR<-p.adjust(diffs$p_val, method="BH")
             diffs$sig<-diffs$FDR<0.05
             diffs$location<-rownames(diffs)
             closest<-ClosestFeature(combined, rownames(diffs))
             diffs<-left_join(diffs, closest, by=c("location"="query_region"))
             des[[celltype]][[name]]<-diffs 
             
        } else{
          message("skipping comparison ", comparison, " for ", celltype, " because there were not enough cells")
        }
    }
}

##create excel file for differential accessibilty
for(celltype in names(des)){
    dat<-des[[celltype]]
    write_xlsx(dat, paste0(celltype, "ATAC_diff_acc.xlsx"))
}


##Code for volcano plot
drp_df1 <- read.csv("ATAC_ipc.csv")
drp_df1 <- subset(drp_df1, type=='utr')
colnames(drp_df1)[colnames(drp_df1) == "avg_log2FC"] <- "avg_log2FC"
drp_df1$gene_label <- drp_df1$gene_name
  drp_df1$group <- with(drp_df1, ifelse(avg_log2FC > 0 & sig == TRUE, "up_sig",
                              ifelse(avg_log2FC > 0 & sig == FALSE, "up_non_sig",
                                     ifelse(avg_log2FC <= 0 & sig == TRUE, "down_sig", "down_non_sig"))))

 drp_df1$neg_log_padj <- -log10(drp_df1$FDR)
replacement_number <- 500
drp_df1$neg_log_padj <- replace(drp_df1$neg_log_padj, is.infinite(drp_df1$neg_log_padj), replacement_number)
volcano1 <- ggplot(drp_df1, aes(x = avg_log2FC, y = neg_log_padj, color = group)) +
  geom_point(size = 1.5) +  # Increase the size of the dots to 3 (adjust as needed)
  theme(legend.position = "none") +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black")
  ) +
  ggtitle("E12_Female vs E12_Male") +
  ylab("-log10(FDR)") +
  xlab("log2 Fold Change") +
  scale_color_manual(values = c("up_sig" = "#1434A4", "down_sig" = "#FF5733", "up_non_sig" = "gray", "down_non_sig" = "gray")) +
  guides(color = FALSE) +
  geom_hline(yintercept = -log10(0.05), alpha = 0.5, linetype = "dashed", color = "gray") +
  geom_vline(xintercept = c(-0.1, 0.1), alpha = 0.5, linetype = "dashed", color = "gray") +
  coord_cartesian(ylim = c(0, 25)) +  # Adjust y-axis range dynamically
  geom_text_repel(
    data = subset(drp_df1, gene_name %in% c("Tnrc18","Gsk3b","Ddx39b","Cirbp","Ccr1","Kat6b","Arl4d","Atp5g2","H3f3a","Scaf11","Lig3","Ackr3","Kif3b","Herc1","Mrpl9","Ugt8a","Zmynd8","Zic3","Zfp26","Bri3","Iqcd","Cnot3","Hexim1","Nbr1","Ube2i","Kat6b","Tnrc6c","Arl4d")), # Highlight specific genes
    aes(label = gene_name),
    color = "black",
    box.padding = unit(.7, "lines"), max.overlaps=Inf,
    hjust= 0.30
  )   

pdf("E12_atac_volcanoplot_IPC1.pdf")
volcano1
dev.off()

drp_df1 <- read.csv("ATAC_TAPs.csv")
drp_df1 <- subset(drp_df1, type=='utr')
colnames(drp_df1)[colnames(drp_df1) == "avg_log2FC"] <- "avg_log2FC"


drp_df1$gene_label <- drp_df1$gene_name
  drp_df1$group <- with(drp_df1, ifelse(avg_log2FC > 0 & sig == TRUE, "up_sig",
                              ifelse(avg_log2FC > 0 & sig == FALSE, "up_non_sig",
                                     ifelse(avg_log2FC <= 0 & sig == TRUE, "down_sig", "down_non_sig"))))

 drp_df1$neg_log_padj <- -log10(drp_df1$FDR)
replacement_number <- 500
drp_df1$neg_log_padj <- replace(drp_df1$neg_log_padj, is.infinite(drp_df1$neg_log_padj), replacement_number)
volcano1 <- ggplot(drp_df1, aes(x = avg_log2FC, y = neg_log_padj, color = group)) +
  geom_point(size = 1.5) +  # Increase the size of the dots to 3 (adjust as needed)
  theme(legend.position = "none") +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black")
  ) +
  ggtitle("E12_Female vs E12_Male") +
  ylab("-log10(FDR)") +
  xlab("log2 Fold Change") +
  scale_color_manual(values = c("up_sig" = "#1434A4", "down_sig" = "#FF5733", "up_non_sig" = "gray", "down_non_sig" = "gray")) +
  guides(color = FALSE) +
  geom_hline(yintercept = -log10(0.05), alpha = 0.5, linetype = "dashed", color = "gray") +
  geom_vline(xintercept = c(-0.1, 0.1), alpha = 0.5, linetype = "dashed", color = "gray") +
  coord_cartesian(ylim = c(0, 25)) +  # Adjust y-axis range dynamically
  geom_text_repel(
    data = subset(drp_df1, gene_name %in% c("Wnt1","Igf2bp3","Id4","Axin2","Rap2a","Pbx1","Cdh26","Ackr3","Creld2","Rybp","Tescl","Veph1","Trmt11","Cdh26","Sox11","Zic4","Otx2","Cdkn1a","Zic1","Far1","Fzd8","Dll4","Azin2","Kdm6a","Tescl","Trmt11","Zim1")), # Highlight specific genes
    aes(label = gene_name),
    color = "black",
    box.padding = unit(.7, "lines"), max.overlaps=Inf,
    hjust= 0.30
  )   

pdf("E12_atac_volcanoplot_TAPs1.pdf")
volcano1
dev.off()


##code for plot for Gprofiler results
# Define colors for each source
source_colors <- c("GO:MF" = "#1f77b4", "GO:BP" = "#ff7f0e", "GO:CC" = "#097969")

# Reorder label within each source category
df$term_id <- factor(df$term_id, levels = unique(df$term_id[order(df$source)]))
df$term_id <- factor(df$term_id, levels = unique(df[df$source %in% c("GO:MF", "GO:BP", "GO:CC"), "term_id"]))
# Create the bar plot
ggplot(df, aes(x = term_id, y = negative_log10_of_adjusted_p_value, fill = source)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = source_colors) +
  labs(x = "GO term", y = "-log10(Adjusted p-value)", title = "GO Enrichment Analysis") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(size = 12),
        axis.ticks = element_line(color = "black"),
        panel.grid.major.y = element_line(color = "lightgray"),
        panel.grid.minor.y = element_blank(),
        axis.line = element_line(color = "black"),
        axis.ticks.length = unit(0.3, "cm"),
        axis.ticks.margin = unit(0.5, "cm"),
        axis.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())



