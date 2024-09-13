library(Signac)
library(Seurat)
library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Mmusculus.UCSC.mm10)
library(patchwork)


#ATAC seq Motif analysis
#Adding motif information to the Seurat object
# Get a list of motif position frequency matrices from the JASPAR database
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)

# add motif information
combined <- AddMotifs(
  object = combined,
  genome = BSgenome.Mmusculus.UCSC.mm10,
  pfm = pfm
)
##Find enriched motifs in TAPs
#Finding overrepresented motifs
da_peaks <- read.csv("E12_TAPS_DEG.csv")

# get top differentially accessible peaks
top.da.peak <- rownames(da_peaks[da_peaks$p_val < 0.005 & da_peaks$pct.1 > 0.2, ])


#enrichment
enriched.motifs <- FindMotifs(
  object = combined,
  features = top.da.peak
)
write.csv(enriched.motifs,"enriched_motifs_TAPs.csv")

##find enriched motifs in IPCs
da_peaks <- read.csv("E12_IPC_DEG.csv")

# get top differentially accessible peaks
top.da.peak <- rownames(da_peaks[da_peaks$p_val < 0.005 & da_peaks$pct.1 > 0.2, ])


#enrichment
enriched.motifs <- FindMotifs(
  object = combined,
  features = top.da.peak
)
write.csv(enriched.motifs,"enriched_motifs_IPC.csv")

##code to create motif bar plots
df <- read.csv("TAPs_top_20_motifs_male_and_female.csv", sep = '\t')

pdf("TAPS_plot.pdf")
ggplot(df1, aes(x = reorder(motif.name, fold.enrichment), y = fold.enrichment, fill = Sex)) +
  geom_bar(stat = "identity", width = 0.5) +
  scale_fill_manual(values = c("Female" = "pink", "Male" = "lightblue")) +
  labs(title = "Top 20 Motifs by Fold Enrichment", x = "Motif", y = "Fold Enrichment") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 12, color = "black"),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 12),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(color = "black")
  ) +
  coord_flip() +
  facet_wrap(~Sex, scales = "free_y")
 
dev.off() 

df <- read.csv("TAPs_top_20_motifs_male_and_female.csv", sep = '\t')

pdf("IPC_plot.pdf")
ggplot(df1, aes(x = reorder(motif.name, fold.enrichment), y = fold.enrichment, fill = Sex)) +
  geom_bar(stat = "identity", width = 0.5) +
  scale_fill_manual(values = c("Female" = "pink", "Male" = "lightblue")) +
  labs(title = "Top 20 Motifs by Fold Enrichment", x = "Motif", y = "Fold Enrichment") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 12, color = "black"),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 12),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(color = "black")
  ) +
  coord_flip() +
  facet_wrap(~Sex, scales = "free_y")
 
dev.off() 
