library(dplyr)
library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Mmusculus.v79)
library(ggplot2)
library(patchwork)
library(GenomeInfoDb)
library(cowplot)


gr <- GRanges(seqnames = "chr1", ranges = IRanges(start = c(20, 70, 300), end = c(120, 200, 400)))
gr
#read peaks file for each sample
peaks.m1 <- read.table(
  file = "/E12-M1/outs/peaks.bed",
  col.names = c("chr", "start", "end")
)
peaks.m2 <- read.table(
  file = "/E12-M2/outs/peaks.bed",
  col.names = c("chr", "start", "end")
)
peaks.m3 <- read.table(
  file = "/E12-M3/outs/peaks.bed",
  col.names = c("chr", "start", "end")
)
peaks.f1 <- read.table(
  file = "/E12-F1/outs/peaks.bed",
  col.names = c("chr", "start", "end")
)
peaks.f2 <- read.table(
  file = "/E12-F2/outs/peaks.bed",
  col.names = c("chr", "start", "end")
)
peaks.f3 <- read.table(
  file = "/E12-F3/outs/peaks.bed",
  col.names = c("chr", "start", "end")
)
peaks.f4 <- read.table(
  file = "/E12-F4/outs/peaks.bed",
  col.names = c("chr", "start", "end")
)

# convert to genomic ranges
gr.m1 <- makeGRangesFromDataFrame(peaks.m1)
gr.m2 <- makeGRangesFromDataFrame(peaks.m2)
gr.m3 <- makeGRangesFromDataFrame(peaks.m3)
gr.f1 <- makeGRangesFromDataFrame(peaks.f1)
gr.f2 <- makeGRangesFromDataFrame(peaks.f2)
gr.f3 <- makeGRangesFromDataFrame(peaks.f3)
gr.f4 <- makeGRangesFromDataFrame(peaks.f4)

# Create a unified set of peaks to quantify in each dataset
combined.peaks <- reduce(x = c(gr.m1, gr.m2, gr.m3, gr.f1,gr.f2,gr.f3, gr.f4))

# Filter out bad peaks based on length
peakwidths <- width(combined.peaks)
combined.peaks <- combined.peaks[peakwidths  < 10000 & peakwidths > 20]
combined.peaks


# load metadata
md.m1 <- read.table(
  file = "/E12-M1/outs/singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ] # remove the first row


md.m2 <- read.table(
  file = "/E12-M2/outs/singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ]

md.m3 <- read.table(
  file = "/E12-M3/outs/singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ]

md.f1 <- read.table(
  file = "/E12-F1/outs/singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ]

md.f2 <- read.table(
  file = "/E12-F2/outs/singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ]

md.f3 <- read.table(
  file = "/E12-F3/outs/singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ]

md.f4 <- read.table(
  file = "/E12-F4/outs/singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ]


# perform an initial filtering of low count cells
md.m1 <- md.m1[md.m1$passed_filters > 500, ]
md.m2 <- md.m2[md.m2$passed_filters > 500, ]
md.m3 <- md.m3[md.m3$passed_filters > 500, ]
md.f1 <- md.f1[md.f1$passed_filters > 500, ] # sequenced deeper so set higher cutoff
md.f2 <- md.f2[md.f2$passed_filters > 500, ]
md.f3 <- md.f3[md.f3$passed_filters > 500, ]
md.f4 <- md.f4[md.f4$passed_filters > 500, ]

# create fragment objects
frags.m1 <- CreateFragmentObject(
  path = "/E12-M1/outs/fragments.tsv.gz",
  cells = rownames(md.m1)
)
## Computing hash
frags.m2 <- CreateFragmentObject(
  path = "/E12-M2/outs/fragments.tsv.gz",
  cells = rownames(md.m2)
)
## Computing hash
frags.m3 <- CreateFragmentObject(
  path = "/E12-M3/outs/fragments.tsv.gz",
  cells = rownames(md.m3)
)
## Computing hash
frags.f1 <- CreateFragmentObject(
  path = "/E12-F1/outs/fragments.tsv.gz",
  cells = rownames(md.f1)
)
## Computing hash
frags.f2 <- CreateFragmentObject(
  path = "/E12-F2/outs/fragments.tsv.gz",
  cells = rownames(md.f2)
)
## Computing hash
frags.f3 <- CreateFragmentObject(
  path = "/E12-F3/outs/fragments.tsv.gz",
  cells = rownames(md.f3)
)
## Computing hash
frags.f4 <- CreateFragmentObject(
  path = "/E12-F4/outs/fragments.tsv.gz",
  cells = rownames(md.f4)
)

#Quantify peaks in each datasset
m1.counts <- FeatureMatrix(
  fragments = frags.m1,
  features = combined.peaks,
  cells = rownames(md.m1)
)

m2.counts <- FeatureMatrix(
  fragments = frags.m2,
  features = combined.peaks,
  cells = rownames(md.m2)
)

m3.counts <- FeatureMatrix(
  fragments = frags.m3,
  features = combined.peaks,
  cells = rownames(md.m3)
)

f1.counts <- FeatureMatrix(
  fragments = frags.f1,
  features = combined.peaks,
  cells = rownames(md.f1)
)

f2.counts <- FeatureMatrix(
  fragments = frags.f2,
  features = combined.peaks,
  cells = rownames(md.f2)
)

f3.counts <- FeatureMatrix(
  fragments = frags.f3,
  features = combined.peaks,
  cells = rownames(md.f3)
)

f4.counts <- FeatureMatrix(
  fragments = frags.f4,
  features = combined.peaks,
  cells = rownames(md.f4)
)

#Create the objects
m1_assay <- CreateChromatinAssay(m1.counts, fragments = frags.m1)
objm1 <- CreateSeuratObject(m1_assay, assay = "ATAC", meta.data=md.m1)

m2_assay <- CreateChromatinAssay(m2.counts, fragments = frags.m2)
objm2 <- CreateSeuratObject(m2_assay, assay = "ATAC", meta.data=md.m2)

m3_assay <- CreateChromatinAssay(m3.counts, fragments = frags.m3)
objm3 <- CreateSeuratObject(m3_assay, assay = "ATAC", meta.data=md.m3)

f1_assay <- CreateChromatinAssay(f1.counts, fragments = frags.f1)
objf1 <- CreateSeuratObject(f1_assay, assay = "ATAC", meta.data=md.f1)

f2_assay <- CreateChromatinAssay(f2.counts, fragments = frags.f2)
objf2 <- CreateSeuratObject(f2_assay, assay = "ATAC", meta.data=md.f2)

f3_assay <- CreateChromatinAssay(f3.counts, fragments = frags.f3)
objf3 <- CreateSeuratObject(f3_assay, assay = "ATAC", meta.data=md.f3)

f4_assay <- CreateChromatinAssay(f4.counts, fragments = frags.f4)
objf4 <- CreateSeuratObject(f4_assay, assay = "ATAC", meta.data=md.f4)

saveRDS(objm1, "M1_objectnopassfilter.rds")
saveRDS(objm2, "M2_objectnopassfilter.rds")
saveRDS(objm3, "M3_objectnopassfilter.rds")
saveRDS(objf1, "F1_objectnopassfilter.rds")
saveRDS(objf2, "F2_objectnopassfilter.rds")
saveRDS(objf3, "F3_objectnopassfilter.rds")
saveRDS(objf4, "F4_objectnopassfilter.rds")

objm1 <- readRDS("M1_object.rds")
objm2 <- readRDS("M2_object.rds")
objm3 <- readRDS("M3_object.rds")
objf1 <- readRDS("F1_object.rds")
objf2 <- readRDS("F2_object.rds")
objf3 <- readRDS("F3_object.rds")
objf4 <- readRDS("F4_object.rds")

#merge objects
objm1$dataset <- 'M1'
objm2$dataset <- 'M2'
objm3$dataset <- 'M3'
objf1$dataset <- 'F1'
objf2$dataset <- 'F2'
objf3$dataset <- 'F3'
objf4$dataset <- 'F4'

# merge all datasets, adding a cell ID to make sure cell names are unique
combined <- merge(
  x = objm1,
  y = list(objm2, objm3, objf1, objf2,objf3, objf4),
  add.cell.ids = c("M1", "M2", "M3", "F1", "F2","F3","F4")
)
combined[["ATAC"]]

saveRDS(combined,"combined.rds")

annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)

# change to UCSC style since the data was mapped to hg19
seqlevels(annotations) <- paste0('chr', seqlevels(annotations))
genome(annotations) <- "mm10"

# add the gene information to the object
Annotation(combined) <- annotations

####Computing QC Metrics
combined <- NucleosomeSignal(object = combined)

combined$nucleosome_group <- ifelse(combined$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
pdf("combined_QCplots1.pdf")
FragmentHistogram(object = combined, group.by = 'nucleosome_group', region = 'chr1-1-10000000')

combined <- TSSEnrichment(combined, fast = FALSE)
combined$high.tss <- ifelse(combined$TSS.enrichment > 2, 'High', 'Low')
pdf("combined_tssplots1.pdf")
TSSPlot(combined, group.by = 'high.tss') + NoLegend()

combined$pct_reads_in_peaks <- combined$peak_region_fragments / combined$passed_filters * 100
combined$blacklist_ratio <- combined$blacklist_region_fragments / combined$peak_region_fragments

VlnPlot(
  object = combined,
  features = c('pct_reads_in_peaks', 'peak_region_fragments',
               'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal'),
  pt.size = 0.1,
  ncol = 5
)

combined <- subset(
  x = combined,
  subset = peak_region_fragments > 3000 &
    peak_region_fragments < 100000 &
    pct_reads_in_peaks > 40 &
    blacklist_ratio < 0.025 &
    nucleosome_signal < 4 &
    TSS.enrichment > 2
)
combined
##Normalization and linear dimensional reduction
pdf("dimplot_combined.pdf")
combined <- RunTFIDF(combined)
combined <- FindTopFeatures(combined, min.cutoff = 20)
combined <- RunSVD(combined)
combined <- RunUMAP(combined, dims = 2:30, reduction = 'lsi')
DimPlot(combined, group.by = 'dataset', pt.size = 0.1)
CoveragePlot(
  object = combined,
  group.by = 'dataset',
  region = "chr14-99700000-99760000"
)
dev.off()
pdf("dimplot_combined2.pdf")
combined <- FindNeighbors(
  object = combined,
  reduction = 'lsi',
  dims = 2:30
)
combined <- FindClusters(
  object = combined,
  algorithm = 3,
  resolution = 1.2,
  verbose = FALSE
)

DimPlot(object = combined, label = TRUE) + NoLegend()

dev.off()

saveRDS(combined, "combined_obj.rds")

##Create a gene activity matrix
gene.activities <- GeneActivity(combined)

# add the gene activity matrix to the Seurat object as a new assay
combined[['RNA']] <- CreateAssayObject(counts = gene.activities)
combined <- NormalizeData(
  object = combined,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(combined$nCount_RNA)
)


####code to make celltype by gene dot plot
allgene <- c("Pax6","Nfix","Sox2","Nfib","Aldoc","Eomes","Gbx2","Dbx1", "Nes", "Zic4","Pax3", "Pax7","Tbr1","Satb2","Neurod1", "Slc17a7", "Slc17a6","Dlx1","Dlx2","Ascl1","Dlx5","Gad2","Nkx2-1","Lhx6","Olig2", "Nkx2-2", 
"Aldh1l1", "Slc1a3", "Pdgfrb", "Cspg4", "Atp13a5","Cx3cr1", "Cd68", "Aif1","Pecam1", "Flt1", "Col4a2")



pdf("dotplot_cc.pdf", width=10, height=8)

# Define the reversed magma color palette
#load library
library(viridis)
reversed_magma_colors <- rev(viridis::viridis(n = 256, option = "magma"))

# Create the DotPlot with the specified order of cell types and genes
dotplot <- DotPlot(combined, features = allgene, group.by = "cell_type") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +  # Rotate x-axis labels
  scale_color_gradientn(
    colors = reversed_magma_colors,  # Apply reversed magma color palette
    limits = c(-2, 3)  # Set the color scale limits from -3 to +3
  ) +
  coord_flip()  # Flip coordinates to have cell types on x-axis and genes on y-axis

# Modify the size of the dots to represent expression level
dotplot <- dotplot +
  geom_point(aes(size = avg.exp), shape = 21, color = "black", show.legend = TRUE, data = dotplot$data) +
  theme(axis.line = element_blank(), axis.ticks = element_blank())

# Print the plot to the PDF device
print(dotplot)

# Close the PDF device
dev.off()


###code to label cell clusters
nucseq_meta<-combined@meta.data
nucseq_meta$cell_type<-NA
nucseq_meta$cell_type[nucseq_meta$seurat_clusters==0]<-"RG_1"
nucseq_meta$cell_type[nucseq_meta$seurat_clusters==1]<-"IPC"
nucseq_meta$cell_type[nucseq_meta$seurat_clusters==2]<-"Interneuron_1"
nucseq_meta$cell_type[nucseq_meta$seurat_clusters==3]<-"Progenitors_DSC1"
nucseq_meta$cell_type[nucseq_meta$seurat_clusters==4]<-"Neuroblast_1_midbrain"
nucseq_meta$cell_type[nucseq_meta$seurat_clusters==5]<-"RG_2"
nucseq_meta$cell_type[nucseq_meta$seurat_clusters==6]<-"Neuroblast_2_Satb2"
nucseq_meta$cell_type[nucseq_meta$seurat_clusters==7]<-"TAPs"
nucseq_meta$cell_type[nucseq_meta$seurat_clusters==8]<-"MGE_Neuroblast"
nucseq_meta$cell_type[nucseq_meta$seurat_clusters==21]<-"Inhibitory_pre_Spinalchord"
nucseq_meta$cell_type[nucseq_meta$seurat_clusters==10]<-"Inhibitory_Neuroblast"
nucseq_meta$cell_type[nucseq_meta$seurat_clusters==11]<-"Interneuron_2_DSC_Pnoc"
nucseq_meta$cell_type[nucseq_meta$seurat_clusters==12]<-"Progenitors_Thalamus_1"
nucseq_meta$cell_type[nucseq_meta$seurat_clusters==13]<-"Excitatory_Neuron_Cerebellum"
nucseq_meta$cell_type[nucseq_meta$seurat_clusters==14]<-"Progenitors_Hindbrain_cerebellum"
nucseq_meta$cell_type[nucseq_meta$seurat_clusters==22]<-"Progenitors_Thalamus_2"
nucseq_meta$cell_type[nucseq_meta$seurat_clusters==16]<-"Astrocytes_precursors"
nucseq_meta$cell_type[nucseq_meta$seurat_clusters==17]<-"TAPs_Hindbrain"
nucseq_meta$cell_type[nucseq_meta$seurat_clusters==23]<-"MGE_progenitors_1"
nucseq_meta$cell_type[nucseq_meta$seurat_clusters==19]<-"Inhibitory_neurons"
nucseq_meta$cell_type[nucseq_meta$seurat_clusters==20]<-"IPC_hindbrain"
nucseq_meta$cell_type[nucseq_meta$seurat_clusters==24]<-"Interneuron_3"
nucseq_meta$cell_type[nucseq_meta$seurat_clusters==26]<-"Endothelial"
nucseq_meta$cell_type[nucseq_meta$seurat_clusters==28]<-"MGE_progenitors_2"
nucseq_meta$cell_type[nucseq_meta$seurat_clusters==29]<-"Exc_Neuron_Purkinje"
nucseq_meta$cell_type[nucseq_meta$seurat_clusters==30]<-"Pericytes"
nucseq_meta$cell_type[nucseq_meta$seurat_clusters==31]<-"Progenitors_DSC2"
nucseq_meta$cell_type[nucseq_meta$seurat_clusters==36]<-"Microglia"


selected_meta<-nucseq_meta[!is.na(nucseq_meta$cell_type), ]

merged_filtered<-combined[, colnames(combined) %in% rownames(selected_meta)]

#code to make umap for labelled clusters
pdf("atac_umap_celltype.pdf", width =15, height =10)
DimPlot(object = merged_filtered, label = TRUE, group.by = 'cell_type') + NoLegend()
DimPlot(object = merged_filtered, label = FALSE, group.by = 'cell_type')
dev.off()