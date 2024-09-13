library(Seurat)
#library(BPCells)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(patchwork)
# set this option when analyzing large datasets
options(future.globals.maxSize = 3e+09)
options(Seurat.object.assay.version = "v5")

# preprocessing of each merfish segmented dataset was performed using following code for each region of MIA and PBS dataset

am6<-ReadVizgen(data.dir="/MIA_2_reg0/", type = c("segmentations", "centroids"), z = 3L)

#add metadata
am6_mtd = read.csv("/MIA_2_reg0/cell_metadata.csv")
rownames(am6_mtd) <- am6_mtd$X; am6_mtd$X <- NULL

#create seurat object

af1.obj <- CreateSeuratObject(counts = af1$transcripts, meta.data = af1_mtd, assay = "Vizgen")
coords <- subset(x = coords, cells = intersect(x = Cells(x = coords[["centroids"]]), y = Cells(x = af1.obj)))

mtd1 <- af1_mtd[,c('center_x', 'center_y')]
coords <- CreateFOV(mtd1, type = "centroids", assay = "Vizgen")
af1.obj[["MIA_2_reg0"]] <- coords

vizgen.obj <- SCTransform(af1.obj, assay = "Vizgen", clip.range = c(-10, 10),vst.flavor="v2")
saveRDS(vizgen.obj, "MIA_2_reg0_SCT.rds")

##create SCT objects for all samples and then integrate them. you can try integrating complete dataset or subset by region and perform integration on subseted data

#read all SCT datasets for MIA and PBS and integrate
mia <- merge(af1.obj, y = c(af2.obj, af3.obj, af4.obj, af5.obj, am6.obj, am7.obj, am8.obj, am9.obj, am10.obj,pf1.obj,pf2.obj, pf3.obj, pf4.obj, pf5.obj, pm6.obj, pm7.obj, pm8.obj, pm9.obj, pm10.obj), project = "Vizgen")


# split assay into 5 layers

# split assay into 20 layers (number of images)

mia[["Vizgen"]] <- split(mia[["Vizgen"]], f = mia$Sample)
mia[["SCT"]] <- split(mia[["SCT"]], f = mia$Sample)
mia <- FindVariableFeatures(mia, verbose = FALSE)

mia <- SketchData(object = mia, ncells = 5000, method = "LeverageScore", sketched.assay = "sketch1")
mia
DefaultAssay(mia) <- "sketch1"
mia <- FindVariableFeatures(mia, verbose = F)
mia <- ScaleData(mia, verbose = F)
mia <- RunPCA(mia, verbose = F)
# integrate the datasets
#mia[["sketch1"]] <- JoinLayers(mia[["sketch1"]])
#mia[["sketch1"]] <- split(mia[["sketch1"]], f = mia$Sample)
mia <- IntegrateLayers(mia, method = RPCAIntegration, orig = "pca", new.reduction = "integrated.rpca1",
   verbose = F)
saveRDS(mia,"mia_int_RPCA.rds")	


# cluster the integrated data
mia <- FindNeighbors(mia, reduction = "integrated.rpca1", dims = 1:30)
mia <- FindClusters(mia, resolution = 0.8)
mia <- FindClusters(mia, resolution = 1)
mia <- RunUMAP(mia, reduction = "integrated.rpca1", dims = 1:30, return.model = T, verbose = F)
mia <- RunTSNE(mia, reduction = "integrated.rpca1", dims = 1:30, return.model = T, verbose = F)

saveRDS(mia,"mia_sketch_umap_tsne_int_rpca.rds")
# you can now rejoin the layers in the sketched assay this is required to perform differential expression
mia[["sketch1"]] <- JoinLayers(mia[["sketch1"]])

# You can now annotate clusters using marker genes.
# results in the 'sketch.celltype' metadata column

plot.s1 <- DimPlot(object, reduction = "umap")
plot.s2 <- DimPlot(object, group.by = "predicted.celltype", reduction = "ref.umap")
plot.s3 <- DimPlot(object, group.by = "predicted.celltype", reduction = "umap")
plot.s1 + plot.s2 + plot_layout(ncol = 1)

#Integrate the full datasets

# resplit the sketched cell assay into layers this is required to project the integration onto
# all cells
mia[["sketch1"]] <- split(mia[["sketch1"]], f = mia$Sample)

#mia <- ProjectIntegration(object = mia, sketched.assay = "sketch1", assay = "Vizgen", reduction = "integrated.rpca1")


#mia <- ProjectData(object = mia, sketched.assay = "sketch1", assay = "SCT", sketched.reduction = "integrated.rpca1.full",
    full.reduction = "integrated.rpca1.full", dims = 1:30, refdata = list(celltype.full = "broad_celltype_est"))
mia <- ProjectData(
  object = mia,
  assay = "Vizgen",
  full.reduction = "integrated.rpca1.full",
  sketched.assay = "sketch1",
  sketched.reduction = "integrated.rpca1",
  dims = 1:30,
  refdata = list(cluster_full = "sketch1_snn_res.0.8")
)
saveRDS(mia,"obj_mia_sketch1_umap_tsne_rpcafull.rds")
mia <- RunUMAP(mia, reduction = "integrated.rpca1.full", dims = 1:30, reduction.name = "umap.full",
    reduction.key = "UMAP_full_")
saveRDS(mia, "obj_mia_sketch_fullumap_clusterfull_rpcaint.rds")
