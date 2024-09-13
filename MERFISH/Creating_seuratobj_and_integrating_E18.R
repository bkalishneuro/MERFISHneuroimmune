# preprocessing of each merfish segmented dataset was performed using following code for each region of MIA, MMD and PBS dataset

af1<-ReadVizgen(data.dir="/E18_MIA/HDF_new/E18M_A3_reg1/", type = c("segmentations", "centroids"), z = 3L)

#add metadata
af1_mtd = read.csv("/E18_MIA/HDF_new/E18M_A3_reg1/cell_metadata.csv")
rownames(af1_mtd) <- af1_mtd$X; af1_mtd$X <- NULL

#create seurat object

af1.obj <- CreateSeuratObject(counts = af1$transcripts, meta.data = af1_mtd, assay = "Vizgen")
coords <- subset(x = coords, cells = intersect(x = Cells(x = coords[["centroids"]]), y = Cells(x = af1.obj)))

mtd1 <- af1_mtd[,c('center_x', 'center_y')]
coords <- CreateFOV(mtd1, type = "centroids", assay = "Vizgen")
af1.obj[["E18M_A3_reg1"]] <- coords

vizgen.obj <- SCTransform(af1.obj, assay = "Vizgen", clip.range = c(-10, 10),vst.flavor="v2")
saveRDS(vizgen.obj, "E18M_A3_reg1_SCT.rds")

##create SCT objects for all samples and then integrate them. you can try integrating complete dataset or subset by region and perform integration on subseted data

#read all SCT datasets and integrate, here showing example of E18 PBS males

E18 <- merge(vz1, y = c(vz2, vz3, vz4,vz5), add.cell.ids = c("PBS_E18M","PBS_E18M","PBS_E18M","PBS_E18M","PBS_E18M"), project = "Vizgen")

# split assay into 5 layers

E18[["Vizgen"]] <- split(E18[["Vizgen"]], f = E18$Sample)
E18[["SCT"]] <- split(E18[["SCT"]], f = E18$Sample)
E18 <- FindVariableFeatures(E18, verbose = FALSE)
#use seurat sketch vignette to integrate big dataset and then transfer clustering to complete dataset
E18 <- SketchData(object = E18, ncells = 20000, method = "LeverageScore", sketched.assay = "sketch1")
E18
DefaultAssay(E18) <- "sketch1"
E18 <- FindVariableFeatures(E18, verbose = F)
E18 <- ScaleData(E18, verbose = F)
E18 <- RunPCA(E18, verbose = F)
# integrate the datasets
#E18[["sketch1"]] <- JoinLayers(E18[["sketch1"]])
#E18[["sketch1"]] <- split(E18[["sketch1"]], f = E18$Sample)
E18 <- IntegrateLayers(E18, method = RPCAIntegration, orig = "pca", new.reduction = "integrated.rpca1",
   verbose = F)
saveRDS(E18,"E18M_int_RPCA.rds")	
E18 <- IntegrateLayers(E18, method = RPCAIntegration, orig = "pca", new.reduction = "integrated.rpca1",,
    verbose = F)	

# cluster the integrated data
E18 <- FindNeighbors(E18, reduction = "integrated.rpca1", dims = 1:30)
E18 <- FindClusters(E18, resolution = 0.8)
E18 <- FindClusters(E18, resolution = 1)
E18 <- RunUMAP(E18, reduction = "integrated.rpca1", dims = 1:30, return.model = T, verbose = F)
E18 <- RunTSNE(E18, reduction = "integrated.rpca1", dims = 1:30, return.model = T, verbose = F)

saveRDS(E18,"E18M_sketch_umap_tsne_int_rpca.rds")
# you can now rejoin the layers in the sketched assay this is required to perform differential
# expression
E18[["sketch1"]] <- JoinLayers(E18[["sketch1"]])

# You can now annotate clusters using marker genes.  We performed this step, and include the
# results in the 'sketch.celltype' metadata column

plot.s1 <- DimPlot(object, reduction = "umap")
plot.s2 <- DimPlot(object, group.by = "predicted.celltype", reduction = "ref.umap")
plot.s3 <- DimPlot(object, group.by = "predicted.celltype", reduction = "umap")
plot.s1 + plot.s2 + plot_layout(ncol = 1)

#Integrate the full datasets

# resplit the sketched cell assay into layers this is required to project the integration onto
# all cells
E18[["sketch1"]] <- split(E18[["sketch1"]], f = E18$Sample)

E18 <- ProjectData(
  object = E18,
  assay = "Vizgen",
  full.reduction = "integrated.rpca1.full",
  sketched.assay = "sketch1",
  sketched.reduction = "integrated.rpca1",
  dims = 1:30,
  refdata = list(cluster_full = "sketch1_snn_res.0.8")
)
saveRDS(E18,"obj_e18M_sketch1_umap_tsne_rpcafull.rds")
E18 <- RunUMAP(E18, reduction = "integrated.rpca1.full", dims = 1:30, reduction.name = "umap.full",
    reduction.key = "UMAP_full_")
saveRDS(E18, "obj_E18M_sketch_fullumap_clusterfull_rpcaint.rds")



