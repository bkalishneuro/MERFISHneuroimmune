#This code was used to create cellchat object and make plots for LR analysis for MIAM vs PBSM, MIAF vs PBSF, Abx vs PBS comparison. This code shows MIM vs PBSM comparison but can be used for other comparison by changing conditions

#load libraries
library(CellChat)
library(patchwork)
library(Seurat)
options(stringsAsFactors = FALSE)
options(future.globals.maxSize = 1024 * 1024^2)  # Set limit to 1 GB
future::plan("multisession", workers = 4) 

#read seurat objects
mia <- readRDS("mia_dataset.rds")
miam <- subset(mia,subset =Condition%in%c("MIAM"))

pbsm <- subset(mia, subset = con%in%c("PBSM"))
seurat <- pbsm

seurat$broad_celltype = droplevels(seurat$broad_celltype, exclude = setdiff(levels(seurat$broad_celltype),unique(seurat$broad_celltype)))
levels(seurat) <- seurat$broad_celltype

replicates <- names(seurat@images)

data.input <- GetAssayData(seurat, slot = "counts", assay = 'Vizgen')
meta <- seurat@meta.data
aggregated_spatial_locs <- data.frame()
all_radius_values <- c()
offset <- 0
group.by <- "broad_celltype_abx1"
scale.distance = 0.2
conditions <- unique(seurat@meta.data[['Condition']])
Idents(seurat) <- seurat$Condition

conversion.factor = 1
spot.size = 10
spatial.factors = data.frame(ratio = conversion.factor, tol = spot.size/2)

#create cellchat object for each condition
cellchat <- createCellChat(object = data.input, meta = meta, group.by = group.by,
                             datatype = "spatial", coordinates = aggregated_spatial_locs, spatial.factors = spatial.factors)
							 
CellChatDB <- CellChatDB.mouse # use CellChatDB.human if running on human data
CellChatDB.use <- CellChatDB # simply use the default CellChatDB
cellchat@DB <- CellChatDB.use
options(future.globals.maxSize = 302400 * 1024^2)
future::plan("multisession", workers = 4)
cellchat <- subsetData(cellchat)
cellchat <- updateCellChat(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0.1,
                              distance.use = TRUE, interaction.range = 250, scale.distance = 0.2,
                              contact.dependent = TRUE, contact.range = 10)
							  
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat_PBSM <- netAnalysis_computeCentrality(cellchat)
saveRDS(cellchat_PBSM, "Cellchat_PBSM_object.rds")


##create cellchat object for MIAM
seurat <- miam

seurat$broad_celltype = droplevels(seurat$broad_celltype, exclude = setdiff(levels(seurat$broad_celltype),unique(seurat$broad_celltype)))
levels(seurat) <- seurat$broad_celltype

replicates <- names(seurat@images)

data.input <- GetAssayData(seurat, slot = "counts", assay = 'Vizgen')
meta <- seurat@meta.data
aggregated_spatial_locs <- data.frame()
all_radius_values <- c()
offset <- 0
group.by <- "broad_celltype_abx1"
scale.distance = 0.2
conditions <- unique(seurat@meta.data[['Condition']])
Idents(seurat) <- seurat$Condition

conversion.factor = 1
spot.size = 10
spatial.factors = data.frame(ratio = conversion.factor, tol = spot.size/2)

#create cellchat object for each condition
cellchat <- createCellChat(object = data.input, meta = meta, group.by = group.by,
                             datatype = "spatial", coordinates = aggregated_spatial_locs, spatial.factors = spatial.factors)
							 
CellChatDB <- CellChatDB.mouse # use CellChatDB.human if running on human data
CellChatDB.use <- CellChatDB # simply use the default CellChatDB
cellchat@DB <- CellChatDB.use
options(future.globals.maxSize = 302400 * 1024^2)
future::plan("multisession", workers = 4)
cellchat <- subsetData(cellchat)
cellchat <- updateCellChat(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0.1,
                              distance.use = TRUE, interaction.range = 250, scale.distance = 0.2,
                              contact.dependent = TRUE, contact.range = 10)
							  
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat_MIAM <- netAnalysis_computeCentrality(cellchat)
saveRDS(cellchat_MIAM, "Cellchat_MIAM_object.rds")


#combine objects for comparison 

object.list <- list(PBSM = cellchat_PBSM, MIAM = cellchat_MIAM)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

#Perform differential comparison

pos.dataset <- 'MIAM'
features.name <- pos.dataset
cellchat <- identifyOverExpressedGenes(cellchat, group.dataset = "datasets",
                                       pos.dataset = pos.dataset,
                                       features.name = features.name,
                                       only.pos = FALSE,
                                       thresh.pc = 0.1,
                                       thresh.fc = 0.1,
                                       thresh.p = 1)

net <- netMappingDEG(cellchat, features.name = features.name)
net.up <- subsetCommunication(cellchat, net = net, datasets = "MIAM",
                                thresh = 0.05,
                                ligand.logFC = 0.1,
                                receptor.logFC = 0.1)
gene.up <- extractGeneSubsetFromPair(net.up, cellchat)

write.csv(net, "MIAM_vs_PBSM_net.csv")
write.csv(net.up, "MIAM_vs_PBSM_netup.csv")
write.csv(gene.up, "MIAM_vs_PBSM_geneup.csv")
net.down <- subsetCommunication(cellchat, net = net, datasets = "PBSF",
                                  thresh = 0.05,
                                  ligand.logFC = - 0.1,
                                  receptor.logFC = - 0.1)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat)
write.csv(net.down, "MIAM_vs_PBSM_netdown.csv")
write.csv(gene.down, "MIAM_vs_PBSM_genedown.csv")

#use following code to create chord diagram for Cxcl12 pathways

pdf("males_ackr3.pdf", width=10,height =10)
ackr_up <- net.up %>% filter(receptor == "Ackr3")					   
 pathways.show <- c("CXCL")
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_chord_cell(object.list[[i]], signaling = pathways.show,slot.name = "netP", pairLR.use = pairLR.use.up,lab.cex = 0.3,small.gap = 2.5,big.gap = 5,color.use = color_vector, remove.isolate = TRUE, annotationTrackHeight = c(0.01),title.name = paste0(pathways.show, " signaling network - ", names(object.list)[i]))
}

par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_chord_cell(object.list[[i]], signaling = pathways.show,slot.name = "netP",pairLR.use = pairLR.use.up,lab.cex = 0.000000000000000000000000000001,remove.isolate = TRUE, color.use = color_vector, small.gap = 2.5,big.gap = 5,annotationTrackHeight = c(0.01),title.name = paste0(pathways.show, " signaling network - ", names(object.list)[i]))
}					   
dev.off()

pdf("males_cxcr4.pdf", width=10,height =10)
ackr_up <- net.up %>% filter(receptor == "Cxcr4")					   
 pathways.show <- c("CXCL")
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_chord_cell(object.list[[i]], signaling = pathways.show,slot.name = "netP", pairLR.use = pairLR.use.up,lab.cex = 0.3,small.gap = 2.5,big.gap = 5,color.use = color_vector, remove.isolate = TRUE, annotationTrackHeight = c(0.01),title.name = paste0(pathways.show, " signaling network - ", names(object.list)[i]))
}

par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_chord_cell(object.list[[i]], signaling = pathways.show,slot.name = "netP",pairLR.use = pairLR.use.up,lab.cex = 0.000000000000000000000000000001,remove.isolate = TRUE, color.use = color_vector, small.gap = 2.5,big.gap = 5,annotationTrackHeight = c(0.01),title.name = paste0(pathways.show, " signaling network - ", names(object.list)[i]))
}					   
dev.off()