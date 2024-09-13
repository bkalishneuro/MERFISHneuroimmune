# Load required libraries
library(sf)        # Simple Features for R
library(dplyr)     # Data manipulation
library(pbapply)   # Progress bar for apply functions
library(reticulate) # Interface to Python
library(lmerTest)  # Linear mixed effects models
library(brms)      # Bayesian regression models
library(nlme)      # Linear and nonlinear mixed effects models
library(Seurat)    # Single-cell RNA-seq analysis
library(ggplot2)   # Plotting
library(transport) # Optimal transport
library(gridExtra) # Arrange multiple grid-based figures

# Source the Python script for cell proximity calculations
source_python("Scripts/Tims_cell_proximity_python_helper.py")

# Function to calculate Cohen's d effect size
cohen_d <- function(x, y) {
  lx <- length(x) - 1
  ly <- length(y) - 1
  md <- mean(x) - mean(y)
  pooled_sd <- sqrt(((lx * var(x)) + (ly * var(y))) / (lx + ly))
  d <- md / pooled_sd
  return(d)
}

# Define input parameters
input <- "Project/MIA_Updated_Object/Abx_PBS_E14_Male.rds" # Path to the Seurat object
output_dir <- "Project/Cell_Proximity_Tool_Dev/Abx_PBS_Results" # Output directory for results
cell_annotation_column <- "broad_celltype_abx1" # Column for cell type annotations
comparison_column <- "Condition" # Column for experimental conditions (must have two conditions)
sample_column <- "Sample" # Column identifying each tissue sample
baseline_condition <- "PBSM" # Baseline condition for comparisons
p_value_cutoff <- 0.01 # P-value cutoff for statistical significance
cut_off_scaler <- 30 # Distance cutoff for filtering outliers

# Read Seurat object
sobject <- readRDS(input) # Load the Seurat object
sobject@meta.data$cell_type <- sobject@meta.data[[cell_annotation_column]] # Annotate cell types
sobject@meta.data$condition <- sobject@meta.data[[comparison_column]] # Annotate conditions
sobject@meta.data$Sample <- sobject@meta.data[[sample_column]] # Annotate samples
all_cell_data <- sobject@meta.data %>% dplyr::select(center_x, center_y, cell_type, Sample, condition) # Select relevant columns
images <- unique(sobject$Sample) # Get unique sample names
distance_dfs <- list() # Initialize list to store distance dataframes

# Iterate over each sample and calculate distances
for (image in images) {
  print(paste0("Now working on sample name: ", image))
  
  sub_csv_dir <- file.path(output_dir, "results_csv")
  dir.create(sub_csv_dir, recursive = TRUE, showWarnings = FALSE) # Create directory for CSV results
  csv_file_name <- file.path(sub_csv_dir, paste0(image, ".csv"))
  
  if (file.exists(csv_file_name)) {
    file <- read.csv(csv_file_name) # Load existing CSV if available
    file$Sample <- image
    cell_data <- all_cell_data[all_cell_data$Sample == image,]
    file$condition <- cell_data$condition
    distance_dfs[[image]] <- file
    print(paste0("File exists, skipping sample: ", image))
    next
  }
  
  print(paste0("Output will be saved to: ", csv_file_name))
  cell_data <- all_cell_data[all_cell_data$Sample == image,]
  cell_data$cell_ID <- rownames(cell_data)
  data <- cell_data %>% dplyr::select(cell_ID, center_x, center_y, cell_type)
  data <- r_to_py(data)
  results <- get_nearest_neighbors(data) # Calculate nearest neighbors using Python
  scale_factor <- calculate_median_nearest_neighbor_distance(data) # Calculate median nearest neighbor distance
  results$condition <- cell_data$condition
  results$scale_factor <- scale_factor
  results$Sample <- image
  stats <- calculate_distance_statistics(data) # Calculate distance statistics
  stats$normalized <- stats$distances_df / stats$median_distance
  
  write.csv(results, csv_file_name) # Save results to CSV
  distance_dfs[[image]] <- results
  
  sub_hist_dir <- file.path(output_dir, "results_hist")
  dir.create(sub_hist_dir, recursive = TRUE) # Create directory for histograms
  
  cellpair_distances_plot_name <- file.path(sub_hist_dir, paste0(image, ".pdf"))
  pdf(cellpair_distances_plot_name)
  hist(stats$distances_df, 
       breaks = seq(min(stats$distances_df), max(stats$distances_df), length.out = 100), 
       main = paste0("Distribution of Cell Pairs Distance of ", image), 
       xlab = "Distance", ylab = "Frequency")
  dev.off()
  
  cellpair_distances_normalized_plot_name <- file.path(sub_hist_dir, paste0(image, "_normalized.pdf"))
  pdf(cellpair_distances_normalized_plot_name)
  hist(stats$normalized, 
       breaks = seq(min(stats$normalized), max(stats$normalized), length.out = 100), 
       main = paste0("Normalized Distribution of Cell Pairs Distance of ", image), 
       xlab = "Normalized Distance", ylab = "Frequency")
  dev.off()
}

# Combine all distance results into one dataframe
combined_df <- bind_rows(distance_dfs)
colnames(combined_df) <- gsub("\\.", "-", colnames(combined_df)) # Replace dots with hyphens in column names

all_conditions <- unique(combined_df$condition) # Get all unique conditions
disease_condition <- all_conditions[all_conditions != baseline_condition] # Identify disease condition

# Set the baseline condition as the first level
combined_df$condition <- factor(combined_df$condition, 
                                levels = c(baseline_condition, disease_condition))

celltypes <- unique(combined_df$cell_type) # Get unique cell types

# Initialize results dataframe
results_df <- tibble(
  CentroidCellType = character(),
  QueryCellType = character(),
  PValue_wilcox = numeric(),
  PValue_lme = numeric(),
  MeanDifference = numeric(),
  MedianDifference = numeric(), 
  FC_Median = numeric(), 
  FC_Mean = numeric(), 
  Effectsize = numeric(),
  Condition1Median = numeric(),
  Condition2Median = numeric()
)

cut_off_scaler <- 50 # Distance cutoff for statistical analysis

# Analyze cell distances and statistical comparisons
for (centroid_cell in celltypes) {
  for (query_cell in celltypes) {
    # Uncomment and modify these lines for specific comparisons
    # centroid_cell <- "Microglia"
    # query_cell <- "Neuroblast"
    
    query_id <- paste0(query_cell, "_id")
    query_dist <- paste0(query_cell, "_dist")
    
    sub_results_df <- combined_df[combined_df$cell_type == centroid_cell,]
    sub_results_pair <- sub_results_df %>% dplyr::select(cell_ID, Sample, cell_type, condition, query_dist, query_id)
    sub_results_df_control <- sub_results_df[sub_results_df$condition == baseline_condition,]
    sub_results_df_disease <- sub_results_df[sub_results_df$condition == disease_condition,]
    sub_results_pair_control <- sub_results_df_control %>% 
      dplyr::select(cell_ID, Sample, cell_type, condition, query_dist, query_id)
    sub_results_pair_disease <- sub_results_df_disease %>% 
      dplyr::select(cell_ID, Sample, cell_type, condition, query_dist, query_id)
    
    distance_control <- sub_results_pair_control[[query_dist]]
    distance_disease <- sub_results_pair_disease[[query_dist]]
    
    distance_control <- distance_control[distance_control < cut_off_scaler]
    distance_disease <- distance_disease[distance_disease < cut_off_scaler]
    
    distance_control <- na.omit(distance_control)
    distance_disease <- na.omit(distance_disease)
    
    hist(distance_control, breaks = 50, main = paste0(query_cell, " Related to ", centroid_cell, " Distance Histogram in ", baseline_condition))
    hist(distance_disease, breaks = 50, main = paste0(query_cell, " Related to ", centroid_cell, " Distance Histogram in ", disease_condition))
    
    median_control <- median(distance_control)
    median_disease <- median(distance_disease)
    
    mean_control <- mean(distance_control)
    mean_disease <- mean(distance_disease)
    
    mean_diff <- mean_disease - mean_control
    median_diff <- median_disease - median_control
    log2fc_median <- log2(median_disease / median_control)
    log2fc_mean <- log2(mean_disease / mean_control)
    
    min_data_points <- 5 # Minimum number of data points required
    
    if (length(distance_control) >= min_data_points && length(distance_disease) >= min_data_points) {
      if (sum(!is.na(distance_control)) >= min_data_points && sum(!is.na(distance_disease)) >= min_data_points) {
        wilcox_result <- wilcox.test(distance_control, distance_disease)
        p_value_wilcox <- wilcox_result$p.value
      } else {
        p_value_wilcox <- NA
        message("Not enough non-missing observations for Wilcox test.")
      }
      
      if (nrow(na.omit(sub_results_pair)) >= min_data_points) {
        full_model <- lmer(get(query_dist) ~ condition + (1 | Sample), data = sub_results_pair)
        null_model <- lmer(get(query_dist) ~ (1 | Sample), data = sub_results_pair)
        lmer_results <- anova(null_model, full_model)
        p_value_lme <- lmer_results$"Pr(>Chisq)"[2]
      } else {
        p_value_lme <- NA
        message("Not enough data points for linear mixed-effects model.")
      }
      
      effect_size <- cohen_d(distance_disease, distance_control)
      
      results_df <- rbind(results_df, tibble(
        CentroidCellType = centroid_cell,
        QueryCellType = query_cell,
        PValue_wilcox = p_value_wilcox,
        PValue_lme = p_value_lme,
        MeanDifference = mean_diff,
        MedianDifference = median_diff, 
        FC_Median = log2fc_median, 
        FC_Mean = log2fc_mean, 
        Effectsize = effect_size,
        ControlMedian = median_control,
        DiseaseMedian = median_disease
      ))
    } else {
      message("Not enough data points after filtration to perform the test.")
      results_df <- rbind(results_df, tibble(
        CentroidCellType = centroid_cell,
        QueryCellType = query_cell,
        PValue_wilcox = NA,
        PValue_lme = NA,
        MeanDifference = NA,
        MedianDifference = NA, 
        FC_Median = NA, 
        FC_Mean = NA, 
        Effectsize = NA,
        ControlMedian = NA,
        DiseaseMedian = NA
      ))
    }
  }
}

# Adjust p-values for multiple comparisons
results_df$PValueAdjusted_wilcox_BH <- p.adjust(results_df$PValue_wilcox, method = "BH")
results_df$PValueAdjusted_lme_BH <- p.adjust(results_df$PValue_lme, method = "BH")

# Annotate significance in the results
results_df <- results_df %>%
  mutate(Comparison_Significance = case_when(PValueAdjusted_wilcox_BH < p_value_cutoff ~ "*", TRUE ~ "")) %>%
  mutate(Model_Significance = case_when(PValueAdjusted_lme_BH < p_value_cutoff ~ "*", TRUE ~ "")) %>%
  mutate(Significance = case_when(Comparison_Significance == "*" & Model_Significance == "*" ~ "*", TRUE ~ ""))

# Save the results to CSV
write.csv(results_df, file.path(output_dir, paste0(disease_condition, "_", baseline_condition, "_results.csv")))

# Plot the results
p <- ggplot(results_df, aes(x = QueryCellType, y = CentroidCellType, fill = Effectsize)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = Significance), vjust = 0.8, color = "black", size = 6) +
  scale_fill_gradient2(low = "#7EA36B", high = "#8C71A7", mid = "white", midpoint = 0, name = "Effect Size") +
  theme_minimal() +
  labs(title = paste0("Cell Distances Changes: ", disease_condition, " vs ", baseline_condition),
       x = "Surrounding Cell Type", y = "Center Cell Type", fill = "Effect Size") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(face = "bold"),
        legend.title.align = 0)

# Save the plot to PDF
pdf(file.path(output_dir, paste0(disease_condition, "_", baseline_condition, "_results.pdf")), width = 5.25, height = 4.47)
print(p)
dev.off()

# Spatial plot generation
dir.create(file.path(output_dir, "Spatial_plot"), recursive = TRUE, showWarnings = FALSE)

centroid_cell <- "Microglia"
query_cell <- "MGE_Int"
query_id <- paste0(query_cell, "_id")
query_dist <- paste0(query_cell, "_dist")

sub_results_df <- combined_df[combined_df$cell_type == centroid_cell,]
sub_results_pair <- sub_results_df %>% dplyr::select(cell_ID, cell_type, condition, query_dist, query_id)
color_values <- setNames(c("lightgray", "blue", "red"), c("other", query_cell, centroid_cell))
all_pair_cells <- c(as.vector(sub_results_pair$cell_ID), as.vector(sub_results_pair[[query_id]]))

coords <- combined_df %>%
  dplyr::select(cell_ID, Sample, cell_type, center_x, center_y, condition) %>%
  dplyr::mutate(cell_type = if_else(cell_ID %in% all_pair_cells, cell_type, "other")) %>%
  dplyr::arrange(cell_type != "other")

# Iterate over each condition to generate spatial plots
for (cond in unique(coords$condition)) {
  condition_data <- coords %>% filter(condition == cond)
  plots <- list()
  plot_names <- list()
  
  for (sample in unique(condition_data$Sample)) {
    sample_data <- condition_data %>% filter(Sample == sample)
    
    p <- ggplot(sample_data, aes(x = center_x, y = center_y, color = cell_type)) +
      geom_point() +  
      labs(title = paste("Spatial Plot ", sample, "- ", centroid_cell, "_vs_", query_cell), 
           x = "Center X", 
           y = "Center Y",
           color = "Cell Type") +
      scale_color_manual(values = color_values) +
      theme_minimal() +
      theme(legend.position = "right")
    
    plots[[sample]] <- p
    plot_names[[sample]] <- sprintf("%s_%s_%s_%s", cond, sample, centroid_cell, query_cell)
  }
  
  g <- do.call(grid.arrange, c(plots, ncol = 2))
  file_name <- file.path(output_dir, "Spatial_plot", sprintf("%s_vs_%s_%s.pdf", centroid_cell, query_cell, cond))
  ggsave(file_name, plot = g, device = "pdf", width = 15, height = 15)
}

# End of script
