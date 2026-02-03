{ # always run this block upon startup
  .libPaths("/oak/stanford/groups/horence/chester/dabs_ref/REnv/seurat_4.4")
  
  suppressPackageStartupMessages({
    # library(svglite)
    library(Matrix)
    library(Seurat)
    library(ggplot2) 
    library(gridExtra)
    library(reshape2)
    library(intrinsicDimension)
    library(tibble)
    library(dplyr)
    library(patchwork)
    library(viridis)
    library(ggforce)
    library(gghalves)
    library(ggnetwork)
    library(ggridges)
    library(ggtree)
    library(scran)
    library(ape)
    library(cluster)
    library(glmGamPoi)
    library(tidyr)
    library(SingleCellExperiment)
    library(presto)
    library(clusterProfiler)
    library(org.Hs.eg.db)
    library(readr)
    library(ggrepel)
    library(scales)
    library(CytoTRACE2)
    library(Scillus)
    library(rlang)
    library(dittoSeq)
    library(ggstatsplot)
    library(ggsignif)
    library(colorspace)
    library(monocle3)
    library(SeuratDisk)
    library(SeuratWrappers)
    library(SingleR)
    # library(scDotPlot)
    library(irlba)
    library(scDblFinder)
    library(metap)
    library(multtest)
    library(mascarade)
    # library(biomaRt)
    library(cluster) # clustering
    library(kernlab) # clustering
    library(apcluster) # clustering
    # library(tricycle) # cell cycle
    library(QuieScore) # compute G0 score
  })
  
  custom_colors <- list()
  
  custom_colors$discrete <- c(
    '#FFC312','#C4E538','#12CBC4','#FDA7DF','#ED4C67',
    '#F79F1F','#A3CB38','#1289A7','#D980FA','#B53471',
    '#EE5A24','#009432','#0652DD','#9980FA','#833471',
    '#EA2027','#006266','#1B1464','#5758BB','#6F1E51',
    '#40407a','#706fd3','#f7f1e3','#34ace0','#33d9b2',
    '#2c2c54','#474787','#aaa69d','#227093','#218c74',
    '#ff5252','#ff793f','#d1ccc0','#ffb142','#ffda79',
    '#b33939','#cd6133','#84817a','#cc8e35','#ccae62'
  )
  
  custom_colors$cell_cycle <- setNames(
    c('#45aaf2', '#f1c40f', '#e74c3c', '#7f8c8d'),
    c('G1',      'S',       'G2M',     '-')
  )
  
  niceFormat <- function(number) {
    formatC(number, format = 'f', big.mark = ',', digits = 0)
  }
  
  working_dir <- "/scratch/users/jiamuyu/proj_botryllus/scRNAseq/250125_00_ciona_stemcell"
  
}

### prepare seurat object
# parent_dir <- "/scratch/groups/ayeletv/scRNA_Ciona_Nov2024/usftp21.novogene.com/03.CountData"
# 
# # List the directories to be merged
# sample_dirs <- c(
#   "TG_CONTROL_CKDL240038025-1A_HNC77DSXC",
#   "TG_STEM_B_1_CKDL240038025-1A_HNC77DSXC",
#   "Undetermined_Undetermined_HNC77DSXC",
#   "TG_ALDH_CKDL240038025-1A_HNC77DSXC",
#   "TG_STEM_A_CKDL240038025-1A_HNC77DSXC",
#   "TG_STEM_B_2_CKDL240038025-1A_HNC77DSXC"
# )

parent_dir <- "/scratch/groups/ayeletv/scRNA_Ciona_Nov2024/usftp21.novogene.com/03.CountData.250521"

# List the directories to be merged
sample_dirs <- c(
  "TG_CONTROL/TG_CONTROL",
  "TG_STEM_B_1/TG_STEM_B_1",
  "Undetermined/Undetermined",
  "TG_ALDH/TG_ALDH",
  "TG_STEM_A/TG_STEM_A",
  "TG_STEM_B_2/TG_STEM_B_2"
)

# Initialize a list to store Seurat objects
seurat_list <- list()

gene_map <- read.csv("/home/groups/ayeletv/CionaGenome/extraInfo/id2name.filtered.nounderscores.csv")

# Loop through each directory, load the data, and create a Seurat object
for (sample_dir in sample_dirs) {
  # Construct the path to the filtered feature-barcode matrix
  data_path <- file.path(parent_dir, sample_dir, "outs", "filtered_feature_bc_matrix")
  
  # Load the data using Read10X
  data <- Read10X(data.dir = data_path)
  
  rownames(data) <- gene_map$gene_name[match(rownames(data), gene_map$gene_id)]
  
  # Create a Seurat object
  seurat_obj <- CreateSeuratObject(
    counts = data,
    project = gsub("^TG_|(_CKDL240038025-1A_HNC77DSXC$|_Undetermined_HNC77DSXC$)", "", sample_dir),
    min.cells = 0,         # Filter out genes detected in fewer than 0 cells
    min.features = 0     # Filter out cells with fewer than 0 genes
  )
  
  # Add the Seurat object to the list
  seurat_list[[sample_dir]] <- seurat_obj
}

# Merge all Seurat objects into one
seurat <- merge(
  x = seurat_list[[1]],
  y = seurat_list[2:length(seurat_list)],
  add.cell.ids = gsub("^TG_|(_CKDL240038025-1A_HNC77DSXC$|_Undetermined_HNC77DSXC$)", "", sample_dirs)  # Add sample IDs as prefixes to cell barcodes
)

print(seurat)
# An object of class Seurat
# 18827 features across 11351 samples within 1 assay
# Active assay: RNA (18827 features, 0 variable features)
# 6 layers present: counts.CONTROL, counts.STEM_B_1, counts.Undetermined, counts.ALDH, counts.STEM_A, counts.STEM_B_2

seurat <- readRDS(file.path(working_dir, "raw.rds"))
# saveRDS(seurat, file.path(working_dir, "raw.rds"))

# ### purge seurat object
# ## merge the layers
# counts_list <- list()
# 
# counts_list[[1]] <- seurat@assays$RNA$"counts.TG_CONTROL_CKDL240038025-1A_HNC77DSXC"
# counts_list[[2]] <- seurat@assays$RNA$"counts.TG_STEM_B_1_CKDL240038025-1A_HNC77DSXC"
# counts_list[[3]] <- seurat@assays$RNA$"counts.Undetermined_Undetermined_HNC77DSXC"
# counts_list[[4]] <- seurat@assays$RNA$"counts.TG_ALDH_CKDL240038025-1A_HNC77DSXC"
# counts_list[[5]] <- seurat@assays$RNA$"counts.TG_STEM_A_CKDL240038025-1A_HNC77DSXC"
# counts_list[[6]] <- seurat@assays$RNA$"counts.TG_STEM_B_2_CKDL240038025-1A_HNC77DSXC"
# 
# merged_counts <- do.call(cbind, counts_list)
# 
# seurat <- CreateSeuratObject(counts = merged_counts)

# ## merge layers finale
# # Extract count matrices from the Seurat object using the specified method
# counts.TG_CONTROL_CKDL240038025_1A_HNC77DSXC <- seurat[["RNA"]]$"counts.CONTROL"
# counts.TG_STEM_B_1_CKDL240038025_1A_HNC77DSXC <- seurat[["RNA"]]$"counts.STEM_B_1"
# counts.Undetermined_Undetermined_HNC77DSXC <- seurat[["RNA"]]$"counts.Undetermined"
# counts.TG_ALDH_CKDL240038025_1A_HNC77DSXC <- seurat[["RNA"]]$"counts.ALDH"
# counts.TG_STEM_A_CKDL240038025_1A_HNC77DSXC <- seurat[["RNA"]]$"counts.STEM_A"
# counts.TG_STEM_B_2_CKDL240038025_1A_HNC77DSXC <- seurat[["RNA"]]$"counts.STEM_B_2"
# 
# # Combine count matrices into a single matrix
# concatenated_counts <- t(rbind(
#   t(counts.TG_CONTROL_CKDL240038025_1A_HNC77DSXC),
#   t(counts.TG_STEM_B_1_CKDL240038025_1A_HNC77DSXC),
#   t(counts.Undetermined_Undetermined_HNC77DSXC),
#   t(counts.TG_ALDH_CKDL240038025_1A_HNC77DSXC),
#   t(counts.TG_STEM_A_CKDL240038025_1A_HNC77DSXC),
#   t(counts.TG_STEM_B_2_CKDL240038025_1A_HNC77DSXC)
# ))

# gene_map <- read.csv("/home/groups/ayeletv/CionaGenome/id2name.filtered.nounderscores.csv")
# 
# rownames(concatenated_counts) <- gene_map$gene_name[match(rownames(concatenated_counts), gene_map$gene_id)]

# seurat[["RNA"]]$concatenated <- concatenated_counts
# seurat[["RNA"]]$`counts.CONTROL` <- NULL
# seurat[["RNA"]]$`counts.STEM_B_1` <- NULL
# seurat[["RNA"]]$`counts.Undetermined` <- NULL
# seurat[["RNA"]]$`counts.ALDH` <- NULL
# seurat[["RNA"]]$`counts.STEM_A` <- NULL
# seurat[["RNA"]]$`counts.STEM_B_2` <- NULL

# seurat <- CreateSeuratObject(counts = concatenated_counts)

# ## add sample of origin
# # seurat@meta.data$sample <- sub("^TG_", "", sub("(_CKDL240038025-1A_HNC77DSXC|_Undetermined_HNC77DSXC)$", "", sub("^(.*)_.+$", "\\1", rownames(seurat@meta.data))))
# seurat@meta.data$sample <- sub("^(.*)_.+$", "\\1", rownames(seurat@meta.data))
# 
# # unique(seurat@meta.data$sample)
# # [1] "CONTROL"      "STEM_B_1"     "Undetermined" "ALDH"         "STEM_A"       "STEM_B_2"
# 
# ## add mitochrondria content
# # seurat <- PercentageFeatureSet(seurat, pattern = "MG", col.name = "percent.mt")
# seurat <- PercentageFeatureSet(seurat, pattern = "Mitochondria", col.name = "percent.mt")

# print(seurat)
# An object of class Seurat 
# 18827 features across 11351 samples within 1 assay 
# Active assay: RNA (18827 features, 0 variable features)
# 1 layer present: counts

add_sample_metadata <- function(seurat) {
  seurat@meta.data$sample <- as.factor(seurat@meta.data$sample)
  desired_order <- c("STEM_B_1", "STEM_B_2", "ALDH", "STEM_A", "CONTROL", "Undetermined")
  seurat@meta.data$sample <- factor(seurat$sample, levels = desired_order)
  
  ### store merged_sample c("STEM_B", "ALDH", "STEM_A", "CONTROL", "Undetermined")
  seurat@meta.data$merged_sample <- as.factor(seurat@meta.data$sample)
  levels(seurat@meta.data$merged_sample) <- c(levels(seurat@meta.data$merged_sample), "STEM_B")
  seurat@meta.data$merged_sample[seurat@meta.data$merged_sample %in% c("STEM_B_1", "STEM_B_2")] <- "STEM_B"
  seurat@meta.data$merged_sample <- droplevels(seurat@meta.data$merged_sample)
  desired_order <- c("STEM_B", "ALDH", "STEM_A", "CONTROL", "Undetermined")
  seurat@meta.data$merged_sample <- factor(seurat$merged_sample, levels = desired_order)
  
  ### store combined_sample c("STEM_B", "ALDH", "CONTROL", "Undetermined")
  seurat@meta.data$combined_sample <- seurat@meta.data$merged_sample
  seurat@meta.data$combined_sample[seurat@meta.data$sample %in% c("STEM_A")] <- "CONTROL"
  seurat@meta.data$combined_sample <- droplevels(seurat@meta.data$combined_sample)
  desired_order <- c("STEM_B", "ALDH", "CONTROL", "Undetermined")
  seurat@meta.data$combined_sample <- factor(seurat$combined_sample, levels = desired_order)
  
  ### store meta_sample c("STEM", "CONTROL", "Undetermined")
  seurat@meta.data$meta_sample <- seurat@meta.data$combined_sample
  levels(seurat@meta.data$meta_sample) <- c(levels(seurat@meta.data$meta_sample), "STEM")
  seurat@meta.data$meta_sample[seurat@meta.data$meta_sample %in% c("STEM_B", "ALDH")] <- "STEM"
  seurat@meta.data$meta_sample <- droplevels(seurat@meta.data$meta_sample)
  desired_order <- c("STEM", "CONTROL", "Undetermined")
  seurat$meta_sample <- factor(seurat$meta_sample, levels = desired_order)
  return(seurat)
}

seurat <- add_sample_metadata(seurat)

seurat <- readRDS("/oak/stanford/groups/horence/chester/proj_botryllus/scRNAseq_old/250125_00_ciona_stemcell/raw.purged.0210.rds")
seurat <- readRDS(file.path(working_dir, "raw.purged.rds"))
# saveRDS(seurat, file.path(working_dir, "raw.purged.rds"))

only3_seurat <- function(seurat) {
  seurat_filtered <- subset(seurat, subset = sample != "STEM_A")
  
  levels(seurat_filtered@meta.data$merged_sample) <- c(levels(seurat_filtered@meta.data$merged_sample), "Stem", "Control")
  seurat_filtered@meta.data$merged_sample[seurat_filtered@meta.data$merged_sample %in% c("STEM_B_1", "STEM_B_2")] <- "Stem"
  seurat_filtered@meta.data$merged_sample[seurat_filtered@meta.data$merged_sample %in% c("CONTROL")] <- "Control"
  seurat_filtered@meta.data$merged_sample <- droplevels(seurat_filtered@meta.data$merged_sample)
  desired_order <- c("Stem", "ALDH", "STEM_A", "Control", "Undetermined")
  seurat_filtered@meta.data$merged_sample <- factor(seurat_filtered$merged_sample, levels = desired_order)
  
  levels(seurat_filtered@meta.data$combined_sample) <- c(levels(seurat_filtered@meta.data$combined_sample), "Stem", "Control")
  seurat_filtered@meta.data$combined_sample[seurat_filtered@meta.data$combined_sample %in% c("STEM_B")] <- "Stem"
  seurat_filtered@meta.data$combined_sample[seurat_filtered@meta.data$combined_sample %in% c("CONTROL")] <- "Control"
  seurat_filtered@meta.data$combined_sample <- droplevels(seurat_filtered@meta.data$combined_sample)
  desired_order <- c("Stem", "ALDH", "Control", "Undetermined")
  seurat_filtered@meta.data$combined_sample <- factor(seurat_filtered$combined_sample, levels = desired_order)
  
  ### store meta_sample c("STEM", "CONTROL", "Undetermined")
  levels(seurat_filtered@meta.data$meta_sample) <- c(levels(seurat_filtered@meta.data$meta_sample), "Stem", "Control")
  seurat_filtered@meta.data$meta_sample[seurat_filtered@meta.data$meta_sample %in% c("STEM")] <- "Stem"
  seurat_filtered@meta.data$meta_sample[seurat_filtered@meta.data$meta_sample %in% c("CONTROL")] <- "Control"
  seurat_filtered@meta.data$meta_sample <- droplevels(seurat_filtered@meta.data$meta_sample)
  desired_order <- c("Stem", "Control", "Undetermined")
  seurat_filtered$meta_sample <- factor(seurat_filtered$meta_sample, levels = desired_order)
  return(seurat_filtered)
}

purge4_seurat <- function(seurat) {
  seurat_filtered <- seurat
  levels(seurat_filtered@meta.data$merged_sample) <- c(levels(seurat_filtered@meta.data$merged_sample), "Stem", "Stem 30 dpt", "Control")
  seurat_filtered@meta.data$merged_sample[seurat_filtered@meta.data$merged_sample %in% c("STEM_B")] <- "Stem"
  seurat_filtered@meta.data$merged_sample[seurat_filtered@meta.data$merged_sample %in% c("STEM_A")] <- "Stem 30 dpt"
  seurat_filtered@meta.data$merged_sample[seurat_filtered@meta.data$merged_sample %in% c("CONTROL")] <- "Control"
  seurat_filtered@meta.data$merged_sample <- droplevels(seurat_filtered@meta.data$merged_sample)
  desired_order <- c("Stem", "ALDH", "Stem 30 dpt", "Control", "Undetermined")
  seurat_filtered@meta.data$merged_sample <- factor(seurat_filtered$merged_sample, levels = desired_order)
  
  levels(seurat_filtered@meta.data$combined_sample) <- c(levels(seurat_filtered@meta.data$combined_sample), "Stem", "Control")
  seurat_filtered@meta.data$combined_sample[seurat_filtered@meta.data$combined_sample %in% c("STEM_B")] <- "Stem"
  seurat_filtered@meta.data$combined_sample[seurat_filtered@meta.data$combined_sample %in% c("STEM_A", "CONTROL")] <- "Control"
  seurat_filtered@meta.data$combined_sample <- droplevels(seurat_filtered@meta.data$combined_sample)
  desired_order <- c("Stem", "ALDH", "Control", "Undetermined")
  seurat_filtered@meta.data$combined_sample <- factor(seurat_filtered$combined_sample, levels = desired_order)
  
  ### store meta_sample c("STEM", "CONTROL", "Undetermined")
  levels(seurat_filtered@meta.data$meta_sample) <- c(levels(seurat_filtered@meta.data$meta_sample), "Stem", "Control")
  seurat_filtered@meta.data$meta_sample[seurat_filtered@meta.data$meta_sample %in% c("STEM")] <- "Stem"
  seurat_filtered@meta.data$meta_sample[seurat_filtered@meta.data$meta_sample %in% c("CONTROL")] <- "Control"
  seurat_filtered@meta.data$meta_sample <- droplevels(seurat_filtered@meta.data$meta_sample)
  desired_order <- c("Stem", "Control", "Undetermined")
  seurat_filtered$meta_sample <- factor(seurat_filtered$meta_sample, levels = desired_order)
  
  seurat_filtered@meta.data$sample <- droplevels(seurat_filtered@meta.data$sample)
  seurat_filtered@meta.data$merged_sample <- droplevels(seurat_filtered@meta.data$merged_sample)
  seurat_filtered@meta.data$combined_sample <- droplevels(seurat_filtered@meta.data$combined_sample)
  seurat_filtered@meta.data$meta_sample <- droplevels(seurat_filtered@meta.data$meta_sample)
  return(seurat_filtered)
}

filtering_seurat <- function(seurat) {
  seurat_filtered <- subset(seurat, subset = sample != "Undetermined")
  seurat_filtered@meta.data$sample <- droplevels(seurat_filtered@meta.data$sample)
  seurat_filtered@meta.data$merged_sample <- droplevels(seurat_filtered@meta.data$merged_sample)
  seurat_filtered@meta.data$combined_sample <- droplevels(seurat_filtered@meta.data$combined_sample)
  seurat_filtered@meta.data$meta_sample <- droplevels(seurat_filtered@meta.data$meta_sample)
  
  # desired_order <- c("STEM_B", "ALDH", "STEM_A", "CONTROL")
  # seurat_filtered@meta.data$merged_sample <- factor(seurat_filtered$merged_sample, levels = desired_order)
  # 
  # desired_order <- c("STEM_B", "ALDH", "CONTROL")
  # seurat_filtered@meta.data$combined_sample <- factor(seurat_filtered$combined_sample, levels = desired_order)
  # 
  # desired_order <- c("STEM", "CONTROL")
  # seurat_filtered$meta_sample <- factor(seurat_filtered$meta_sample, levels = desired_order)
  return(seurat_filtered)
}

seurat <- readRDS("/scratch/users/jiamuyu/proj_botryllus/scRNAseq/250125_00_ciona_stemcell/raw.purged.filtered.rds")
# saveRDS(seurat, "/scratch/users/jiamuyu/proj_botryllus/scRNAseq/250125_00_ciona_stemcell/raw.purged.filtered.rds")

### organise the metrics used for filtering
cells <- tibble(
  cell = colnames(seurat),
  merged_sample = as.factor(seurat@meta.data$merged_sample),
  nCount = seurat@meta.data$nCount_RNA,
  nFeature = seurat@meta.data$nFeature_RNA,
  percent_MT = seurat@meta.data$percent.mt
)

averagenCountnFeature <- function(cells) {
  output <- tibble(
    'sample' = character(),
    'mean(nCount)' = numeric(),
    'median(nCount)' = numeric(),
    'mean(nFeature)' = numeric(),
    'median(nFeature)' = numeric()
  )
  for ( i in levels(cells$sample) ) {
    tmp <- tibble(
      'sample' = i,
      'mean(nCount)' = cells %>% filter(sample == i) %>% pull(nCount) %>% mean(),
      'median(nCount)' = cells %>% filter(sample == i) %>% pull(nCount) %>% median(),
      'mean(nFeature)' = cells %>% filter(sample == i) %>% pull(nFeature) %>% mean(),
      'median(nFeature)' = cells %>% filter(sample == i) %>% pull(nFeature) %>% median()
    )
    output <- bind_rows(output, tmp)
  }
  return(output)
}

averagenCountnFeature(cells) %>% knitr::kable()

# |sample       | mean(nCount)| median(nCount)| mean(nFeature)| median(nFeature)|
# |:------------|------------:|--------------:|--------------:|----------------:|
# |ALDH         |     2860.726|         1524.5|       615.7336|            443.0|
# |CONTROL      |     2231.769|         1158.5|       624.8353|            465.0|
# |STEM_A       |     2335.527|         1186.5|       648.4033|            473.5|
# |STEM_B_1     |     1973.400|         1221.0|       507.3091|            371.0|
# |STEM_B_2     |     2327.187|         1359.0|       541.5212|            383.0|
# |Undetermined |     2015.537|         1294.0|       544.0640|            408.0|

# single_filtered
# |sample  | mean(nCount)| median(nCount)| mean(nFeature)| median(nFeature)|
# |:-------|------------:|--------------:|--------------:|----------------:|
# |STEM_B  |     2312.870|         1932.0|       552.1753|            412.0|
# |ALDH    |     1905.490|         1429.0|       522.9364|            413.0|
# |STEM_A  |     1613.825|         1158.5|       588.7235|            469.0|
# |CONTROL |     1577.288|         1107.0|       579.5664|            465.5|

### filter doublets
sce <- SingleCellExperiment(assays = list(counts = seurat@assays$RNA$counts))
sce <- scDblFinder(sce)

cells$multiplet_class <- colData(sce)$scDblFinder.class

# Calculate the number of cells per sample
cell_counts <- cells %>%
  group_by(merged_sample) %>%
  dplyr::summarise(cell_count = n(), .groups = 'drop')

# custom_order <- c("ALDH (n = 1036)", "STEM_B_1 (n = 1359)", "STEM_B_2 (n = 1606)", "CONTROL (n = 1506)", "STEM_A (n = 1644)", "Undetermined (n = 4200)") 

# Create a new label with sample names and counts
cell_counts <- cell_counts %>%
  mutate(sample_label = paste0(merged_sample, " (n = ", cell_count, ")"))

# Modify the original data to use the new labels
cells <- cells %>%
  left_join(cell_counts %>% dplyr::select(merged_sample, sample_label), by = "merged_sample") %>% 
  mutate(sample_label = factor(sample_label))

p <- cells %>%
  filter(multiplet_class != 'singlet') %>%
  group_by(sample_label) %>%
  dplyr::summarize(count = n()) %>%
  ggplot(aes(x = sample_label, y = count, fill = sample_label)) +
  geom_col(color = 'black') +
  theme_bw() +
  scale_fill_manual(values = custom_colors$discrete) +
  scale_x_discrete(limits = rev(levels(cells$sample_label))) +
  scale_y_continuous(name = 'Number of doublets', labels = scales::comma) +
  theme(
    axis.title.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    legend.position = 'none'
  ) +
  coord_flip()

ggsave(
  file.path(working_dir, 'qc_number_of_doublets_by_sample.png'),
  # file.path(working_dir, 'qc_number_of_doublets_by_sample.filtered.png'),
  p, height = 3, width = 6
)

p1 <- ggplot(cells, aes(x = multiplet_class, y = nCount, fill = multiplet_class)) +
  geom_violin(draw_quantiles = c(0.5), scale = 'area', trim = TRUE) +
  theme_bw() +
  scale_fill_manual(values = custom_colors$discrete) +
  scale_x_discrete(limits = rev(levels(cells$multiplet_class))) +
  scale_y_continuous(labels = scales::comma) +
  labs(title = 'Number of transcripts', subtitle = 'linear scale') +
  theme(
    axis.title = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    legend.position = 'none'
  ) +
  coord_flip()

p2 <- ggplot(cells, aes(x = multiplet_class, y = nCount, fill = multiplet_class)) +
  geom_violin(draw_quantiles = c(0.5), scale = 'area', trim = TRUE) +
  theme_bw() +
  scale_fill_manual(values = custom_colors$discrete) +
  scale_x_discrete(limits = rev(levels(cells$multiplet_class))) +
  scale_y_log10(labels = scales::comma) +
  labs(title = 'Number of transcripts', subtitle = 'log-scale') +
  theme(
    axis.title = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    legend.position = 'none'
  ) +
  coord_flip()

p3 <- ggplot(cells, aes(x = multiplet_class, y = nFeature, fill = multiplet_class)) +
  geom_violin(draw_quantiles = c(0.5), scale = 'area', trim = TRUE) +
  theme_bw() +
  scale_fill_manual(values = custom_colors$discrete) +
  scale_x_discrete(limits = rev(levels(cells$multiplet_class))) +
  scale_y_continuous(labels = scales::comma) +
  labs(title = 'Number of expressed genes', subtitle = 'linear scale') +
  theme(
    axis.title = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    legend.position = 'none'
  ) +
  coord_flip()

p4 <- ggplot(cells, aes(x = multiplet_class, y = nFeature, fill = multiplet_class)) +
  geom_violin(draw_quantiles = c(0.5), scale = 'area', trim = TRUE) +
  theme_bw() +
  scale_fill_manual(values = custom_colors$discrete) +
  scale_x_discrete(limits = rev(levels(cells$multiplet_class))) +
  scale_y_log10(labels = scales::comma) +
  labs(title = 'Number of expressed genes', subtitle = 'log-scale') +
  theme(
    axis.title = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    legend.position = 'none'
  ) +
  coord_flip()

ggsave(
  file.path(working_dir, 'qc_ncount_nfeature_by_multiplet.png'),
  # file.path(working_dir, 'qc_ncount_nfeature_by_multiplet.filtered.png'),
  p1 + p3 +
    p2 + p4 + plot_layout(ncol = 2),
  height = 7, width = 7
)

cells <- filter(cells, multiplet_class == 'singlet')

averagenCountnFeature(cells) %>% knitr::kable()

# |sample       | mean(nCount)| median(nCount)| mean(nFeature)| median(nFeature)|
# |:------------|------------:|--------------:|--------------:|----------------:|
# |ALDH         |     2745.892|         1436.0|       579.4779|            427.0|
# |CONTROL      |     2231.769|         1158.5|       624.8353|            465.0|
# |STEM_A       |     2335.527|         1186.5|       648.4033|            473.5|
# |STEM_B_1     |     1942.926|         1170.0|       491.4336|            361.0|
# |STEM_B_2     |     2236.114|         1256.0|       507.2828|            369.0|
# |Undetermined |     1985.134|         1259.0|       534.7416|            402.0|

# Calculate the number of cells per sample
cells <- subset(cells, select = -sample_label)

cell_counts <- cells %>%
  group_by(merged_sample) %>%
  dplyr::summarise(cell_count = n(), .groups = 'drop')

# Create a new label with sample names and counts
cell_counts <- cell_counts %>%
  mutate(sample_label = paste0(merged_sample, " (n = ", cell_count, ")"))

cell_counts$sample_label
# custom_order <- c("ALDH (n = 1036)", "STEM_B_1 (n = 1359)", "STEM_B_2 (n = 1606)", "CONTROL (n = 1506)", "STEM_A (n = 1644)", "Undetermined (n = 4200)") 
# custom_order <- c("STEM_B_1 (n = 1273)", "STEM_B_2 (n = 1471)", "ALDH (n = 973)", "CONTROL (n = 1506)", "STEM_A (n = 1644)", "Undetermined (n = 4094)") 
custom_order <- c("STEM_B (n = 1797)", "ALDH (n = 723)", "STEM_A (n = 1226)", "CONTROL (n = 1144)")

# Modify the original data to use the new labels
cells <- cells %>%
  left_join(cell_counts %>% dplyr::select(merged_sample, sample_label), by = "merged_sample") %>% 
  mutate(sample_label = factor(sample_label, levels = custom_order))

p1 <- ggplot(cells, aes(x = sample_label, y = nCount, fill = sample_label)) +
  geom_violin(draw_quantiles = c(0.5), scale = 'area', trim = TRUE) +
  theme_bw() +
  scale_fill_manual(values = custom_colors$discrete) +
  scale_x_discrete(limits = rev(levels(cells$sample_label))) +
  scale_y_continuous(labels = scales::comma) +
  labs(title = 'Number of transcripts', subtitle = 'linear scale') +
  theme(
    axis.title = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    legend.position = 'none'
  ) +
  coord_flip()

p2 <- ggplot(cells, aes(x = sample_label, y = nCount, fill = sample_label)) +
  geom_violin(draw_quantiles = c(0.5), scale = 'area', trim = TRUE) +
  theme_bw() +
  scale_fill_manual(values = custom_colors$discrete) +
  scale_x_discrete(limits = rev(levels(cells$sample_label))) +
  scale_y_log10(labels = scales::comma) +
  labs(title = 'Number of transcripts', subtitle = 'log-scale') +
  theme(
    axis.title = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    legend.position = 'none'
  ) +
  coord_flip()

p3 <- ggplot(cells, aes(x = sample_label, y = nFeature, fill = sample_label)) +
  geom_violin(draw_quantiles = c(0.5), scale = 'area', trim = TRUE) +
  theme_bw() +
  scale_fill_manual(values = custom_colors$discrete) +
  scale_x_discrete(limits = rev(levels(cells$sample_label))) +
  scale_y_continuous(labels = scales::comma) +
  labs(title = 'Number of expressed genes', subtitle = 'linear scale') +
  theme(
    axis.title = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    legend.position = 'none'
  ) +
  coord_flip()

p4 <- ggplot(cells, aes(x = sample_label, y = nFeature, fill = sample_label)) +
  geom_violin(draw_quantiles = c(0.5), scale = 'area', trim = TRUE) +
  theme_bw() +
  scale_fill_manual(values = custom_colors$discrete) +
  scale_x_discrete(limits = rev(levels(cells$sample_label))) +
  scale_y_log10(labels = scales::comma) +
  labs(title = 'Number of expressed genes', subtitle = 'log-scale') +
  theme(
    axis.title = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    legend.position = 'none'
  ) +
  coord_flip()

p5 <- ggplot(cells, aes(x = sample_label, y = percent_MT / 100, fill = sample_label)) +
  geom_violin(draw_quantiles = c(0.5), scale = 'area', trim = TRUE) +
  theme_bw() +
  scale_fill_manual(values = custom_colors$discrete) +
  scale_x_discrete(limits = rev(levels(cells$sample_label))) +
  scale_y_continuous(labels = scales::percent) +
  labs(title = 'Mitochondrial transcripts %', subtitle = 'linear scale') +
  theme(
    axis.title = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    legend.position = 'none'
  ) +
  coord_flip()

ggsave(
  # file.path(working_dir, 'qc_violin.png'),
  file.path(working_dir, 'qc_violin.single.png'),
  # file.path(working_dir, 'qc_violin.singlet.png'),
  p1 + p3 + p5 +
    p2 + p4 + plot_layout(ncol = 3),
  height = 8, width = 12
)

p1 <- ggplot(cells, aes(x = nCount, y = nFeature, color = sample_label)) +  # Remove color = percent_MT from aes()
  geom_point(size = 0.5, alpha = 0.5) +
  scale_color_manual(values = custom_colors$discrete) +
  scale_x_continuous(name = 'Number of transcripts', labels = scales::comma) +
  scale_y_continuous(name = 'Number of expressed genes', labels = scales::comma) +
  theme_bw() +
  geom_density_2d(color = 'black', size = 0.5, alpha = 0.6) +  # Density contours
  labs(title = 'Number of expressed genes vs. transcripts', subtitle = 'linear scale') +
  theme(legend.position = "none")

p2 <- ggplot(cells, aes(x = nCount, y = nFeature, color = sample_label)) +  # Remove color = percent_MT from aes()
  geom_point(size = 0.5, alpha = 0.5) +
  scale_color_manual(values = custom_colors$discrete) +
  scale_x_log10(name = 'Number of transcripts', labels = scales::comma) +
  scale_y_log10(name = 'Number of expressed genes', labels = scales::comma) +
  theme_bw() +
  geom_density_2d(color = 'black', size = 0.5, alpha = 0.6) +  # Density contours
  labs(title = 'Number of expressed genes vs. transcripts', subtitle = 'log-scale', color = 'Samples')

p3 <- ggplot(cells, aes(x = nCount, y = nFeature, color = percent_MT)) +
  geom_point(size = 0.5, alpha = 0.5) +
  scale_x_continuous(name = 'Number of transcripts', labels = scales::comma) +
  scale_y_continuous(name = 'Number of expressed genes', labels = scales::comma) +
  theme_bw() +
  scale_color_viridis(
    name = 'Mitochondrial\ntranscripts',
    limits = c(0, 100),
    labels = scales::percent_format(scale = 1),
    guide = guide_colorbar(frame.colour = 'black', ticks.colour = 'black')
  ) +
  labs(title = 'Number of expressed genes vs. transcripts', subtitle = 'linear scale') +
  theme(legend.position = "none")

p4 <- ggplot(cells, aes(x = nCount, y = nFeature, color = percent_MT)) +
  geom_point(size = 0.5, alpha = 0.5) +
  scale_x_log10(name = 'Number of transcripts', labels = scales::comma) +
  scale_y_log10(name = 'Number of expressed genes', labels = scales::comma) +
  theme_bw() +
  scale_color_viridis(
    name = 'Mitochondrial\ntranscripts',
    limits = c(0, 100),
    labels = scales::percent_format(scale = 1),
    guide = guide_colorbar(frame.colour = 'black', ticks.colour = 'black')
  ) +
  labs(title = 'Number of expressed genes vs. transcripts', subtitle = 'log-scale')

ggsave(
  # file.path(working_dir, 'qc_scatter.png'), 
  file.path(working_dir, 'qc_scatter.singlet.png'), 
  p1 + p2 + 
    p3 + p4 + plot_layout(ncol = 2), 
  height = 8, width = 10
)

# # Split the data frame by the 'sample' column
# sample_list <- split(cells, cells$sample)
# 
# # Create a list to store the plots
# plots <- list()
# 
# # Loop through each sample and create a plot
# for (sample_name in names(sample_list)) {
#   # Subset the data for the current sample
#   sample_data <- sample_list[[sample_name]]
#   
#   # Create the plot
#   p <- ggplot(sample_data, aes(x = nCount, y = nFeature, color = sample)) +
#     geom_point(size = 0.5) +
#     scale_color_manual(values = custom_colors$discrete) +
#     scale_x_log10(name = 'Number of transcripts', labels = scales::comma) +
#     scale_y_log10(name = 'Number of expressed genes', labels = scales::comma) +
#     theme_bw() +
#     labs(title = paste('Number of expressed genes vs. transcripts'), subtitle = sample_name) +
#     theme(legend.position = "none")
# 
#   # Add the plot to the list
#   plots[[sample_name]] <- p
# }
# 
# # Save the combined plot to a file
# ggsave(file.path(working_dir, 'qc_nfeature_over_ncount_thresholds_combined.png'), wrap_plots(plots, nrow = 2, ncol = 3), height = 8, width = 12)

## with contour
# Split the data frame by the 'sample' column
sample_list <- split(cells, cells$sample)

# Create a list to store the plots
plots <- list()

# Loop through each sample and create a plot
for (sample_name in names(sample_list)) {
  # Subset the data for the current sample
  sample_data <- sample_list[[sample_name]]
  
  # Create the plot with scatter points and density contours
  p <- ggplot(sample_data, aes(x = nCount, y = nFeature, color = percent_MT)) +
    geom_point(size = 0.5, alpha = 0.5) +  # Scatter plot with some transparency
    geom_density_2d(color = custom_colors$discrete[38], size = 0.5) +  # Density contours
    scale_color_viridis(
      name = 'Percent MT\ntranscripts',
      limits = c(0, 100),
      labels = scales::percent_format(scale = 1),
      guide = guide_colorbar(frame.colour = 'black', ticks.colour = 'black')
    ) +    
    scale_x_log10(name = 'Number of transcripts', labels = scales::comma) +
    scale_y_log10(name = 'Number of expressed genes', labels = scales::comma) +
    theme_bw() +
    labs(title = paste('Number of expressed genes vs. transcripts'), subtitle = sample_name) +
    labs(title = 'Number of expressed genes vs. transcripts', subtitle = paste0(merged_sample_name, ' (n = ', nrow(sample_data), ')')) +
    theme(legend.position = ifelse(sample_name == tail(names(sample_list), 1), "right", "none"))  
  # Add the plot to the list
  plots[[sample_name]] <- p
}

# Save the combined plot to a file
ggsave(file.path(working_dir, 'qc_scatter_by_sample.png'), wrap_plots(plots, nrow = 2, ncol = 3) + theme(legend.position = "right"), height = 8, width = 14)

### Filtering
## slope filter
# # Calculate the line values
# cells$line_value <- -0.2 + 0.85 * log10(cells$nCount)
# 
# # Create a new column to categorize points
# cells$point_color <- ifelse(log10(cells$nFeature) > cells$line_value, "Above", "Below")
# 
# # Define the colors for the points
# point_colors <- c("Above" = custom_colors$discrete[1], "Below" = "black")
# 
# # Calculate the percentage of points not being filtered
# percent_not_filtered <- round(100 * sum(cells$point_color == "Above") / nrow(cells), 2)
# 
# point_labels <- c("Above" = paste0("Unfiltered (", percent_not_filtered, "%)"), "Below" = paste0("Filtered (", 100 - percent_not_filtered, "%)"))
# 
# p1 <- ggplot(cells, aes(x = nCount, y = nFeature, color = point_color)) +  # Remove color = percent_MT from aes()
#   geom_point(size = 0.5, alpha = 0.5) +
#   scale_color_manual(values = point_colors, labels = point_labels) +
#   scale_x_continuous(name = 'Number of transcripts', labels = scales::comma) +
#   scale_y_continuous(name = 'Number of expressed genes', labels = scales::comma) +
#   theme_bw() +
#   geom_density_2d(color = custom_colors$discrete[38], size = 0.5) +  # Density contours
#   labs(title = 'Number of expressed genes vs. transcripts', subtitle = paste0('linear scale (n = ', nrow(cells), ')')) +
#   theme(legend.position = "none")
# 
# p2 <- ggplot(cells, aes(x = nCount, y = nFeature, color = point_color)) +
#   geom_point(size = 0.5, alpha = 0.5) +
#   scale_color_manual(values = point_colors, labels = point_labels) +
#   scale_x_log10(name = 'Number of transcripts', labels = scales::comma) +
#   scale_y_log10(name = 'Number of expressed genes', labels = scales::comma) +
#   theme_bw() +
#   geom_density_2d(color = custom_colors$discrete[38], size = 0.5) +  # Density contours
#   labs(title = 'Number of expressed genes vs. transcripts', subtitle = paste0('log-scale (n = ', nrow(cells), ')')) +
#   theme(legend.position = "right") +
#   guides(color = guide_legend(title = NULL))
# 
# ggsave(file.path(working_dir, 'qc_scatter_thresholds.png'), p1 + p2, height = 4, width = 10)
# 
# cells <- cells %>%
#   filter(log10(nFeature) > line_value)

## naive filter
median_nCount <- median(cells$nCount)
# 1089 after slope filter
# 1271 before slope filter
mad_nCount <- mad(cells$nCount)
# 711 after slope filter
# 944 before slope filter

median_nFeature <- median(cells$nFeature)
# 458 after slope filter
# 418 before slope filter
mad_nFeature <- mad(cells$nFeature)
# 260 after slope filter
# 243 before slope filter

median_percent_MT <- median(cells$percent_MT)
# 9.948542 after slope filter
# 9.948542 before slope filter
mad_percent_MT <- mad(cells$percent_MT)
# 9.812097 after slope filter
# 9.812097 before slope filter

thresholds_nCount <- c(0, median_nCount + 5*mad_nCount)
thresholds_nFeature <- c(0, median_nFeature + 5*mad_nFeature)
thresholds_percent_MT <- c(0, median_percent_MT + 5*mad_percent_MT)

p1 <- ggplot(cells, aes(x = sample_label, y = nCount, fill = sample_label)) +
  geom_violin(draw_quantiles = c(0.5), scale = 'area', trim = TRUE) +
  geom_hline(yintercept = median_nCount, color = 'black') +
  geom_hline(yintercept = thresholds_nCount, color = 'red') +
  theme_bw() +
  scale_fill_manual(values = custom_colors$discrete) +
  scale_x_discrete(limits = rev(levels(cells$sample_label))) +
  scale_y_continuous(labels = scales::comma) +
  labs(title = 'Number of transcripts', subtitle = 'linear scale') +
  theme(
    axis.title = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    legend.position = 'none'
  ) +
  coord_flip()

p2 <- ggplot(cells, aes(x = sample_label, y = nCount, fill = sample_label)) +
  geom_violin(draw_quantiles = c(0.5), scale = 'area', trim = TRUE) +
  geom_hline(yintercept = median_nCount, color = 'black') +
  geom_hline(yintercept = pmax(thresholds_nCount), color = 'red') +
  theme_bw() +
  scale_fill_manual(values = custom_colors$discrete) +
  scale_x_discrete(limits = rev(levels(cells$sample_label))) +
  scale_y_log10(labels = scales::comma) +
  labs(title = 'Number of transcripts', subtitle = 'log-scale') +
  theme(
    axis.title = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    legend.position = 'none'
  ) +
  coord_flip()

p3 <- ggplot(cells, aes(x = sample_label, y = nFeature, fill = sample_label)) +
  geom_violin(draw_quantiles = c(0.5), scale = 'area', trim = TRUE) +
  geom_hline(yintercept = median_nFeature, color = 'black') +
  geom_hline(yintercept = thresholds_nFeature, color = 'red') +
  theme_bw() +
  scale_fill_manual(values = custom_colors$discrete) +
  scale_x_discrete(limits = rev(levels(cells$sample_label))) +
  scale_y_continuous(labels = scales::comma) +
  labs(title = 'Number of expressed genes', subtitle = 'linear scale') +
  theme(
    axis.title = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    legend.position = 'none'
  ) +
  coord_flip()

p4 <- ggplot(cells, aes(x = sample_label, y = nFeature, fill = sample_label)) +
  geom_violin(draw_quantiles = c(0.5), scale = 'area', trim = TRUE) +
  geom_hline(yintercept = median_nFeature, color = 'black') +
  geom_hline(yintercept = pmax(thresholds_nFeature), color = 'red') +
  theme_bw() +
  scale_fill_manual(values = custom_colors$discrete) +
  scale_x_discrete(limits = rev(levels(cells$sample_label))) +
  scale_y_log10(labels = scales::comma) +
  labs(title = 'Number of expressed genes', subtitle = 'log-scale') +
  theme(
    axis.title = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    legend.position = 'none'
  ) +
  coord_flip()

p5 <- ggplot(cells, aes(x = sample_label, y = percent_MT, fill = sample_label)) +
  geom_violin(draw_quantiles = c(0.5), scale = 'area', trim = TRUE) +
  geom_hline(yintercept = median_percent_MT, color = 'black') +
  geom_hline(yintercept = thresholds_percent_MT, color = 'red') +
  theme_bw() +
  scale_fill_manual(values = custom_colors$discrete) +
  scale_x_discrete(limits = rev(levels(cells$sample_label))) +
  scale_y_continuous(labels = scales::percent) +
  labs(title = 'Mitochondrial transcripts %', subtitle = 'linear scale') +
  theme(
    axis.title = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    legend.position = 'none'
  ) +
  coord_flip()

ggsave(
  # file.path(working_dir, 'qc_violin_thresholds.single.png'),
  # file.path(working_dir, 'qc_violin_thresholds.singlet.png'),
  file.path(working_dir, 'qc_violin_thresholds.5.singlet.png'),
  p1 + p3 + p5 +
    p2 + p4 + plot_layout(ncol = 3),
  height = 8, width = 12
)

filtered_cells <- cells %>%
  dplyr::filter(
    nCount >= thresholds_nCount[1],
    nCount <= thresholds_nCount[2],
    nFeature >= thresholds_nFeature[1],
    nFeature <= thresholds_nFeature[2],
    percent_MT >= thresholds_percent_MT[1],
    percent_MT <= thresholds_percent_MT[2]
  )

p1 <- ggplot(filtered_cells, aes(x = merged_sample, y = nCount, fill = sample)) +
  geom_violin(draw_quantiles = c(0.5), scale = 'area', trim = TRUE) +
  geom_hline(yintercept = median_nCount, color = 'black') +
  geom_hline(yintercept = thresholds_nCount, color = 'red') +
  theme_bw() +
  scale_fill_manual(values = custom_colors$discrete) +
  scale_x_discrete(limits = rev(levels(filtered_cells$sample))) +
  scale_y_continuous(labels = scales::comma) +
  labs(title = 'Number of transcripts', subtitle = 'linear scale') +
  theme(
    axis.title = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    legend.position = 'none'
  ) +
  coord_flip()

p2 <- ggplot(filtered_cells, aes(x = merged_sample, y = nCount, fill = sample)) +
  geom_violin(draw_quantiles = c(0.5), scale = 'area', trim = TRUE) +
  geom_hline(yintercept = median_nCount, color = 'black') +
  geom_hline(yintercept = thresholds_nCount, color = 'red') +
  theme_bw() +
  scale_fill_manual(values = custom_colors$discrete) +
  scale_x_discrete(limits = rev(levels(filtered_cells$sample))) +
  scale_y_log10(labels = scales::comma) +
  labs(title = 'Number of transcripts', subtitle = 'log-scale') +
  theme(
    axis.title = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    legend.position = 'none'
  ) +
  coord_flip()

p3 <- ggplot(filtered_cells, aes(x = merged_sample, y = nFeature, fill = sample)) +
  geom_violin(draw_quantiles = c(0.5), scale = 'area', trim = TRUE) +
  geom_hline(yintercept = median_nFeature, color = 'black') +
  geom_hline(yintercept = thresholds_nFeature, color = 'red') +
  theme_bw() +
  scale_fill_manual(values = custom_colors$discrete) +
  scale_x_discrete(limits = rev(levels(filtered_cells$sample))) +
  scale_y_continuous(labels = scales::comma) +
  labs(title = 'Number of expressed genes', subtitle = 'linear scale') +
  theme(
    axis.title = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    legend.position = 'none'
  ) +
  coord_flip()

p4 <- ggplot(filtered_cells, aes(x = merged_sample, y = nFeature, fill = sample)) +
  geom_violin(draw_quantiles = c(0.5), scale = 'area', trim = TRUE) +
  geom_hline(yintercept = median_nFeature, color = 'black') +
  geom_hline(yintercept = thresholds_nFeature, color = 'red') +
  theme_bw() +
  scale_fill_manual(values = custom_colors$discrete) +
  scale_x_discrete(limits = rev(levels(filtered_cells$sample))) +
  scale_y_log10(labels = scales::comma) +
  labs(title = 'Number of expressed genes', subtitle = 'log-scale') +
  theme(
    axis.title = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    legend.position = 'none'
  ) +
  coord_flip()

p5 <- ggplot(filtered_cells, aes(x = merged_sample, y = percent_MT, fill = sample)) +
  geom_violin(draw_quantiles = c(0.5), scale = 'area', trim = TRUE) +
  geom_hline(yintercept = median_percent_MT, color = 'black') +
  geom_hline(yintercept = thresholds_percent_MT, color = 'red') +
  theme_bw() +
  scale_fill_manual(values = custom_colors$discrete) +
  scale_x_discrete(limits = rev(levels(cells$sample))) +
  scale_y_continuous(labels = scales::percent) +
  labs(title = 'Mitochondrial transcripts %', subtitle = 'linear scale') +
  theme(
    axis.title = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    legend.position = 'none'
  ) +
  coord_flip()

ggsave(
  # file.path(working_dir, 'qc_violin_filtered.single.png'),
  # file.path(working_dir, 'qc_violin_filtered.singlet.png'),
  file.path(working_dir, 'qc_violin_filtered.5.singlet.png'),
  p1 + p3 + p5 +
    p2 + p4 + plot_layout(ncol = 3),
  height = 8, width = 12
)

p1 <- ggplot(filtered_cells, aes(x = nCount, y = nFeature, color = percent_MT)) +  # Remove color = percent_MT from aes()
  geom_point(size = 0.5, alpha = 0.5) +
  geom_hline(yintercept = thresholds_nFeature, color = 'red') +
  geom_vline(xintercept = thresholds_nCount, color = 'red') +
  scale_x_continuous(name = 'Number of transcripts', labels = scales::comma) +
  scale_y_continuous(name = 'Number of expressed genes', labels = scales::comma) +
  theme_bw() +
  geom_density_2d(color = custom_colors$discrete[38], size = 0.5) +  # Density contours
  scale_color_viridis(
    name = 'Percent MT\ntranscripts',
    limits = c(0, 100),
    labels = scales::percent_format(scale = 1),
    guide = guide_colorbar(frame.colour = 'black', ticks.colour = 'black')
  ) +
  labs(title = 'Number of expressed genes vs. transcripts', subtitle = paste0('linear scale (n = ', nrow(filtered_cells), ')')) +
  theme(legend.position = "none")

p2 <- ggplot(filtered_cells, aes(x = nCount, y = nFeature, color = percent_MT)) +  # Remove color = percent_MT from aes()
  geom_point(size = 0.5, alpha = 0.5) +
  geom_hline(yintercept = thresholds_nFeature, color = 'red') +
  geom_vline(xintercept = thresholds_nCount, color = 'red') +
  scale_x_log10(name = 'Number of transcripts', labels = scales::comma) +
  scale_y_log10(name = 'Number of expressed genes', labels = scales::comma) +
  theme_bw() +
  geom_density_2d(color = custom_colors$discrete[38], size = 0.5) +  # Density contours
  scale_color_viridis(
    name = 'Percent MT\ntranscripts',
    limits = c(0, 100),
    labels = scales::percent_format(scale = 1),
    guide = guide_colorbar(frame.colour = 'black', ticks.colour = 'black')
  ) +
  labs(title = 'Number of expressed genes vs. transcripts', subtitle = paste0('linear scale (n = ', nrow(filtered_cells), ')'))

ggsave(
  # file.path(working_dir, 'qc_scatter_filtered.single.png'), 
  file.path(working_dir, 'qc_scatter_filtered.5.singlet.png'), 
  p1 + p2, 
  height = 6, width = 10
)

cells_to_keep <- filtered_cells$cell
length(cells_to_keep)
# 6588 after slope filter
# 8465 before slope filter
# 8119 singlet single filter
# 10046 5.singlet single filter
ncol(seurat)
# 11351

### if the layers haven't been merged
# layers <- seurat@assays$RNA@layers
# 
# # Initialize a list to store binary matrices for each layer
# binary_matrices <- vector("list", length(layers))
# 
# # Convert each layer to a binary matrix and store it in the list
# for (i in seq_along(layers)) {
#   binary_matrices[[i]] <- as.matrix(layers[[i]])
#   binary_matrices[[i]][binary_matrices[[i]] > 0] <- 1  # Convert to binary
#   binary_matrices[[i]] <- rowSums(binary_matrices[[i]])
# }
# 
# # Sum the binary matrices across all layers
# gene_expression_counts <- Reduce("+", binary_matrices)
# names(gene_expression_counts) <- rownames(seurat)

### if the layers have been merged
gene_expression_counts <- as.matrix(seurat[['RNA']]$counts)
gene_expression_counts[gene_expression_counts > 0] <- 1
gene_expression_counts <- rowSums(gene_expression_counts)

filtered_genes <- gene_expression_counts[gene_expression_counts > 5]

genes_to_keep <- names(filtered_genes)

length(genes_to_keep)
# 14144 singlet single filter
# 14144 5.singlet single filter
nrow(seurat)
# 18827

sorted_counts <- sort(gene_expression_counts, decreasing = TRUE)

# Create an index for the sorted genes
gene_indices <- seq_along(sorted_counts)

# Calculate cumulative percentage
cumulative_percentage <- cumsum(sorted_counts) / sum(sorted_counts) * 100

counts_df <- data.frame(
  Gene_Index = gene_indices,
  Expression_Counts = sorted_counts,
  Cumulative_Percentage = cumulative_percentage,
  Gene_Names = names(sorted_counts)
)

p <- ggplot(counts_df, aes(x = Gene_Index, y = Expression_Counts)) +
  geom_bar(stat = "identity", fill = custom_colors$discrete[24], color = NA) +
  geom_line(aes(y = Cumulative_Percentage * max(Expression_Counts) / 100), 
            color = custom_colors$discrete[1], size = 1) +
  scale_y_continuous(sec.axis = sec_axis(~ . * 100 / max(counts_df$Expression_Counts), name = "Cumulative Percentage (%)")) +
  labs(title = "Gene Expression Counts",
       x = "Sorted Gene Index",
       y = "Number of Cells Expressing the Gene") +
  theme_minimal() +
  scale_x_continuous(breaks = seq(0, length(sorted_counts), by = 3000))

ggsave(
  # file.path(working_dir, "gene_elbow_plot.png"), 
  # file.path(working_dir, "gene_elbow_plot.singlet.png"), 
  file.path(working_dir, "gene_elbow_plot.5.singlet.png"), 
  p, width = 4, height = 6
)

# Create a tibble to track cells per sample before and after filtering
cells_per_sample_after_filtering <- tibble(
  Sample = levels(as.factor(seurat@meta.data$sample)),
  before = numeric(length(levels(as.factor(seurat@meta.data$sample)))),
  after = numeric(length(levels(as.factor(seurat@meta.data$sample))))
)

# Count the number of cells for each sample before filtering
for (i in seq_along(cells_per_sample_after_filtering$Sample)) {
  sample_name <- cells_per_sample_after_filtering$Sample[i]
  
  # Count cells for the current sample before filtering
  cells_per_sample_after_filtering$before[i] <- sum(seurat@meta.data$sample == sample_name)
}

# Filter cells and genes in the Seurat object
seurat <- seurat[genes_to_keep, cells_to_keep]

# Count the number of cells for each sample after filtering
for (i in seq_along(cells_per_sample_after_filtering$Sample)) {
  sample_name <- cells_per_sample_after_filtering$Sample[i]
  
  # Count cells for the current sample after filtering
  cells_per_sample_after_filtering$after[i] <- sum(seurat@meta.data$sample == sample_name)
}

# Display the table using knitr::kable
knitr::kable(cells_per_sample_after_filtering)
### double filtering
# |Sample       | before| after|
# |:------------|------:|-----:|
# |ALDH         |   1036|   485|
# |CONTROL      |   1506|  1096|
# |STEM_A       |   1644|  1168|
# |STEM_B_1     |   1359|   554|
# |STEM_B_2     |   1606|   557|
# |Undetermined |   4200|  2728|

### single slope filtering
# |Sample       | before| after|
# |:------------|------:|-----:|
# |ALDH         |   1036|   742|
# |CONTROL      |   1506|  1407|
# |STEM_A       |   1644|  1538|
# |STEM_B_1     |   1359|  1071|
# |STEM_B_2     |   1606|  1199|
# |Undetermined |   4200|  3328|

### 1.5 singlet, without slope filtering
# |Sample       | before| after|
# |:------------|------:|-----:|
# |ALDH         |   1036|   723|
# |CONTROL      |   1506|  1144|
# |STEM_A       |   1644|  1226|
# |STEM_B_1     |   1359|   844|
# |STEM_B_2     |   1606|   953|
# |Undetermined |   4200|  3575|

### 5 singlet, without slope filtering
# |Sample       | before| after|
# |:------------|------:|-----:|
# |ALDH         |   1036|   823|
# |CONTROL      |   1506|  1371|
# |STEM_A       |   1644|  1478|
# |STEM_B_1     |   1359|  1182|
# |STEM_B_2     |   1606|  1310|
# |Undetermined |   4200|  3882|

# saveRDS(seurat, file.path(working_dir, "filtered.filtered.rds"))
# seurat <- readRDS(file.path(working_dir, "filtered.filtered.rds"))

# saveRDS(seurat, file.path(working_dir, "filtered.single.rds"))
seurat <- readRDS(file.path(working_dir, "filtered.single.rds"))

# saveRDS(seurat, file.path(working_dir, "gfp.filtered.single.rds"))
seurat <- readRDS(file.path(working_dir, "gfp.filtered.single.rds"))

# saveRDS(seurat, file.path(working_dir, "filtered.1.5.single.rds"))
seurat <- readRDS(file.path(working_dir, "filtered.1.5.singlet.rds"))

# saveRDS(seurat, file.path(working_dir, "filtered.5.single.rds"))
seurat <- readRDS(file.path(working_dir, "filtered.5.singlet.rds"))

### Processing
## integration
{
  seurat_list <- SplitObject(seurat, split.by = 'sample')
  
  # seurat_list <- lapply(
  #   X = seurat_list,
  #   FUN = function(x) {
  #     x <- SCTransform(x)
  # return(x)
  #   }
  # )
  
  # seurat_list <- lapply(
  #   X = seurat_list,
  #   FUN = function(x) {
  #     x <- SCTransform(
  #       x,
  #       assay = 'RNA',
  #       new.assay.name = 'SCT',
  #       vars.to.regress = c('percent.mt', 'nFeature_RNA', 'nCount_RNA')
  #     )
  #     return(x)
  #   }
  # )
  
  seurat_list <- lapply(
    X = seurat_list,
    FUN = function(x) {
      x <- NormalizeData(
        x,
        assay = 'RNA',
        normalization.method = 'LogNormalize',
        scale.factor = median(x@meta.data$nCount_RNA)
      )
      return(x)
    }
  )
  
  # new_s_genes <- read.csv(file.path(working_dir, "s_genes_mapping.csv"))$mapped
  # new_g2m_genes <- read.csv(file.path(working_dir, "g2m_genes_mapping.csv"))$mapped
  # 
  # seurat_list <- lapply(
  #   X = seurat_list,
  #   FUN = function(x) {
  #     x <- CellCycleScoring(
  #       x,
  #       assay = 'RNA',
  #       s.features = new_s_genes,
  #       g2m.features = new_g2m_genes
  #     )
  #     
  #     x$CC.Difference <- x$S.Score - x$G2M.Score
  #     
  #     return(x)  
  #   }
  # )
  
  seurat_list <- lapply(
    X = seurat_list,
    FUN = function(x) {
      x <- SCTransform(
        x,
        assay = 'RNA',
        new.assay.name = 'SCT',
        # vars.to.regress = c('percent.mt', 'nFeature_RNA', 'nCount_RNA', 'S.Score', 'G2M.Score')
        # vars.to.regress = c('percent.mt', 'nFeature_RNA', 'nCount_RNA', 'CC.Difference')
        vars.to.regress = c('percent.mt', 'nFeature_RNA', 'nCount_RNA')
      )
      return(x)
    }
  )
  
  seurat_features <- SelectIntegrationFeatures(
    seurat_list,
    nfeatures = 3000
  )
  
  seurat_list <- PrepSCTIntegration(
    seurat_list,
    anchor.features = seurat_features
  )
  
  seurat_list <- lapply(
    X = seurat_list,
    FUN = RunPCA,
    verbose = FALSE,
    features = seurat_features
  )
  
  seurat_anchors <- FindIntegrationAnchors(
    object.list = seurat_list,
    anchor.features = seurat_features,
    normalization.method = 'SCT',
    reduction = 'rpca',
    reference = which(grepl(names(seurat_list), pattern = 'AML') == FALSE)
  )
  
  rm(seurat_list, seurat_features); gc(verbose=FALSE)
  
  seurat <- IntegrateData(
    anchorset = seurat_anchors,
    normalization.method = 'SCT',
    preserve.order = TRUE
  )
  
  rm(seurat_anchors); gc(verbose=FALSE)
  
  # make a concatenated column
  counts.1 <- seurat[["RNA"]]$counts.1
  counts.2 <- seurat[["RNA"]]$counts.2
  counts.3 <- seurat[["RNA"]]$counts.3
  counts.4 <- seurat[["RNA"]]$counts.4
  counts.5 <- seurat[["RNA"]]$counts.5
  counts.6 <- seurat[["RNA"]]$counts.6
  
  concatenated <- t(rbind(t(counts.1), t(counts.2), t(counts.3), t(counts.4), t(counts.5), t(counts.6)))
  
  seurat[["RNA"]]$counts <- concatenated
  seurat[["RNA"]]$counts.1 <- NULL
  seurat[["RNA"]]$counts.2 <- NULL
  seurat[["RNA"]]$counts.3 <- NULL
  seurat[["RNA"]]$counts.4 <- NULL
  seurat[["RNA"]]$counts.5 <- NULL
  seurat[["RNA"]]$counts.6 <- NULL
  
  rm(counts.1, counts.2, counts.3, counts.4, counts.5, counts.6, concatenated); gc(verbose = FALSE)
  
  # seurat@meta.data$cell_cycle_seurat <- seurat@meta.data$Phase
  # seurat@meta.data$Phase <- NULL
  # seurat@meta.data$cell_cycle_seurat <- factor(
  #   seurat@meta.data$cell_cycle_seurat, levels = c('G1', 'S', 'G2M')
  # )
}

# 
# {
#   seurat <- SCTransform(
#     seurat,
#     assay = 'RNA',
#     new.assay.name = 'SCT',
#     vars.to.regress = c('percent.mt', 'nFeature_RNA', 'nCount_RNA')
#   )
#   
#   new_s_genes <- read.csv(file.path(working_dir, "s_genes_mapping.csv"))$mapped
#   new_g2m_genes <- read.csv(file.path(working_dir, "g2m_genes_mapping.csv"))$mapped
#   
#   seurat <- CellCycleScoring(
#     seurat,
#     s.features = new_s_genes,
#     g2m.features = new_g2m_genes,
#     assay = 'SCT',
#     set.ident = TRUE
#   )
#   
#   seurat$CC.Difference <- seurat$S.Score - seurat$G2M.Score
#   
#   seurat <- SCTransform(
#     seurat,
#     assay = 'RNA',
#     new.assay.name = 'SCT',
#     vars.to.regress = c('percent.mt', 'nFeature_RNA', 'nCount_RNA', 'S.Score', 'G2M.Score')
#     # vars.to.regress = c('percent.mt', 'nFeature_RNA', 'nCount_RNA', 'CC.Difference')
#   )
# }

# saveRDS(seurat, file.path(working_dir, "filtered.filtered.preoptimisation.single.rds"))
seurat <- readRDS(file.path(working_dir, "filtered.filtered.preoptimisation.single.rds"))

# saveRDS(seurat, file.path(working_dir, "filtered.filtered.preoptimisation.5.singlet.rds"))
seurat <- readRDS(file.path(working_dir, "filtered.preoptimisation.5.singlet.rds"))

# saveRDS(seurat, file.path(working_dir, "filtered.filtered.preoptimisation.5.singlet.ln+.rds"))
seurat <- readRDS(file.path(working_dir, "filtered.filtered.preoptimisation.5.singlet.ln+.rds"))

# saveRDS(seurat, file.path(working_dir, "filtered.filtered.preoptimisation.5.singlet.ln-.rds"))
seurat <- readRDS(file.path(working_dir, "filtered.filtered.preoptimisation.5.singlet.ln-.rds"))

# saveRDS(seurat, file.path(working_dir, "filtered.filtered.preoptimisation.5.singlet.na+.rds"))
seurat <- readRDS(file.path(working_dir, "filtered.preoptimisation.5.singlet.na+.rds"))

# saveRDS(seurat, file.path(working_dir, "filtered.gfp.preoptimisation.single.s-g.rds"))
seurat <- readRDS(file.path(working_dir, "filtered.gfp.preoptimisation.single.s-g.rds"))

# saveRDS(seurat, file.path(working_dir, "filtered.preoptimisation.gfp.single.rds"))
seurat <- readRDS(file.path(working_dir, "filtered.preoptimisation.gfp.single.rds"))

# saveRDS(seurat, file.path(working_dir, "filtered.preoptimisation.gfp.single.nocc.rds"))
seurat <- readRDS(file.path(working_dir, "filtered.preoptimisation.gfp.single.nocc.rds"))

# # fast PCA, clustering, UMAP, tSNE
# seurat <- RunPCA(seurat, assay = 'integrated', npcs = 50)
# intrinsicDimension::maxLikGlobalDimEst(seurat@reductions$pca@cell.embeddings, k = 10)
# # 7.080667 5.singlet
# intrinsicDimension::maxLikGlobalDimEst(seurat@reductions$pca@cell.embeddings, k = 20)
# # 7.883863
# pca_cutoff <- 12
# seurat <- FindNeighbors(seurat, reduction = 'pca', dims = 1:pca_cutoff)
# seurat <- FindClusters(seurat, resolution = 0.3)
# distance_matrix <- dist(Embeddings(seurat[['pca']])[, 1:pca_cutoff])
# clusters <- seurat@meta.data$seurat_clusters
# silhouette <- silhouette(as.numeric(clusters), dist = distance_matrix)
# seurat@meta.data$silhouette_score <- silhouette[, 3]
# mean(seurat@meta.data$silhouette_score)
# seurat <- RunUMAP(seurat, dims = 1:pca_cutoff, reduction = "pca", n.components = 2, seed.use = 42)
# DimPlot(seurat) + DimPlot(seurat, group.by = "merged_sample")

# Scaling the data before PCA
# seurat <- ScaleData(seurat, assay = 'integrated')
# seurat <- ScaleData(seurat, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(seurat), assay = 'integrated')
# seurat$CC.Difference <- seurat$S.Score - seurat$G2M.Score
# seurat <- ScaleData(seurat, vars.to.regress = "CC.Difference", features = rownames(seurat), assay = 'integrated')

# prudent PCA, clustering, UMAP, tSNE
{
  seurat <- RunPCA(seurat, assay = 'integrated', npcs = 50)
  
  intrinsicDimension::maxLikGlobalDimEst(seurat@reductions$pca@cell.embeddings, k = 10)
  # 6.367669
  # 6.574318 5.singlet.s+g
  # 6.513917 5.singlet.s-g
  # 6.498403 5.singlet.s-g.int
  # 5.736749 gfp.single.s-g
  # 5.609283 gfp.single
  # 5.609283 gfp.single.nocc
  
  intrinsicDimension::maxLikGlobalDimEst(seurat@reductions$pca@cell.embeddings, k = 20)
  # 6.970843
  # 7.6299 5.singlet,s+g
  # 7.597453 5.singlet.s-g
  # 7.561695 5.singlet.s-g.int
  # 6.874308 gfp.single.s-g
  # 6.71184 gfp.single
  # 6.71184 gfp.single.nocc
  
  pca_results <- seurat@reductions$pca@stdev
  cumulative_variance <- cumsum(pca_results) / sum(pca_results)
  
  p <- tibble(
    PC = 1:50,
    stdev = pca_results,
    cumulative_variance = cumulative_variance
  ) %>%
    ggplot(aes(PC)) +
    geom_point(aes(y = seurat@reductions$pca@stdev)) +
    geom_line(aes(y = cumulative_variance * max(seurat@reductions$pca@stdev))) +  # Scale cumulative variance
    geom_vline(xintercept = 10, color = 'blue') +
    geom_vline(xintercept = 20, color = 'red') +
    scale_y_continuous(sec.axis = sec_axis(~ . / max(seurat@reductions$pca@stdev), name = "Cumulative Variance Explained")) +
    theme_bw() +
    labs(x = 'Principal components', y = 'Standard deviation') +
    ggtitle("PCA Elbow Plot with Cumulative Variance Explained")
  
  # Save the PCA elbow plot
  ggsave(
    # file.path(working_dir, "pca_elbow_plot_integrated.single.png"), 
    # file.path(working_dir, "pca_elbow_plot_integrated.5.singlet.s+g.png"), 
    # file.path(working_dir, "pca_elbow_plot_integrated.5.singlet.s-g.png"), 
    # file.path(working_dir, "pca_elbow_plot_integrated.5.singlet.s-g.int.png"), 
    # file.path(working_dir, "pca_elbow_plot_integrated.5.singlet.ln+.png"), 
    file.path(working_dir, paste0("pca_elbow.", runid, ".pdf")), 
    p, 
    height = 6, width = 8
  )
}

### Perform clustering
print(runid)
{
  pca_cutoff <- 10
  
  resolution <- 0.2
  
  clustering_methods <- list(
    Louvain = 1,
    Hierarchical = 3,
    DBSCAN = 4
  )
  method_name <- "Louvain"
  algorithm <- clustering_methods[[method_name]]
  
  seurat <- FindNeighbors(seurat, reduction = 'pca', dims = 1:pca_cutoff) 
  seurat <- FindClusters(seurat, resolution = resolution, algorithm = algorithm, random.seed = 42) # 0.8 by default
  
  ### Silhouette plot
  temp_labels <- seurat@meta.data %>%
    group_by(seurat_clusters) %>%
    tally()
  
  distance_matrix <- dist(Embeddings(seurat[['pca']])[, 1:pca_cutoff])
  clusters <- seurat@meta.data$seurat_clusters
  silhouette <- silhouette(as.numeric(clusters), dist = distance_matrix)
  seurat@meta.data$silhouette_score <- silhouette[, 3]
  mean_silhouette_score <- mean(seurat@meta.data$silhouette_score)
  
  p <- seurat@meta.data %>%
    mutate(barcode = rownames(.)) %>%
    arrange(seurat_clusters, -silhouette_score) %>%
    mutate(barcode = factor(barcode, levels = barcode)) %>%
    ggplot() +
    geom_col(aes(barcode, silhouette_score, fill = as.factor(seurat_clusters)), show.legend = TRUE) +
    geom_hline(yintercept = mean_silhouette_score, color = 'red', linetype = 'dashed') +
    scale_x_discrete(name = 'Cells') +
    scale_y_continuous(name = 'Silhouette score') +
    scale_fill_manual(values = custom_colors$discrete) +
    labs(fill = "Cluster") +
    theme_bw() +
    ggtitle(paste("Silhouette Plot (PCA Cutoff =", pca_cutoff, ", Resolution =", resolution, ", Method =", method_name, ", Mean Silhouette score =", round(mean_silhouette_score, 2), ")")) +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
  
  ggsave(
    # file.path(working_dir, 'silhouette_plot.single.png'), 
    # file.path(working_dir, 'silhouette_plot.5.singlet.s+g.png'), 
    # file.path(working_dir, 'silhouette_plot.5.singlet.s-g.int.png'), 
    # file.path(working_dir, 'silhouette_plot.5.singlet.ln+.png'), 
    # file.path(working_dir, 'silhouette_plot.single.gfp.s-g.pdf'), 
    # file.path(working_dir, 'silhouette_plot.single.gfp.png'), 
    file.path(working_dir, paste0('silhouette.', runid, '.png')), 
    p, height = 8, width = 12
  )
}

## integrated parameter optimisation
print(runid)
{
  # runid <- "single.gfp.nocc"
  
  clustering_methods <- list(
    Louvain = 1,
    Hierarchical = 3,
    DBSCAN = 4
  )
  
  for (pca_cutoff in seq(8, 14, 2)) {
    for (resolution in seq(0.2, 0.6, by = 0.2)) {
      for (method_name in names(clustering_methods)) {
        algorithm <- clustering_methods[[method_name]]
        
        seurat <- FindNeighbors(seurat, reduction = 'pca', dims = 1:pca_cutoff)
        seurat <- FindClusters(seurat, resolution = resolution, algorithm = algorithm)
        
        distance_matrix <- dist(Embeddings(seurat[['pca']])[, 1:pca_cutoff])
        clusters <- seurat@meta.data$seurat_clusters
        silhouette <- silhouette(as.numeric(clusters), dist = distance_matrix)
        seurat@meta.data$silhouette_score <- silhouette[, 3]
        mean_silhouette_score <- mean(seurat@meta.data$silhouette_score)
        
        p <- seurat@meta.data %>%
          mutate(barcode = rownames(.)) %>%
          arrange(seurat_clusters, -silhouette_score) %>%
          mutate(barcode = factor(barcode, levels = barcode)) %>%
          ggplot() +
          geom_col(aes(barcode, silhouette_score, fill = as.factor(seurat_clusters)), show.legend = TRUE) +
          geom_hline(yintercept = mean_silhouette_score, color = 'red', linetype = 'dashed') +
          scale_x_discrete(name = 'Cells') +
          scale_y_continuous(name = 'Silhouette score') +
          scale_fill_manual(values = custom_colors$discrete) +
          labs(fill = "Cluster") +
          theme_bw() +
          ggtitle(paste("Silhouette Plot (PCA Cutoff =", pca_cutoff,
                        ", Resolution =", resolution,
                        ", Method =", method_name,
                        ", Mean Silhouette score =", round(mean_silhouette_score, 2), ")")) +
          theme(
            axis.title.x = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()
          )
        
        file_name <- file.path(working_dir, paste0('silhouette_', runid, '/silhouette_plot_',
                                                   method_name,
                                                   '_pca_cutoff_', pca_cutoff,
                                                   '_resolution_', resolution,
                                                   '_mean_silhouette_', round(mean_silhouette_score, 2),
                                                   '.png'))
        ggsave(file_name, p, height = 6, width = 10)
      }
    }
  }
}

### Expression profile per cluster
temp_labels <- seurat@meta.data %>%
  group_by(seurat_clusters) %>%
  tally()

p1 <- ggplot() +
  geom_half_violin(
    data = seurat@meta.data, aes(seurat_clusters, nCount_RNA, fill = seurat_clusters),
    side = 'l', show.legend = FALSE, trim = TRUE
  ) +
  geom_half_boxplot(
    data = seurat@meta.data, aes(seurat_clusters, nCount_RNA, fill = seurat_clusters),
    side = 'r', outlier.color = NA, width = 0.4, show.legend = FALSE
  ) +
  geom_text(
    data = temp_labels,
    aes(x = seurat_clusters, y = -Inf, label = paste0('n = ', format(n, big.mark = ',', trim = TRUE)), vjust = -1),
    color = 'black', size = 2.8
  ) +
  scale_color_manual(values = custom_colors$discrete) +
  scale_fill_manual(values = custom_colors$discrete) +
  scale_y_continuous(name = 'Number of transcripts', labels = scales::comma, expand = c(0.08, 0)) +
  theme_bw() +
  theme(
    panel.grid.major.x = element_blank(),
    axis.title.x = element_blank()
  )

p2 <- ggplot() +
  geom_half_violin(
    data = seurat@meta.data, aes(seurat_clusters, nFeature_RNA, fill = seurat_clusters),
    side = 'l', show.legend = FALSE, trim = TRUE
  ) +
  geom_half_boxplot(
    data = seurat@meta.data, aes(seurat_clusters, nFeature_RNA, fill = seurat_clusters),
    side = 'r', outlier.color = NA, width = 0.4, show.legend = FALSE
  ) +
  geom_text(
    data = temp_labels,
    aes(x = seurat_clusters, y = -Inf, label = paste0('n = ', format(n, big.mark = ',', trim = TRUE)), vjust = -1),
    color = 'black', size = 2.8
  ) +
  scale_color_manual(values = custom_colors$discrete) +
  scale_fill_manual(values = custom_colors$discrete) +
  scale_y_continuous(name = 'Number of expressed genes', labels = scales::comma, expand = c(0.08, 0)) +
  theme_bw() +
  theme(
    panel.grid.major.x = element_blank(),
    axis.title.x = element_blank()
  )

ggsave(
  # file.path(working_dir, 'ncount_nfeature_by_cluster.single.png'),
  file.path(working_dir, 'ncount_nfeature_by_cluster.5.singlet.s-g.png'),
  p1 + p2 + plot_layout(ncol = 1),
  height = 7, width = 14
)

### Cluster stability failed due to R thus scran version. Need R/4.0 and scran/3.16
sce <- as.SingleCellExperiment(seurat)
reducedDim(sce, 'PCA_sub') <- reducedDim(sce, 'PCA')[,1:pca_cutoff, drop = FALSE]

ass_prob <- bluster::bootstrapStability(sce, FUN = function(x) {
  g <- buildSNNGraph(x, use.dimred = 'PCA_sub')
  igraph::cluster_walktrap(g)$membership
},
clusters = sce$seurat_clusters
)

p <- ass_prob %>%
  as_tibble() %>%
  rownames_to_column(var = 'cluster_1') %>%
  pivot_longer(
    cols = 2:ncol(.),
    names_to = 'cluster_2',
    values_to = 'probability'
  ) %>%
  mutate(
    cluster_1 = as.character(as.numeric(cluster_1) - 1),
    cluster_1 = factor(cluster_1, levels = rev(unique(cluster_1))),
    cluster_2 = factor(cluster_2, levels = unique(cluster_2))
  ) %>%
  ggplot(aes(cluster_2, cluster_1, fill = probability)) +
  geom_tile(color = 'white') +
  geom_text(aes(label = round(probability, digits = 2)), size = 2.5) +
  scale_x_discrete(name = 'Cluster', position = 'top') +
  scale_y_discrete(name = 'Cluster') +
  scale_fill_gradient(
    name = 'Probability', low = 'white', high = '#c0392b', na.value = '#bdc3c7',
    limits = c(0, max(as.matrix(ass_prob), na.rm = TRUE)),
    guide = guide_colorbar(
      frame.colour = 'black', ticks.colour = 'black', title.position = 'left',
      title.theme = element_text(hjust = 1, angle = 90),
      barwidth = 0.75, barheight = 10
    )
  ) +
  coord_fixed() +
  theme_bw() +
  theme(
    legend.position = 'right',
    panel.grid.major = element_blank()
  )
# write.csv(as.matrix(ass_prob), file = "ass_prob.csv")

ggsave(file.path(working_dir, 'cluster_stability.png'), p, height = 6, width = 7)

### cluster similarity

sce <- as.SingleCellExperiment(seurat)

reducedDim(sce, 'PCA_sub') <- reducedDim(sce, 'PCA')[,1:pca_cutoff, drop = FALSE]

g <- scran::buildSNNGraph(sce, use.dimred = 'PCA_sub')

ratio <- bluster::pairwiseModularity(g, seurat@meta.data$seurat_clusters, as.ratio = TRUE)

ratio_to_plot <- log10(ratio+1)

p <- ratio_to_plot %>%
  as_tibble() %>%
  rownames_to_column(var = 'cluster_1') %>%
  pivot_longer(
    cols = 2:ncol(.),
    names_to = 'cluster_2',
    values_to = 'probability'
  ) %>%
  mutate(
    cluster_1 = as.character(as.numeric(cluster_1) - 1),
    cluster_1 = factor(cluster_1, levels = rev(unique(cluster_1))),
    cluster_2 = factor(cluster_2, levels = unique(cluster_2))
  ) %>%
  ggplot(aes(cluster_2, cluster_1, fill = probability)) +
  geom_tile(color = 'white') +
  geom_text(aes(label = round(probability, digits = 2)), size = 2.5) +
  scale_x_discrete(name = 'Cluster', position = 'top') +
  scale_y_discrete(name = 'Cluster') +
  scale_fill_gradient(
    name = 'log10(ratio)', low = 'white', high = '#c0392b', na.value = '#bdc3c7',
    guide = guide_colorbar(
      frame.colour = 'black', ticks.colour = 'black', title.position = 'left',
      title.theme = element_text(hjust = 1, angle = 90),
      barwidth = 0.75, barheight = 10
    )
  ) +
  coord_fixed() +
  theme_bw() +
  theme(
    legend.position = 'right',
    panel.grid.major = element_blank()
  )

ggsave(file.path(working_dir, 'cluster_similarity.png'), p, height = 6, width = 7)

### Cluster tree
{
  seurat <- BuildClusterTree(
    seurat,
    dims = 1:pca_cutoff,
    reorder = FALSE,
    reorder.numeric = FALSE
  )
  
  tree <- seurat@tools$BuildClusterTree
  tree$tip.label <- paste0("Cluster ", tree$tip.label)
  
  p <- ggtree::ggtree(tree, aes(x, y)) +
    scale_y_reverse() +
    ggtree::geom_tree() +
    ggtree::theme_tree() +
    ggtree::geom_tiplab(offset = 1) +
    ggtree::geom_tippoint(color = custom_colors$discrete[1:length(tree$tip.label)], shape = 16, size = 5) +
    coord_cartesian(clip = 'off') +
    theme(plot.margin = unit(c(0,2.5,0,0), 'cm'))
  
  ggsave(
    # file.path(working_dir, 'cluster_tree.single.png'), 
    # file.path(working_dir, 'cluster_tree.single.gfp.s-g.pdf'), 
    file.path(working_dir, paste0('cluster_tree.', runid, '.png')), 
    p, 
    height = 8, 
    width = 6
  )
}

### composition of samples and clusters
seurat <- filtering_seurat(seurat)
seurat <- only3_seurat(seurat)
seurat <- purge4_seurat(seurat)
levels(seurat@meta.data$merged_sample)
seurat <- subset(seurat, subset = sample != "Undetermined")

cluster_composition <- function(
    seurat, 
    runid, 
    resolution = "merged_sample", # sample, merged_sample, combined_sample, meta_sample
    plot_width = 16 # 18 for 23 columns, 16 for 16 columns, 14 for 15 columns
) {
  
  seurat@meta.data$resolution <- seurat@meta.data[[resolution]]
  
  table_samples_by_clusters <- seurat@meta.data %>%
    group_by(resolution, seurat_clusters) %>%
    dplyr::summarise(count = n()) %>%
    spread(seurat_clusters, count, fill = 0) %>%
    ungroup() %>%
    mutate(total_cell_count = rowSums(.[c(2:ncol(.))])) %>%
    dplyr::select(c('resolution', 'total_cell_count', everything())) %>%
    arrange(factor(resolution, levels = levels(seurat@meta.data$resolution)))
  
  knitr::kable(table_samples_by_clusters)
  
  table_clusters_by_samples <- seurat@meta.data %>%
    dplyr::rename('cluster' = 'seurat_clusters') %>%
    group_by(cluster, resolution) %>%
    dplyr::summarize(count = n()) %>%
    spread(resolution, count, fill = 0) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(total_cell_count = rowSums(.[c(2:ncol(.))])) %>%
    dplyr::select(c('cluster', 'total_cell_count', everything())) %>%
    arrange(factor(cluster, levels = levels(seurat@meta.data$seurat_clusters)))
  
  knitr::kable(table_clusters_by_samples)
  
  write.csv(rbind(cbind(table_type = "Samples by Clusters - Absolute Counts", table_samples_by_clusters), setNames(data.frame(table_type = "", resolution = "", total_cell_count = "", t(rep("", ncol(table_samples_by_clusters) - 2))), names(cbind(table_type = "Samples by Clusters - Absolute Counts", table_samples_by_clusters))), cbind(table_type = "Samples by Clusters - Percentages", resolution = table_samples_by_clusters$resolution, total_cell_count = table_samples_by_clusters$total_cell_count, round(table_samples_by_clusters[,3:ncol(table_samples_by_clusters)] / table_samples_by_clusters$total_cell_count * 100, 2))), file.path(working_dir, "temp_samples_by_clusters.csv"), row.names = FALSE, quote = FALSE)
  
  write.csv(rbind(cbind(table_type = "Clusters by Samples - Absolute Counts", table_clusters_by_samples), setNames(data.frame(table_type = "", cluster = "", total_cell_count = "", t(rep("", ncol(table_clusters_by_samples) - 2))), names(cbind(table_type = "Clusters by Samples - Absolute Counts", table_clusters_by_samples))), cbind(table_type = "Clusters by Samples - Percentages", cluster = table_clusters_by_samples$cluster, total_cell_count = table_clusters_by_samples$total_cell_count, round(table_clusters_by_samples[,3:ncol(table_clusters_by_samples)] / table_clusters_by_samples$total_cell_count * 100, 2))), file.path(working_dir, "temp_clusters_by_samples.csv"), row.names = FALSE, quote = FALSE)
  
  system(paste0("cat '", file.path(working_dir, "temp_samples_by_clusters.csv"), "' > '", file.path(working_dir, paste0("cluster_sample_summary.", runid, ".csv")), "' && echo '' >> '", file.path(working_dir, paste0("cluster_sample_summary.", runid, ".csv")), "' && tail -n +2 '", file.path(working_dir, "temp_clusters_by_samples.csv"), "' >> '", file.path(working_dir, "cluster_sample_summary.single.csv"), "' && rm '", file.path(working_dir, "temp_samples_by_clusters.csv"), "' '", file.path(working_dir, "temp_clusters_by_samples.csv"), "'"))
  
  # absolute count
  temp_labels <- seurat@meta.data %>%
    group_by(resolution) %>%
    tally()
  
  p1 <- table_samples_by_clusters %>%
    dplyr::select(-c('total_cell_count')) %>%
    reshape2::melt(id.vars = 'resolution') %>%
    ggplot(aes(resolution, value)) +
    geom_bar(aes(fill = variable), position = 'stack', stat = 'identity') +
    geom_text(
      data = temp_labels,
      aes(x = resolution, y = Inf, label = paste0('n = ', format(n, big.mark = ',', trim = TRUE)), vjust = -1),
      color = 'black', size = 2.8
    ) +
    scale_fill_manual(name = 'Clusters', values = custom_colors$discrete) +
    scale_y_continuous(name = 'Number of cells', labels = scales::comma, expand = c(0.01, 0)) +
    coord_cartesian(clip = 'off') +
    theme_bw() +
    theme(
      legend.position = 'left',
      plot.title = element_text(hjust = 0.5),
      text = element_text(size = 16),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.title.x = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      plot.margin = margin(t = 20, r = 0, b = 0, l = 0, unit = 'pt')
    )
  
  temp_labels <- seurat@meta.data %>%
    group_by(seurat_clusters) %>%
    tally() %>%
    dplyr::rename('cluster' = seurat_clusters)
  
  p2 <- table_clusters_by_samples %>%
    dplyr::select(-c('total_cell_count')) %>%
    reshape2::melt(id.vars = 'cluster') %>%
    dplyr::mutate(cluster = factor(cluster, levels = levels(seurat@meta.data$seurat_clusters))) %>%
    ggplot(aes(cluster, value)) +
    geom_bar(aes(fill = variable), position = 'stack', stat = 'identity') +
    geom_text(
      data = temp_labels,
      aes(x = cluster, y = Inf, label = paste0('n = ', format(n, big.mark = ',', trim = TRUE)), vjust = -1),
      color = 'black', size = 2.8
    ) +
    scale_fill_manual(name = 'Samples', values = custom_colors$discrete) +
    scale_y_continuous(labels = scales::comma, expand = c(0.01, 0)) +
    coord_cartesian(clip = 'off') +
    theme_bw() +
    theme(
      legend.position = 'right',
      plot.title = element_text(hjust = 0.5),
      text = element_text(size = 16),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.title = element_blank(),
      plot.margin = margin(t = 20, r = 0, b = 0, l = 10, unit = 'pt')
    )
  
  ggsave(
    # file.path(working_dir, 'composition_samples_clusters_by_number.single.png'),
    # file.path(working_dir, 'final/composition.absolute.single.purge4.png'),
    # file.path(working_dir, 'final/composition.absolute.gfp.single.s-g.png'),
    file.path(working_dir, paste0('final/composition_cluster.absolute.', runid, '.pdf')),
    # file.path(working_dir, 'composition_samples_clusters_by_number.5.singlet.s-g.png'),
    p1 + p2 +
      plot_layout(ncol = 2, widths = c(
        seurat@meta.data$resolution %>% unique() %>% length() + 1,
        seurat@meta.data$seurat_clusters %>% unique() %>% length()
      )),
    width = plot_width, 
    height = 8
  )
  
  # composition in percentage
  temp_labels <- seurat@meta.data %>%
    group_by(resolution) %>%
    tally()
  
  p1 <- table_samples_by_clusters %>%
    dplyr::select(-c('total_cell_count')) %>%
    reshape2::melt(id.vars = 'resolution') %>%
    ggplot(aes(resolution, value)) +
    geom_bar(aes(fill = variable), position = 'fill', stat = 'identity') +
    geom_text(
      data = temp_labels,
      aes(x = resolution, y = Inf, label = paste0('n = ', format(n, big.mark = ',', trim = TRUE)), vjust = -1),
      color = 'black', size = 2.8
    ) +
    scale_fill_manual(name = 'Cluster', values = custom_colors$discrete) +
    scale_y_continuous(name = 'Percentage of cells (%)', labels = scales::percent_format(), expand = c(0.01, 0)) +
    coord_cartesian(clip = 'off') +
    theme_bw() +
    theme(
      legend.position = 'left',
      plot.title = element_text(hjust = 0.5),
      text = element_text(size = 16),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.title.x = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      plot.margin = margin(t = 20, r = 0, b = 0, l = 0, unit = 'pt')
    )
  
  temp_labels <- seurat@meta.data %>%
    group_by(seurat_clusters) %>%
    tally() %>%
    dplyr::rename('cluster' = seurat_clusters)
  
  p2 <- table_clusters_by_samples %>%
    dplyr::select(-c('total_cell_count')) %>%
    reshape2::melt(id.vars = 'cluster') %>%
    dplyr::mutate(cluster = factor(cluster, levels = levels(seurat@meta.data$seurat_clusters))) %>%
    ggplot(aes(cluster, value)) +
    geom_bar(aes(fill = variable), position = 'fill', stat = 'identity') +
    geom_text(
      data = temp_labels,
      aes(x = cluster, y = Inf, label = paste0('n = ', format(n, big.mark = ',', trim = TRUE)), vjust = -1),
      color = 'black', size = 2.8
    ) +
    scale_fill_manual(name = 'Sample', values = custom_colors$discrete) +
    scale_y_continuous(labels = scales::percent_format(), expand = c(0.01, 0)) +
    coord_cartesian(clip = 'off') +
    theme_bw() +
    theme(
      legend.position = 'right',
      plot.title = element_text(hjust = 0.5),
      text = element_text(size = 16),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.title = element_blank(),
      plot.margin = margin(t = 20, r = 0, b = 0, l = 10, unit = 'pt')
    )
  
  ggsave(
    # file.path(working_dir, 'final/composition.percentage.single.purge4.svg'),
    file.path(working_dir, paste0('final/composition_cluster.percentage.', runid, '.pdf')),
    # file.path(working_dir, 'composition_samples_clusters_by_number.5.singlet.s-g.png'),  
    p1 + p2 +
      plot_layout(ncol = 2, widths = c(
        seurat@meta.data$resolution %>% unique() %>% length() + 1,
        seurat@meta.data$seurat_clusters %>% unique() %>% length()
      )),
    width = plot_width, 
    height = 8
  )
  
  seurat@meta.data$resolution <- NULL
}

print(runid)
cluster_composition(seurat, runid, resolution = "merged_sample", plot_width = 14)

### cell cycle
# prepare gene lists
# {
#   # Get human and mouse gene databases
#   human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
#   mouse <- useMart("ensembl", dataset = "musculus_gene_ensembl")
#   
#   # Get mouse orthologs for cell cycle genes
#   s.genes.mouse <- getLDS(attributes = c("hgnc_symbol"),
#                           filters = "hgnc_symbol",
#                           values = cc.genes$s.genes,
#                           mart = human,
#                           attributesL = c("mgi_symbol"),
#                           martL = mouse,
#                           uniqueRows = TRUE)
#   
#   g2m.genes.mouse <- getLDS(attributes = c("hgnc_symbol"),
#                             filters = "hgnc_symbol",
#                             values = cc.genes$g2m.genes,
#                             mart = human,
#                             attributesL = c("mgi_symbol"),
#                             martL = mouse,
#                             uniqueRows = TRUE)
#   # Function to find potential matches for a gene in your Seurat object
#   find_gene_variations <- function(gene, all_genes, max_results = 5) {
#     # Direct match
#     if(gene %in% all_genes) {
#       return(data.frame(original = gene, match = gene, match_type = "exact", stringsAsFactors = FALSE))
#     }
#     
#     # Search specifically for genes that start with our gene symbol followed by a dot
#     dot_suffix_matches <- grep(paste0("^", gene, "\\."), all_genes, value = TRUE)
#     
#     # Create results data frame
#     results <- data.frame(
#       original = character(),
#       match = character(),
#       match_type = character(),
#       stringsAsFactors = FALSE
#     )
#     
#     # Add dot suffix matches
#     if(length(dot_suffix_matches) > 0) {
#       dot_suffix_df <- data.frame(
#         original = rep(gene, length(dot_suffix_matches)),
#         match = dot_suffix_matches,
#         match_type = "dot_suffix_match",
#         stringsAsFactors = FALSE
#       )
#       results <- rbind(results, dot_suffix_df)
#     }
#     
#     # If no matches found, return NA
#     if(nrow(results) == 0) {
#       return(data.frame(original = gene, match = NA, match_type = "no_match", stringsAsFactors = FALSE))
#     }
#     
#     # Limit results if there are too many
#     if(nrow(results) > max_results) {
#       results <- results[1:max_results,]
#     }
#     
#     return(results)
#   }
#   
#   # Get all genes in your Seurat object
#   all_genes <- rownames(seurat)
#   
#   # Search for matches for S phase genes
#   s_genes_matches <- lapply(cc.genes$s.genes, find_gene_variations, all_genes = all_genes)
#   s_genes_df <- do.call(rbind, s_genes_matches)
#   
#   # Search for matches for G2M phase genes
#   g2m_genes_matches <- lapply(cc.genes$g2m.genes, find_gene_variations, all_genes = all_genes)
#   g2m_genes_df <- do.call(rbind, g2m_genes_matches)
#   
#   # Function to print results in a readable format
#   print_matches <- function(matches_df, title) {
#     cat("\n", title, "\n", sep = "")
#     cat(paste0(rep("-", nchar(title)), collapse = ""), "\n\n")
#     
#     # Group by original gene and then by match type
#     genes_by_match <- split(matches_df, matches_df$original)
#     
#     for(gene in names(genes_by_match)) {
#       gene_df <- genes_by_match[[gene]]
#       cat("Gene:", gene, "\n")
#       
#       if(any(!is.na(gene_df$match))) {
#         exact_matches <- gene_df[gene_df$match_type == "exact", "match"]
#         prefix_matches <- gene_df[gene_df$match_type == "prefix_match", "match"]
#         contains_matches <- gene_df[gene_df$match_type == "contains_match", "match"]
#         
#         if(length(exact_matches) > 0) {
#           cat("  Exact match:", paste(exact_matches, collapse = ", "), "\n")
#         }
#         if(length(prefix_matches) > 0) {
#           cat("  Prefix matches:", paste(prefix_matches, collapse = ", "), "\n")
#         }
#         if(length(contains_matches) > 0) {
#           cat("  Contains matches:", paste(contains_matches, collapse = ", "), "\n")
#         }
#       } else {
#         cat("  No matches found\n")
#       }
#       cat("\n")
#     }
#   }
#   
#   # Print results
#   print_matches(s_genes_df, "S PHASE GENES MATCHES")
#   print_matches(g2m_genes_df, "G2M PHASE GENES MATCHES")
#   
#   # Create a mapping dictionary for use in CellCycleScoring
#   create_gene_mapping <- function(matches_df) {
#     # Start with exact matches
#     exact_matches <- matches_df[matches_df$match_type == "exact",]
#     mapping <- setNames(exact_matches$match, exact_matches$original)
#     
#     # For genes without exact matches, look for dot suffix matches
#     missing_genes <- setdiff(unique(matches_df$original), names(mapping))
#     
#     for(gene in missing_genes) {
#       gene_df <- matches_df[matches_df$original == gene,]
#       
#       # Look for dot suffix matches
#       dot_suffix_matches <- gene_df[gene_df$match_type == "dot_suffix_match", "match"]
#       if(length(dot_suffix_matches) > 0) {
#         # If multiple matches, use the simplest one (usually the one with shortest suffix)
#         suffix_lengths <- nchar(dot_suffix_matches) - nchar(gene) - 1  # -1 for the dot
#         simplest_match <- dot_suffix_matches[which.min(suffix_lengths)]
#         mapping[gene] <- simplest_match
#       }
#     }
#     
#     return(mapping)
#   }
#   
#   # Create mapping dictionaries
#   s_genes_mapping <- create_gene_mapping(s_genes_df)
#   g2m_genes_mapping <- create_gene_mapping(g2m_genes_df)
#   
#   # Replace NA values with original gene names (these will be filtered out later)
#   s_genes_mapping[is.na(s_genes_mapping)] <- names(s_genes_mapping)[is.na(s_genes_mapping)]
#   g2m_genes_mapping[is.na(g2m_genes_mapping)] <- names(g2m_genes_mapping)[is.na(g2m_genes_mapping)]
#   
#   # Create new gene lists with mappings
#   new_s_genes <- s_genes_mapping[s_genes_mapping %in% all_genes]
#   new_g2m_genes <- g2m_genes_mapping[g2m_genes_mapping %in% all_genes]
#   
#   cat("\nFinal gene counts for cell cycle scoring:\n")
#   cat("S phase genes:", length(new_s_genes), "out of", length(cc.genes$s.genes), "\n")
#   cat("G2M phase genes:", length(new_g2m_genes), "out of", length(cc.genes$g2m.genes), "\n\n")
#   
#   # You can also save the mapping for future reference
#   cell_cycle_mappings <- list(
#     s_genes = data.frame(
#       original = names(s_genes_mapping),
#       mapped = as.character(s_genes_mapping),
#       in_dataset = s_genes_mapping %in% all_genes,
#       stringsAsFactors = FALSE
#     ),
#     g2m_genes = data.frame(
#       original = names(g2m_genes_mapping),
#       mapped = as.character(g2m_genes_mapping),
#       in_dataset = g2m_genes_mapping %in% all_genes,
#       stringsAsFactors = FALSE
#     )
#   )
#   
#   # Save to CSV for reference
#   write.csv(cell_cycle_mappings$s_genes, file.path(working_dir, "s_genes_mapping.csv"), row.names = FALSE)
#   write.csv(cell_cycle_mappings$g2m_genes, file.path(working_dir, "g2m_genes_mapping.csv"), row.names = FALSE)
# }

# cellcyclescore
new_s_genes <- read.csv(file.path(working_dir, "s_genes_mapping.csv"))$mapped
new_g2m_genes <- read.csv(file.path(working_dir, "g2m_genes_mapping.csv"))$mapped

seurat <- CellCycleScoring(
  seurat,
  # assay = 'integrated',
  assay = 'SCT',
  # assay = 'RNA',
  s.features = new_s_genes,
  g2m.features = new_g2m_genes
  #   s.features = cc.genes$s.genes,
  #   g2m.features = cc.genes$g2m.genes
)

seurat@meta.data$cell_cycle_seurat <- seurat@meta.data$Phase
seurat@meta.data$Phase <- NULL
seurat@meta.data$cell_cycle_seurat <- factor(
  seurat@meta.data$cell_cycle_seurat, levels = c('G1', 'S', 'G2M')
)

### add g0 cellcycle modules
rownames(seurat[['integrated']])[grep("(?i)KLF4", rownames(seurat[['integrated']]))]

new_s_genes <- read.csv(file.path(working_dir, "s_genes_mapping.csv"))$mapped
new_g2m_genes <- read.csv(file.path(working_dir, "g2m_genes_mapping.csv"))$mapped
# g0_genes <- c("RB1CC1.2", "TP53INP1", "TP53I11",  "TP53I3.2", "TP53I13", "Foxo (FOXO4)", "SLC2A1.9", "SLC2A1", "MTOR", "Myc (MYCN)")
# 
# # Step 1: Run standard Seurat cell cycle analysis for S and G2M
# seurat <- CellCycleScoring(seurat, 
#                            assay = 'SCT',
#                            s.features = new_s_genes, 
#                            g2m.features = new_g2m_genes, 
#                            set.ident = FALSE)
# 
# # Step 2: Add G0 quiescence score
# # Separate G0 genes by expected expression direction
# g0_high <- c("RB1CC1.2", "TP53INP1", "TP53I11", "TP53I3.2", "TP53I13", "Foxo (FOXO4)")
# g0_low <- c("MTOR", "Myc (MYCN)", "SLC2A1.9", "SLC2A1")
# 
# # Check which genes are available
# available_g0_high <- intersect(g0_high, rownames(seurat))
# available_g0_low <- intersect(g0_low, rownames(seurat))
# 
# # Calculate G0 scores
# seurat <- AddModuleScore(seurat, features = list(available_g0_high), name = "G0_High")
# seurat <- AddModuleScore(seurat, features = list(available_g0_low), name = "G0_Low")
# 
# # Combined G0 score
# seurat$G0_Score <- seurat$G0_High1 - seurat$G0_Low1
# 
# # Step 3: Classify cells into three phases
# classify_cell_cycle_phases <- function(seurat_obj) {
#   # Get scores
#   s_score <- seurat_obj$S.Score
#   g2m_score <- seurat_obj$G2M.Score
#   g0_score <- seurat_obj$G0_Score
#   
#   # Normalize scores to 0-1 range for comparison
#   s_norm <- (s_score - min(s_score)) / (max(s_score) - min(s_score))
#   g2m_norm <- (g2m_score - min(g2m_score)) / (max(g2m_score) - min(g2m_score))
#   g0_norm <- (g0_score - min(g0_score)) / (max(g0_score) - min(g0_score))
#   
#   # Classification logic
#   phases <- rep("G1", length(s_score))  # Default to G1
#   
#   # High G0 score = G0/Quiescent
#   phases[g0_norm > quantile(g0_norm, 0.7)] <- "G0"
#   
#   # High S score (and not G0) = S phase
#   phases[s_norm > quantile(s_norm, 0.6) & phases != "G0"] <- "S"
#   
#   # High G2M score (and not G0) = G2M phase
#   phases[g2m_norm > quantile(g2m_norm, 0.6) & phases != "G0"] <- "G2M"
#   
#   return(phases)
# }
# 
# # Apply classification
# seurat$custom_cell_cycle_seurat <- classify_cell_cycle_phases(seurat)
# 
# seurat@meta.data$cell_cycle_seurat <- seurat@meta.data$Phase
# seurat@meta.data$Phase <- NULL
# seurat@meta.data$cell_cycle_seurat <- factor(
#   seurat@meta.data$cell_cycle_seurat, levels = c('G1', 'S', 'G2M')
# )
# seurat@meta.data$custom_cell_cycle_seurat <- factor(
#   seurat@meta.data$custom_cell_cycle_seurat, levels = c('G0', 'G1', 'S', 'G2M')
# )

cell_cycle_composition <- function(
    seurat, 
    runid,
    cycle = 3, # 3 for G1/S/G2M, 4 for G0/G1/S/G2M
    resolution = "merged_sample", # sample, merged_sample, combined_sample, meta_sample
    plot_width = 16 # 18 for 23 columns, 16 for 16 columns, 14 for 15 columns
) {
  
  # Determine cell cycle column and stages based on cycle parameter
  if (cycle == 4) {
    cell_cycle_col <- "custom_cell_cycle_seurat"
    cycle_stages <- c("G0", "G1", "S", "G2M")
    cycle_suffix <- ".4cycle"
  } else {
    cell_cycle_col = "cell_cycle_seurat"
    cycle_stages <- c("G1", "S", "G2M")
    cycle_suffix <- ".3cycle"
  }
  
  # Create temporary columns for dynamic referencing
  seurat@meta.data$resolution <- seurat@meta.data[[resolution]]
  seurat@meta.data$cell_cycle <- seurat@meta.data[[cell_cycle_col]]
  
  table_samples_by_cell_cycle <- seurat@meta.data %>%
    group_by(resolution, cell_cycle) %>%
    dplyr::summarize(count = n(), .groups = 'drop') %>%
    spread(cell_cycle, count, fill = 0) %>%
    ungroup() %>%
    mutate(total_cell_count = rowSums(.[c(2:ncol(.))])) %>%
    dplyr::select(c('resolution', 'total_cell_count', everything())) %>%
    arrange(factor(resolution, levels = levels(seurat@meta.data$resolution)))
  
  knitr::kable(table_samples_by_cell_cycle)
  
  table_clusters_by_cell_cycle <- seurat@meta.data %>%
    dplyr::rename(cluster = seurat_clusters) %>%
    group_by(cluster, cell_cycle) %>%
    dplyr::summarize(count = n(), .groups = 'drop') %>%
    spread(cell_cycle, count, fill = 0) %>%
    ungroup() %>%
    mutate(total_cell_count = rowSums(.[c(2:ncol(.))])) %>%
    dplyr::select(c('cluster', 'total_cell_count', everything())) %>%
    arrange(factor(cluster, levels = levels(seurat@meta.data$seurat_clusters)))
  
  knitr::kable(table_clusters_by_cell_cycle)
  
  # Write CSV files
  write.csv(rbind(
    cbind(table_type = "Samples by Cell Cycle - Absolute Counts", table_samples_by_cell_cycle), 
    setNames(data.frame(table_type = "", resolution = "", total_cell_count = "", t(rep("", ncol(table_samples_by_cell_cycle) - 2))), names(cbind(table_type = "Samples by Cell Cycle - Absolute Counts", table_samples_by_cell_cycle))), 
    cbind(table_type = "Samples by Cell Cycle - Percentages", resolution = table_samples_by_cell_cycle$resolution, total_cell_count = table_samples_by_cell_cycle$total_cell_count, round(table_samples_by_cell_cycle[,3:ncol(table_samples_by_cell_cycle)] / table_samples_by_cell_cycle$total_cell_count * 100, 2))
  ), file.path(working_dir, "temp_samples_by_cell_cycle.csv"), row.names = FALSE, quote = FALSE)
  
  write.csv(rbind(
    cbind(table_type = "Clusters by Cell Cycle - Absolute Counts", table_clusters_by_cell_cycle), 
    setNames(data.frame(table_type = "", cluster = "", total_cell_count = "", t(rep("", ncol(table_clusters_by_cell_cycle) - 2))), names(cbind(table_type = "Clusters by Cell Cycle - Absolute Counts", table_clusters_by_cell_cycle))), 
    cbind(table_type = "Clusters by Cell Cycle - Percentages", cluster = table_clusters_by_cell_cycle$cluster, total_cell_count = table_clusters_by_cell_cycle$total_cell_count, round(table_clusters_by_cell_cycle[,3:ncol(table_clusters_by_cell_cycle)] / table_clusters_by_cell_cycle$total_cell_count * 100, 2))
  ), file.path(working_dir, "temp_clusters_by_cell_cycle.csv"), row.names = FALSE, quote = FALSE)
  
  system(paste0("cat '", file.path(working_dir, "temp_samples_by_cell_cycle.csv"), "' > '", file.path(working_dir, paste0("cluster_sample_cellcycle_summary.", runid, ".csv")), "' && echo '' >> '", file.path(working_dir, paste0("cluster_sample_cellcycle_summary.", runid, ".csv")), "' && tail -n +2 '", file.path(working_dir, "temp_clusters_by_cell_cycle.csv"), "' >> '", file.path(working_dir, paste0("cluster_sample_cellcycle_summary.", runid, ".csv")), "' && rm '", file.path(working_dir, "temp_samples_by_cell_cycle.csv"), "' '", file.path(working_dir, "temp_clusters_by_cell_cycle.csv"), "'"))
  
  # Prepare data for visualization - dynamically handle cycle stages
  if (cycle == 4) {
    # Ensure all 4 cycle stages exist in the data (add missing columns with 0s)
    for (stage in cycle_stages) {
      if (!stage %in% names(table_samples_by_cell_cycle)) {
        table_samples_by_cell_cycle[[stage]] <- 0
      }
      if (!stage %in% names(table_clusters_by_cell_cycle)) {
        table_clusters_by_cell_cycle[[stage]] <- 0
      }
    }
    
    table_samples_for_viz <- table_samples_by_cell_cycle %>%
      mutate(non_na_total = rowSums(dplyr::select(., all_of(cycle_stages)), na.rm = TRUE))
    
    table_clusters_for_viz <- table_clusters_by_cell_cycle %>%
      mutate(non_na_total = rowSums(dplyr::select(., all_of(cycle_stages)), na.rm = TRUE))
    
  } else {
    # Original 3-cycle logic
    table_samples_for_viz <- table_samples_by_cell_cycle %>%
      mutate(non_na_total = rowSums(dplyr::select(., G1, S, G2M), na.rm = TRUE))  
    
    table_clusters_for_viz <- table_clusters_by_cell_cycle %>%
      mutate(non_na_total = rowSums(dplyr::select(., G1, S, G2M), na.rm = TRUE))  
  }
  
  # Absolute composition of cell cycle
  temp_labels <- table_samples_for_viz %>%
    dplyr::select(resolution, total_cell_count) %>%
    mutate(
      non_na_count = table_samples_for_viz$non_na_total,
      label = paste0('n = ', format(total_cell_count, big.mark = ',', trim = TRUE))
    )
  
  # Define colors based on cycle type
  if (cycle == 4) {
    # Use first 4 colors from your custom palette
    cycle_colors <- custom_colors$discrete[1:4]
    names(cycle_colors) <- cycle_stages
  } else {
    # Use original 3 colors (assuming they match G1, S, G2M order)
    cycle_colors <- custom_colors$discrete[1:3]
    names(cycle_colors) <- cycle_stages
  }
  
  # Absolute count plot for samples
  p1 <- table_samples_for_viz %>%
    dplyr::select(resolution, all_of(cycle_stages)) %>%  
    reshape2::melt(id.vars = 'resolution') %>%
    mutate(
      sample = factor(resolution, levels = levels(seurat@meta.data$resolution)),
      variable = factor(variable, levels = cycle_stages)
    ) %>%
    ggplot(aes(resolution, value)) +
    geom_bar(aes(fill = variable), position = 'stack', stat = 'identity') +
    geom_text(
      data = temp_labels,
      aes(x = resolution, y = Inf, label = label, vjust = -1),
      color = 'black', size = 2.8
    ) +
    scale_fill_manual(name = 'Cell Cycle', values = cycle_colors) +
    scale_y_continuous(name = 'Number of cells', labels = scales::comma, expand = c(0.01, 0)) +
    coord_cartesian(clip = 'off') +
    theme_bw() +
    theme(
      legend.position = 'left',
      plot.title = element_text(hjust = 0.5),
      text = element_text(size = 16),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.title.x = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      plot.margin = margin(t = 20, r = 0, b = 0, l = 0, unit = 'pt')
    )
  
  temp_labels_clusters <- seurat@meta.data %>%
    group_by(seurat_clusters) %>%
    tally() %>%
    dplyr::rename('cluster' = seurat_clusters) %>%
    left_join(
      table_clusters_for_viz %>% dplyr::select(cluster, non_na_total),
      by = "cluster"
    ) %>%
    mutate(label = paste0('n = ', format(n, big.mark = ',', trim = TRUE)))
  
  p2 <- table_clusters_for_viz %>%
    dplyr::select(cluster, all_of(cycle_stages)) %>%  
    reshape2::melt(id.vars = 'cluster') %>%
    dplyr::mutate(
      cluster = factor(cluster, levels = levels(seurat@meta.data$seurat_clusters)),
      variable = factor(variable, levels = cycle_stages)
    ) %>%
    ggplot(aes(cluster, value)) +
    geom_bar(aes(fill = variable), position = 'stack', stat = 'identity') +
    geom_text(
      data = temp_labels_clusters,
      aes(x = cluster, y = Inf, label = label, vjust = -1),
      color = 'black', size = 2.8
    ) +
    scale_fill_manual(name = 'Cell Cycle', values = cycle_colors) +
    scale_y_continuous(labels = scales::comma, expand = c(0.01, 0)) +
    coord_cartesian(clip = 'off') +
    theme_bw() +
    theme(
      legend.position = 'right',
      plot.title = element_text(hjust = 0.5),
      text = element_text(size = 16),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.title = element_blank(),
      plot.margin = margin(t = 20, r = 0, b = 0, l = 10, unit = 'pt')
    )
  
  ggsave(
    file.path(working_dir, paste0('final/composition_cell_cycle.', runid, cycle_suffix, '.absolute.pdf')),
    p1 + p2 +
      plot_layout(ncol = 2, widths = c(
        seurat@meta.data$resolution %>% unique() %>% length() + 1,
        seurat@meta.data$seurat_clusters %>% unique() %>% length()
      )),
    width = plot_width,
    height = 8
  )
  
  message("Plot saved to: ", file.path(working_dir, paste0('final/composition_cell_cycle.', runid, cycle_suffix, '.absolute.pdf')))
  
  
  # Composition of cell cycle in percentage
  p1 <- table_samples_for_viz %>%
    dplyr::select(resolution, all_of(cycle_stages)) %>%  
    reshape2::melt(id.vars = 'resolution') %>%
    mutate(
      sample = factor(resolution, levels = levels(seurat@meta.data$resolution)),
      variable = factor(variable, levels = cycle_stages)
    ) %>%
    ggplot(aes(resolution, value)) +
    geom_bar(aes(fill = variable), position = 'fill', stat = 'identity') +
    geom_text(
      data = temp_labels,
      aes(x = resolution, y = Inf, label = label, vjust = -1),
      color = 'black', size = 2.8
    ) +
    scale_fill_manual(name = 'Cell Cycle', values = cycle_colors) +
    scale_y_continuous(name = 'Percentage of cells (%)', labels = scales::percent_format(), expand = c(0.01, 0)) +
    coord_cartesian(clip = 'off') +
    theme_bw() +
    theme(
      legend.position = 'left',
      plot.title = element_text(hjust = 0.5),
      text = element_text(size = 16),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.title.x = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      plot.margin = margin(t = 20, r = 0, b = 0, l = 0, unit = 'pt')
    )
  
  p2 <- table_clusters_for_viz %>%
    dplyr::select(cluster, all_of(cycle_stages)) %>%  
    reshape2::melt(id.vars = 'cluster') %>%
    dplyr::mutate(
      cluster = factor(cluster, levels = levels(seurat@meta.data$seurat_clusters)),
      variable = factor(variable, levels = cycle_stages)
    ) %>%
    ggplot(aes(cluster, value)) +
    geom_bar(aes(fill = variable), position = 'fill', stat = 'identity') +
    geom_text(
      data = temp_labels_clusters,
      aes(x = cluster, y = Inf, label = label, vjust = -1),
      color = 'black', size = 2.8
    ) +
    scale_fill_manual(name = 'Cell Cycle', values = cycle_colors) +
    scale_y_continuous(labels = scales::percent_format(), expand = c(0.01, 0)) +
    coord_cartesian(clip = 'off') +
    theme_bw() +
    theme(
      legend.position = 'right',
      plot.title = element_text(hjust = 0.5),
      text = element_text(size = 16),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.title = element_blank(),
      plot.margin = margin(t = 20, r = 0, b = 0, l = 10, unit = 'pt')
    )
  
  ggsave(
    file.path(working_dir, paste0('final/composition_cell_cycle.', runid, cycle_suffix, '.percentage.pdf')),
    p1 + p2 +
      plot_layout(ncol = 2, widths = c(
        seurat@meta.data$resolution %>% unique() %>% length() + 1,
        seurat@meta.data$seurat_clusters %>% unique() %>% length()
      )),
    width = plot_width, 
    height = 8
  )
  
  message("Plot saved to: ", file.path(working_dir, paste0('final/composition_cell_cycle.', runid, cycle_suffix, '.percentage.pdf')))
  
  # Clean up temporary columns
  seurat@meta.data$resolution <- NULL
  seurat@meta.data$cell_cycle <- NULL
}

cell_cycle_composition(seurat, runid, cycle = 3)
cell_cycle_composition(seurat, runid, cycle = 4)

### G0 scoring
# genelist_full_name <- list(
#   # [[1]] - Genes upregulated in quiescence
#   c("CFLAR", "CALCOCO1", "YPEL3", "CST3", "SERINC1", "CLIP4", "PCYOX1", "TMEM59", "RGS2", "YPEL5", "CD63", "KIAA1109", "CDH13", "GSN", "MR1", "CYB5R1", "AZGP1", "ZFYVE1", "DMXL1", "EPS8L2", "PTTG1IP", "MIR22HG", "PSAP", "GOLGA8B", "NEAT1", "TXNIP", "MTRNR2L12"),
#   
#   # [[2]] - Genes downregulated in quiescence (proliferation genes)
#   c("NCAPD2", "PTBP1", "MPHOSPH9", "NUCKS1", "TCOF1", "SART3", "SNRPA", "KIF22", "HSP90AA1", "WBP11", "CAD", "SF3B2", "KHSRP", "WDR76", "NUP188", "HSP90AB1", "HNRNPM", "SMARCB1", "PNN", "RBBP7", "NPRL3", "USP10", "SGTA", "MRPL4", "PSMD3", "KPNB1", "CBX1", "LRRC59", "TMEM97", "NSD2", "PRPF19", "PTGES3", "CPSF6", "SRSF3", "TCERG1", "SMC4", "EIF4G1", "ZNF142", "MSH6", "MRPL37", "SFPQ", "STMN1", "ARID1A", "PROSER1", "DDX39A", "EXOSC9", "USP22", "DEK", "DUT", "ILF3", "DNMT1", "NASP", "HMGB1P5", "SRRM1", "GNL2", "RNF138", "SRSF1", "TRA2B", "SMPD4", "ANP32B", "HMGA1", "MDC1", "HADH", "ARHGDIA", "PRCC", "HDGF", "SF3B4", "UBAP2L", "ILF2", "PARP1", "LBR", "CNOT9", "PPRC1", "SSRP1", "CCT5", "DLAT", "HNRNPU", "LARP1", "SCAF4", "RRP1B", "RRP1", "CHCHD4", "GMPS", "RFC4", "SLBP", "PSIP1", "HNRNPK", "SKA3", "DIS3L", "USP39", "GPS1", "PA2G4", "HCFC1", "SLC19A1", "ETV4", "RAD23A", "DCTPP1", "RCC1", "EWSR1", "ALYREF", "PTMA", "HMGB1", "POM121", "MCMBP", "TEAD4", "CHAMP1", "TOP1", "PRRC2A", "RBM14", "HMGB1P6", "POM121C", "UHRF1")
# )

query <- c("NCAPD2", "PTBP1", "MPHOSPH9", "NUCKS1", "TCOF1", "SART3", "SNRPA", "KIF22", "HSP90AA1", "WBP11", "CAD", "SF3B2", "KHSRP", "WDR76", "NUP188", "HSP90AB1", "HNRNPM", "SMARCB1", "PNN", "RBBP7", "NPRL3", "USP10", "SGTA", "MRPL4", "PSMD3", "KPNB1", "CBX1", "LRRC59", "TMEM97", "NSD2", "PRPF19", "PTGES3", "CPSF6", "SRSF3", "TCERG1", "SMC4", "EIF4G1", "ZNF142", "MSH6", "MRPL37", "SFPQ", "STMN1", "ARID1A", "PROSER1", "DDX39A", "EXOSC9", "USP22", "DEK", "DUT", "ILF3", "DNMT1", "NASP", "HMGB1P5", "SRRM1", "GNL2", "RNF138", "SRSF1", "TRA2B", "SMPD4", "ANP32B", "HMGA1", "MDC1", "HADH", "ARHGDIA", "PRCC", "HDGF", "SF3B4", "UBAP2L", "ILF2", "PARP1", "LBR", "CNOT9", "PPRC1", "SSRP1", "CCT5", "DLAT", "HNRNPU", "LARP1", "SCAF4", "RRP1B", "RRP1", "CHCHD4", "GMPS", "RFC4", "SLBP", "PSIP1", "HNRNPK", "SKA3", "DIS3L", "USP39", "GPS1", "PA2G4", "HCFC1", "SLC19A1", "ETV4", "RAD23A", "DCTPP1", "RCC1", "EWSR1", "ALYREF", "PTMA", "HMGB1", "POM121", "MCMBP", "TEAD4", "CHAMP1", "TOP1", "PRRC2A", "RBM14", "HMGB1P6", "POM121C", "UHRF1")

rownames(seurat[['integrated']])[grep("(?i)<query>", rownames(seurat[['integrated']]))]

for (gene in query) {
  result <- rownames(seurat[['integrated']])[grep(paste0("(?i)", gene), rownames(seurat[['integrated']]))]
  if (length(result) > 0) {
    quoted_results <- paste0('"', result, '"')
    cat('Gene: "', gene, '" -> Found: ', paste(quoted_results, collapse = ", "), "\n", sep = "")
  } else {
    cat('Gene: "', gene, '" -> Not found\n', sep = "")
  }
}

ciona_genelist_full <- list(
  # [[1]] = proliferation genes (downregulated in G0)
  c("CFLAR.2", "SERINC1.2", "SERINC1", "TMEM59.2", "CD63.3", "CD63.5", "CD63.2", "KIAA1109.3", "GSN.2", "GSN", "PTTG1IP.2", "PSAP"),
  
  # [[2]] = quiescence genes (upregulated in G0)
  c("NCAPD2", "PTBP1.2", "SART3", "SNRPA1", "HSP90AA1", "WBP11", "CAD", "SF3B2", "HNRNPM", "SMARCB1", "USP10", "KPNB1", "LRRC59.2", "Ci-ZF307 (NSD2)", "PTGES3", "TCERG1.2", "SMC4", "DEK", "DUT", "NASP", "SRRM1", "GNL2", "SRSF1", "HADH", "ARHGDIA", "ARHGDIA.2", "PARP1.2",  "PARP1", "LBR", "PPRC1", "Ssrp (SSRP1)", "DLAT", "HNRNPUL1", "RRP1B.2", "CHCHD4", "SLBP", "HCFC1", "PRRC2A", "UHRF1")
)

assignInNamespace("genelist_full_name", ciona_genelist_full, ns = "QuieScore")
# assignInNamespace("genelist_core_name", ciona_genelist_core, ns = "QuieScore")

print("New full gene list:")
print(QuieScore:::genelist_full_name)
# print("New core gene list:")
# print(QuieScore:::genelist_core_name)

counts <- GetAssayData(seurat, assay = "SCT", slot = "data")
counts <- as.data.frame(as.matrix(counts))

processedData <- processInput(counts, cancer_type = "ESCA", 
                              gene_naming = "name", log_transformed = TRUE)

G0scores <- QuiescenceScore(processedData)

seurat$g0_score_raw <- G0scores$q_score_raw
seurat$g0_score_normalized <- G0scores$q_score_normalized
seurat$g0_percentile <- G0scores$q_percentile

g0_score_violin_plots <- function(
    seurat, 
    runid,
    resolution = "merged_sample", # sample, merged_sample, combined_sample, meta_sample
    score_column = "g0_score_normalized", # g0_score_raw, g0_score_normalized, g0_percentile
    plot_width = 16
) {
  
  # Create temporary columns for dynamic referencing
  seurat@meta.data$resolution <- seurat@meta.data[[resolution]]
  seurat@meta.data$score <- seurat@meta.data[[score_column]]
  
  # Prepare data for plotting
  plot_data_samples <- seurat@meta.data %>%
    dplyr::select(resolution, score) %>%
    dplyr::filter(!is.na(score)) %>%
    dplyr::mutate(resolution = factor(resolution, levels = levels(seurat@meta.data$resolution)))
  
  plot_data_clusters <- seurat@meta.data %>%
    dplyr::select(seurat_clusters, score) %>%
    dplyr::filter(!is.na(score)) %>%
    dplyr::mutate(seurat_clusters = factor(seurat_clusters, levels = levels(seurat@meta.data$seurat_clusters)))
  
  n_samples <- length(unique(plot_data_samples$resolution))
  n_clusters <- length(unique(plot_data_clusters$seurat_clusters))
  
  sample_colors <- custom_colors$discrete[1:n_samples]
  cluster_colors <- custom_colors$discrete[1:n_clusters]
  
  p1 <- ggbetweenstats(
    data = plot_data_samples,
    x = resolution,
    y = score,
    type = "nonparametric",
    pairwise.comparisons = TRUE,
    pairwise.display = "significant",
    centrality.plotting = TRUE,
    centrality.point.args = list(size = 0, color = "darkred"),
    centrality.label.args = list(alpha = 0),
    p.adjust.method = "bonferroni",
    violin.args = list(width = 0.8, alpha = 0.7),
    point.args = list(position = ggplot2::position_jitterdodge(jitter.width = 0.8), 
                      alpha = 0.4, size = 1.5, stroke = 0, na.rm = TRUE),
    title = "G0 Score by Sample",
    results.subtitle = FALSE,
    xlab = "",
    ylab = paste(score_column, "Score"),
    ggplot.component = list(
      scale_fill_manual(values = sample_colors),
      scale_color_manual(values = sample_colors),
      theme(axis.title.y.right = element_blank(), 
            axis.text.y.right = element_blank(), 
            axis.ticks.y.right = element_blank())
    )
  )
  
  # Create violin plot for clusters
  p2 <- ggbetweenstats(
    data = plot_data_clusters,
    x = seurat_clusters,
    y = score,
    type = "nonparametric",
    pairwise.comparisons = TRUE,
    pairwise.display = "none",
    centrality.plotting = TRUE,
    centrality.point.args = list(size = 0, color = "darkred"),
    centrality.label.args = list(alpha = 0),
    p.adjust.method = "bonferroni",
    violin.args = list(width = 0.8, alpha = 0.7),
    point.args = list(position = ggplot2::position_jitterdodge(jitter.width = 0.8), 
                      alpha = 0.4, size = 1.5, stroke = 0, na.rm = TRUE),
    title = "G0 Score by Cluster",
    results.subtitle = FALSE,
    xlab = "Cluster",
    ylab = paste(score_column, "Score"),
    ggplot.component = list(
      scale_fill_manual(values = cluster_colors),
      scale_color_manual(values = cluster_colors),
      theme(axis.title.y.right = element_blank(), 
            axis.text.y.right = element_blank(), 
            axis.ticks.y.right = element_blank())
    )
  )
  
  # # Save individual plots
  # ggsave(
  #   file.path(working_dir, paste0('final/g0_score_violin_samples.', runid, '.', score_column, '.pdf')),
  #   p1,
  #   width = max(8, length(unique(plot_data_samples$resolution)) * 1.5),
  #   height = 8
  # )
  # 
  # ggsave(
  #   file.path(working_dir, paste0('final/g0_score_violin_clusters.', runid, '.', score_column, '.pdf')),
  #   p2,
  #   width = max(8, length(unique(plot_data_clusters$seurat_clusters)) * 0.8),
  #   height = 8
  # )
  
  ggsave(
    file.path(working_dir, paste0('final/g0_score_violin.', runid, ".", score_column, '.pdf')),
    p1 + p2 +
      plot_layout(ncol = 2, widths = c(
        seurat@meta.data$resolution %>% unique() %>% length() + 1,
        seurat@meta.data$seurat_clusters %>% unique() %>% length()
      )),
    width = plot_width,
    height = 8
  )
  
  cat("Summary statistics for", score_column, "by", resolution, ":\n")
  summary_samples <- plot_data_samples %>%
    group_by(resolution) %>%
    dplyr::summarise(
      n = n(),
      mean = round(mean(score, na.rm = TRUE), 3),
      median = round(median(score, na.rm = TRUE), 3),
      sd = round(sd(score, na.rm = TRUE), 3),
      .groups = 'drop'
    )
  print(summary_samples)
  
  cat("\nSummary statistics for", score_column, "by cluster:\n")
  summary_clusters <- plot_data_clusters %>%
    group_by(seurat_clusters) %>%
    dplyr::summarise(
      n = n(),
      mean = round(mean(score, na.rm = TRUE), 3),
      median = round(median(score, na.rm = TRUE), 3),
      sd = round(sd(score, na.rm = TRUE), 3),
      .groups = 'drop'
    )
  print(summary_clusters)
  
  # Prepare data for CSV with separator
  # Add section headers and separator
  samples_with_header <- rbind(
    data.frame(resolution = paste("Summary statistics for", score_column, "by", resolution), 
               n = "", mean = "", median = "", sd = "", stringsAsFactors = FALSE),
    setNames(data.frame(summary_samples, stringsAsFactors = FALSE), c("resolution", "n", "mean", "median", "sd"))
  )
  
  separator_row <- data.frame(resolution = "", n = "", mean = "", median = "", sd = "", stringsAsFactors = FALSE)
  
  clusters_with_header <- rbind(
    data.frame(resolution = paste("Summary statistics for", score_column, "by cluster"), 
               n = "", mean = "", median = "", sd = "", stringsAsFactors = FALSE),
    setNames(data.frame(summary_clusters, stringsAsFactors = FALSE), c("resolution", "n", "mean", "median", "sd"))
  )
  
  # Combine all parts
  final_table <- rbind(
    samples_with_header,
    separator_row,
    clusters_with_header
  )
  
  # Write to CSV
  write.csv(
    final_table,
    file.path(working_dir, paste0('g0_score_summary.', runid, ".", score_column, '.csv')),
    row.names = FALSE
  )
  
  message("Plots saved to: ", file.path(working_dir, paste0('final/g0_score_violin.', runid, ".", score_column, '.pdf')))
  message("Summary statistics saved to: ", file.path(working_dir, paste0('final/g0_score_violin.', runid, ".", score_column, ".csv")))
  
  # Clean up temporary columns
  seurat@meta.data$resolution <- NULL
  seurat@meta.data$score <- NULL

}  

g0_score_violin_plots(seurat, runid = "single.4", 
                      resolution = "merged_sample",
                      score_column = "g0_score_raw")

g0_score_violin_plots(seurat, runid = "single.3", 
                      resolution = "combined_sample",
                      score_column = "g0_score_raw")

g0_score_violin_plots(seurat, runid = "single.4", 
                      resolution = "merged_sample",
                      score_column = "g0_percentile")

g0_score_violin_plots(seurat, runid = "single.3", 
                      resolution = "merged_sample",
                      score_column = "g0_percentile")

g0_table <- seurat@meta.data %>%
  group_by(seurat_clusters, merged_sample) %>%
  dplyr::summarise(avg_g0_score = round(mean(g0_score_raw, na.rm = TRUE), 3),
            .groups = 'drop') %>%
  tidyr::pivot_wider(names_from = merged_sample, 
                     values_from = avg_g0_score, 
                     values_fill = NA) %>%
  arrange(as.numeric(as.character(seurat_clusters)))

g0_table <- seurat@meta.data %>%
  group_by(cell_type, merged_sample) %>%
  dplyr::summarise(avg_g0_score = round(mean(g0_score_raw, na.rm = TRUE), 3),
                   .groups = 'drop') %>%
  tidyr::pivot_wider(names_from = merged_sample, 
                     values_from = avg_g0_score, 
                     values_fill = NA) %>%
  arrange(as.numeric(as.character(cell_type)))

knitr::kable(g0_table,
             caption = "Average G0 Score by Cell Type and Sample",
             align = c("l", rep("c", ncol(g0_table)-1)),
             digits = 3,
             na = "-")

g0_score_by_celltype_violin <- function(
    seurat, 
    runid,
    resolution = "merged_sample", # sample, merged_sample, combined_sample, meta_sample
    score_column = "g0_score_raw", # g0_score_raw, g0_score_normalized, g0_percentile
    plot_width = 18
) {
  
  # Create temporary columns for dynamic referencing
  seurat@meta.data$resolution <- seurat@meta.data[[resolution]]
  seurat@meta.data$score <- seurat@meta.data[[score_column]]
  
  # Get unique cell types
  cell_types <- unique(seurat@meta.data$cell_type)
  cell_types <- cell_types[!is.na(cell_types)]
  
  # Determine colors
  n_samples <- length(unique(seurat@meta.data$resolution))
  sample_colors <- custom_colors$discrete[1:n_samples]
  
  # Create list to store plots
  plot_list <- list()
  
  # Create a plot for each cell type
  for(i in 1:length(cell_types)) {
    ct <- cell_types[i]
    
    # Filter data for this cell type
    plot_data <- seurat@meta.data %>%
      dplyr::filter(cell_type == ct) %>%
      dplyr::select(resolution, score) %>%
      dplyr::filter(!is.na(score)) %>%
      dplyr::mutate(resolution = factor(resolution, levels = levels(seurat@meta.data$resolution)))
    
    # Skip if no data for this cell type
    if(nrow(plot_data) == 0) next
    
    # Create violin plot for this cell type
    p <- ggbetweenstats(
      data = plot_data,
      x = resolution,
      y = score,
      type = "nonparametric",
      pairwise.comparisons = FALSE,  # Set to FALSE to avoid overcrowding
      centrality.plotting = TRUE,
      centrality.point.args = list(size = 0, color = "darkred"),
      centrality.label.args = list(alpha = 0),
      p.adjust.method = "bonferroni",
      violin.args = list(width = 0.8, alpha = 0.7),
      point.args = list(position = ggplot2::position_jitterdodge(jitter.width = 0.8), 
                        alpha = 0.4, size = 1.0, stroke = 0, na.rm = TRUE),
      title = paste("G0 Score:", ct),
      results.subtitle = FALSE,
      xlab = "",
      ylab = ifelse(i == 1, paste(score_column, "Score"), ""),  # Only show y-label on first plot
      ggplot.component = list(
        scale_fill_manual(values = sample_colors),
        scale_color_manual(values = sample_colors),
        theme(
          axis.title.y.right = element_blank(), 
          axis.text.y.right = element_blank(), 
          axis.ticks.y.right = element_blank(),
          plot.title = element_text(size = 12, hjust = 0.5)
        )
      )
    )
    
    plot_list[[ct]] <- p
  }
  
  # Combine plots using patchwork
  if(length(plot_list) == 3) {
    combined_plot <- plot_list[[1]] | plot_list[[2]] | plot_list[[3]]
  } else if(length(plot_list) == 2) {
    combined_plot <- plot_list[[1]] | plot_list[[2]]
  } else if(length(plot_list) == 1) {
    combined_plot <- plot_list[[1]]
  } else {
    stop("No valid cell types found")
  }
  
  # Save plot
  ggsave(
    file.path(working_dir, paste0('final/g0_score_by_celltype_violin.', runid, '.', score_column, '.pdf')),
    combined_plot,
    width = plot_width,
    height = 8
  )
  
  # Clean up temporary columns
  seurat@meta.data$resolution <- NULL
  seurat@meta.data$score <- NULL
  
  message("Plot saved to: ", file.path(working_dir, paste0('final/g0_score_by_celltype_violin.', runid, '.', score_column, '.pdf')))
  
}

g0_score_by_celltype_violin(seurat, runid, 
  resolution = "merged_sample",
  score_column = "g0_score_raw",
  plot_width = 10
)

g0_score_by_resolution_violin <- function(
    seurat, 
    runid,
    resolution = "merged_sample", # X-axis variable (sample-describing parameter)
    spread = "seurat_clusters", # What to create separate plots for
    score_column = "g0_score_raw", # g0_score_raw, g0_score_normalized, g0_percentile
    individual_width = 3.5,
    individual_height = 8
) {
  
  # Create temporary columns for dynamic referencing
  seurat@meta.data$x_var <- seurat@meta.data[[resolution]]  # X-axis variable
  seurat@meta.data$spread_var <- seurat@meta.data[[spread]]  # What to split plots by
  seurat@meta.data$score <- seurat@meta.data[[score_column]]
  
  spread_groups <- unique(seurat@meta.data$spread_var)
  spread_groups <- spread_groups[!is.na(spread_groups)]
  
  if(is.factor(seurat@meta.data$spread_var)) {
    spread_groups <- levels(seurat@meta.data$spread_var)[levels(seurat@meta.data$spread_var) %in% spread_groups]
  }
  
  # Determine colors for x-axis variable
  n_x_groups <- length(unique(seurat@meta.data$x_var))
  x_colors <- custom_colors$discrete[1:n_x_groups]
  
  # Create list to store all plots
  all_plots <- list()
  
  # Create a plot for each spread group
  for(i in 1:length(spread_groups)) {
    group <- spread_groups[i]
    
    # Filter data for this spread group
    plot_data <- seurat@meta.data %>%
      dplyr::filter(spread_var == group) %>%
      dplyr::select(x_var, score) %>%
      dplyr::filter(!is.na(score)) %>%
      dplyr::mutate(x_var = factor(x_var, levels = levels(seurat@meta.data$x_var)))
    
    # Skip if no data for this group
    if(nrow(plot_data) == 0) next
    
    # Create title based on spread variable
    if(spread == "seurat_clusters") {
      plot_title <- paste("Cluster:", group)
    } else if(spread == "cell_type") {
      plot_title <- paste("Cell Type:", group)
    } else {
      plot_title <- paste(gsub("_", " ", tools::toTitleCase(spread)), ":", group)
    }
    
    if (resolution == "merged_sample" | resolution == "combined_sample" | resolution == "sample" | resolution == "meta_sample") {
      x_label <- "Sample"
    } else {
      x_label <- gsub("_", " ", tools::toTitleCase(resolution))
    }
    
    # Create violin plot for this group
    p <- ggbetweenstats(
      data = plot_data,
      x = x_var,
      y = score,
      type = "nonparametric",
      pairwise.comparisons = FALSE,
      centrality.plotting = TRUE,
      centrality.point.args = list(size = 0, color = "darkred"),
      centrality.label.args = list(alpha = 0),
      p.adjust.method = "bonferroni",
      violin.args = list(width = 0.8, alpha = 0.3),
      point.args = list(position = ggplot2::position_jitterdodge(jitter.width = 0.8), 
                        alpha = 0.4, size = 1.0, stroke = 0, na.rm = TRUE),
      title = plot_title,
      results.subtitle = FALSE,
      xlab = x_label,
      ylab = paste(score_column),
      ggplot.component = list(
        scale_fill_manual(values = x_colors),
        scale_color_manual(values = x_colors),
        theme(
          axis.title.y.right = element_blank(), 
          axis.text.y.right = element_blank(), 
          axis.ticks.y.right = element_blank(),
          plot.title = element_text(size = 12, hjust = 0.5)
        )
      )
    )
    
    all_plots[[paste0(spread, "_", group)]] <- p
  }
  
  # Calculate layout dimensions
  total_plots <- length(all_plots)
  n_cols <- ceiling(sqrt(total_plots))  # Roughly square layout
  n_rows <- ceiling(total_plots / n_cols)
  
  # Alternative: arrange in a single row if not too many plots
  if(total_plots <= 6) {
    n_cols <- total_plots
    n_rows <- 1
  }
  
  # Combine all plots using patchwork
  if(total_plots > 0) {
    combined_plot <- wrap_plots(all_plots, ncol = n_cols, nrow = n_rows)
    
    # Add overall title
    spread_label <- if(spread == "seurat_clusters") {
      "Clusters"
    } else if(spread == "cell_type") {
      "Cell Types"
    } else {
      gsub("_", " ", tools::toTitleCase(spread))
    }
    
    # Calculate total plot dimensions
    total_width <- n_cols * individual_width
    total_height <- n_rows * individual_height + 1  # +1 for title space
    
    # Save plot
    plot_suffix <- paste0("by_", resolution, "_across_", spread)
    ggsave(
      file.path(working_dir, paste0('final/g0_score_violin_', plot_suffix, '.', runid, '.', score_column, '.pdf')),
      combined_plot,
      width = total_width,
      height = total_height
    )
    
    message("Plot saved to: ", file.path(working_dir, paste0('final/g0_score_violin_', plot_suffix, '.', runid, '.', score_column, '.pdf')))
    
  } else {
    stop("No valid plots could be created")
  }
  
  seurat@meta.data$x_var <- NULL
  seurat@meta.data$spread_var <- NULL
  seurat@meta.data$score <- NULL
}

g0_score_by_resolution_violin(seurat, runid, 
  resolution = "merged_sample", 
  spread = "cell_type",
  score_column = "g0_score_raw"
)

g0_score_by_resolution_violin(seurat, runid, 
  resolution = "merged_sample", 
  spread = "seurat_clusters",
  score_column = "g0_score_raw"
)


# ## reclustering based on cell cycle regression
# test <- ScaleData(seurat, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(seurat))
# test <- FindVariableFeatures(test)
# test <- RunPCA(test, features = VariableFeatures(test), npcs = 50)
# DimPlot(test, reduction = "pca") + DimPlot(seurat, reduction = "pca")
# 
# intrinsicDimension::maxLikGlobalDimEst(test@reductions$pca@cell.embeddings, k = 10)
# # 7.194341 c.f. 6.367669
# intrinsicDimension::maxLikGlobalDimEst(test@reductions$pca@cell.embeddings, k = 20)
# # 7.591661 c.f. 6.970843
# 
# pca_cutoff <- 10
# 
# resolution <- 0.4
# 
# clustering_methods <- list(
#   Louvain = 1,
#   Hierarchical = 3,
#   DBSCAN = 4
# )
# method_name <- "Hierarchical"
# algorithm <- clustering_methods[[method_name]]
# 
# test <- FindNeighbors(test, reduction = 'pca', dims = 1:pca_cutoff)
# test <- FindClusters(test, resolution = resolution, algorithm = algorithm, random.seed = 42) # 0.8 by default
# 
# test <- RunUMAP(
#   test,
#   reduction.name = 'UMAP',
#   reduction = 'pca',
#   dims = 1:pca_cutoff,
#   n.components = 2,
#   seed.use = 100
# )
# 
# DimPlot(test, reduction = "UMAP") + DimPlot(seurat, reduction = "UMAP")
# 
# test <- RunTSNE(
#   test,
#   reduction.name = 'tSNE',
#   reduction = 'pca',
#   dims = 1:pca_cutoff,
#   seed.use = 100
# )
# 
# DimPlot(test, reduction = "tSNE") + DimPlot(seurat, reduction = "tSNE")
# 
# test <- readRDS(file.path(working_dir, "s+g2m.rds"))
# # saveRDS(test, file.path(working_dir, "s+g2m.rds"))
# 
# ## regress on s - g2m
# test <- seurat
# test$CC.Difference <- test$S.Score - test$G2M.Score
# test <- ScaleData(test, vars.to.regress = "CC.Difference")
# 
# test <- FindVariableFeatures(test)
# test <- RunPCA(test, features = VariableFeatures(test), npcs = 50)
# DimPlot(test, reduction = "pca") + DimPlot(seurat, reduction = "pca")
# 
# intrinsicDimension::maxLikGlobalDimEst(test@reductions$pca@cell.embeddings, k = 10)
# # 7.119778 c.f. 6.367669
# intrinsicDimension::maxLikGlobalDimEst(test@reductions$pca@cell.embeddings, k = 20)
# # 7.449233 c.f. 6.970843
# 
# pca_cutoff <- 10
# 
# resolution <- 0.4
# 
# clustering_methods <- list(
#   Louvain = 1,
#   Hierarchical = 3,
#   DBSCAN = 4
# )
# method_name <- "Hierarchical"
# algorithm <- clustering_methods[[method_name]]
# 
# test <- FindNeighbors(test, reduction = 'pca', dims = 1:pca_cutoff)
# test <- FindClusters(test, resolution = resolution, algorithm = algorithm, random.seed = 42) # 0.8 by default
# 
# test <- RunUMAP(
#   test,
#   reduction.name = 'UMAP',
#   reduction = 'pca',
#   dims = 1:pca_cutoff,
#   n.components = 2,
#   seed.use = 100
# )
# 
# DimPlot(test, reduction = "UMAP") + DimPlot(seurat, reduction = "UMAP")
# 
# test <- RunTSNE(
#   test,
#   reduction.name = 'tSNE',
#   reduction = 'pca',
#   dims = 1:pca_cutoff,
#   seed.use = 100
# )
# 
# DimPlot(test, reduction = "tSNE") + DimPlot(seurat, reduction = "tSNE")
# 
# test <- readRDS(file.path(working_dir, "s-g2m.rds"))
# # saveRDS(test, file.path(working_dir, "s-g2m.rds"))
# 
# test <- readRDS(file.path(working_dir, "5.singlet.s-g2m.rds"))
# # saveRDS(test, file.path(working_dir, "t.singlet.s-g2m.rds"))
# 
# ## visualisation of cluster identity assignment among different clustering settings
# cluster_df <- data.frame(
#   Original = seurat$seurat_clusters,
#   Regressed = test$seurat_clusters
# )
# 
# # Create a contingency table
# cluster_matrix <- table(cluster_df)
# prop_matrix <- prop.table(cluster_matrix, margin = 1)
# 
# # Convert the table to a data frame for ggplot
# matrix_df <- as.data.frame.table(prop_matrix)
# colnames(matrix_df) <- c("Original_Cluster", "Regressed_Cluster", "Proportion")
# 
# # Create a heatmap
# p <- ggplot(matrix_df, aes(x = Regressed_Cluster, y = Original_Cluster, fill = Proportion)) +
#   geom_tile(color = 'white') +
#   geom_text(aes(label = round(Proportion, digits = 2)), size = 4) +
#   scale_x_discrete(name = 'Cell Cycle Regressed Clusters', position = 'top') +
#   scale_y_discrete(
#     name = 'Clusters',
#     limits = rev(levels(matrix_df$Original_Cluster))  
#   ) +  scale_fill_gradient(
#     low = 'white', high = '#c0392b', na.value = '#bdc3c7',
#     guide = guide_colorbar(
#       frame.colour = 'black', ticks.colour = 'black', title.position = 'left',
#       title.theme = element_text(hjust = 1, angle = 90),
#       barwidth = 0.75, barheight = 10
#     )
#   ) +  
#   theme_bw() +
#   coord_fixed() +
#   labs(title = "Cluster Assignment Matrix",
#        # x = "Cell Cycle Regressed Clusters",
#        # y = "Original Clusters",
#        fill = "Proportion") +
#   theme(
#     legend.position = 'right',
#     panel.grid.major = element_blank()
#   )
# 
# ggsave(
#   file.path(working_dir, 'cluster_asgn_s+g.png'),
#   # file.path(working_dir, 'cluster_asgn_s-g.png'),
#   p,
#   height = 5, width = 11
# )
# 
# sce <- as.SingleCellExperiment(test)
# 
# reducedDim(sce, 'PCA_sub') <- reducedDim(sce, 'PCA')[,1:pca_cutoff, drop = FALSE]
# 
# g <- scran::buildSNNGraph(sce, use.dimred = 'PCA_sub')
# 
# ratio <- bluster::pairwiseModularity(g, test@meta.data$seurat_clusters, as.ratio = TRUE)
# 
# ratio_to_plot <- log10(ratio+1)
# 
# # library(reshape2)
# reshaped_df <- melt(as.matrix(ratio_to_plot))
# colnames(reshaped_df) <- c("cluster_1", "cluster_2", "probability")
# 
# reshaped_df$cluster_1 <- factor(reshaped_df$cluster_1, levels = rownames(ratio_to_plot))
# reshaped_df$cluster_2 <- factor(reshaped_df$cluster_2, levels = colnames(ratio_to_plot))
# 
# p <- ggplot(reshaped_df, aes(cluster_2, cluster_1, fill = probability)) +
#   geom_tile(color = 'white') +
#   geom_text(aes(label = round(probability, digits = 2)), size = 2.5) +
#   scale_x_discrete(name = 'Cluster', position = 'top') +
#   scale_y_discrete(name = 'Cluster', limits = rev(levels(reshaped_df$cluster_1))) +
#   scale_fill_gradient(
#     name = 'log10(ratio)', low = 'white', high = '#c0392b', na.value = '#bdc3c7',
#     guide = guide_colorbar(
#       frame.colour = 'black', ticks.colour = 'black', title.position = 'left',
#       title.theme = element_text(hjust = 1, angle = 90),
#       barwidth = 0.75, barheight = 10
#     )
#   ) +
#   coord_fixed() +
#   theme_bw() +
#   theme(
#     legend.position = 'right',
#     panel.grid.major = element_blank()
#   )
# 
# ggsave(
#   # file.path(working_dir, 'cluster_similarity_s+g.png'), 
#   file.path(working_dir, 'cluster_similarity_s-g.png'),
#   p, 
#   height = 6, 
#   width = 7
# )

## single r annotation
library(SingleR)
singler_ref <- HumanPrimaryCellAtlasData()

singler_results_blueprintencode_main <- SingleR::SingleR(
  test = GetAssayData(test, assay = 'integrated', slot = 'scale.data'),
  ref = singler_ref,
  labels = singler_ref@colData@listData$label.main
)

singler_results_blueprintencode_fine <- SingleR::SingleR(
  test = GetAssayData(test, assay = 'integrated', slot = 'scale.data'),
  ref = singler_ref,
  labels = singler_ref@colData@listData$label.fine
)

saveRDS(singler_results_blueprintencode_main, file.path(working_dir, "singler_main_s-g.rds"))
saveRDS(singler_results_blueprintencode_fine, file.path(working_dir, "singler_fine_s-g.rds"))
singler_results_blueprintencode_main <- readRDS(file.path(working_dir, "singler_main_s-g.rds"))

# p <- plotScoreHeatmap(
#   singler_results_blueprintencode_main,
#   clusters = test@meta.data$seurat_clusters,
#   show.labels = TRUE,
#   cluster_rows = FALSE,
#   # cells.use = sorted_cell_order,
#   # cells.order = sorted_cell_order,
# )
# 
# # Create a dataframe with both labels and clusters, properly matched by rownames
# merged_df <- data.frame(
#   cell_id = rownames(singler_results_blueprintencode_main),
#   label = singler_results_blueprintencode_main$labels,
#   stringsAsFactors = FALSE
# )
# 
# # Ensure the cell IDs exist in the test metadata before merging
# valid_cells <- merged_df$cell_id[merged_df$cell_id %in% rownames(test@meta.data)]
# merged_df <- merged_df[merged_df$cell_id %in% valid_cells, ]
# 
# # Add the cluster information
# merged_df$cluster <- test@meta.data[merged_df$cell_id, "seurat_clusters"]
# 
# # Now sort by label first, then by cluster
# merged_df <- merged_df[order(merged_df$label, merged_df$cluster), ]
# 
# # Get the ordered cell IDs
# sorted_cell_order <- merged_df$cell_id
# 
# # Extract the scores matrix
# scores_matrix <- as.matrix(singler_results_blueprintencode_main$scores)
# rownames(scores_matrix) <- rownames(singler_results_blueprintencode_main)
# 
# # Reorder the scores matrix using the sorted cell order
# ordered_scores <- scores_matrix[sorted_cell_order, ]
# 
# # Create annotation dataframe for the heatmap
# annotation_df <- data.frame(
#   SingleR_label = merged_df$label,
#   Seurat_cluster = merged_df$cluster,
#   row.names = merged_df$cell_id
# )
# 
# # Extract unique clusters and labels for color assignment
# clusters <- sort(unique(annotation_df$Seurat_cluster))
# 
# labels <- sort(unique(annotation_df$SingleR_label))
# 
# # Create the heatmap
# library(pheatmap)
# p <- pheatmap(
#   ordered_scores,
#   cluster_rows = FALSE,
#   cluster_cols = TRUE,
#   annotation_row = annotation_df,
#   # annotation_colors = custom_colors$discrete,
#   show_rownames = FALSE,
#   color = viridis::viridis(10, option = "E"),
#   main = "SingleR Cell Type Scores"
# )
# 
# ggsave(
#   file.path(working_dir, 'singler_dimheatmap_main.png'), p,
#   height = 12, width = 12
# )
# 
# p <- plotScoreHeatmap(
#   singler_results_blueprintencode_fine,
#   show.labels = TRUE,
#   cluster_rows = FALSE,
#   annotation_col = data.frame(
#     donor = test@meta.data$sample,
#     row.names = rownames(singler_results_blueprintencode_fine)
#   )
# )
# 
# ggsave(
#   file.path(working_dir, 'singler_dimheatmap_fine.png'), p,
#   height = 12, width = 12
# )
# 
# test@meta.data$cell_type_singler_blueprintencode_main <- singler_results_blueprintencode_main@listData$labels
# 
# singler_scores <- singler_results_blueprintencode_main@listData$scores %>%
#   as_tibble() %>%
#   dplyr::mutate(assigned_score = NA)
# 
# for (i in seq_len(nrow(singler_scores))) {
#   singler_scores$assigned_score[i] <- singler_scores[[singler_results_blueprintencode_main@listData$labels[i]]][i]
# }
# 
# test@meta.data$cell_type_singler_blueprintencode_main_score <- singler_scores$assigned_score
# 
# test@meta.data$cell_type_singler_blueprintencode_fine <- singler_results_blueprintencode_fine@listData$labels
# 
# singler_scores <- singler_results_blueprintencode_fine@listData$scores %>%
#   as_tibble() %>%
#   dplyr::mutate(assigned_score = NA)
# 
# for (i in seq_len(nrow(singler_scores))) {
#   singler_scores$assigned_score[i] <- singler_scores[[singler_results_blueprintencode_fine@listData$labels[i]]][i]
# }
# 
# test@meta.data$cell_type_singler_blueprintencode_fine_score <- singler_scores$assigned_score
# 
# temp_labels <- test@meta.data %>%
#   group_by(merged_sample) %>%
#   tally()
# 
# p1 <- ggplot() +
#   geom_half_violin(
#     data = test@meta.data,
#     aes(
#       x = merged_sample,
#       y = cell_type_singler_blueprintencode_main_score,
#       fill = sample
#     ),
#     side = 'l', show.legend = FALSE, trim = FALSE
#   ) +
#   geom_half_boxplot(
#     data = test@meta.data,
#     aes(
#       x = merged_sample,
#       y = cell_type_singler_blueprintencode_main_score,
#       fill = sample
#     ),
#     side = 'r', outlier.color = NA, width = 0.4, show.legend = FALSE
#   ) +
#   geom_text(
#     data = temp_labels,
#     aes(
#       x = merged_sample,
#       y = -Inf,
#       label = paste0('n = ', format(n, big.mark = ',', trim = TRUE)),
#       vjust = -1
#     ),
#     color = 'black', size = 2.8
#   ) +
#   scale_color_manual(values = custom_colors$discrete) +
#   scale_fill_manual(values = custom_colors$discrete) +
#   scale_y_continuous(
#     name = 'Assignment score',
#     labels = scales::comma,
#     limits = c(0,1)
#   ) +
#   theme_bw() +
#   labs(title = 'Samples') +
#   theme(
#     axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
#     axis.title.x = element_blank(),
#     panel.grid.major.x = element_blank()
#   )
# 
# temp_labels <- test@meta.data %>%
#   group_by(seurat_clusters) %>%
#   tally()
# 
# p2 <- ggplot() +
#   geom_half_violin(
#     data = test@meta.data,
#     aes(
#       x = seurat_clusters,
#       y = cell_type_singler_blueprintencode_main_score,
#       fill = seurat_clusters
#     ),
#     side = 'l', show.legend = FALSE, trim = FALSE
#   ) +
#   geom_half_boxplot(
#     data = test@meta.data,
#     aes(
#       x = seurat_clusters,
#       y = cell_type_singler_blueprintencode_main_score,
#       fill = seurat_clusters
#     ),
#     side = 'r', outlier.color = NA, width = 0.4, show.legend = FALSE
#   ) +
#   geom_text(
#     data = temp_labels,
#     aes(
#       x = seurat_clusters,
#       y = -Inf,
#       label = paste0('n = ', format(n, big.mark = ',', trim = TRUE)),
#       vjust = -1
#     ),
#     color = 'black', size = 2.8
#   ) +
#   scale_color_manual(values = custom_colors$discrete) +
#   scale_fill_manual(values = custom_colors$discrete) +
#   scale_y_continuous(
#     name = 'Assignment score',
#     labels = scales::comma,
#     limits = c(0,1)
#   ) +
#   labs(title = 'Clusters') +
#   theme_bw() +
#   theme(
#     axis.title.x = element_blank(),
#     panel.grid.major.x = element_blank()
#   )
# 
# cell_counts <- table(test@meta.data$cell_type_singler_blueprintencode_main)
# 
# # Identify cell types with more than one observation
# valid_cell_types <- names(cell_counts[cell_counts > 1])
# 
# # Filter the metadata to include only these cell types
# filtered_meta <- test@meta.data %>%
#   filter(cell_type_singler_blueprintencode_main %in% valid_cell_types)
# 
# # Update temp_labels to match filtered data
# temp_labels <- filtered_meta %>%
#   group_by(cell_type_singler_blueprintencode_main) %>%
#   tally()
# 
# p3 <- ggplot() +
#   gghalves::geom_half_violin(
#     data = filtered_meta,
#     mapping = aes(
#       x = cell_type_singler_blueprintencode_main,
#       y = cell_type_singler_blueprintencode_main_score,
#       fill = cell_type_singler_blueprintencode_main
#     ),
#     side = "left",
#     show.legend = FALSE,
#     trim = FALSE
#   ) +
#   gghalves::geom_half_boxplot(
#     data = filtered_meta,
#     mapping = aes(
#       x = cell_type_singler_blueprintencode_main,
#       y = cell_type_singler_blueprintencode_main_score,
#       fill = cell_type_singler_blueprintencode_main
#     ),
#     side = "right",
#     outlier.color = NA, 
#     width = 0.4, 
#     show.legend = FALSE
#   ) +
#   geom_text(
#     data = temp_labels,
#     aes(
#       x = cell_type_singler_blueprintencode_main,
#       y = 0,  # Using 0 instead of -Inf for better positioning
#       label = paste0('n = ', format(n, big.mark = ',', trim = TRUE))
#     ),
#     vjust = -1,
#     color = 'black', 
#     size = 2.8
#   ) +
#   scale_color_manual(values = custom_colors$discrete) +
#   scale_fill_manual(values = custom_colors$discrete) +
#   scale_y_continuous(
#     name = 'Assignment score',
#     labels = scales::comma,
#     limits = c(0, 1)
#   ) +
#   labs(title = 'Cell types') +
#   theme_bw() +
#   theme(
#     axis.title.x = element_blank(),
#     panel.grid.major.x = element_blank(),
#     axis.text.x = element_text(angle = 45, hjust = 1)  # Angled labels for readability
#   )
# 
# ggsave(
#   file.path(working_dir, 'singler_scores_s-g.png'),
#   p1 + p2 + p3 + plot_layout(ncol = 1), height = 13, width = 16
# )

# table_samples_by_cell_type <- test@meta.data %>%
#   dplyr::group_by(seurat_clusters, cell_type_singler_blueprintencode_main) %>%
#   dplyr::summarize(count = n()) %>%
#   tidyr::spread(cell_type_singler_blueprintencode_main, count, fill = 0) %>%
#   dplyr::ungroup() %>%
#   dplyr::mutate(total_cell_count = rowSums(.[c(2:ncol(.))])) %>%
#   dplyr::select(c("seurat_clusters", "total_cell_count", dplyr::everything()))
# 
# # Prepare the data for the heatmap
# heatmap_data <- table_samples_by_cell_type[, -which(names(table_samples_by_cell_type) == "total_cell_count")]
# 
# # Reshape to long format
# heatmap_data_long <- melt(heatmap_data, id.vars = "seurat_clusters", variable.name = "cell_type", value.name = "count")
# 
# # Normalize by row max
# heatmap_data_long <- heatmap_data_long %>%
#   group_by(seurat_clusters) %>%
#   mutate(count = count / sum(count)) %>%
#   ungroup()
# 
# # Create a matrix for clustering
# heatmap_matrix <- acast(heatmap_data_long, seurat_clusters ~ cell_type, value.var = "count")
# heatmap_matrix <- sweep(heatmap_matrix, 1, rowSums(heatmap_matrix), FUN = "/")
# 
# # Perform hierarchical clustering
# dist_matrix <- dist(t(heatmap_matrix))  # Distance on columns
# hc <- hclust(dist_matrix)
# 
# # Reorder columns based on clustering
# ordered_cols <- heatmap_matrix[, hc$order]
# 
# # Create a new long format data frame for the reordered heatmap
# heatmap_data_ordered <- melt(ordered_cols)
# colnames(heatmap_data_ordered) <- c("seurat_clusters", "cell_type", "likelihood")
# 
# # Convert seurat_clusters to numeric for proper ordering
# heatmap_data_ordered <- heatmap_data_ordered %>%
#   mutate(seurat_clusters = as.numeric(as.character(seurat_clusters)))
# 
# # Create the heatmap with ordered columns
# p <- ggplot(heatmap_data_ordered, aes(x = cell_type, y = factor(seurat_clusters), fill = likelihood)) +  # Use factor for proper ordering
#   geom_tile(color = "white") +
#   scale_fill_gradient(low = "white", high = "blue") +
#   theme_minimal() +
#   labs(title = "Per-cluster Cell Types Prediction by SingleR",
#        x = "Cell Type",
#        y = "Seurat Clusters") +
#   scale_y_discrete(expand = c(0, 0)) +  # Ensure y-axis ticks for each row
#   theme(axis.text.x = element_text(angle = 45, hjust = 1),
#         axis.ticks.x = element_line(),
#         axis.ticks.y = element_line(),
#         panel.background = element_rect(fill = "white"),  
#         plot.background = element_rect(fill = "white")
#   )  
# 
# ggsave(
#   file.path(working_dir, 'singler_per_cluster_heatmap_s-g.png'),
#   p,
#   height = 6,  
#   width = 8
# )
# 
# # Create a matrix for clustering
# heatmap_matrix <- acast(heatmap_data_long, seurat_clusters ~ cell_type, value.var = "count")
# heatmap_matrix <- sweep(heatmap_matrix, 2, colSums(heatmap_matrix), FUN = "/")
# 
# # Perform hierarchical clustering
# dist_matrix <- dist(t(heatmap_matrix))  # Distance on columns
# hc <- hclust(dist_matrix)
# 
# # Reorder columns based on clustering
# ordered_cols <- heatmap_matrix[, hc$order]
# 
# # Create a new long format data frame for the reordered heatmap
# heatmap_data_ordered <- melt(ordered_cols)
# colnames(heatmap_data_ordered) <- c("seurat_clusters", "cell_type", "likelihood")
# 
# # Convert seurat_clusters to numeric for proper ordering
# heatmap_data_ordered <- heatmap_data_ordered %>%
#   mutate(seurat_clusters = as.numeric(as.character(seurat_clusters)))
# 
# # Create the heatmap with ordered columns
# p <- ggplot(heatmap_data_ordered, aes(x = cell_type, y = factor(seurat_clusters), fill = likelihood)) +  # Use factor for proper ordering
#   geom_tile(color = "white") +
#   scale_fill_gradient(low = "white", high = "blue") +
#   theme_minimal() +
#   labs(title = "Per-celltype Cell Types Prediction by SingleR",
#        x = "Cell Type",
#        y = "Seurat Clusters") +
#   scale_y_discrete(expand = c(0, 0)) +  # Ensure y-axis ticks for each row
#   theme(axis.text.x = element_text(angle = 45, hjust = 1),
#         axis.ticks.x = element_line(),
#         axis.ticks.y = element_line(),
#         panel.background = element_rect(fill = "white"),  
#         plot.background = element_rect(fill = "white")
#   )  
# 
# ggsave(
#   file.path(working_dir, 'singler_per_celltype_heatmap_s-g.png'),
#   p,
#   height = 6,  
#   width = 8
# )
# 
# # Prepare the data for the heatmap
# heatmap_data <- table_samples_by_cell_type[, -which(names(table_samples_by_cell_type) == "total_cell_count")]
# 
# # Reshape to long format
# heatmap_data_long <- melt(heatmap_data, id.vars = "seurat_clusters", variable.name = "cell_type", value.name = "count")
# 
# # Normalize by row total
# heatmap_data_long <- heatmap_data_long %>%
#   group_by(seurat_clusters) %>%
#   mutate(count = count / sum(count)) %>%  # Normalize by total for each cluster
#   ungroup()
# 
# # Create a matrix for clustering
# heatmap_matrix <- acast(heatmap_data_long, seurat_clusters ~ cell_type, value.var = "count")
# 
# # Perform hierarchical clustering
# dist_matrix <- dist(t(heatmap_matrix))  # Distance on columns
# hc <- hclust(dist_matrix)
# 
# # Reorder columns based on clustering
# ordered_cols <- heatmap_matrix[, hc$order]
# 
# # Create a new long format data frame for the reordered heatmap
# heatmap_data_ordered <- melt(ordered_cols)
# colnames(heatmap_data_ordered) <- c("seurat_clusters", "cell_type", "likelihood")
# 
# # Convert seurat_clusters to numeric for proper ordering
# heatmap_data_ordered <- heatmap_data_ordered %>%
#   mutate(seurat_clusters = as.numeric(as.character(seurat_clusters)))
# 
# # Create the heatmap with ordered columns
# p <- ggplot(heatmap_data_ordered, aes(x = cell_type, y = factor(seurat_clusters), fill = likelihood)) +
#   geom_tile(color = "white") +
#   scale_fill_gradient(low = "white", high = "blue") +
#   theme_minimal() +
#   labs(title = "Per-cluster Cell Types Prediction by SingleR",
#        x = "Cell Type",
#        y = "Seurat Clusters") +
#   scale_y_discrete(expand = c(0, 0)) +  # Ensure y-axis ticks for each row
#   theme(axis.text.x = element_text(angle = 45, hjust = 1),
#         axis.ticks.x = element_line(),
#         axis.ticks.y = element_line(),
#         panel.background = element_rect(fill = "white"),  
#         plot.background = element_rect(fill = "white")
#   )  
# 
# # Save the heatmap
# ggsave(
#   file.path(working_dir, 'singler_heatmap.png'),
#   p,
#   height = 6,  
#   width = 8
# )
# 
# table_samples_by_cell_type %>% knitr::kable()
# 
# # |sample | total_cell_count| Astrocyte| B_cell| Chondrocytes| DC| Endothelial_cells| Epithelial_cells| Fibroblasts| Gametocytes| Hepatocytes| HSC_-G-CSF| Keratinocytes| Macrophage| Monocyte| Neuroepithelial_cell| Neurons| Neutrophils| NK_cell| Osteoblasts| Platelets| Smooth_muscle_cells| T_cells| Tissue_stem_cells|
# # |:------|----------------:|---------:|------:|------------:|--:|-----------------:|----------------:|-----------:|-----------:|-----------:|----------:|-------------:|----------:|--------:|--------------------:|-------:|-----------:|-------:|-----------:|---------:|-------------------:|-------:|-----------------:|
# # |Brain  |              470|        14|      1|           10| 44|                 6|               22|           1|           9|           2|          1|            50|          3|       23|                    2|       6|         155|      25|          76|         1|                   6|       2|                11|
# 
# 
# table_clusters_by_cell_type <- test@meta.data %>%
#   dplyr::group_by(seurat_clusters, cell_type_singler_blueprintencode_main) %>%
#   dplyr::rename(cluster = seurat_clusters) %>%
#   dplyr::summarize(count = n()) %>%
#   tidyr::spread(cell_type_singler_blueprintencode_main, count, fill = 0) %>%
#   dplyr::ungroup() %>%
#   dplyr::mutate(total_cell_count = rowSums(.[c(2:ncol(.))])) %>%
#   dplyr::select(c("cluster", "total_cell_count", dplyr::everything()))
# 
# table_clusters_by_cell_type %>% knitr::kable()
# 
# # |cluster | total_cell_count| Astrocyte| B_cell| Chondrocytes| DC| Endothelial_cells| Epithelial_cells| Fibroblasts| Gametocytes| Hepatocytes| HSC_-G-CSF| Keratinocytes| Macrophage| Monocyte| Neuroepithelial_cell| Neurons| Neutrophils| NK_cell| Osteoblasts| Platelets| Smooth_muscle_cells| T_cells| Tissue_stem_cells|
# # |:-------|----------------:|---------:|------:|------------:|--:|-----------------:|----------------:|-----------:|-----------:|-----------:|----------:|-------------:|----------:|--------:|--------------------:|-------:|-----------:|-------:|-----------:|---------:|-------------------:|-------:|-----------------:|
# # |0       |              154|         8|      0|            7|  3|                 1|               18|           1|           6|           1|          0|            12|          1|        4|                    2|       0|          40|       5|          35|         0|                   2|       0|                 8|
# # |1       |              119|         0|      1|            1| 22|                 3|                3|           0|           2|           1|          0|            11|          1|       10|                    0|       0|          36|      14|           9|         1|                   1|       2|                 1|
# # |2       |               75|         5|      0|            1|  5|                 0|                0|           0|           0|           0|          0|             0|          0|        1|                    0|       6|          25|       0|          28|         0|                   3|       0|                 1|
# # |3       |               62|         1|      0|            1|  5|                 1|                0|           0|           0|           0|          0|            19|          0|        4|                    0|       0|          23|       5|           2|         0|                   0|       0|                 1|
# # |4       |               36|         0|      0|            0|  5|                 1|                0|           0|           0|           0|          0|             0|          0|        2|                    0|       0|          26|       1|           1|         0|                   0|       0|                 0|
# # |5       |               24|         0|      0|            0|  4|                 0|                1|           0|           1|           0|          1|             8|          1|        2|                    0|       0|           5|       0|           1|         0|                   0|       0|                 0|
# # 
# 
# UMAP_centers_cluster <- tibble(
#   UMAP_1 = as.data.frame(test@reductions$UMAP@cell.embeddings)$UMAP_1,
#   UMAP_2 = as.data.frame(test@reductions$UMAP@cell.embeddings)$UMAP_2,
#   cluster = test@meta.data$cell_type_singler_blueprintencode_main
# ) %>%
#   group_by(cluster) %>%
#   dplyr::summarize(x = median(UMAP_1), y = median(UMAP_2))
# 
# # Create the UMAP plot with repelled labels
# plot_umap_by_cluster <- bind_cols(test@meta.data, as.data.frame(test@reductions$UMAP@cell.embeddings)) %>%
#   ggplot(aes(UMAP_1, UMAP_2, color = seurat_clusters)) +
#   geom_point(size = 0.2) +
#   geom_label_repel(  # Use geom_label_repel instead of geom_label
#     data = UMAP_centers_cluster,
#     mapping = aes(x, y, label = cluster),
#     size = 4.5,
#     fill = 'white',
#     color = 'black',
#     fontface = 'bold',
#     alpha = 0.7,
#     label.size = 0,
#     show.legend = FALSE,
#     box.padding = 0.5,  # Adjust padding around labels
#     point.padding = 0.5,  # Adjust padding between labels and points
#     segment.color = NA  # Color of the connecting line
#   ) +
#   theme_bw() +
#   scale_color_manual(
#     name = 'Cluster', values = custom_colors$discrete,
#     guide = guide_legend(override.aes = list(size = 2))
#   ) +
#   theme(legend.position = 'right') +
#   coord_fixed() +
#   annotate(
#     geom = 'text', x = Inf, y = -Inf,
#     label = paste0('n = ', format(nrow(test@meta.data), big.mark = ',', trim = TRUE)),
#     vjust = -1.5, hjust = 1.25, color = 'black', size = 2.5
#   )
# 
# tSNE_centers_cluster <- tibble(
#   tSNE_1 = as.data.frame(test@reductions$tSNE@cell.embeddings)$tSNE_1,
#   tSNE_2 = as.data.frame(test@reductions$tSNE@cell.embeddings)$tSNE_2,
#   cluster = test@meta.data$cell_type_singler_blueprintencode_main
# ) %>%
#   group_by(cluster) %>%
#   dplyr::summarize(x = median(tSNE_1), y = median(tSNE_2))
# 
# plot_tsne_by_cluster <- bind_cols(test@meta.data, as.data.frame(test@reductions$tSNE@cell.embeddings)) %>%
#   ggplot(aes(tSNE_1, tSNE_2, color = seurat_clusters)) +
#   geom_point(size = 0.2) +
#   geom_label_repel(  # Use geom_label_repel instead of geom_label
#     data = tSNE_centers_cluster,
#     mapping = aes(x, y, label = cluster),
#     size = 4.5,
#     fill = 'white',
#     color = 'black',
#     fontface = 'bold',
#     alpha = 0.7,
#     label.size = 0,
#     show.legend = FALSE,
#     box.padding = 0.5,  # Adjust padding around labels
#     point.padding = 0.5,  # Adjust padding between labels and points
#     segment.color = NA  # Color of the connecting line
#   ) +
#   theme_bw() +
#   scale_color_manual(
#     name = 'Cluster', values = custom_colors$discrete,
#     guide = guide_legend(override.aes = list(size = 2))
#   ) +
#   theme(legend.position = 'right') +
#   coord_fixed() +
#   annotate(
#     geom = 'text', x = Inf, y = -Inf,
#     label = paste0('n = ', format(nrow(test@meta.data), big.mark = ',', trim = TRUE)),
#     vjust = -1.5, hjust = 1.25, color = 'black', size = 2.5
#   )
# 
# ggsave(
#   file.path(working_dir, 'singler_annotation_main.brain.png'),
#   plot_umap_by_cluster + plot_tsne_by_cluster,
#   height = 6,  
#   width = 12
# )

### SNN graph
RNA_snn <- seurat@graphs$RNA_snn %>%
  as.matrix() %>%
  ggnetwork() %>%
  left_join(seurat@meta.data %>% mutate(vertex.names = rownames(.)), by = 'vertex.names')

# Create the second plot colored by 'seurat_clusters'
p1 <- ggplot(RNA_snn, aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_edges(color = 'grey50', alpha = 0.05) +
  geom_nodes(aes(color = seurat_clusters), size = 0.5) +
  scale_color_manual(
    name = 'Cluster', values = custom_colors$discrete,
    guide = guide_legend(ncol = 1, override.aes = list(size = 2))
  ) +
  theme_void() +
  theme(legend.position = 'right') +
  annotate(
    geom = 'text', x = Inf, y = -Inf,
    label = paste0('n = ', format(nrow(seurat@meta.data), big.mark = ',', trim = TRUE)),
    vjust = -1.5, hjust = 1.25, color = 'black', size = 2.5
  )

p2 <- ggplot(RNA_snn, aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_edges(color = 'grey50', alpha = 0.05) +
  geom_nodes(aes(color = Sample), size = 0.5) +
  scale_color_manual(
    name = 'Sample', values = custom_colors$discrete,
    guide = guide_legend(ncol = 1, override.aes = list(size = 2))
  ) +
  theme_blank() +
  theme(legend.position = 'left') +
  annotate(
    geom = 'text', x = Inf, y = -Inf,
    label = paste0('n = ', format(nrow(seurat@meta.data), big.mark = ',', trim = TRUE)),
    vjust = -1.5, hjust = 1.25, color = 'black', size = 2.5
  )

# Save the combined plot
# ggsave(
#  file.path(working_dir, 'snn_graph_by_sample_cluster.png'),
#   p1 + p2 + plot_layout(ncol = 2),
#   height = 5, width = 11
# )

ggsave(
  file.path(working_dir, 'snn.png'),
  p1 + p2 + plot_layout(ncol = 2),
  height = 5, width = 11
)

### UMAP
{
  seurat <- RunUMAP(
    seurat,
    reduction.name = 'UMAP',
    reduction = 'pca',
    dims = 1:pca_cutoff,
    n.components = 2,
    seed.use = 100
  )
  
  ## all UMAP
  plot_umap_by_nCount <- bind_cols(seurat@meta.data, as.data.frame(seurat@reductions$UMAP@cell.embeddings)) %>%
    ggplot(aes(UMAP_1, UMAP_2, color = 10 * log(1+nCount_RNA))) +
    geom_point(size = 0.2) +
    theme_bw() +
    scale_color_viridis(
      guide = guide_colorbar(frame.colour = 'black', ticks.colour = 'black'),
      labels = scales::comma,
    ) +
    labs(color = 'Number of\ntranscripts') +
    theme(legend.position = 'right') +
    coord_fixed() +
    annotate(
      geom = 'text', x = Inf, y = -Inf,
      label = paste0('n = ', format(nrow(seurat@meta.data), big.mark = ',', trim = TRUE)),
      vjust = -1.5, hjust = 1.25, color = 'black', size = 2.5
    )
  
  plot_umap_by_nFeature <- bind_cols(seurat@meta.data, as.data.frame(seurat@reductions$UMAP@cell.embeddings)) %>%
    ggplot(aes(UMAP_1, UMAP_2, color = 10 * log(1+nFeature_RNA))) +
    geom_point(size = 0.2) +
    theme_bw() +
    scale_color_viridis(
      guide = guide_colorbar(frame.colour = 'black', ticks.colour = 'black'),
      labels = scales::comma,
    ) +
    labs(color = 'Number of\ngenes') +
    theme(legend.position = 'right') +
    coord_fixed() +
    annotate(
      geom = 'text', x = Inf, y = -Inf,
      label = paste0('n = ', format(nrow(seurat@meta.data), big.mark = ',', trim = TRUE)),
      vjust = -1.5, hjust = 1.25, color = 'black', size = 2.5
    )
  
  UMAP_centers_cluster <- tibble(
    UMAP_1 = as.data.frame(seurat@reductions$UMAP@cell.embeddings)$UMAP_1,
    UMAP_2 = as.data.frame(seurat@reductions$UMAP@cell.embeddings)$UMAP_2,
    cluster = seurat@meta.data$seurat_clusters
  ) %>%
    group_by(cluster) %>%
    dplyr::summarize(x = median(UMAP_1), y = median(UMAP_2))
  
  plot_umap_by_cluster <- bind_cols(seurat@meta.data, as.data.frame(seurat@reductions$UMAP@cell.embeddings)) %>%
    ggplot(aes(UMAP_1, UMAP_2, color = seurat_clusters)) +
    geom_point(size = 0.2) +
    geom_label(
      data = UMAP_centers_cluster,
      mapping = aes(x, y, label = cluster),
      size = 4.5,
      fill = 'white',
      color = 'black',
      fontface = 'bold',
      alpha = 0.5,
      label.size = 0,
      show.legend = FALSE
    ) +
    theme_bw() +
    scale_color_manual(
      name = 'Cluster', values = custom_colors$discrete,
      guide = guide_legend(override.aes = list(size = 2))
    ) +
    theme(legend.position = 'right') +
    coord_fixed() +
    annotate(
      geom = 'text', x = Inf, y = -Inf,
      label = paste0('n = ', format(nrow(seurat@meta.data), big.mark = ',', trim = TRUE)),
      vjust = -1.5, hjust = 1.25, color = 'black', size = 2.5
    )
  
  plot_umap_by_sample <- bind_cols(seurat@meta.data, as.data.frame(seurat@reductions$UMAP@cell.embeddings)) %>%
    ggplot(aes(UMAP_1, UMAP_2, color = merged_sample)) +
    geom_point(size = 0.2) +
    theme_bw() +
    scale_color_manual(
      name = 'Sample', values = custom_colors$discrete[c(1, 4, 5, 2, 3, 14)],
      guide = guide_legend(ncol = 1, override.aes = list(size = 2))
    ) +
    theme(legend.position = 'right') +
    coord_fixed() +
    annotate(
      geom = 'text', x = Inf, y = -Inf,
      label = paste0('n = ', format(nrow(seurat@meta.data), big.mark = ',', trim = TRUE)),
      vjust = -1.5, hjust = 1.25, color = 'black', size = 2.5
    )
  
  ggsave(
    # file.path(working_dir, 'UMAP.single.png'),
    # file.path(working_dir, 'UMAP.5.singlet.s+g.png'),
    # file.path(working_dir, 'UMAP.5.singlet.s-g.int.png'),
    file.path(working_dir, paste0('UMAP.', runid, '.pdf')),
    plot_umap_by_nCount + plot_umap_by_cluster + 
      plot_umap_by_nFeature + plot_umap_by_sample + plot_layout(ncol = 2),
    height = 12,
    width = 12
  )
  
  # UMAP by samples
  temp_labels <- seurat@meta.data %>%
    group_by(merged_sample) %>%
    tally()
  
  ncount_umap_plot <- bind_cols(seurat@meta.data, as.data.frame(seurat@reductions$UMAP@cell.embeddings)) %>%
    ggplot(aes(UMAP_1, UMAP_2, color = 10 * log(1+nCount_RNA))) +
    geom_point(size = 0.2) +
    geom_text(
      data = temp_labels,
      aes(x = Inf, y = -Inf, label = paste0('n = ', format(n, big.mark = ',', trim = TRUE)), vjust = -1.5, hjust = 1.25),
      color = 'black', size = 2.8
    ) +
    theme_bw() +
    scale_color_viridis(
      guide = guide_colorbar(frame.colour = 'black', ticks.colour = 'black'),
      labels = scales::comma,
    ) +
    labs(color = 'Number of\ntranscripts') +
    theme(legend.position = 'right') +
    coord_fixed() +
    facet_wrap(~merged_sample, ncol = 6)
  
  nfeature_umap_plot <- bind_cols(seurat@meta.data, as.data.frame(seurat@reductions$UMAP@cell.embeddings)) %>%
    ggplot(aes(UMAP_1, UMAP_2, color = 10 * log(1+nFeature_RNA))) +
    geom_text(
      data = temp_labels,
      aes(x = Inf, y = -Inf, label = paste0('n = ', format(n, big.mark = ',', trim = TRUE)), vjust = -1.5, hjust = 1.25),
      color = 'black', size = 2.8
    ) +
    geom_point(size = 0.2) +
    theme_bw() +
    scale_color_viridis(
      guide = guide_colorbar(frame.colour = 'black', ticks.colour = 'black'),
      labels = scales::comma,
    ) +
    labs(color = 'Number of\ngenes') +
    theme(legend.position = 'right') +
    coord_fixed() +
    facet_wrap(~merged_sample, ncol = 6)  
  
  umap_plot <- bind_cols(seurat@meta.data, as.data.frame(seurat@reductions$UMAP@cell.embeddings)) %>%
    ggplot(aes(UMAP_1, UMAP_2, color = merged_sample)) +
    geom_point(size = 0.2, show.legend = FALSE) +
    geom_text(
      data = temp_labels,
      aes(x = Inf, y = -Inf, label = paste0('n = ', format(n, big.mark = ',', trim = TRUE)), vjust = -1.5, hjust = 1.25),
      color = 'black', size = 2.8
    ) +
    theme_bw() +
    scale_color_manual(values = custom_colors$discrete) +
    coord_fixed() +
    theme(
      legend.position = 'none',
      strip.text = element_text(face = 'bold', size = 10)
    ) +
    facet_wrap(~merged_sample, ncol = 6)  # 6 columns for UMAP
  
  ggsave(
    # file.path(working_dir, 'UMAP_by_sample.single.png'),
    # file.path(working_dir, 'UMAP_by_sample.5.singlet.s+g.png'),
    file.path(working_dir, paste0('UMAP_by_sample.', runid, '.pdf')),
    ncount_umap_plot / nfeature_umap_plot / umap_plot,
    height = 12,  
    width = 18
  )
}

### tSNE
{
  seurat <- RunTSNE(
    seurat,
    reduction.name = 'tSNE',
    reduction = 'pca',
    dims = 1:pca_cutoff,
    seed.use = 100
  )
  
  ## all tSNE
  plot_tsne_by_nCount <- bind_cols(seurat@meta.data, as.data.frame(seurat@reductions$tSNE@cell.embeddings)) %>%
    ggplot(aes(tSNE_1, tSNE_2, color = 10 * log(1 + nCount_RNA))) +
    geom_point(size = 0.2) +
    theme_bw() +
    scale_color_viridis(
      guide = guide_colorbar(frame.colour = 'black', ticks.colour = 'black'),
      labels = scales::comma,
    ) +
    labs(color = 'Number of\ntranscripts') +
    theme(legend.position = 'right') +
    coord_fixed() +
    annotate(
      geom = 'text', x = Inf, y = -Inf,
      label = paste0('n = ', format(nrow(seurat@meta.data), big.mark = ',', trim = TRUE)),
      vjust = -1.5, hjust = 1.25, color = 'black', size = 2.5
    )
  
  plot_tsne_by_nFeature <- bind_cols(seurat@meta.data, as.data.frame(seurat@reductions$tSNE@cell.embeddings)) %>%
    ggplot(aes(tSNE_1, tSNE_2, color = 10 * log(1 + nFeature_RNA))) +
    geom_point(size = 0.2) +
    theme_bw() +
    scale_color_viridis(
      guide = guide_colorbar(frame.colour = 'black', ticks.colour = 'black'),
      labels = scales::comma,
    ) +
    labs(color = 'Number of\ngenes') +
    theme(legend.position = 'right') +
    coord_fixed() +
    annotate(
      geom = 'text', x = Inf, y = -Inf,
      label = paste0('n = ', format(nrow(seurat@meta.data), big.mark = ',', trim = TRUE)),
      vjust = -1.5, hjust = 1.25, color = 'black', size = 2.5
    )
  
  tSNE_centers_cluster <- tibble(
    tSNE_1 = as.data.frame(seurat@reductions$tSNE@cell.embeddings)$tSNE_1,
    tSNE_2 = as.data.frame(seurat@reductions$tSNE@cell.embeddings)$tSNE_2,
    cluster = seurat@meta.data$seurat_clusters
  ) %>%
    group_by(cluster) %>%
    dplyr::summarize(x = median(tSNE_1), y = median(tSNE_2))
  
  plot_tsne_by_cluster <- bind_cols(seurat@meta.data, as.data.frame(seurat@reductions$tSNE@cell.embeddings)) %>%
    ggplot(aes(tSNE_1, tSNE_2, color = seurat_clusters)) +
    geom_point(size = 0.2) +
    geom_label(
      data = tSNE_centers_cluster,
      mapping = aes(x, y, label = cluster),
      size = 4.5,
      fill = 'white',
      color = 'black',
      fontface = 'bold',
      alpha = 0.5,
      label.size = 0,
      show.legend = FALSE
    ) +
    theme_bw() +
    scale_color_manual(
      name = 'Cluster', values = custom_colors$discrete,
      guide = guide_legend(override.aes = list(size = 2))
    ) +
    theme(legend.position = 'right') +
    coord_fixed() +
    annotate(
      geom = 'text', x = Inf, y = -Inf,
      label = paste0('n = ', format(nrow(seurat@meta.data), big.mark = ',', trim = TRUE)),
      vjust = -1.5, hjust = 1.25, color = 'black', size = 2.5
    )
  
  plot_tsne_by_sample <- bind_cols(seurat@meta.data, as.data.frame(seurat@reductions$tSNE@cell.embeddings)) %>%
    ggplot(aes(tSNE_1, tSNE_2, color = merged_sample)) +
    geom_point(size = 0.2) +
    theme_bw() +
    scale_color_manual(
      name = 'Sample', values = custom_colors$discrete[c(1, 4, 5, 2, 3, 14)],
      guide = guide_legend(ncol = 1, override.aes = list(size = 2))
    ) +
    theme(legend.position = 'right') +
    coord_fixed() +
    annotate(
      geom = 'text', x = Inf, y = -Inf,
      label = paste0('n = ', format(nrow(seurat@meta.data), big.mark = ',', trim = TRUE)),
      vjust = -1.5, hjust = 1.25, color = 'black', size = 2.5
    )
  
  ggsave(
    # file.path(working_dir, 'tSNE.single.png'),
    # file.path(working_dir, 'tSNE.5.singlet.s+g.png'),
    # file.path(working_dir, 'tSNE.5.singlet.s-g.int.png'),
    file.path(working_dir, paste0('tSNE.', runid, '.pdf')),
    plot_tsne_by_nCount + plot_tsne_by_cluster + 
      plot_tsne_by_nFeature + plot_tsne_by_sample + plot_layout(ncol = 2),
    height = 12,
    width = 12
  )
  
  # tSNE by samples
  temp_labels <- seurat@meta.data %>%
    group_by(merged_sample) %>%
    tally()
  
  ncount_tsne_plot <- bind_cols(seurat@meta.data, as.data.frame(seurat@reductions$tSNE@cell.embeddings)) %>%
    ggplot(aes(tSNE_1, tSNE_2, color = 10 * log(1 + nCount_RNA))) +
    geom_point(size = 0.2) +
    geom_text(
      data = temp_labels,
      aes(x = Inf, y = -Inf, label = paste0('n = ', format(n, big.mark = ',', trim = TRUE)), vjust = -1.5, hjust = 1.25),
      color = 'black', size = 2.8
    ) +
    theme_bw() +
    scale_color_viridis(
      guide = guide_colorbar(frame.colour = 'black', ticks.colour = 'black'),
      labels = scales::comma,
    ) +
    labs(color = 'Number of\ntranscripts') +
    theme(legend.position = 'right') +
    coord_fixed() +
    facet_wrap(~merged_sample, ncol = 6)
  
  nfeature_tsne_plot <- bind_cols(seurat@meta.data, as.data.frame(seurat@reductions$tSNE@cell.embeddings)) %>%
    ggplot(aes(tSNE_1, tSNE_2, color = 10 * log(1 + nFeature_RNA))) +
    geom_text(
      data = temp_labels,
      aes(x = Inf, y = -Inf, label = paste0('n = ', format(n, big.mark = ',', trim = TRUE)), vjust = -1.5, hjust = 1.25),
      color = 'black', size = 2.8
    ) +
    geom_point(size = 0.2) +
    theme_bw() +
    scale_color_viridis(
      guide = guide_colorbar(frame.colour = 'black', ticks.colour = 'black'),
      labels = scales::comma,
    ) +
    labs(color = 'Number of\ngenes') +
    theme(legend.position = 'right') +
    coord_fixed() +
    facet_wrap(~merged_sample, ncol = 6)
  
  tsne_plot <- bind_cols(seurat@meta.data, as.data.frame(seurat@reductions$tSNE@cell.embeddings)) %>%
    ggplot(aes(tSNE_1, tSNE_2, color = merged_sample)) +
    geom_point(size = 0.2, show.legend = FALSE) +
    geom_text(
      data = temp_labels,
      aes(x = Inf, y = -Inf, label = paste0('n = ', format(n, big.mark = ',', trim = TRUE)), vjust = -1.5, hjust = 1.25),
      color = 'black', size = 2.8
    ) +
    theme_bw() +
    scale_color_manual(values = custom_colors$discrete) +
    coord_fixed() +
    theme(
      legend.position = 'none',
      strip.text = element_text(face = 'bold', size = 10)
    ) +
    facet_wrap(~merged_sample, ncol = 6)  # 6 columns for tSNE
  
  ggsave(
    # file.path(working_dir, 'tSNE_by_sample.single.png'),
    # file.path(working_dir, 'tSNE_by_sample.5.singlet.s-g.int.png'),
    file.path(working_dir, paste0('tSNE_by_sample.', runid, '.pdf')),
    ncount_tsne_plot / nfeature_tsne_plot / tsne_plot,
    height = 12,  
    width = 18
  )
}

seurat <- readRDS(file.path(working_dir, "processed.single.rds"))
# saveRDS(seurat, file.path(working_dir, "processed.single.rds"))

seurat <- readRDS(file.path(working_dir, "processed.gfp.single.s-g.rds"))
# saveRDS(seurat, file.path(working_dir, "processed.gfp.single.s-g.rds"))

seurat <- readRDS(file.path(working_dir, "processed.gfp.single.rds"))
# saveRDS(seurat, file.path(working_dir, "processed.gfp.single.rds"))

seurat <- readRDS(file.path(working_dir, "processed.gfp.single.nocc.rds"))
# saveRDS(seurat, file.path(working_dir, "processed.gfp.single.nocc.rds"))

seurat <- readRDS(file.path(working_dir, "processed.5.singlet.s+g.rds"))
# saveRDS(seurat, file.path(working_dir, "processed.5.singlet.s+g.rds"))

seurat <- readRDS(file.path(working_dir, "processed.5.singlet.s-g.int.rds"))
# saveRDS(seurat, file.path(working_dir, "processed.5.singlet.s-g.int.rds"))

### Dim plots
plots <- list()

for (i in 1:pca_cutoff) {
  # Create the heatmap for each dimension
  heatmap_plot <- DimHeatmap(seurat, dims = i, nfeatures = 30, fast = FALSE) + 
    scale_fill_viridis(option = "B", name = 'Expression') +
    ggtitle(paste("PC", i, sep = "_")) +  # Use paste to create the title
    theme(plot.title = element_text(hjust = 0.5, size = 14))
  
  # Remove legend for all but the last plot
  if (i < pca_cutoff) {
    heatmap_plot <- heatmap_plot + NoLegend()
  }
  
  # Store the plot in the list
  plots[[i]] <- heatmap_plot
}

ggsave(file.path(working_dir, 'dimheatmap.single.png'), wrap_plots(plots, nrow = 2, ncol = 5), height = 8, width = 16)

### FeaturePlots
rownames(seurat[['RNA']])[grep("(?i)flk", rownames(seurat[['RNA']]))]
rownames(seurat[['integrated']])[grep("(?i)Myc", rownames(seurat[['integrated']]))]

sapply(rownames(seurat[['RNA']])[grep("MYO", rownames(seurat[['RNA']]))], function(g) mean(seurat[["RNA"]]$counts[g,][seurat[["RNA"]]$counts[g,] > 0]))

# individual violins
{
  # naming
  query <- "4" # lowercases will suffice
  DefaultAssay(object = seurat) <- "RNA"
  runid <- "0718"
  
  # query 
  # genes <- rownames(seurat)[grep(paste0("(?i)", query), rownames(seurat))]
  
  # hard-coding
  # genes <- c("CD53", "FLT3.3", "ALPL.4", "DDX56", "MYO9A", "EFHD2", "Sox14/15/21 (SOX14)", "FABP2")
  # genes <- c("CD53", "FLT3.3", "ALPL.4", "DDX56")
  # genes <- c("Hes. (HES1)", "Hes. (HES1).2", "Hes. (HES1).3")
  # genes <- c("KIT", "CD34")
  # genes <- c(
  #   "Pou-r (POU6F1)", "DDX56", "MEX3C", "Myc (MYCN)", "ALPL.4", "NOP56", "NOP58", "HSPE1", # 0
  #   "CD53", "NFIA", "SPON1", "SLC12A3", "SHANK2", # 1
  #   "NRG1", "THBD", "CD34", # 2
  #   "SH3D19", "SEMA3E", "EFHD2", # 3
  #   "EYA1", "FBLN7", # 4
  #   "FLT3.3", "CD163", "CLBA1", # 5
  #   # 6
  #   "CYBB", "SHISAL2B", # 8
  #   "PIGK", "Serum (SRF)", "SPHK1", # 9
  #   "CTSA", "MALRD1", "ITGA2", # 10
  #   # 11
  #   # 12
  #   "CCDC18", "GBP2", "MIF", "LGALS4.2", "CTSL", "VIL1", # 13
  #   "ACSL1.2", # 14
  #   "FTCD", # 15
  #   "SYCP2", # 16
  #   "PTPN23", # 17
  #   "TMSB4X", # 18
  #   "SLC30A3", "SNAP25" # 19
  # )
  # genes <- c("NOP56", "NOP58", "HSPE1", "ITGA2")
  # genes <- c(
  #   "Pou-r (POU6F1)", "DDX56", "MEX3C", "Myc (MYCN)", "ALPL.4", "NOP56", "NOP58", "HSPE1", # 0
  #   "CD53", "NFIA", "SPON1", "SHANK2", # 1
  #   "NRG1", "THBD", "CD34", # 2
  #   "SEMA3E", "EFHD2", # 3
  #   "EYA1", # 4
  #   "FLT3.3", "CD163", # 5
  #   # 6
  #   "CYBB", "SHISAL2B", # 8
  #   "Serum (SRF)", "SPHK1", # 9
  #   "CTSA", "MALRD1", "ITGA2", # 10
  #   # 11
  #   # 12
  #   "GBP2", "MIF", "CTSL", "VIL1", # 13
  #   "ACSL1.2", # 14
  #   "PAX6", # 15
  #   # 16
  #   "PTPN23", # 17
  #   "TMSB4X", # 18
  #   "SLC30A3", "SNAP25" # 19
  # )
  # 0718 markers
  genes <- c(
    "Pou-r (POU6F1)", "DDX56", "MEX3C", "Myc (MYCN)", "ALPL.4", "NOP56", "NOP58", "HSPE1", # 0
    "CD53", "NFIA", "SPON1", "SHANK2", # 1
    "NRG1", "THBD", "CD34", # 2
    "SEMA3E", "EFHD2", # 3
    "EYA1", # 4
    "FLT3.3", "CD163", # 5
    # 6
    "CYBB", "SHISAL2B", # 8
    "Serum (SRF)", "SPHK1", # 9
    "CTSA", "MALRD1", "ITGA2", # 10
    # 11
    # 12
    "GBP2", "MIF", "CTSL", "VIL1", # 13
    "ACSL1.2", # 14
    # 15
    # 16
    "PTPN23", # 17
    "TMSB4X", # 18
    "SLC30A3", "SNAP25", # 19
    "PAX6", "Hes. (HES1).3", "CISH.2", "RHBDL2"
  )
  
  if (query == "3") {
    expr_data <- FetchData(seurat, vars = c(genes, "combined_sample"), layer = 'counts')
  } else if (query == "4") {
    expr_data <- FetchData(seurat, vars = c(genes, "merged_sample"), layer = 'counts')
  } 
  
  DefaultAssay(object = seurat) <- "integrated"
  
  expr_long <- expr_data %>%
    pivot_longer(
      cols = all_of(genes),
      names_to = "gene",
      values_to = "expr_value"
    ) 
  
  if (query == "3") {
    for (gene in genes) {
      gene_data <- expr_long %>% filter(gene == !!gene)
      
      plot <- ggstatsplot::ggbetweenstats(
        data = gene_data,
        x = combined_sample,
        y = expr_value,
        type = "nonparametric",
        pairwise.comparisons = TRUE,
        pairwise.display = "significant",
        centrality.plotting = TRUE,
        centrality.point.args = list(size = 0, color = "darkred"),
        centrality.label.args = list(alpha = 0),
        p.adjust.method  = "bonferroni",
        palette = "Set1",
        package = "RColorBrewer",
        violin.args = list(width = 0.8, alpha = 0.7),
        point.args = list(position = ggplot2::position_jitterdodge(jitter.width = 0.8), alpha = 0.4, size = 1.5, stroke = 0, na.rm = TRUE),
        title = gene,
        results.subtitle = FALSE,
        xlab = "",
        ylab = "Expression",
        ggplot.component = list(theme(axis.title.y.right = element_blank(), 
                                      axis.text.y.right = element_blank(), 
                                      axis.ticks.y.right = element_blank()))
      )
      
      # Clean gene name for filename
      clean_gene_name <- gsub(" ", "_", gene)
      clean_gene_name <- gsub("/", ".", clean_gene_name)
      clean_gene_name <- gsub("\\(", "_", clean_gene_name)
      clean_gene_name <- gsub("\\)", "_", clean_gene_name)
      clean_gene_name <- gsub("__", "_", clean_gene_name)
      clean_gene_name <- gsub("_\\.", ".", clean_gene_name)  
      clean_gene_name <- gsub("\\._", ".", clean_gene_name) 
      clean_gene_name <- gsub("_+$", "", clean_gene_name)
      
      ggsave(
        file.path(working_dir, paste0(runid, "/violin.", query, ".", clean_gene_name, ".pdf")), 
        plot, 
        width = 3, # 3 for 3 samples; 4 for 4 samples?
        height = 6, # 8 for more than 4 pairs of comparison; 6 is good for 5
        dpi = 300, 
        limitsize = FALSE
      )
    }
  } else if (query == "4") {
    for (gene in genes) {
      gene_data <- expr_long %>% filter(gene == !!gene)
      
      plot <- ggstatsplot::ggbetweenstats(
        data = gene_data,
        x = merged_sample,
        y = expr_value,
        type = "nonparametric",
        pairwise.comparisons = TRUE,
        pairwise.display = "significant",
        centrality.plotting = TRUE,
        centrality.point.args = list(size = 0, color = "darkred"),
        centrality.label.args = list(alpha = 0),
        p.adjust.method  = "bonferroni",
        palette = "Set1",
        package = "RColorBrewer",
        violin.args = list(width = 0.8, alpha = 0.7),
        point.args = list(position = ggplot2::position_jitterdodge(jitter.width = 0.8), alpha = 0.4, size = 1.5, stroke = 0, na.rm = TRUE),
        title = gene,
        results.subtitle = FALSE,
        xlab = "",
        ylab = "Expression",
        ggplot.component = list(theme(axis.title.y.right = element_blank(), 
                                      axis.text.y.right = element_blank(), 
                                      axis.ticks.y.right = element_blank()))
      )
      
      # Clean gene name for filename
      clean_gene_name <- gsub(" ", "_", gene)
      clean_gene_name <- gsub("/", ".", clean_gene_name)
      clean_gene_name <- gsub("\\(", "_", clean_gene_name)
      clean_gene_name <- gsub("\\)", "_", clean_gene_name)
      clean_gene_name <- gsub("__", "_", clean_gene_name)
      clean_gene_name <- gsub("_\\.", ".", clean_gene_name)  
      clean_gene_name <- gsub("\\._", ".", clean_gene_name) 
      clean_gene_name <- gsub("_+$", "", clean_gene_name)
      
      ggsave(
        file.path(working_dir, paste0(runid, "/violin.", query, ".", clean_gene_name, ".pdf")), 
        plot, 
        width = 4, # 3 for 3 samples; 4 for 4 samples?
        height = 6, # 8 for more than 4 pairs of comparison; 6 is good for 5?
        dpi = 300, 
        limitsize = FALSE
      )
    }
  }
}

# violins in a combined plot
{
  # naming
  query <- "0716_GFP" # lowercases will suffice
  DefaultAssay(object = seurat) <- "RNA"
  
  # query 
  # genes <- rownames(seurat)[grep(paste0("(?i)", query), rownames(seurat))]
  
  # hard-coding
  # genes <- c("CD53", "FLT3.3", "ALPL.4", "DDX56", "MYO9A", "EFHD2", "Sox14/15/21 (SOX14)", "FABP2")
  # genes <- c("CD53", "FLT3.3", "ALPL.4", "DDX56")
  # genes <- c("Hes. (HES1)", "Hes. (HES1).2", "Hes. (HES1).3")
  # genes <- c("KIT", "CD34")
  # genes <- c(
  #   "Pou-r (POU6F1)", "DDX56", "MEX3C", "Myc (MYCN)", "ALPL.4", "NOP56", "NOP58", "HSPE1", # 0
  #   "CD53", "NFIA", "SPON1", "SLC12A3", "SHANK2", # 1
  #   "NRG1", "THBD", "CD34", # 2
  #   "SH3D19", "SEMA3E", "EFHD2", # 3
  #   "EYA1", "FBLN7", # 4
  #   "FLT3.3", "CD163", "CLBA1", # 5
  #   # 6
  #   "CYBB", "SHISAL2B", # 8
  #   "PIGK", "Serum (SRF)", "SPHK1", # 9
  #   "CTSA", "MALRD1", "ITGA2", # 10
  #   # 11
  #   # 12
  #   "CCDC18", "GBP2", "MIF", "LGALS4.2", "CTSL", "VIL1", # 13
  #   "ACSL1.2", # 14
  #   "FTCD", # 15
  #   "SYCP2", # 16
  #   "PTPN23", # 17
  #   "TMSB4X", # 18
  #   "SLC30A3", "SNAP25" # 19
  # )
  genes <- markers <- c(
    "Pou-r (POU6F1)", "DDX56", "MEX3C", "Myc (MYCN)", "ALPL.4", "NOP56", "NOP58", "HSPE1", # 0
    "CD53", "NFIA", "SPON1", "SHANK2", # 1
    "NRG1", "THBD", "CD34", # 2
    "SEMA3E", "EFHD2", # 3
    "EYA1", # 4
    "FLT3.3", "CD163", # 5
    # 6
    "CYBB", "SHISAL2B", # 8
    "Serum (SRF)", "SPHK1", # 9
    "CTSA", "MALRD1", "ITGA2", # 10
    # 11
    # 12
    "GBP2", "MIF", "CTSL", "VIL1", # 13
    "ACSL1.2", # 14
    "PAX6", # 15
    # 16
    "PTPN23", # 17
    "TMSB4X", # 18
    "SLC30A3", "SNAP25" # 19
  )
  
  expr_data <- FetchData(seurat, vars = c(genes, "combined_sample"), layer = 'counts')
  # expr_data <- FetchData(seurat, vars = c(genes, "merged_sample"), layer = 'counts')
  DefaultAssay(object = seurat) <- "integrated"
  
  expr_long <- expr_data %>%
    pivot_longer(
      cols = all_of(genes),
      names_to = "gene",
      values_to = "expr_value"
    ) 
  
  # Create individual plots with ggstatsplot and combine
  plot_list <- lapply(genes, function(gene) {
    gene_data <- expr_long %>% filter(gene == !!gene)
    
    ggstatsplot::ggbetweenstats(
      data = gene_data,
      x = combined_sample,
      # x = merged_sample,
      y = expr_value,
      type = "nonparametric",
      pairwise.comparisons = TRUE,
      pairwise.display = "significant",
      centrality.plotting = TRUE,
      centrality.point.args = list(size = 0, color = "darkred"),
      centrality.label.args = list(alpha = 0),
      p.adjust.method  = "bonferroni",
      palette = "Set1",
      package = "RColorBrewer",
      violin.args = list(width = 0.8, alpha = 0.7),
      point.args = list(position = ggplot2::position_jitterdodge(jitter.width = 0.8), alpha = 0.4, size = 1.5, stroke = 0, na.rm = TRUE),
      title = gene,
      results.subtitle = FALSE,
      xlab = "",
      ylab = "Expression",
      ggplot.component = list(theme(axis.title.y.right = element_blank(), 
                                    axis.text.y.right = element_blank(), 
                                    axis.ticks.y.right = element_blank()))
    )
  })
  
  # num_plots <- length(plot_list)
  # num_cols <- ceiling(sqrt(num_plots))
  # num_rows <- ceiling(num_plots / num_cols)
  # 
  num_cols <- 4
  num_rows <- 1
  
  combined_plot <- grid.arrange(grobs = plot_list, ncol = num_cols, nrow = num_rows)
  
  combined_plot
}

ggsave(
  file.path(working_dir, paste0("final/violin.", query, ".png")), 
  combined_plot, 
  width = 3*num_cols, # 3 for 3 samples; 4 for 4 samples?
  height = 5*num_rows, # 8 for more than 4 pairs of comparison; 5 is good for 5?
  dpi = 300, limitsize = FALSE
)

normalize_to_centered_median <- function(data_matrix) {
  # Convert to vector for easier manipulation
  data_vec <- as.vector(data_matrix)
  
  # Find min, max, and median
  min_val <- min(data_vec, na.rm = TRUE)
  max_val <- max(data_vec, na.rm = TRUE)
  median_val <- median(data_vec, na.rm = TRUE)
  
  # Calculate the target median (halfway between min and max)
  target_median <- (min_val + max_val) / 2
  
  # Create separate scaling for values below and above median
  normalized <- data_matrix
  
  # For values below median
  below_median <- data_matrix <= median_val
  if (any(below_median)) {
    # Scale to range from min_val to target_median
    normalized[below_median] <- min_val + 
      (data_matrix[below_median] - min_val) * 
      (target_median - min_val) / (median_val - min_val)
  }
  
  # For values above median
  above_median <- data_matrix > median_val
  if (any(above_median)) {
    # Scale to range from target_median to max_val
    normalized[above_median] <- target_median + 
      (data_matrix[above_median] - median_val) * 
      (max_val - target_median) / (max_val - median_val)
  }
  
  return(normalized)
}

# heatmap combined
{
  query <- "all" # lowercases will suffice
  DefaultAssay(object = seurat) <- "RNA"
  # genes <- rownames(seurat)[grep(paste0("(?i)", query), rownames(seurat))]
  # genes <- c("CD53", "FLT3.3", "ALPL.4", "DDX56", "MYO9A", "EFHD2", "Sox14/15/21 (SOX14)", "FABP2")
  # genes <- c("CD53", "FLT3.3", "ALPL.4", "DDX56", "MYO9A", "EFHD2", "Sox14/15/21 (SOX14)", "FABP2", "CASKIN1", "ABR", "STX7", "SNAP25", "TPM4", "CD34", "CDK6", "CDK6.2", "Hes. (HES1)", "Hes. (HES1).2", "Hes. (HES1).3")
  # genes <- c("CD53", "FLT3.3", "ALPL.4", "DDX56", "CD34", "TPM4", "EFHD2", "FABP2", "MYO9A", "Sox14/15/21 (SOX14)", "CASKIN1", "ABR", "STX7", "SNAP25", "CDK6", "Hes. (HES1).3")
  genes <- c("KIT", "CD34")
  
  expr_data <- FetchData(seurat, vars = c(genes, "meta_sample"), layer = 'counts')
  DefaultAssay(object = seurat) <- "integrated"
  
  expr_long <- expr_data %>%
    pivot_longer(
      cols = all_of(genes),
      names_to = "gene",
      values_to = "expr_value"
    ) 
  
  result_matrix <- tapply(expr_long$expr_value, 
                          list(expr_long$gene, expr_long$meta_sample), 
                          mean, na.rm = TRUE)
  
  # result_matrix <- log1p(result_matrix)
  
  result_matrix <- normalize_to_centered_median(result_matrix)
  
  print(c(max(result_matrix), median(result_matrix), min(result_matrix)))
  
  heatmap_data <- as.data.frame(result_matrix) %>%
    tibble::rownames_to_column(var = "gene") %>%
    reshape2::melt(id.vars = "gene", variable.name = "meta_sample", value.name = "mean_expression")
  
  heatmap_data$gene <- factor(heatmap_data$gene, levels = (genes))
  heatmap_data$meta_sample <- factor(heatmap_data$meta_sample, levels = rev(c("STEM", "CONTROL")))
  
  heatmap_plot <- ggplot(heatmap_data, aes(y = meta_sample, x = gene, fill = mean_expression)) +
    geom_tile() +
    scale_fill_gradientn(
      colors = c(sequential_hcl(50, h = 260, c = 65, l = c(37, 100), power = 1), rev(sequential_hcl(50, h = 13, c = 126, l = c(50, 100), power = 1))),
      limits = quantile(heatmap_data$mean_expression, c(0, 1), na.rm = TRUE),
    ) +
    theme_minimal() +
    labs(
      title = "Normalised Mean Gene Expression by Sample Group",
      y = "Sample Group",
      x = "Gene",
      fill = "Normalised\nMean Expression"
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid = element_blank(),
      plot.title = element_text(hjust = 0.5)
    )
  
  heatmap_plot
}

ggsave(file.path(working_dir, paste0("heatmap.", query, ".pdf")), plot = heatmap_plot, width = 8, height = 6, dpi = 300, limitsize = FALSE)

seurat <- only3_seurat(seurat)

# heatmap individual
{
  query <- "stemb" # lowercases will suffice
  DefaultAssay(object = seurat) <- "RNA"
  # genes <- rownames(seurat)[grep(paste0("(?i)", query), rownames(seurat))]
  genes <- c("CD53", "FLT3.3", "ALPL.4", "DDX56", "MYO9A", "EFHD2", "Sox14/15/21 (SOX14)", "FABP2")
  # genes <- c("CD53", "FLT3.3", "ALPL.4", "DDX56", "MYO9A", "EFHD2", "Sox14/15/21 (SOX14)", "FABP2", "CASKIN1", "ABR", "STX7", "SNAP25", "TPM4", "CD34", "CDK6", "CDK6.2", "Hes. (HES1)", "Hes. (HES1).2", "Hes. (HES1).3")
  # genes <- c("CD53", "FLT3.3", "ALPL.4", "DDX56", "CD34", "TPM4", "EFHD2", "FABP2", "MYO9A", "Sox14/15/21 (SOX14)", "CASKIN1", "ABR", "STX7", "SNAP25", "CDK6")
  # genes <- c("CD53", "FLT3.3", "ALPL.4", "DDX56", "CD34", "TPM4", "EFHD2", "FABP2", "MYO9A", "Sox14/15/21 (SOX14)", "CASKIN1", "ABR", "STX7", "SNAP25", "CDK6", "Hes. (HES1).3")
  # genes <- c("KIT", "CD34")
  
  expr_data <- FetchData(seurat, vars = c(genes, "combined_sample"), layer = 'counts')
  DefaultAssay(object = seurat) <- "integrated"
  
  expr_long <- expr_data %>%
    pivot_longer(
      cols = all_of(genes),
      names_to = "gene",
      values_to = "expr_value"
    ) 
  
  result_matrix <- tapply(expr_long$expr_value, 
                          list(expr_long$gene, expr_long$combined_sample), 
                          mean, na.rm = TRUE)
  
  # result_matrix <- log1p(result_matrix)
  
  result_matrix <- normalize_to_centered_median(result_matrix)
  
  print(c(max(result_matrix), median(result_matrix), min(result_matrix)))
  
  heatmap_data <- as.data.frame(result_matrix) %>%
    tibble::rownames_to_column(var = "gene") %>%
    reshape2::melt(id.vars = "gene", variable.name = "combined_sample", value.name = "mean_expression")
  
  # heatmap_data$gene <- factor(heatmap_data$gene, levels = genes)
  # heatmap_data$combined_sample <- factor(heatmap_data$combined_sample, levels = rev(c("STEM_B", "CONTROL")))
  heatmap_data$combined_sample <- factor(heatmap_data$combined_sample, levels = c("Stem", "ALDH", "Control"))
  
  heatmap_plot <- ggplot(heatmap_data, aes(y = combined_sample, x = gene, fill = mean_expression)) +
    geom_tile() +
    scale_fill_gradientn(
      colors = c(sequential_hcl(50, h = 260, c = 65, l = c(37, 100), power = 1), rev(sequential_hcl(50, h = 13, c = 126, l = c(50, 100), power = 1))),
      limits = quantile(heatmap_data$mean_expression, c(0, 1), na.rm = TRUE),
    ) +
    theme_minimal() +
    labs(
      title = "Normalised Mean Gene Expression by Sample Group",
      y = "Sample Group",
      x = "Gene",
      fill = "Normalised\nMean Expression"
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid = element_blank(),
      plot.title = element_text(hjust = 0.5)
    )
  
  heatmap_plot
}

# hclplot(sequential_hcl(3, h = 260, c = 65, l = c(37, 100), power = 1))
# hclplot(sequential_hcl(3, h = 13, c = 126, l = c(50, 100), power = 1))
# print(c(sequential_hcl(3, h = 260, c = 65, l = c(37, 100), power = 1), rev(sequential_hcl(50, h = 13, c = 126, l = c(50, 100), power = 1))))
# print(sequential_hcl(50, h = 13, c = 126, l = c(50, 100), power = 1))
# [1] "#3E53A0" "#FFFFFF" "#D93A33"
# (260, 65, 37), (13, 126, 50)
# [1] "#DA3B36" "#FBA09F" "#FFFFFF"
# [1] "#3C529F" "#9FA5CB" "#FFFFFF"

ggsave(file.path(working_dir, paste0("heatmap.", query, ".png")), plot = heatmap_plot, width = 8, height = 6, dpi = 300, limitsize = FALSE)

# featureplot combined
{
  # gene_ids <- c("integrated_KY21.Chr1.1900","sct_KY21.Chr3.1389","integrated_KY21.Chr9.235")
  # gene_ids <- c("FSTL5", "EFHD2", "KY21.Chr1.408", "FTH1.2", "ARHGAP15", "TP53INP1", "OR4B1", "BCL2L12")
  # gene_ids <- c("MEX3C", "KIFC1", "DDX56", "PRKD3", "HSPE1", "IMP4", "EXOSC6", "CD53", "FLT3.3", "ALPL.4", "PIWIL1", "Ci-ZF118 (DDX4)")
  # gene_ids <- c("CD53", "FLT3.3", "ALPL.4", "DDX56", "MYO9A", "EFHD2", "Sox14/15/21 (SOX14)", "FABP2")
  gene_ids <- c("CD53", "FLT3.3", "ALPL.4", "DDX56")
  # gene_ids <- c("MYO3B", "MYO9A", "MYO1C.2", "MYO10.2", "MYO1E", "MYO7A", "MYO18A", "MYO18B.3")
  # gene_ids <- c("SYAP1")
  # gene_ids <- c("Sox14/15/21 (SOX14)", "Sox4/11/12 (SOX4)" )
  # gene_ids <- c("IFITM3")
  # gene_ids <- c("FABP2")
  # gene_ids <- c("PRSS8.3")
  # gene_ids <- c("IFI30", "IFI30.2", "IFI30.3")
  # gene_ids <- c("VWF", "VWF.5")
  # gene_ids <- c( "FTH1.2",  "FTH1P19", "FTH1.3", "EFHD2")
  # gene_ids <- c("KIT", "CD34")
  gene_ids <- c(
    "Myc (MYCN)", "ALPL.4", "DDX56", "HSPA9", "OGT.2", # 0 - Stem
    "SHANK2", "IRX1", "SPON1", "PDIA3", # 1 - Stem
    "CD34", "NRG1", "THBD", # 2 - Stem
    "SNAP25", "SYAP1", "SEMA3E", "EFHD2", # 3 - Neural
    "Hes. (HES1).3", "GBP2", # 4 - Hematopoietic
    "PTPN23", "LGALS4.2", "Sfrp1/5 (SFRP5)", # 5 - Epithelial
    "MZB1", "CISH.2", # 6 - Hematopoietic progenitor
    "PAX6", "Atf4/5 (ATF5)", # 7 - Neural progenitor
    "Nr6a (NR6A1)", "CDC42", # 8 - Germ
    # "SEMA3E", # 9 - Neural
    # 9 - Neural
    "CNTNAP4", "GLUL", # 10 - Neural
    # 11 - Neural
    # "SYAP1", # 12 - Neural
    "SLC30A3", # 12 - Neural
    "CTSL", "MIF", "TRIB2", # 13 - Hematopoietic
    # "GBP2", # 13 - Hematopoietic
    "Sox4/11/12 (SOX4)", "DYRK1A", "ACSL1.2", # 14 - Neural progenitor
    "RHBDL2", "CYB5D2", # 15 - Epithelial
    "Hhex (HHEX)", "Pax2/5/8.b (PAX5)", # 16 - Hematopoietic
    "FHL2", "ACTN2", # 17 - Muscle
    "TSPAN6", "VIL1", # 18 - Epithelial
    "TPM1.2", "MYF5" # 19 - Muscle
  )
  gene_ids <- c(
    "Pou-r (POU6F1)", "DDX56", "MEX3C", "Myc (MYCN)", "ALPL.4", "NOP56", "NOP58", "HSPE1", # 0
    "CD53", "NFIA", "SPON1", "SLC12A3", "SHANK2", # 1
    "NRG1", "THBD", "CD34", # 2
    "SH3D19", "SEMA3E", "EFHD2", # 3
    "EYA1", "FBLN7", # 4
    "FLT3.3", "CD163", "CLBA1", # 5
    # 6
    "CYBB", "SHISAL2B", # 8
    "PIGK", "Serum (SRF)", "SPHK1", # 9
    "CTSA", "MALRD1", "ITGA2", # 10
    # 11
    # 12
    "CCDC18", "GBP2", "MIF", "LGALS4.2", "CTSL", "VIL1", # 13
    "ACSL1.2", # 14
    "FTCD", # 15
    "SYCP2", # 16
    "PTPN23", # 17
    "TMSB4X", # 18
    "SLC30A3", "SNAP25" # 19
  )
  
  gene_ids <- c("NOP56", "NOP58", "HSPE1", "ITGA2")
  
  plot_list <- list()
  DefaultAssay(seurat) <- 'RNA'
  
  for (gene_id in gene_ids) {
    plot <- FeaturePlot(object = seurat, features = gene_id, order = TRUE, slot = 'scale.data', cols = colorRampPalette(c("#FFF2F2", "#990000"))(2), pt.size = 1.5)
    # plot <- FeaturePlot(object = seurat, features = gene_id, order = TRUE, cols = colorRampPalette(c("#FFF2F2", "#990000"))(2), pt.size = 1.5)
    plot_list[[gene_id]] <- plot
  }
  
  # num_plots <- length(plot_list)
  # num_cols <- ceiling(sqrt(num_plots))
  # num_rows <- ceiling(num_plots / num_cols)
  
  # num_cols <- 4
  # num_rows <- 1
  
  # Combine all the individual plots into a single figure without empty rows or columns
  combined_plot <- grid.arrange(grobs = plot_list, ncol = num_cols, nrow = num_rows)
  combined_plot
  DefaultAssay(seurat) <- 'integrated'
  
  ggsave(
    file.path(working_dir, paste0("0716/featureplot.", query, ".pdf")), 
    combined_plot, 
    width = 4.5*num_cols, 
    height = 4.5*num_rows, 
    dpi = 300, limitsize = FALSE
  )
}

# featureplot individual
{
  # runid <- "0718"
  
  # 0718 markers
  genes <- c(
    "Pou-r (POU6F1)", "DDX56", "MEX3C", "Myc (MYCN)", "ALPL.4", "NOP56", "NOP58", "HSPE1", # 0
    "CD53", "NFIA", "SPON1", "SHANK2", # 1
    "NRG1", "THBD", "CD34", # 2
    "SEMA3E", "EFHD2", # 3
    "EYA1", # 4
    "FLT3.3", "CD163", # 5
    # 6
    "CYBB", "SHISAL2B", # 8
    "Serum (SRF)", "SPHK1", # 9
    "CTSA", "MALRD1", "ITGA2", # 10
    # 11
    # 12
    "GBP2", "MIF", "CTSL", "VIL1", # 13
    "ACSL1.2", # 14
    # 15
    # 16
    "PTPN23", # 17
    "TMSB4X", # 18
    "SLC30A3", "SNAP25", # 19
    "PAX6", "Hes. (HES1).3", "CISH.2", "RHBDL2"
  )
  
  plot_list <- list()
  DefaultAssay(seurat) <- 'RNA'
  
  for (gene in genes) {
    plot <- FeaturePlot(object = seurat, features = gene, order = TRUE, slot = 'scale.data', cols = colorRampPalette(c("#FFF2F2", "#990000"))(2), pt.size = 1.5)
    # plot <- FeaturePlot(object = seurat, features = gene_id, order = TRUE, cols = colorRampPalette(c("#FFF2F2", "#990000"))(2), pt.size = 1.5)
    plot_list[[gene]] <- plot
    
    # Replace spaces with underscores in genes for filename
    clean_gene_name <- gsub(" ", "_", gene)
    clean_gene_name <- gsub("/", ".", clean_gene_name)
    clean_gene_name <- gsub("\\(", "_", clean_gene_name)
    clean_gene_name <- gsub("\\)", "_", clean_gene_name)
    clean_gene_name <- gsub("__", "_", clean_gene_name)
    clean_gene_name <- gsub("_\\.", ".", clean_gene_name)  
    clean_gene_name <- gsub("\\._", ".", clean_gene_name) 
    clean_gene_name <- gsub("_+$", "", clean_gene_name)
    
    # Save individual plot
    ggsave(
      file.path(working_dir, paste0(runid, "/featureplot.", clean_gene_name, ".pdf")), 
      plot, 
      width = 4.5, 
      height = 4.5, 
      dpi = 300, 
      limitsize = FALSE
    )
  }
  
  DefaultAssay(seurat) <- 'integrated'
}

## batch plotting of marker cells for each cluster
# Load the cluster markers
cluster_markers <- readRDS(file.path(working_dir, "cluster_markers.rds"))

# Define the number of top genes to select
top_gene_counts <- 10

# Read the gene name to gene ID mapping
id_name_mapping <- read.csv("/home/groups/ayeletv/CionaGenome/id2name.filtered.nounderscores.csv")

# Loop through each cluster's markers to extract top genes and generate plots
for (cluster_id in names(cluster_markers)) {
  markers <- cluster_markers[[cluster_id]]  # Access the markers for the current cluster
  
  if (!is.null(markers) && "avg_log2FC" %in% colnames(markers)) {
    # Sort markers by average log fold change
    sorted_markers <- markers[order(-markers$avg_log2FC), ]
    
    # Select the top n genes
    top_genes <- head(rownames(sorted_markers), top_gene_counts)
    
    # Initialize an empty named vector for gene_id to display_name mapping
    gene_id_names <- setNames(c(), c())
    
    # Loop through each top gene to find its display name
    for (gene in top_genes) {
      # Find the display name from the mapping
      display_name <- id_name_mapping$gene_name[id_name_mapping$gene_id == gene]
      
      # If a display name is found, add it to the named vector
      if (length(display_name) > 0) {
        gene_id_names[display_name] <- gene
      }
    }
    
    # Initialize an empty list for plots
    plot_list <- list()
    
    # Loop through the named vector to create FeaturePlots
    for (display_name in names(gene_id_names)) {
      gene_id <- gene_id_names[[display_name]]
      
      # Create the FeaturePlot
      plot <- FeaturePlot(object = seurat, 
                          features = gene_id, 
                          order = TRUE, 
                          cols = colorRampPalette(c("#FFF2F2", "#990000"))(2), 
                          pt.size = 1.5) +
        ggtitle(display_name)  # Add custom title
      
      # Store the plot in the list
      plot_list[[display_name]] <- plot
    }
    
    # Calculate the number of rows and columns for the grid (2 rows and 5 columns)
    num_plots <- length(plot_list)
    num_cols <- 5  # Fixed to 5 columns
    num_rows <- ceiling(num_plots / num_cols)
    
    # Combine all the individual plots into a single figure
    combined_plot <- grid.arrange(grobs = plot_list, ncol = num_cols, nrow = num_rows)
    
    # Save the combined plot for the current cluster
    ggsave(file.path(working_dir, "markers", paste0("cluster_", cluster_id, "_top_marker_genes_featureplots.png")), 
           plot = combined_plot, 
           width = 18, height = 8, dpi = 300, limitsize = FALSE)
  }
}

write.table(all_top_genes_info, 
            file.path(working_dir, "top_marker_genes_all_clusters.txt"), 
            row.names = FALSE, 
            col.names = TRUE, 
            sep = "\t", 
            quote = FALSE)

### save cellnames
gene_lists <- list(
  "Ig_3" = c("HMCN2_1"),
  "Ig_C17orf99" = c("CCDC14", "COL7A1_1", "DPH1", "FUN_008790", "FUN_008793", "FUN_008795", "FUN_008799", "FUN_008801", "FUN_008802", "FUN_008803", "FUN_008804", "FUN_008805", "FUN_008808", "FUN_008809", "FUN_008810", "FUN_008811", "FUN_008812", "MCPH1", "OSBPL1A", "RIOK1", "STMN1", "WDR16", "WDR27"),
  "ig" = c("HMCN2_1"),
  "LRR_10" = c("FUN_032691", "FUN_032692", "FUN_032693", "FUN_032694", "FUN_032695", "FUN_032696", "FUN_032697"),
  "LRR_12" = c("FUN_020775"),
  "LRR_2" = c("ABT1", "ATG13", "BIRC2_73", "CDAN1", "CEP78", "CHD1", "CLEC19A_3", "DDX58_1", "EHMT2_2", "ERN2_2", "FLT1_2", "FUN_015787", "FUN_015791", "FUN_015792", "FUN_015793", "FUN_015795", "FUN_015797", "FUN_015798", "FUN_015799", "FUN_015800", "FUN_016711", "FUN_016712", "FUN_016714", "FUN_016715", "FUN_016722", "FUN_016724", "FUN_016726", "FUN_016727", "FUN_016729", "FUN_016730", "FUN_016731", "FUN_016737", "FUN_016739", "FUN_016740", "FUN_016742", "FUN_016743", "FUN_016746", "FUN_016747", "FUN_016748", "FUN_016750", "FUN_016752", "FUN_016753", "FUN_016754", "FUN_016757", "FUN_016759", "FUN_016762", "FUN_016764", "FUN_016766", "FUN_016767", "FUN_016768", "FUN_016769", "FUN_016770", "FUN_016771", "FUN_016772", "FUN_016773", "FUN_025575", "FUN_025576", "FUN_025577", "FUN_025579", "FUN_025580", "FUN_025581", "FUN_025582", "FUN_043962", "FUN_043963", "FUN_043964", "GJA8_2", "GJA8_3", "GJA8_4", "GJD2_2", "GPR64_3", "HOMER2", "IFIH1_1", "IFIH1_2", "IQANK1", "LRRCC1", "M6PR", "MCM2", "MDM1", "NAT10", "OTUD5", "POLD4", "RFC4", "SALL1", "SETD8", "TBX1_1", "TBX1_2", "TLK2_2", "TLR8", "TMEM98", "TRMT10C"),
  "LRR_6" = c("FUN_028373", "MICALL1"),
  "M157" = c("ATP2A2_2", "CAMSAP1_1", "FUN_014141", "FUN_014142", "FUN_014147", "mchr2", "PALD1", "PLBD1", "POLR3A", "SSH3"),
  "S-antigen" = c("CCDC14", "COL7A1_1", "DPH1", "FUN_008790", "FUN_008793", "FUN_008795", "FUN_008799", "FUN_008801", "FUN_008802", "FUN_008803", "FUN_008804", "FUN_008805", "FUN_008808", "FUN_008809", "FUN_008810", "FUN_008811", "FUN_008812", "MCPH1", "OSBPL1A", "RIOK1", "STMN1", "WDR16", "WDR27"),
  "Tcell_CD4_C" = c("FUN_040913", "FUN_040914", "FUN_040915"),
  "zf-RAG1" = c("ALDH9A1_1", "CD36_2", "CD36_3", "CYP4V2_15", "DEF8", "FUN_028448", "FUN_028449", "FUN_028452", "FUN_028453", "FUN_028454", "FUN_028456", "FUN_028458", "FUN_028459", "FUN_028460", "FUN_028461", "FUN_028462", "FUN_028463", "FUN_028464", "FUN_028465", "FUN_028466", "FUN_028469", "FUN_028470", "FUN_028471", "FUN_028474", "FUN_028475", "FUN_028477", "FUN_028480", "FUN_044423", "FUN_044424", "Heatr9", "PLB1_1", "PRSS22", "PSMD9", "RRM2B_2", "STX1A_2", "STX7", "TCEA1_1")
)

gene_lists <- list(
  "LRR_10.cluster" = c("FUN_032691", "FUN_032692", "FUN_032693", "FUN_032694", "FUN_032695", "FUN_032696", "FUN_032697"),
  "LRR_6.cluster" = c("FUN_028373", "MICALL1"),
  "Tcell_CD4_C.cluster" = c("FUN_040913", "FUN_040914", "FUN_040915")
) # VLR cluster genes

gene_lists <- list(
  "LRR_10.cluster" = c("FUN_032691", "FUN_032692", "FUN_032693", "FUN_032694", "FUN_032695", "FUN_032696", "FUN_032697"),
  "Tcell_CD4_C.cluster" = c("FUN_040913", "FUN_040914", "FUN_040915")
) # VLR cluster genes except for LRR6 whose expression pattern is  a bit weird

gene_lists <- list(
  "LRR_10.cluster" = c("FUN_032691", "FUN_032692", "FUN_032693", "FUN_032694", "FUN_032695", "FUN_032696", "FUN_032697")
) 

gene_lists <- list(
  "Tcell_CD4_C.cluster" = c("FUN_040913", "FUN_040914", "FUN_040915")
) 

gene_lists <- list(
  "pos.marker" = c("FUN_026261")
) 


### Obtain cell names expressing certain genes 
# # Convert RNA data to a sparse matrix (assuming you have the Seurat object)
# sparse_rna_data <- as(seurat@assays$RNA$data, "dgCMatrix")
# 
# sample_sheet <- data.frame(read.table("/home/groups/ayeletv/CionaGenome/id2name.filtered.nounderscores.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE))
# # Create a mapping from gene IDs to gene names
# gene_id_to_name <- setNames(sample_sheet$gene_name, sample_sheet$gene_id)
# 
# # Get the current row names of the sparse RNA data
# rna_row_names <- rownames(sparse_rna_data)
# 
# # Map gene IDs to gene names for row names
# new_row_names <- sapply(rna_row_names, function(gene_id) {
#   if (gene_id %in% names(gene_id_to_name)) {
#     return(gene_id_to_name[gene_id])
#   } else {
#     return(gene_id)  # Keep the original ID if no mapping is found
#   }
# })
# 
# # Update the row names of the sparse RNA data with gene names
# rownames(sparse_rna_data) <- new_row_names
# 
# # Check the updated row names
# head(rownames(sparse_rna_data))
# 
# # Create a list to store the results
# result_list <- list()
# 
# setwd("/scratch/users/jiamuyu/proj_botryllus/scRNAseq/240804_00_botryllus_ss3_new")
# # Find expressing cells for each gene list
# for (pfam_entry in names(gene_lists)) {
#     genenames <- gene_lists[[pfam_entry]]
#     
#     # Open a file to write the results in TSV format
#     file_path <- paste0(pfam_entry, ".tsv")
#     file_conn <- file(file_path, "w")
#     
#     # Clear the content of the file
#     truncate(file_conn, 0)
#     
#     # Find expressing cells for each gene
#     for (i in seq_along(genenames)) {
#         gene_index <- which(rownames(sparse_rna_data) == genenames[i])
#         if (length(gene_index) > 0) {
#             expressing_cells <- colnames(sparse_rna_data)[sparse_rna_data[gene_index, ] > 0]
#             # Write gene name and expressing cells to the TSV file
#             for (cell in expressing_cells) {
#                 writeLines(paste(genenames[i], cell, sep = "\t"), file_conn)
#             }
#         } else {
#             cat("No such gene found in the expression matrix:", genenames[i], "\n")
#         }
#     }
#     
#     close(file_conn) # Close the TSV file
# }
# for (pfam_entry in names(gene_lists)) {
#     genenames <- gene_lists[[pfam_entry]]

#     # Initialize a list to store expressing cells for each gene
#     expressing_cells_list <- vector("list", length = length(genenames))

#     # Find expressing cells for each gene
#     for (i in seq_along(genenames)) {
#         gene_index <- which(rownames(sparse_rna_data) == genenames[i])
#         if (length(gene_index) > 0) {
#             expressing_cells <- colnames(sparse_rna_data)[sparse_rna_data[gene_index, ] > 0]
#         } else {
#             expressing_cells <- character(0)
#             cat("No such gene found in the expression matrix:", genenames[i], "\n")
#         }
#         expressing_cells_list[[i]] <- expressing_cells
#     }

#     # Store the results in the list
#     for (i in seq_along(genenames)) {
#         result_list[[genenames[i]]] <- expressing_cells_list[[i]]
#     }
#     writeLines(expressing_cells, paste0(pfam_entry, ".txt"))
# }
# for (pfam_entry in names(gene_lists)) {
#     genenames <- gene_lists[[pfam_entry]]
#     expressing_cells <- colnames(sparse_rna_data)[colSums(sparse_rna_data[genenames, ]) > 0]

#     # Store the results in the list
#     for (gene in genenames) {
#         result_list[[gene]] <- expressing_cells
#     }
# }

# # Write the results to a TSV file
# write.table(data.frame(gene = unlist(lapply(result_list, names)), cell = unlist(result_list), row.names = FALSE, sep = "\t", file = "gene_cell_mapping.tsv"))


### Check cluster/sample distribution of cells expressing certain genes 
clusters <- Idents(seurat)

# Create a list to store expressing cells and their cluster IDs
expressing_cells_with_clusters <- list()

for (pfam_entry in names(gene_lists)) {
  genenames <- gene_lists[[pfam_entry]]
  
  for (gene in genenames) {
    gene_index <- which(rownames(sparse_rna_data) == gene)
    if (length(gene_index) > 0) {
      expressing_cells <- colnames(sparse_rna_data)[sparse_rna_data[gene_index, ] > 0]
      for (cell in expressing_cells) {
        cluster_id <- clusters[cell]
        expressing_cells_with_clusters[[cell]] <- cluster_id
      }
    }
  }
}

# Convert list to a dataframe
expressing_cells_df <- data.frame(
  cell = names(expressing_cells_with_clusters),
  cluster = unlist(expressing_cells_with_clusters),
  stringsAsFactors = FALSE
)

# Filter the Seurat object to include only the selected cells
selected_cells <- expressing_cells_df$cell
seurat_filtered <- subset(seurat, cells = selected_cells)

### composition of samples and clusters
table_samples_by_clusters <- seurat_filtered@meta.data %>%
  group_by(merged_sample, seurat_clusters) %>%
  summarize(count = n()) %>%
  spread(seurat_clusters, count, fill = 0) %>%
  ungroup() %>%
  mutate(total_cell_count = rowSums(.[c(2:ncol(.))])) %>%
  dplyr::select(c('merged_sample', 'total_cell_count', everything())) %>%
  arrange(factor(Sample, levels = levels(seurat_filtered@meta.data$sample)))

# knitr::kable(table_samples_by_clusters)

table_clusters_by_samples <- seurat_filtered@meta.data %>%
  dplyr::rename('cluster' = 'seurat_clusters') %>%
  group_by(cluster, merged_sample) %>%
  summarize(count = n()) %>%
  spread(merged_sample, count, fill = 0) %>%
  ungroup() %>%
  mutate(total_cell_count = rowSums(.[c(2:ncol(.))])) %>%
  select(c('cluster', 'total_cell_count', everything())) %>%
  arrange(factor(cluster, levels = levels(seurat_filtered@meta.data$seurat_clusters)))

# knitr::kable(table_clusters_by_samples)

# Prepare labels for the sample distribution plot
temp_labels_samples <- seurat_filtered@meta.data %>%
  group_by(merged_sample) %>%
  tally()

# Plot for the samples by clusters
p1 <- table_samples_by_clusters %>%
  select(-c('total_cell_count')) %>%
  melt(id.vars = 'merged_sample') %>%
  mutate(Sample = factor(Sample, levels = unique(seurat_filtered@meta.data$sample))) %>%
  ggplot(aes(merged_sample, value)) +
  geom_bar(aes(fill = variable), position = 'stack', stat = 'identity') +
  geom_text(
    data = temp_labels_samples,
    aes(x = merged_sample, y = Inf, label = paste0('n = ', format(n, big.mark = ',', trim = TRUE)), vjust = -1),
    color = 'black', size = 2.8
  ) +
  scale_fill_manual(name = 'Cluster', values = custom_colors$discrete) +
  scale_y_continuous(name = 'Number of cells', labels = scales::comma, expand = c(0.01, 0)) +
  coord_cartesian(clip = 'off') +
  theme_bw() +
  theme(
    legend.position = 'left',
    plot.title = element_text(hjust = 0.5),
    text = element_text(size = 16),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    plot.margin = margin(t = 20, r = 0, b = 0, l = 0, unit = 'pt')
  )

# Prepare labels for the cluster distribution plot
temp_labels_clusters <- seurat_filtered@meta.data %>%
  group_by(seurat_clusters) %>%
  tally() %>%
  dplyr::rename('cluster' = seurat_clusters)

# Plot for the clusters by samples
p2 <- table_clusters_by_samples %>%
  select(-c('total_cell_count')) %>%
  melt(id.vars = 'cluster') %>%
  mutate(cluster = factor(cluster, levels = levels(seurat_filtered@meta.data$seurat_clusters))) %>%
  ggplot(aes(cluster, value)) +
  geom_bar(aes(fill = variable), position = 'stack', stat = 'identity') +
  geom_text(
    data = temp_labels_clusters,
    aes(x = cluster, y = Inf, label = paste0('n = ', format(n, big.mark = ',', trim = TRUE)), vjust = -1),
    color = 'black', size = 2.8
  ) +
  scale_fill_manual(name = 'Sample', values = custom_colors$discrete) +
  scale_y_continuous(labels = scales::comma, expand = c(0.01, 0)) +
  coord_cartesian(clip = 'off') +
  theme_bw() +
  theme(
    legend.position = 'right',
    plot.title = element_text(hjust = 0.5),
    text = element_text(size = 16),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title = element_blank(),
    plot.margin = margin(t = 20, r = 0, b = 0, l = 10, unit = 'pt')
  )

# Save the combined plots
ggsave(
  file.path(working_dir, 'composition_selected_cells_samples_clusters_by_number.pos.marker.png'),
  p1 + p2 +
    plot_layout(ncol = 2, widths = c(
      seurat_filtered@meta.data$merged_sample %>% unique() %>% length(),
      seurat_filtered@meta.data$seurat_clusters %>% unique() %>% length()
    )),
  width = 18, height = 8
)

vlr_transp <- rownames(seurat_filtered@meta.data)[seurat_filtered@meta.data$seurat_clusters
                                                  %in% c(3, 4)]

de_vlr_transp <- FindMarkers(seurat, ident.1 = c(3, 4), cells = vlr_transp, logfc.threshold = log(2), min.diff.pct = 0.25)

de_vlr_transp <- subset(de_vlr_transp, p_val_adj < 0.05)

output_file <- paste0("/scratch/users/jiamuyu/proj_botryllus/scRNAseq/240804_00_botryllus_ss3_new/", "vlr_transp", "_markers.csv")
# row.names(markers) <- markers$gene
write.csv(de_vlr_transp, file = output_file, row.names = TRUE)

vlr_brain_cellislands <- rownames(seurat_filtered@meta.data)[seurat_filtered@meta.data$seurat_clusters
                                                             %in% c(7, 9)]

de_vlr_brain_cellislands <- FindMarkers(seurat, ident.1 = c(7, 9), cells = vlr_brain_cellislands, logfc.threshold = log(2), min.diff.pct = 0.25)

de_vlr_brain_cellislands <- subset(de_vlr_brain_cellislands, p_val_adj < 0.05)

output_file <- paste0("/scratch/users/jiamuyu/proj_botryllus/scRNAseq/240804_00_botryllus_ss3_new/", "vlr_brain_cellislands", "_markers.csv")
# row.names(markers) <- markers$gene
write.csv(de_vlr_brain_cellislands, file = output_file, row.names = TRUE)

### marker gene analysis
# ## by sample, integrated
markers <- FindAllMarkers(seurat, assay = "RNA", group.by = "merged_sample", logfc.threshold = log(2))
markers <- subset(markers, p_val_adj < 0.05)
colnames(markers)[colnames(markers) == "cluster"] <- "sample"
markers_sorted <- markers %>%
  arrange(sample, desc(avg_log2FC))

write.csv(markers, file.path(working_dir, paste0("markers/markers_by_4sample_pval.csv")), row.names = FALSE, quote = FALSE)
write.csv(markers_sorted, file.path(working_dir, paste0("markers/markers_by_4sample_fc.csv")), row.names = FALSE, quote = FALSE)

write.csv(markers, file.path(working_dir, paste0("markers/markers_by_mergedsample_pval.csv")), row.names = FALSE, quote = FALSE)
write.csv(markers_sorted, file.path(working_dir, paste0("markers/markers_by_mergedsample_fc.csv")), row.names = FALSE, quote = FALSE)

markers <- read.csv(file.path(working_dir, "markers/markers_by_mergedsample_pval.csv"),stringsAsFactors = FALSE)
markers <- read.csv(file.path(working_dir, "markers/markers_by_mergedsample_fc.csv"), stringsAsFactors = FALSE)

top_genes_heatmap <- markers %>%
  group_by(sample) %>%
  top_n(n = 10, wt = avg_log2FC) %>%
  pull(gene) %>%
  unique()

p <- DoHeatmap(
  seurat_filtered, 
  features = top_genes_heatmap, 
  group.by = "merged_sample", 
  assay = "RNA", 
  angle = 0, group.colors = custom_colors$discrete, 
  size = 3, group.bar.height = 0.01, hjust = 0.5
) +
  ggtitle("Top Markers per Sample") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_viridis(option = "B")

ggsave(
  file.path(working_dir, 'final/doheatmap.sample.single.svg'), 
  p, 
  height = 8, 
  width = 8
)

# 4_conserved markers for each cluster
conserved_markers_list <- list()
for (cluster in unique(Idents(seurat_filtered))) {
  # Find markers conserved in this cluster across all samples
  markers_conserved <- FindConservedMarkers(
    seurat_filtered,
    ident.1 = cluster,          # Which cluster to analyze
    grouping.var = "merged_sample",  # Group by sample
    assay = "RNA",
    logfc.threshold = log(2)
  )
  
  if (nrow(markers_conserved) > 0) {
    markers_conserved$gene <- rownames(markers_conserved)
    markers_conserved$cluster <- cluster
    conserved_markers_list[[as.character(cluster)]] <- markers_conserved
  }
}

all_conserved_markers <- bind_rows(conserved_markers_list)
all_conserved_markers <- all_conserved_markers %>%
  rowwise() %>%
  dplyr::mutate(max_pval = max(c_across(ends_with("_p_val_adj")))) %>%
  filter(max_pval < 0.05)

write.csv(all_conserved_markers, 
          file.path(working_dir, "markers/4_conserved_markers_by_cluster_pval.csv"), 
          row.names = FALSE, quote = FALSE)
all_conserved_markers <- read.csv(file.path(working_dir, "markers/4_conserved_markers_by_cluster_pval.csv"),stringsAsFactors = FALSE)

all_conserved_markers <- all_conserved_markers %>%
  rowwise() %>%
  dplyr::mutate(min_logfc = min(c_across(ends_with("_avg_log2FC")))) %>%
  arrange(cluster, desc(min_logfc))

write.csv(all_conserved_markers, 
          file.path(working_dir, "markers/4_conserved_markers_by_cluster_fc.csv"), 
          row.names = FALSE, quote = FALSE)
all_conserved_markers <- read.csv(file.path(working_dir, "markers/4_conserved_markers_by_cluster_fc.csv"),stringsAsFactors = FALSE)

top_conserved_genes <- all_conserved_markers %>%
  dplyr::mutate(cluster = as.numeric(cluster)) %>% 
  arrange(cluster) %>%                      
  group_by(cluster) %>%
  top_n(n = 5, wt = min_logfc) %>%
  pull(gene) %>%
  unique()

p <- DoHeatmap(
  seurat_filtered, 
  features = top_conserved_genes, 
  group.by = "seurat_clusters",  
  assay = "RNA", 
  angle = 0, 
  group.colors = custom_colors$discrete, 
  size = 3, 
  group.bar.height = 0.01, 
  hjust = 0.5
) +
  ggtitle("Conserved Markers Across Samples") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_viridis(option = "B")

ggsave(
  file.path(working_dir, 'final/doheatmap.4_conserved_markers.png'), 
  p, 
  height = 10, 
  width = 12
)

# 4_conserved markers for each cluster
conserved_markers_list <- list()
for (cluster in unique(Idents(seurat_filtered))) {
  # Find markers conserved in this cluster across all samples
  markers_conserved <- FindConservedMarkers(
    seurat_filtered,
    ident.1 = cluster,          
    grouping.var = "combined_sample",  
    assay = "RNA",
    logfc.threshold = log(2)
  )
  
  if (nrow(markers_conserved) > 0) {
    markers_conserved$gene <- rownames(markers_conserved)
    markers_conserved$cluster <- cluster
    conserved_markers_list[[as.character(cluster)]] <- markers_conserved
  }
}

all_conserved_markers <- bind_rows(conserved_markers_list)
all_conserved_markers <- all_conserved_markers %>%
  rowwise() %>%
  dplyr::mutate(min_pval = min(c_across(ends_with("_p_val_adj")))) %>%
  filter(min_pval < 0.05)

write.csv(all_conserved_markers, 
          file.path(working_dir, "markers/4_lenient_conserved_markers_by_cluster_pval.csv"), 
          row.names = FALSE, quote = FALSE)
all_conserved_markers <- read.csv(file.path(working_dir, "markers/4_lenient_conserved_markers_by_cluster_pval.csv"),stringsAsFactors = FALSE)

all_conserved_markers <- all_conserved_markers %>%
  rowwise() %>%
  dplyr::mutate(max_logfc = max(c_across(ends_with("_avg_log2FC")))) %>%
  arrange(cluster, desc(max_logfc))

write.csv(all_conserved_markers, 
          file.path(working_dir, "markers/4_lenient_conserved_markers_by_cluster_fc.csv"), 
          row.names = FALSE, quote = FALSE)
all_conserved_markers <- read.csv(file.path(working_dir, "markers/4_lenient_conserved_markers_by_cluster_fc.csv"),stringsAsFactors = FALSE)

top_conserved_genes <- all_conserved_markers %>%
  dplyr::mutate(cluster = as.numeric(cluster)) %>% 
  arrange(cluster) %>%                      
  group_by(cluster) %>%
  top_n(n = 5, wt = max_logfc) %>%
  pull(gene) %>%
  unique()

p <- DoHeatmap(
  # seurat_filtered, 
  seurat, 
  features = top_conserved_genes, 
  group.by = "seurat_clusters",  
  assay = "RNA", 
  angle = 0, 
  group.colors = custom_colors$discrete, 
  size = 3, 
  group.bar.height = 0.01, 
  hjust = 0.5
) +
  ggtitle("Conserved Markers Across Samples") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_viridis(option = "B")

ggsave(
  file.path(working_dir, 'final/doheatmap.4_conserved_markers.png'), 
  p, 
  height = 10, 
  width = 12
)

# merged 
sample_specific_clusters <- c(3, 6, 9, 10, 12, 14, 17, 19)

all_markers <- FindAllMarkers(
  seurat_filtered,
  assay = "RNA",
  slot = "data",
  logfc.threshold = log(2),
  test.use = "wilcox",
  only.pos = TRUE,
  min.pct = 0.25
)

all_markers <- all_markers %>%
  filter(p_val_adj < 0.05) %>%
  mutate(gene = rownames(.))

conserved_markers_list <- list()
for (cluster in setdiff(unique(Idents(seurat_filtered)), sample_specific_clusters)) {
  markers_conserved <- FindConservedMarkers(
    seurat_filtered,
    ident.1 = cluster,
    grouping.var = "merged_sample",
    assay = "RNA",
    slot = "data",
    logfc.threshold = log(2),
    min.pct = 0.25
  )
  
  if (nrow(markers_conserved) > 0) {
    markers_conserved$gene <- rownames(markers_conserved)
    markers_conserved$cluster <- cluster
    conserved_markers_list[[as.character(cluster)]] <- markers_conserved
  }
}

all_conserved_markers <- bind_rows(conserved_markers_list) %>%
  rowwise() %>%
  dplyr::mutate(
    max_pval = max(c_across(ends_with("_p_val_adj"))),
    min_logfc = min(c_across(ends_with("_avg_log2FC")))
  ) %>%
  filter(max_pval < 0.05) %>%
  arrange(cluster, desc(min_logfc))

sample_specific_markers <- all_markers %>%
  filter(cluster %in% sample_specific_clusters) %>%
  dplyr::mutate(min_logfc = avg_log2FC)  

merged_markers <- bind_rows(
  as.data.frame(all_conserved_markers) %>%    
    ungroup() %>%
    dplyr::select(gene, cluster, min_logfc) %>%
    mutate(marker_type = "conserved"),
  as.data.frame(sample_specific_markers) %>%  
    ungroup() %>%
    dplyr::select(gene, cluster, min_logfc) %>%
    mutate(marker_type = "sample_specific")
)

write.csv(merged_markers,
          file.path(working_dir, "markers/merged_markers.csv"),
          row.names = FALSE, quote = FALSE
)
merged_markers <- read.csv(file.path(working_dir, "markers/merged_markers.csv"),stringsAsFactors = FALSE)

top_genes <- merged_markers %>%
  dplyr::mutate(cluster = as.numeric(cluster)) %>%
  arrange(cluster) %>%
  group_by(cluster) %>%
  # filter(!cluster %in% c(11)) %>%  # Add this line to exclude clusters 6 and 11
  top_n(n = 5, wt = min_logfc) %>%
  pull(gene) %>%
  unique()

p <- DoHeatmap(
  seurat_filtered,
  features = top_genes,
  group.by = "seurat_clusters",
  assay = "RNA",
  angle = 0,
  group.colors = custom_colors$discrete,
  size = 3,
  group.bar.height = 0.01,
  hjust = 0.5
) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_viridis(option = "B")

ggsave(
  file.path(working_dir, 'final/doheatmap.merged_markers.svg'),
  p,
  height = 12,
  width = 12
)

# customised markers
query <- "FDZ"
rownames(seurat)[grep(paste0("(?i)", query), rownames(seurat))]

nonspecific <- c("NPM1", "OSBPL9", "HLA", "DDX11", "IRF-like-7 (IRF8)", "TLR8", "PRDX6", "CSNK1G1", "ZEB1", "GCSH", "SYTL4", "KDM8", "ASPM", "VEGFR (FLT1)", "XRN1", "PDE4B.2", "PDE3B", "ASIC2", "TDRD1", "NOB1", "SLIT1.2", "LYG2", "CUBN", "PXDN", "SEMA5B", "CUZD1", "UMOD", "NOTCH1", "ACKR3", "SLC23A1", "AQP8", "MMP16", "ZEB1", "MRC1", "C9orf72", "ASIC5", "Hmgtun1 (TFAM)", "FN1", "GLIPR1", )
suboptimal <- c("SPON1", "UNC13C", "POGLUT2", "ZEB1", "PRKCB", "PDE6D", "PDE4B.2", "PDE3B", "DAXX", "GLIPR2", "Hox4 (HOXC4)", )

# my markers
markers <- c(
  "DDX56", "MEX3C", "Pou-r (POU6F1)", "ALPL.4", # 0
  "CD53", "NFIA", "SLC12A3", "SHANK2", # 1
  "CD34", "CRYAB", # 2
  "SH3D19", # 3
  "FBLN7", "EYA1", # 4
  "CD163", "FLT3.3", "CLBA1", # 5
  # 6
  "FUT6", # 7
  "CYBB", "SHISAL2B", # 8
  "PIGK", "Serum (SRF)", "SPHK1", # 9
  "CTSA", "MALRD1", # 10
  "SLC32A1", # 12
  "CCDC18", # 13
  "SYCP2", # 16
  "TMSB4X", # 18
  "SNAP25" # 19
)

# TG raw
markers <- c(
  "Myc (MYCN)", "ALPL.4", "DDX56", "HSPA9", "OGT.2", # 0 - Stem
  "SHANK2", "IRX1", # 1 - Stem
  "CD34", "NRG1", # 2 - Stem
  "SYAP1", "SNAP25", "SEMA3E", # 3 - Neural
  "Hes. (HES1).3", "GBP2", # 4 - Hematopoietic
  "PTPN23", "LGALS4.2", "Sfrp1/5 (SFRP5)", # 5 - Epithelial
  "MZB1", "CISH.2", # 6 - Hematopoietic progenitor
  "PAX6", "Atf4/5 (ATF5)", # 7 - Neural progenitor
  "Nr6a (NR6A1)", "CDC42", # 8 - Germ
  # "SEMA3E", # 9 - Neural
  "EFHD2", # 9 - Neural
  "CNTNAP4", "GLUL", # 10 - Neural
  "SPON1", "PDIA3", # 11 - Neural
  # "SYAP1", # 12 - Neural
  "SLC30A3", # 12 - Neural
  "CTSL", "MIF", "TRIB2", # 13 - Hematopoietic
  # "GBP2", # 13 - Hematopoietic
  "Sox4/11/12 (SOX4)", "DYRK1A", "ACSL1.2", # 14 - Neural progenitor
  "RHBDL2", "CYB5D2", # 15 - Epithelial
  "Hhex (HHEX)", "Pax2/5/8.b (PAX5)", "THBD", # 16 - Hematopoietic
  "FHL2", "ACTN2", # 17 - Muscle
  "TSPAN6", "VIL1", # 18 - Epithelial
  "TPM1.2", "MYF5" # 19 - Muscle
)

# TG refined 
markers <- c(
  "Myc (MYCN)", "ALPL.4", "DDX56", "HSPA9", "OGT.2", # 0 - Stem
  "SHANK2", "IRX1", "SPON1", "PDIA3", # 1 - Stem
  "CD34", "NRG1", "THBD", # 2 - Stem
  "SNAP25", "SYAP1", "SEMA3E", "EFHD2", # 3 - Neural
  "Hes. (HES1).3", "GBP2", # 4 - Hematopoietic
  "PTPN23", "LGALS4.2", "Sfrp1/5 (SFRP5)", # 5 - Epithelial
  "MZB1", "CISH.2", # 6 - Hematopoietic progenitor
  "PAX6", "Atf4/5 (ATF5)", # 7 - Neural progenitor
  "Nr6a (NR6A1)", "CDC42", # 8 - Germ
  # "SEMA3E", # 9 - Neural
  # 9 - Neural
  "CNTNAP4", "GLUL", # 10 - Neural
  # 11 - Neural
  # "SYAP1", # 12 - Neural
  "SLC30A3", # 12 - Neural
  "CTSL", "MIF", "TRIB2", # 13 - Hematopoietic
  # "GBP2", # 13 - Hematopoietic
  "Sox4/11/12 (SOX4)", "DYRK1A", "ACSL1.2", # 14 - Neural progenitor
  "RHBDL2", "CYB5D2", # 15 - Epithelial
  "Hhex (HHEX)", "Pax2/5/8.b (PAX5)", # 16 - Hematopoietic
  "FHL2", "ACTN2", # 17 - Muscle
  "TSPAN6", "VIL1", # 18 - Epithelial
  "TPM1.2", "MYF5" # 19 - Muscle
)

# merged with my markers
# markers <- c(
#   "Pou-r (POU6F1)", "DDX56", "MEX3C", "Myc (MYCN)", "ALPL.4", "NOP56", "NOP58", "HSPE1", # 0
#   "CD53", "NFIA", "SPON1", "SLC12A3", "SHANK2", # 1
#   "NRG1", "THBD", "CD34", # 2
#   "SH3D19", "SEMA3E", "EFHD2", # 3
#   "EYA1", "FBLN7", # 4
#   "FLT3.3", "CD163", "CLBA1", # 5
#   # 6
#   "CYBB", "SHISAL2B", # 8
#   "PIGK", "Serum (SRF)", "SPHK1", # 9
#   "CTSA", "MALRD1", "ITGA2", # 10
#   # 11
#   # 12
#   "CCDC18", "GBP2", "MIF", "LGALS4.2", "CTSL", "VIL1", # 13
#   "ACSL1.2", # 14
#   "FTCD", # 15
#   "SYCP2", # 16
#   "PTPN23", # 17
#   "TMSB4X", # 18
#   "SLC30A3", "SNAP25" # 19
# )

# # 0718 markers
# markers <- c(
#   "Pou-r (POU6F1)", "DDX56", "MEX3C", "Myc (MYCN)", "ALPL.4", "NOP56", "NOP58", "HSPE1", # 0
#   "CD53", "NFIA", "SPON1", "SHANK2", # 1
#   "NRG1", "THBD", "CD34", # 2
#   "SEMA3E", "EFHD2", # 3
#   "EYA1", # 4
#   "FLT3.3", "CD163", # 5
#   # 6
#   "CYBB", "SHISAL2B", # 8
#   "Serum (SRF)", "SPHK1", # 9
#   "CTSA", "MALRD1", "ITGA2", # 10
#   # 11
#   # 12
#   "GBP2", "MIF", "CTSL", "VIL1", # 13
#   "ACSL1.2", # 14
#   # 15
#   # 16
#   "PTPN23", # 17
#   "TMSB4X", # 18
#   "SLC30A3", "SNAP25", # 19
#   "PAX6", "Hes. (HES1).3", "CISH.2", "RHBDL2"
# )

# 0720 markers
markers <- c(
  "Pou-r (POU6F1)", "DDX56", "MEX3C", "Myc (MYCN)", "ALPL.4", "NOP56", "NOP58", "HSPE1", # 0
  "CD53", "NFIA", "SPON1", "SHANK2", # 1
  "NRG1", "THBD", "CD34", # 2
  "SEMA3E", "EFHD2", # 3
  "EYA1", # 4
  "FLT3.3", "CD163", # 5
  # 6
  "CYBB", "SHISAL2B", # 8
  "Serum (SRF)", "SPHK1", # 9
  "CTSA", "MALRD1", "ITGA2", # 10
  # 11
  # 12
  "GBP2", "MIF", "CTSL", "VIL1", # 13
  "ACSL1.2", # 14
  "PAX6", # 15
  # 16
  "PTPN23", # 17
  "TMSB4X", # 18
  "SLC30A3", "SNAP25" # 19
)

# rownames(seurat[['integrated']])[grep("(?i)RHBDL2", rownames(seurat[['integrated']]))]

p <-DotPlot(
  seurat,
  # seurat_filtered,
  markers,
  group.by = "seurat_clusters",
  assay = "RNA",
  scale = TRUE
) + 
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y = element_text(size = 8),
    plot.margin = margin(t = 5, r = 5, b = 5, l = 50, unit = "pt")
  ) +
  labs(x = "Genes", y = "Clusters")

ggsave(
  file.path(working_dir, 'final/dotplot.custom_markers.0720.pdf'),
  p,
  height = 8,
  width = 16
)

top_genes <- merged_markers %>%
  dplyr::mutate(cluster = as.numeric(cluster)) %>%
  arrange(cluster) %>%
  group_by(cluster) %>%
  # filter(!cluster %in% c(6, 11)) %>%  # Add this line to exclude clusters 6 and 11
  top_n(n = 3, wt = min_logfc) %>%
  pull(gene) %>%
  unique()

p <- DotPlot(
  seurat_filtered,
  features = top_genes,
  group.by = "seurat_clusters",
  assay = "RNA",
  scale = TRUE
) + 
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y = element_text(size = 8),
    plot.margin = margin(t = 5, r = 5, b = 5, l = 50, unit = "pt")
  ) +
  labs(x = "Genes", y = "Cluster")

ggsave(
  file.path(working_dir, 'final/dotplot.merged_markers.svg'),
  p,
  height = 8,
  width = 16
)

# conserved markers for each cluster
conserved_markers_list <- list()
for (cluster in unique(Idents(seurat_filtered))) {
  # Find markers conserved in this cluster across all samples
  markers_conserved <- FindConservedMarkers(
    seurat_filtered,
    ident.1 = cluster,          # Which cluster to analyze
    grouping.var = "meta_sample",  # Group by sample
    assay = "RNA",
    logfc.threshold = log(2)
  )
  
  if (nrow(markers_conserved) > 0) {
    markers_conserved$gene <- rownames(markers_conserved)
    markers_conserved$cluster <- cluster
    conserved_markers_list[[as.character(cluster)]] <- markers_conserved
  }
}

all_conserved_markers <- bind_rows(conserved_markers_list)
all_conserved_markers <- all_conserved_markers %>%
  rowwise() %>%
  dplyr::mutate(max_pval = max(c_across(ends_with("_p_val_adj")))) %>%
  filter(max_pval < 0.05)

write.csv(all_conserved_markers, 
          file.path(working_dir, "markers/2_conserved_markers_by_cluster_pval.csv"), 
          row.names = FALSE, quote = FALSE)
all_conserved_markers <- read.csv(file.path(working_dir, "markers/2_conserved_markers_by_cluster_pval.csv"),stringsAsFactors = FALSE)

all_conserved_markers <- all_conserved_markers %>%
  rowwise() %>%
  dplyr::mutate(min_logfc = min(c_across(ends_with("_avg_log2FC")))) %>%
  arrange(cluster, desc(min_logfc))

write.csv(all_conserved_markers, 
          file.path(working_dir, "markers/2_conserved_markers_by_cluster_fc.csv"), 
          row.names = FALSE, quote = FALSE)
all_conserved_markers <- read.csv(file.path(working_dir, "markers/2_conserved_markers_by_cluster_fc.csv"),stringsAsFactors = FALSE)

top_conserved_genes <- all_conserved_markers %>%
  dplyr::mutate(cluster = as.numeric(cluster)) %>% 
  arrange(cluster) %>%                      
  group_by(cluster) %>%
  top_n(n = 5, wt = min_logfc) %>%
  pull(gene) %>%
  unique()

p <- DoHeatmap(
  # seurat_filtered, 
  seurat,
  # features = top_conserved_genes,
  features = markers,
  group.by = "seurat_clusters",  
  assay = "RNA", 
  angle = 0, 
  group.colors = custom_colors$discrete, 
  size = 3, 
  group.bar.height = 0.01, 
  hjust = 0.5
) +
  # ggtitle("Conserved Markers Across Samples") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_viridis(option = "B")

ggsave(
  # file.path(working_dir, 'final/doheatmap.2_conserved_markers.svg'),
  file.path(working_dir, 'final/doheatmap.custom_markers.0720.pdf'),
  p, 
  height = 10, 
  width = 12
)

# conserved markers for each cluster
conserved_markers_list <- list()
for (cluster in unique(Idents(seurat_filtered))) {
  # Find markers conserved in this cluster across all samples
  markers_conserved <- FindConservedMarkers(
    seurat_filtered,
    ident.1 = cluster,          # Which cluster to analyze
    grouping.var = "meta_sample",  # Group by sample
    assay = "RNA",
    logfc.threshold = log(2)
  )
  
  if (nrow(markers_conserved) > 0) {
    markers_conserved$gene <- rownames(markers_conserved)
    markers_conserved$cluster <- cluster
    conserved_markers_list[[as.character(cluster)]] <- markers_conserved
  }
}

all_conserved_markers <- bind_rows(conserved_markers_list)
all_conserved_markers <- all_conserved_markers %>%
  rowwise() %>%
  dplyr::mutate(min_pval = min(across(ends_with("_p_val_adj"), .fns = identity)[[1]])) %>%
  filter(min_pval < 0.05)

write.csv(all_conserved_markers, 
          file.path(working_dir, "markers/2_lenient_conserved_markers_by_cluster_pval.csv"), 
          row.names = FALSE, quote = FALSE)
all_conserved_markers <- read.csv(file.path(working_dir, "markers/2_lenient_conserved_markers_by_cluster_pval.csv"),stringsAsFactors = FALSE)

all_conserved_markers <- all_conserved_markers %>%
  rowwise() %>%
  dplyr::mutate(max_logfc = max(c_across(ends_with("_avg_log2FC")))) %>%
  arrange(cluster, desc(max_logfc))

write.csv(all_conserved_markers, 
          file.path(working_dir, "markers/2_lenient_conserved_markers_by_cluster_fc.csv"), 
          row.names = FALSE, quote = FALSE)
all_conserved_markers <- read.csv(file.path(working_dir, "markers/2_lenient_conserved_markers_by_cluster_fc.csv"),stringsAsFactors = FALSE)

top_conserved_genes <- all_conserved_markers %>%
  dplyr::mutate(cluster = as.numeric(cluster)) %>% 
  arrange(cluster) %>%                      
  group_by(cluster) %>%
  # top_n(n = 5, wt = max_logfc) %>%
  top_n(n = 5, wt = min_pval) %>%
  pull(gene) %>%
  unique()

p <- DoHeatmap(
  seurat_filtered, 
  features = top_conserved_genes, 
  group.by = "seurat_clusters",  
  assay = "RNA", 
  angle = 0, 
  group.colors = custom_colors$discrete, 
  size = 3, 
  group.bar.height = 0.01, 
  hjust = 0.5
) +
  ggtitle("Conserved Markers Across Samples") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_viridis(option = "B")

ggsave(
  file.path(working_dir, 'final/doheatmap.2_lenient_conserved_markers.png'),
  p, 
  height = 10, 
  width = 12
)

# conserved markers for each cluster
conserved_markers_list <- list()
for (cluster in unique(Idents(seurat_filtered))) {
  # Find markers conserved in this cluster across all samples
  markers_conserved <- FindConservedMarkers(
    seurat_filtered,
    ident.1 = cluster,          # Which cluster to analyze
    grouping.var = "meta_sample",  # Group by sample
    assay = "RNA",
    logfc.threshold = log(2)
  )
  
  if (nrow(markers_conserved) > 0) {
    markers_conserved$gene <- rownames(markers_conserved)
    markers_conserved$cluster <- cluster
    conserved_markers_list[[as.character(cluster)]] <- markers_conserved
  }
}

all_conserved_markers <- bind_rows(conserved_markers_list)
all_conserved_markers <- all_conserved_markers %>%
  rowwise() %>%
  dplyr::mutate(max_pval = max(c_across(ends_with("_p_val_adj")))) %>%
  filter(max_pval < 0.05)

write.csv(all_conserved_markers, 
          file.path(working_dir, "markers/2_conserved_markers_by_cluster_pval.csv"), 
          row.names = FALSE, quote = FALSE)
all_conserved_markers <- read.csv(file.path(working_dir, "markers/2_conserved_markers_by_cluster_pval.csv"),stringsAsFactors = FALSE)

all_conserved_markers <- all_conserved_markers %>%
  rowwise() %>%
  dplyr::mutate(min_logfc = min(c_across(ends_with("_avg_log2FC")))) %>%
  arrange(cluster, desc(min_logfc))

write.csv(all_conserved_markers, 
          file.path(working_dir, "markers/2_conserved_markers_by_cluster_fc.csv"), 
          row.names = FALSE, quote = FALSE)
all_conserved_markers <- read.csv(file.path(working_dir, "markers/2_conserved_markers_by_cluster_fc.csv"),stringsAsFactors = FALSE)

top_conserved_genes <- all_conserved_markers %>%
  dplyr::mutate(cluster = as.numeric(cluster)) %>% 
  arrange(cluster) %>%                      
  group_by(cluster) %>%
  top_n(n = 5, wt = min_logfc) %>%
  pull(gene) %>%
  unique()

p <- DoHeatmap(
  seurat_filtered, 
  features = top_conserved_genes, 
  group.by = "seurat_clusters",  
  assay = "integrated", 
  angle = 0, 
  group.colors = custom_colors$discrete, 
  size = 3, 
  group.bar.height = 0.01, 
  hjust = 0.5
) +
  ggtitle("Conserved Markers Across Samples") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_viridis(option = "B")

ggsave(
  file.path(working_dir, 'final/doheatmap.2_conserved_markers.int.png'),
  p, 
  height = 10, 
  width = 12
)

markers <- FindConservedMarkers(seurat_filtered, assay = "RNA", grouping.var = "merged_sample", logfc.threshold = log(2))
markers <- subset(markers, p_val_adj < 0.05)
colnames(markers)[colnames(markers) == "cluster"] <- "sample"
markers_sorted <- markers %>%
  arrange(sample, desc(avg_log2FC))

write.csv(markers, file.path(working_dir, paste0("markers/markers_by_mergedsample_pval.csv")), row.names = FALSE, quote = FALSE)
write.csv(markers_sorted, file.path(working_dir, paste0("markers/markers_by_mergedsample_fc.csv")), row.names = FALSE, quote = FALSE)

markers <- read.csv(file.path(working_dir, "markers/markers_by_mergedsample_pval.csv"),stringsAsFactors = FALSE)
markers <- read.csv(file.path(working_dir, "markers/markers_by_mergedsample_fc.csv"), stringsAsFactors = FALSE)


n_values <- c(1, 3, 5, 10)

for (n in n_values) {
  top_genes <- markers %>%
    group_by(sample) %>%
    top_n(n = n, wt = avg_log2FC) %>%
    pull(gene)
  
  dot_plot <- DotPlot(seurat_filtered, features = unique(top_genes)) +
    ggtitle(paste("Top", n, "Genes by Cluster")) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    theme(plot.background = element_rect(fill = "white", color = NA),
          panel.grid.major = element_line(color = "grey90"))
  
  ggsave(file.path(working_dir, paste("markers/dot_plot_by_mergedsample_fc_top_", n, ".png", sep = "")), plot = dot_plot, width = max(16, min(length(unique(Idents(seurat))) * n / 4, 40)), height = 8, dpi = 300)
}

# 
top_genes_heatmap <- markers %>%
  group_by(sample) %>%
  top_n(n = 5, wt = avg_log2FC) %>%
  pull(gene) %>%
  unique()

p <- DoHeatmap(seurat_filtered, features = top_genes_heatmap, group.by = "combined_sample", angle = 0, group.colors = custom_colors$discrete, size = 3, group.bar.height = 0.01, hjust = 0.5) +
  ggtitle("Top 5 Markers per Sample") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_viridis(option = "B")

ggsave(file.path(working_dir, 'doheatmap.sample.single.png'), p, height = 16, width = 8)

## by meta sample, integrated
markers <- FindAllMarkers(seurat, assay = "integrated", group.by = "meta_sample", slot = "scale.data", logfc.threshold = log(2))
markers <- subset(markers, p_val_adj < 0.05)
colnames(markers)[colnames(markers) == "cluster"] <- "meta_sample"

markers_sorted <- markers %>%
  arrange(meta_sample, desc(avg_log2FC))

write.csv(markers, file.path(working_dir, paste0("markers/markers_by_metasample_pval.csv")), row.names = FALSE, quote = FALSE)
write.csv(markers_sorted, file.path(working_dir, paste0("markers/markers_by_metasample_fc.csv")), row.names = FALSE, quote = FALSE)

n_values <- c(1, 3, 5, 10)

for (n in n_values) {
  top_genes <- markers %>%
    group_by(cluster) %>%
    top_n(n = n, wt = avg_log2FC) %>%
    pull(gene)
  
  dot_plot <- DotPlot(seurat_filtered, features = unique(top_genes)) +
    ggtitle(paste("Top", n, "Genes by Meta Sample")) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    theme(plot.background = element_rect(fill = "white", color = NA),
          panel.grid.major = element_line(color = "grey90"))
  
  ggsave(file.path(working_dir, paste("markers/dot_plot_integrated_by_metasample_top_", n, ".png", sep = "")), plot = dot_plot, width = max(16, min(length(unique(Idents(seurat))) * n / 4, 40)), height = 8, dpi = 300)
}

top_genes_heatmap <- markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) %>%
  pull(gene) %>%
  unique()

p <- DoHeatmap(seurat_filtered, features = top_genes_heatmap, group.by = "meta_sample", angle = 0, group.colors = custom_colors$discrete, size = 3, group.bar.height = 0.01, hjust = 0.5) +
  ggtitle("Top 5 Markers per Meta Sample") +
  theme(plot.title = element_text(hjust = 0.5)) + 
  scale_fill_viridis(option = "B")

ggsave(file.path(working_dir, 'doheatmap.metasample.single.png'), p, height = 16, width = 8)

## by supersample
markers <- FindAllMarkers(seurat_filtered, assay = "integrated", group.by = "combined_sample", slot = "scale.data", logfc.threshold = log(2))
markers <- subset(markers, p_val_adj < 0.05)
colnames(markers)[colnames(markers) == "cluster"] <- "super_sample"

markers_sorted <- markers %>%
  arrange(super_sample, desc(avg_log2FC))

write.csv(markers, file.path(working_dir, paste0("markers/markers_by_supersample_pval.csv")), row.names = FALSE, quote = FALSE)
write.csv(markers_sorted, file.path(working_dir, paste0("markers/markers_by_supersample_fc.csv")), row.names = FALSE, quote = FALSE)

## by cluster
markers <- FindAllMarkers(seurat, assay = "integrated", min.diff.pct = 0.1, logfc.threshold = log(2))
markers <- subset(markers, p_val_adj < 0.05)

markers_sorted <- markers %>%
  arrange(cluster, desc(avg_log2FC))

write.csv(markers, file.path(working_dir, paste0("markers/markers_by_cluster_pval_all.csv")), row.names = FALSE, quote = FALSE)
write.csv(markers_sorted, file.path(working_dir, paste0("markers/markers_by_cluster_fc_all.csv")), row.names = FALSE, quote = FALSE)

c(
  "Pou-r (POU6F1)", "DDX56", "MEX3C", "Myc (MYCN)", "ALPL.4", "NOP56", "NOP58", "HSPE1", # 0
  "CD53", "NFIA", "SPON1", "SHANK2", # 1
  "NRG1", "THBD", "CD34", # 2
  "SEMA3E", "EFHD2", # 3
  "EYA1", # 4
  "FLT3.3", "CD163", # 5
  # 6
  "CYBB", "SHISAL2B", # 8
  "Serum (SRF)", "SPHK1", # 9
  "CTSA", "MALRD1", "ITGA2", # 10
  # 11
  # 12
  "GBP2", "MIF", "CTSL", "VIL1", # 13
  "ACSL1.2", # 14
  "PAX6", # 15
  # 16
  "PTPN23", # 17
  "TMSB4X", # 18
  "SLC30A3", "SNAP25" # 19
)

markers <- read.csv(file.path(working_dir, "markers/integrated_by_cluster.csv"))

n_values <- c(1, 3, 5, 10)

for (n in n_values) {
  top_genes <- markers %>%
    group_by(cluster) %>%
    top_n(n = n, wt = avg_log2FC) %>%
    pull(gene)
  
  dot_plot <- DotPlot(seurat, features = unique(top_genes)) +
    ggtitle(paste("Top", n, "Genes by Cluster")) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    theme(plot.background = element_rect(fill = "white", color = NA),
          panel.grid.major = element_line(color = "grey90"))
  
  ggsave(file.path(working_dir, paste("markers/dot_plot_top_", n, ".png", sep = "")), plot = dot_plot, width = max(16, min(length(unique(Idents(seurat))) * n / 4, 40)), height = 8, dpi = 300)
}

top_genes_heatmap <- markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC) %>%
  pull(gene) %>%
  unique()

seurat <- NormalizeData(seurat, assay = "RNA")
seurat <- ScaleData(seurat, assay = "RNA")

p <- DoHeatmap(seurat, features = top_genes_heatmap, angle = 0, group.colors = custom_colors$discrete, size = 3, group.bar.height = 0.01, label = FALSE) +
  ggtitle("Top 10 Markers per Cluster") +
  theme(plot.title = element_text(hjust = 0.5)) + 
  scale_fill_viridis(option = "E")

p_rna <- DoHeatmap(seurat, features = top_genes_heatmap, angle = 0, group.colors = custom_colors$discrete, size = 3, group.bar.height = 0.01, label = FALSE, assay = 'RNA') +
  ggtitle("Top 10 Markers per Cluster") +
  theme(plot.title = element_text(hjust = 0.5)) + 
  scale_fill_viridis(option = "E")

p + p_rna + plot_layout(ncol = 2)

ggsave(file.path(working_dir, 'doheatmap.cluster.single.png'), p, height = 16, width = 8)

# p <- dittoHeatmap(seurat, top_genes_heatmap,
#                   annot.by = c("seurat_clusters"),
#                   order.by = c("seurat_clusters"),
#                   heatmap.colors = viridis::viridis(20, option = "E"),
#                   annot.colors = c(custom_colors$discrete[1:20], custom_colors$discrete[31:40]), 
#                   breaks = seq(min(seurat_filtered[["integrated"]]$scale.data[top_genes_heatmap, , drop = FALSE]), max(seurat_filtered[["integrated"]]$scale.data[top_genes_heatmap, , drop = FALSE]), length.out = 20))
# 
# ggsave(file.path(working_dir, 'dittoheatmap.cluster.single.png'), p, height = 16, width = 8)

## multisampling using scillus
## data preparation
seurat_filtered <- subset(seurat, subset = sample != "Undetermined")
seurat_filtered@meta.data$sample <- droplevels(seurat_filtered@meta.data$sample)
seurat_filtered@meta.data$combined_sample <- droplevels(seurat_filtered@meta.data$combined_sample)
seurat_filtered@meta.data$meta_sample <- droplevels(seurat_filtered@meta.data$meta_sample)

# ordered_clusters <- table_clusters_by_samples %>%
#   arrange(desc(STEM)) %>%
#   pull(cluster)

ordered_clusters <- factor(c(0, 1, 2, 5, 11, 10, 7, 6, 4, 16, 15, 8, 13, 18, 17, 19, 12, 14, 3, 9), levels = 0:19)

seurat_filtered@meta.data$seurat_clusters <- factor(
  seurat_filtered@meta.data$seurat_clusters,
  levels = ordered_clusters
)

by_sample_markers <- read.csv(file.path(working_dir, "markers/integrated_by_metasample.csv"))
by_sample_markers <- read.csv(file.path(working_dir, "markers/markers_by_metasample_fc.csv"))

top_genes_by_sample <- by_sample_markers %>%
  group_by(meta_sample) %>%
  top_n(n = 50, wt = avg_log2FC) %>%
  pull(gene) %>%
  unique()

by_cluster_markers <- read.csv(file.path(working_dir, "markers/integrated_by_cluster.csv"))
by_cluster_markers <- read.csv(file.path(working_dir, "markers/markers_by_cluster_fc.csv"))

by_sample_markers_renamed <- by_sample_markers %>%
  rename_with(~ ifelse(grepl("sample$", .), ., paste0(., "_sample")), -gene)

# Rename overlapping columns in by_cluster_markers
by_cluster_markers_renamed <- by_cluster_markers %>%
  rename_with(~ ifelse(grepl("cluster$", .), ., paste0(., "_cluster")), -gene)

merged_markers <- full_join(by_sample_markers_renamed, by_cluster_markers_renamed, by = "gene", relationship = "many-to-many")

merged_markers_sorted <- merged_markers %>%
  arrange(meta_sample, cluster, desc(avg_log2FC_sample))

write.csv(markers, file.path(working_dir, paste0("markers/markers_by_cluster+metasample_pval.csv")), row.names = FALSE, quote = FALSE)
write.csv(markers_sorted, file.path(working_dir, paste0("markers/markers_by_cluster+metasample_fc.csv")), row.names = FALSE, quote = FALSE)

merged_markers <- merged_markers %>%
  mutate(
    avg_log2FC_max = pmax(avg_log2FC_sample, avg_log2FC_cluster, na.rm = TRUE),
    avg_log2FC_min = pmin(avg_log2FC_sample, avg_log2FC_cluster, na.rm = TRUE),
    avg_log2FC_mean = (avg_log2FC_sample + avg_log2FC_cluster) / 2
  )

write.csv(merged_markers, file.path(working_dir, paste0("markers/merged_by_cluster_metasample_.csv")), row.names = FALSE, quote = FALSE)

### supersample+cluster
by_sample_markers <- read.csv(file.path(working_dir, "markers/markers_by_supersample_fc.csv"))
by_cluster_markers <- read.csv(file.path(working_dir, "markers/markers_by_cluster_fc.csv"))

by_sample_markers_renamed <- by_sample_markers %>%
  rename_with(~ ifelse(grepl("sample$", .), ., paste0(., "_sample")), -gene)

# Rename overlapping columns in by_cluster_markers
by_cluster_markers_renamed <- by_cluster_markers %>%
  rename_with(~ ifelse(grepl("cluster$", .), ., paste0(., "_cluster")), -gene)

merged_markers <- full_join(by_sample_markers_renamed, by_cluster_markers_renamed, by = "gene", relationship = "many-to-many")

merged_markers <- merged_markers[, c("gene", "super_sample", "cluster", setdiff(names(merged_markers), c("gene", "super_sample", "cluster")))]

# Rename the column
names(merged_markers)[names(merged_markers) == "super_sample"] <- "sample"

merged_markers <- merged_markers %>%
  arrange(sample, cluster, desc(avg_log2FC_sample))

# merged_markers_sorted <- merged_markers %>%
#   arrange(cluster, sample, desc(avg_log2FC_sample))

write.csv(merged_markers, file.path(working_dir, paste0("markers/markers_by_cluster+sample_fc.csv")), row.names = FALSE, quote = FALSE)
# write.csv(merged_markers_sorted, file.path(working_dir, paste0("markers/markers_by_cluster+sample_fc_2.csv")), row.names = FALSE, quote = FALSE)

# merged_markers <- merged_markers %>%
#   mutate(
#     avg_log2FC_max = pmax(avg_log2FC_sample, avg_log2FC_cluster, na.rm = TRUE),
#     avg_log2FC_min = pmin(avg_log2FC_sample, avg_log2FC_cluster, na.rm = TRUE),
#     avg_log2FC_mean = (avg_log2FC_sample + avg_log2FC_cluster) / 2
#   )
# 
# write.csv(merged_markers, file.path(working_dir, paste0("markers/merged_by_cluster_metasample_.csv")), row.names = FALSE, quote = FALSE)

top_genes_heatmap <- merged_markers %>%
  filter(!is.na(meta_sample)) %>% 
  mutate(meta_sample = factor(meta_sample, levels = c("STEM", "CONTROL"))) %>% 
  mutate(cluster = factor(cluster, levels = ordered_clusters)) %>% 
  # filter(cluster %in% ordered_clusters[c(1, 4)]) %>%
  arrange(meta_sample, cluster) %>%
  # arrange(cluster) %>% 
  group_by(cluster) %>% 
  # top_n(n = 5, wt = avg_log2FC_min) %>%
  top_n(n = 10, wt = avg_log2FC_mean) %>%
  # top_n(n = 5, wt = avg_log2FC_max) %>%
  pull(gene) %>% 
  unique()

p <- dittoHeatmap(seurat_filtered, top_genes_heatmap,
                  annot.by = c("cluster", "meta_sample"),
                  order.by = c("meta_sample", "cluster"),
                  main = "Top Dual-Markers for Stem-Samples per Cluster",
                  heatmap.colors = viridis::viridis(10, option = "E"),
                  annot.colors = c(custom_colors$discrete[as.numeric(ordered_clusters)], custom_colors$discrete[37:38]), 
                  breaks = seq(min(seurat_filtered[["integrated"]]$scale.data[top_genes_heatmap, , drop = FALSE]), max(seurat_filtered[["integrated"]]$scale.data[top_genes_heatmap, , drop = FALSE]), length.out = 10))

ggsave(file.path(working_dir, 'dittoheatmap.single.filtered.mean.png'), p, height = 12, width = 12)
ggsave(file.path(working_dir, 'dittoheatmap.single.min.png'), p, height = 12, width = 12)
ggsave(file.path(working_dir, 'dittoheatmap.single.mean.png'), p, height = 20, width = 12)
ggsave(file.path(working_dir, 'dittoheatmap.single.max.png'), p, height = 12, width = 12)

p <- DoHeatmap(seurat_filtered, 
               features = top_genes_heatmap,
               # features = top_genes_by_sample,
               group.by = "cluster", 
               angle = 0, 
               group.colors = custom_colors$discrete[as.numeric(ordered_clusters)], 
               size = 3, 
               group.bar.height = 0.01, 
               hjust = 0.5, 
               label = FALSE) +
  # ggtitle("Top Markers for the Stem Cluster 0, 5") + # there could be a better name for it
  ggtitle("Top Dual-Markers for Stem-Samples per Cluster") + # there could be a better name for it
  # ggtitle("Top Markers for Stem-Samples") +
  theme(plot.title = element_text(hjust = 0.5)) + 
  scale_fill_viridis(option = "E")

ggsave(file.path(working_dir, 'doheatmap.single.filtered.mean.png'), p, height = 12, width = 12) # features = top_genes_heatmap, filter(cluster %in% ordered_clusters[1:7]) %>% top_n(n = 15, wt = avg_log2FC_mean)
ggsave(file.path(working_dir, 'doheatmap.single.mean.10.png'), p, height = 16, width = 12) # features = top_genes_heatmap, top_n(n = 10, wt = avg_log2FC_mean)
ggsave(file.path(working_dir, 'doheatmap.single.mean.7.png'), p, height = 16, width = 12) # features = top_genes_heatmap, top_n(n = 7, wt = avg_log2FC_mean)
ggsave(file.path(working_dir, 'doheatmap.metasample.single.png'), p, height = 12, width = 12) # features = top_genes_by_sample, top_n(n = 50, wt = avg_log2FC_mean)
ggsave(file.path(working_dir, 'doheatmap.cytotrace.single.png'), p, height = 12, width = 12) # features = top_genes_heatmap, filter(cluster %in% ordered_clusters[1:1]) %>% top_n(n = 100, wt = avg_log2FC_max)
ggsave(file.path(working_dir, 'doheatmap.cytotrace.0,2.single.png'), p, height = 12, width = 12) # features = top_genes_heatmap, filter(cluster %in% ordered_clusters[c(1, 4)]) %>% top_n(n = 50, wt = avg_log2FC_max)

# # png(file.path(working_dir, 'scillus.single.min.png'), height = 5400, width = 2400, res = 300)
png(file.path(working_dir, 'scillus.single.mean.png'), height = 8000, width = 2400, res = 300)
# # png(file.path(working_dir, 'scillus.single.mean.filtered.png'), height = 5400, width = 2400, res = 300)
# # png(file.path(working_dir, 'scillus.single.max.png'), height = 5400, width = 2400, res = 300)
plot_heatmap(dataset = seurat_filtered,
             top_genes_heatmap,
             sort_var = c("meta_sample", "cluster"),
             # sort_var = c("seurat_clusters", "meta_sample"),
             anno_var = c("meta_sample", "cluster"),
             anno_colors = list(custom_colors$discrete[21:22], custom_colors$discrete[1:20]),
             hm_colors = viridis::viridis(3, option = "E"))
dev.off()

### pseudobulk correlation
## stemb technical replicate
# layers_to_compare <- c("counts.2", "counts.6")
layers_to_compare <- c("STEM_B_1", "STEM_B_2")
# layers_to_compare <- c("STEM_A", "CONTROL")
pseudobulk1 <- rowSums(seurat@assays$RNA$counts[, seurat@meta.data$sample == layers_to_compare[1]]) # STEM_B_1
pseudobulk2 <- rowSums(seurat@assays$RNA$counts[, seurat@meta.data$sample == layers_to_compare[2]]) # STEM_B_2

# Calculate Pearson and Spearman correlation values
pearson_value <- cor(pseudobulk1, pseudobulk2, method = "pearson", use = "complete.obs")
spearman_value <- cor(pseudobulk1, pseudobulk2, method = "spearman", use = "complete.obs")

pseudobulk_summary <- data.frame(
  Gene = rownames(seurat@assays$RNA$counts),  # Assuming these are the gene names
  setNames(data.frame(pseudobulk1, pseudobulk2), layers_to_compare)
)

# Create the scatter plot
p <- ggplot(pseudobulk_summary, aes(x = (get(layers_to_compare[1]) + 1), y = (get(layers_to_compare[2]) + 1))) +
  geom_point(size = 0.5, alpha = 0.5, color = "grey30") +
  geom_smooth(method = "lm", color = adjustcolor(custom_colors$discrete[1], alpha.f = 0.7), se = FALSE) +
  scale_x_log10(name = paste(layers_to_compare[1], 'Gene count'), labels = scales::comma) +
  scale_y_log10(name = paste(layers_to_compare[2], 'Gene count'), labels = scales::comma) +
  labs(title = paste0("Pseudobulk Correlation: ", layers_to_compare[1], ' (n = ', ncol(seurat@assays$RNA$counts[, seurat@meta.data$sample == layers_to_compare[1]]), '), ', layers_to_compare[2], ' (n = ', ncol(seurat@assays$RNA$counts[, seurat@meta.data$sample == layers_to_compare[2]]), ')'), subtitle = paste("Pearson:", round(pearson_value, 2), "| Spearman:", round(spearman_value, 2))) +
  theme_bw()

ggsave(file.path(working_dir, 'pseudobulk.stema.png'), p, height = 6, width = 8)

## Pre-trial on proper values of k and b
layers_to_compare <- c("ALDH", "STEM_B_1", "STEM_B_2", "CONTROL", "STEM_A", "Undetermined")

k_values <- seq(0.75, 0.95, length.out = 11)
b_values <- 10^(seq(0, 2, length.out = 11))

# Initialize a list to store heatmap plots for each comparison
heatmap_plots <- list()

for (i in 1:length(layers_to_compare)) {
  for (j in 1:i) {
    if (i != j) {
      # Initialize a matrix to store proportions for the current comparison
      proportions_matrix <- matrix(NA, nrow = length(k_values), ncol = length(b_values))
      
      # Get counts for each sample
      pseudobulk1 <- rowSums(seurat[['RNA']]$counts[, seurat@meta.data$sample == layers_to_compare[i]]) 
      pseudobulk2 <- rowSums(seurat[['RNA']]$counts[, seurat@meta.data$sample == layers_to_compare[j]])
      
      total_genes <- length(pseudobulk1)
      
      # Loop over k and b values to calculate proportions
      for (k_index in 1:length(k_values)) {
        for (b_index in 1:length(b_values)) {
          k <- k_values[k_index]
          b <- b_values[b_index]
          
          # Define the inside region for the current k and b
          outside_region <- pseudobulk2 > (k * (pseudobulk1 - b)) & pseudobulk2 < (1 / k * pseudobulk1 + b)
          
          # Count the proportion of genes outside the region
          outside_count <- sum(!outside_region, na.rm = TRUE)
          proportions_matrix[k_index, b_index] <- outside_count / total_genes
        }
      }
      
      # Create a heatmap for the current comparison
      heatmap_df <- melt(proportions_matrix)
      colnames(heatmap_df) <- c("k_index", "b_index", "proportion")
      
      # Create a column for marking with distinct labels
      heatmap_df$mark <- ifelse(heatmap_df$proportion < 0.05, "Under 5%",
                                ifelse(heatmap_df$proportion < 0.10, "Under 10%", NA))
      
      # Plotting the heatmap with Viridis color scale
      heatmap_plot <- ggplot(heatmap_df, aes(x = b_values[b_index], y = k_values[k_index], fill = proportion)) +
        geom_tile(color = "white") +  # Add white borders to tiles
        scale_fill_viridis(option = "C", name = "Proportion Outside", na.value = "grey50") +
        scale_x_log10(labels = scales::comma) +  # Log10 scale for x-axis
        labs(x = "b values", y = "k values", title = paste(layers_to_compare[i], "vs", layers_to_compare[j])) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
      
      heatmap_plot <- heatmap_plot +
        geom_text(data = subset(heatmap_df, !is.na(mark) & mark == "Under 5%"),
                  aes(label = "⁑"), color = "white", size = 3.5, vjust = 0.7) +
        geom_text(data = subset(heatmap_df, !is.na(mark) & mark == "Under 10%"),
                  aes(label = "∗"), color = "white", size = 2.5, vjust = 0.5)
      
      # Store the plot in the list
      heatmap_plots[[paste0(layers_to_compare[i], "_vs_", layers_to_compare[j])]] <- heatmap_plot
    }
  }
}

# Combine all heatmaps into one large graph
combined_heatmap <- grid.arrange(grobs = heatmap_plots, ncol = 4)

# Save the combined heatmap as a single image
ggsave(file.path(working_dir, "pseudobulk.heatmaps.single.png"), plot = combined_heatmap, width = 16, height = 16)

rownames(seurat[['integrated']])[grep("(?i)DDX56", rownames(seurat[['integrated']]))]

## mutual comparison pseudobulk
{
  runid <- "0729_stem"
  query <- "4"
  
  if (query == "3") {
    layers_to_compare <- levels(seurat$combined_sample)
  } else if (query == "4") {
    layers_to_compare <- levels(seurat$merged_sample)
  }
  
  highlight_genes <- c("ALPL.4", "CD53", "FLT3.3", "DDX56", "Pou-r (POU6F1)", "MEX3C", "Myc (MYCN)", "NOP56", "NOP58", "HSPE1", "CD53", "NFIA", "SPON1", "SHANK2", "NRG1", "THBD", "CD34")
  
  # Debug: Check what genes are in your data
  all_genes <- rownames(seurat[['RNA']])
  message("Total genes in dataset: ", length(all_genes))
  message("First 10 genes in dataset: ", paste(head(all_genes, 10), collapse = ", "))
  
  # Check which highlight genes exist in the data
  found_genes <- highlight_genes[highlight_genes %in% all_genes]
  missing_genes <- highlight_genes[!highlight_genes %in% all_genes]
  
  message("Highlight genes found in data: ", length(found_genes))
  message("Found genes: ", paste(found_genes, collapse = ", "))
  message("Missing genes: ", paste(missing_genes, collapse = ", "))
  
  k <- 0.85
  b <- 15
  
  plots <- list()  # Initialize the list to store plots
  outside_info <- list()  # Initialize the list to store outside information
  
  plotname <- paste0('pseudobulk.RNA.merged.single.', runid, '.pdf')
  title <- 'Pseudobulk gene count comparison normalised by #cells in each sample'
  
  for (i in 1:length(layers_to_compare)) {
    for (j in 1:i) {
      if (i != j) {
        if (query == "3") {
          mask1 <- seurat@meta.data$combined_sample == layers_to_compare[i]
          mask2 <- seurat@meta.data$combined_sample == layers_to_compare[j]
        } else if (query == "4") {
          mask1 <- seurat@meta.data$merged_sample == layers_to_compare[i]
          mask2 <- seurat@meta.data$merged_sample == layers_to_compare[j]
        }
        pseudobulk1 <- rowMeans(seurat[['RNA']]$counts[, mask1, drop = FALSE]) * median(rowSums(seurat[['RNA']]$counts > 0))
        pseudobulk2 <- rowMeans(seurat[['RNA']]$counts[, mask2, drop = FALSE]) * median(rowSums(seurat[['RNA']]$counts > 0))
        Genes <- rownames(seurat[['RNA']])
        
        # Calculate Pearson and Spearman correlation values
        pearson_value <- cor(pseudobulk1, pseudobulk2, method = "pearson", use = "complete.obs")
        spearman_value <- cor(pseudobulk1, pseudobulk2, method = "spearman", use = "complete.obs")
        
        # Create the summary data frame and add outside_region directly
        pseudobulk_summary <- data.frame(
          Gene = Genes,
          pseudobulk1 = pseudobulk1,
          pseudobulk2 = pseudobulk2,
          outside_region = !(pseudobulk2 > (k * (pseudobulk1 - b)) & pseudobulk2 < (1 / k * pseudobulk1 + b)) & !(pseudobulk2 > (1 / k * pseudobulk1 - b) & pseudobulk2 < (k * (pseudobulk1 + b)))
        )
        
        # Sort by pseudobulk2 in decreasing order
        pseudobulk_summary <- pseudobulk_summary[order(-pseudobulk_summary$pseudobulk2), ]
        
        # Write the sorted data frame to a CSV file
        write.csv(pseudobulk_summary, file = file.path(working_dir, paste0("pseudobulk/", gsub("\\.pdf", "_", plotname), layers_to_compare[i], "_vs_", layers_to_compare[j], ".csv")), row.names = FALSE, quote = FALSE)
        
        highlight_data <- pseudobulk_summary[pseudobulk_summary$Gene %in% highlight_genes, ]
        
        # Debug information for each comparison
        message("\n--- Comparison: ", layers_to_compare[j], " vs ", layers_to_compare[i], " ---")
        message("Total genes in pseudobulk_summary: ", nrow(pseudobulk_summary))
        message("Highlight genes found in this comparison: ", nrow(highlight_data))
        
        if (nrow(highlight_data) > 0) {
          message("Highlight genes in this comparison: ", paste(highlight_data$Gene, collapse = ", "))
          message("Their expression values (x,y): ")
          for (k in 1:nrow(highlight_data)) {
            message("  ", highlight_data$Gene[k], ": (", round(highlight_data$pseudobulk1[k], 2), ", ", round(highlight_data$pseudobulk2[k], 2), ")")
          }
        }
        
        # Calculate proportions
        total_genes <- length(pseudobulk1)
        outside_count <- sum(pseudobulk_summary$outside_region, na.rm = TRUE)  # Count TRUE values for outside
        inside_count <- total_genes - outside_count
        outside_proportion <- (outside_count / total_genes) * 100
        inside_proportion <- (inside_count / total_genes) * 100
        
        if (grepl("RNA", plotname)) {
          p <- ggplot(pseudobulk_summary, aes(x = (pseudobulk1 + 1), y = (pseudobulk2 + 1), color = factor(outside_region))) +
            geom_point(size = 0.5, alpha = 0.5) +
            scale_x_log10(name = paste0(layers_to_compare[i], " Expression"), labels = scales::comma) +
            scale_y_log10(name = paste0(layers_to_compare[j], " Expression"), labels = scales::comma) +
            labs(title = paste(layers_to_compare[j], "vs", layers_to_compare[i]), 
                 subtitle = paste("Pearson:", round(pearson_value, 2), " | Spearman:", round(spearman_value, 2))) +
            theme_bw() +
            scale_color_manual(values = c("FALSE" = custom_colors$discrete[23],  # Inside
                                          "TRUE" = custom_colors$discrete[5]),  # Outside
                               labels = c(paste0("Unchanged (", round(inside_proportion, 1), "%)"),
                                          paste0("Deviated (", round(outside_proportion, 1), "%)"))) +
            guides(color = guide_legend(title = NULL)) +
            theme(legend.position = "bottom",
                  plot.margin = margin(5, 5, 5, 5, "mm"))
          
          # Debug: Always try to add labels if there are highlight genes
          if (nrow(highlight_data) > 0) {
            message("Adding ", nrow(highlight_data), " gene labels to plot")
            p <- p + geom_text(data = highlight_data, 
                               aes(x = (pseudobulk1 + 1), 
                                   y = (pseudobulk2 + 1), 
                                   label = Gene), 
                               color = "red", 
                               size = 4, 
                               fontface = 'bold',
                               vjust = -0.5, 
                               hjust = 0.5,
                               check_overlap = FALSE)  # Changed to FALSE to see all labels
          } else {
            message("No highlight genes found for labeling")
          }
        } else if (grepl("integrated", plotname)) {
          p <- ggplot(pseudobulk_summary, aes(x = (pseudobulk1 + 1), y = (pseudobulk2 + 1), color = factor(outside_region))) +
            geom_point(size = 0.5, alpha = 0.5) +
            scale_x_continuous(name = paste0(layers_to_compare[i], " Expression"), labels = scales::comma) +
            scale_y_continuous(name = paste0(layers_to_compare[j], " Expression"), labels = scales::comma) +
            labs(title = paste(layers_to_compare[j], "vs", layers_to_compare[i]), 
                 subtitle = paste("Pearson:", round(pearson_value, 2), " | Spearman:", round(spearman_value, 2))) +
            theme_bw() +
            scale_color_manual(values = c("FALSE" = custom_colors$discrete[23],  # Inside
                                          "TRUE" = custom_colors$discrete[5]),  # Outside
                               labels = c(paste0("Unchanged (", round(inside_proportion, 1), "%)"),
                                          paste0("Deviated (", round(outside_proportion, 1), "%)"))) +
            guides(color = guide_legend(title = NULL)) +
            theme(legend.position = "bottom",
                  plot.margin = margin(5, 5, 5, 5, "mm"))
          
          if (nrow(highlight_data) > 0) {
            message("Adding ", nrow(highlight_data), " gene labels to plot")
            p <- p + geom_text(data = highlight_data, 
                               aes(x = (pseudobulk1 + 1), 
                                   y = (pseudobulk2 + 1), 
                                   label = Gene), 
                               color = "red", 
                               size = 4, 
                               fontface = 'bold',
                               vjust = -0.5, 
                               hjust = 0.5,
                               check_overlap = FALSE)
          } else {
            message("No highlight genes found for labeling")
          }
        }
        
        # Store the plot in the list
        plots[[paste0("plot_", i, "_", j)]] <- p
        
        # Store the count and proportion in the list
        outside_info[[paste0(layers_to_compare[i], "_vs_", layers_to_compare[j])]] <- list(
          count = outside_count,
          proportion = outside_proportion / 100
        )
      }
    }
  }
  
  # Rest of the code remains the same...
  final_plots <- list()
  
  for (i in 2:length(layers_to_compare)) {
    for (j in 1:(i-1)) {
      plot_key <- paste0("plot_", i, "_", j)
      if (!is.null(plots[[plot_key]])) {
        final_plots <- append(final_plots, list(plots[[plot_key]]))
      }
    }
  }
  
  if (length(final_plots) > 0) {
    n_plots <- length(final_plots)
    n_layers <- length(layers_to_compare)
    
    if (n_layers > 1) {
      ncol_grid <- n_layers - 1
      nrow_grid <- ceiling(n_plots / ncol_grid)
    } else {
      ncol_grid <- 1
      nrow_grid <- 1
    }
    
    message(paste("Creating grid with", n_plots, "plots,", nrow_grid, "rows and", ncol_grid, "columns"))
    
    library(cowplot)
    grid_plot <- plot_grid(plotlist = final_plots, ncol = ncol_grid, nrow = nrow_grid)
    
    final_plot <- plot_grid(grid_plot, 
                            ggdraw() + draw_label(title), 
                            ncol = 1, 
                            rel_heights = c(1, 0.05))
    
    ggsave(
      file.path(working_dir, plotname), 
      final_plot, 
      height = 6 * nrow_grid, 
      width = 6 * ncol_grid
    )
    message("Plot saved as PDF: ", file.path(working_dir, plotname))
  } else {
    message("No plots to arrange.")
  }
}

## mutual comparison pseudobulk, functional but with padding and display issues. 
{
  runid <- "0729_stem"
  query <- "4"
  
  if (query == "3") {
    layers_to_compare <- levels(seurat$combined_sample)
  } else if (query == "4") {
    layers_to_compare <- levels(seurat$merged_sample)
  }
  
  highlight_genes <- c("ALPL.4", "CD53", "FLT3.3", "DDX56", "Pou-r (POU6F1)", "MEX3C", "Myc (MYCN)", "NOP56", "NOP58", "HSPE1", "CD53", "NFIA", "SPON1", "SHANK2", "NRG1", "THBD", "CD34")
  
  # Debug: Check what genes are in your data
  all_genes <- rownames(seurat[['RNA']])
  message("Total genes in dataset: ", length(all_genes))
  
  # Check which highlight genes exist in the data
  found_genes <- highlight_genes[highlight_genes %in% all_genes]
  missing_genes <- highlight_genes[!highlight_genes %in% all_genes]
  
  message("Highlight genes found in data: ", length(found_genes))
  message("Found genes: ", paste(found_genes, collapse = ", "))
  message("Missing genes: ", paste(missing_genes, collapse = ", "))
  
  k <- 0.85
  b <- 15
  
  plots <- list()  # Initialize the list to store plots
  outside_info <- list()  # Initialize the list to store outside information
  
  plotname <- paste0('pseudobulk.RNA.merged.single.', runid, '.pdf')
  title <- 'Pseudobulk gene count comparison normalised by #cells in each sample'
  
  for (i in 1:length(layers_to_compare)) {
    for (j in 1:i) {
      if (i != j) {
        if (query == "3") {
          mask1 <- seurat@meta.data$combined_sample == layers_to_compare[i]
          mask2 <- seurat@meta.data$combined_sample == layers_to_compare[j]
        } else if (query == "4") {
          mask1 <- seurat@meta.data$merged_sample == layers_to_compare[i]
          mask2 <- seurat@meta.data$merged_sample == layers_to_compare[j]
        }
        pseudobulk1 <- rowMeans(seurat[['RNA']]$counts[, mask1, drop = FALSE]) * median(rowSums(seurat[['RNA']]$counts > 0))
        pseudobulk2 <- rowMeans(seurat[['RNA']]$counts[, mask2, drop = FALSE]) * median(rowSums(seurat[['RNA']]$counts > 0))
        Genes <- rownames(seurat[['RNA']])
        
        # Calculate Pearson and Spearman correlation values
        pearson_value <- cor(pseudobulk1, pseudobulk2, method = "pearson", use = "complete.obs")
        spearman_value <- cor(pseudobulk1, pseudobulk2, method = "spearman", use = "complete.obs")
        
        # Create the summary data frame and add outside_region directly
        pseudobulk_summary <- data.frame(
          Gene = Genes,
          pseudobulk1 = pseudobulk1,
          pseudobulk2 = pseudobulk2,
          outside_region = !(pseudobulk2 > (k * (pseudobulk1 - b)) & pseudobulk2 < (1 / k * pseudobulk1 + b)) & !(pseudobulk2 > (1 / k * pseudobulk1 - b) & pseudobulk2 < (k * (pseudobulk1 + b)))
        )
        
        # Sort by pseudobulk2 in decreasing order
        pseudobulk_summary <- pseudobulk_summary[order(-pseudobulk_summary$pseudobulk2), ]
        
        # Write the sorted data frame to a CSV file
        write.csv(pseudobulk_summary, file = file.path(working_dir, paste0("pseudobulk/", gsub("\\.pdf", "_", plotname), layers_to_compare[i], "_vs_", layers_to_compare[j], ".csv")), row.names = FALSE, quote = FALSE)
        
        highlight_data <- pseudobulk_summary[pseudobulk_summary$Gene %in% highlight_genes, ]
        
        # Debug information for each comparison
        message("\n--- Comparison: ", layers_to_compare[j], " vs ", layers_to_compare[i], " ---")
        message("Highlight genes found in this comparison: ", nrow(highlight_data))
        
        # Calculate proportions
        total_genes <- length(pseudobulk1)
        outside_count <- sum(pseudobulk_summary$outside_region, na.rm = TRUE)  # Count TRUE values for outside
        inside_count <- total_genes - outside_count
        outside_proportion <- (outside_count / total_genes) * 100
        inside_proportion <- (inside_count / total_genes) * 100
        
        if (grepl("RNA", plotname)) {
          p <- ggplot(pseudobulk_summary, aes(x = (pseudobulk1 + 1), y = (pseudobulk2 + 1), color = factor(outside_region))) +
            geom_point(size = 0.8, alpha = 0.6) +
            scale_x_log10(name = paste0(layers_to_compare[i], " Expression"), labels = scales::comma) +
            scale_y_log10(name = paste0(layers_to_compare[j], " Expression"), labels = scales::comma) +
            labs(title = paste(layers_to_compare[j], "vs", layers_to_compare[i]), 
                 subtitle = paste("Pearson:", round(pearson_value, 2), " | Spearman:", round(spearman_value, 2))) +
            theme_bw() +
            scale_color_manual(values = c("FALSE" = "#7F7F7F",  # Softer gray for unchanged
                                          "TRUE" = "#D62728"),   # Softer red for deviated
                               labels = c(paste0("Unchanged (", round(inside_proportion, 1), "%)"),
                                          paste0("Deviated (", round(outside_proportion, 1), "%)"))) +
            guides(color = guide_legend(title = NULL)) +
            theme(legend.position = "bottom",
                  plot.margin = margin(15, 15, 15, 15, "mm"),
                  text = element_text(size = 10),
                  axis.text = element_text(size = 9),
                  legend.text = element_text(size = 9))
          
          # Add highlighted points with a different color
          if (nrow(highlight_data) > 0) {
            p <- p + geom_point(data = highlight_data, 
                                aes(x = (pseudobulk1 + 1), 
                                    y = (pseudobulk2 + 1)), 
                                color = "#1F77B4",  # Softer blue for highlighted points
                                size = 2.5, 
                                alpha = 0.9,
                                inherit.aes = FALSE)
          }
          
          # Use geom_text with nudging instead of geom_label_repel for grid compatibility
          if (nrow(highlight_data) > 0 && nrow(highlight_data) <= 30) {
            message("Adding ", nrow(highlight_data), " gene labels to plot")
            
            # Create nudged positions to avoid overlap
            highlight_data$nudge_x <- 0.2
            highlight_data$nudge_y <- 0.2
            
            p <- p + 
              # Add background boxes manually
              geom_label(
                data = highlight_data, 
                aes(x = (pseudobulk1 + 1) * 1.2, 
                    y = (pseudobulk2 + 1) * 1.2, 
                    label = Gene), 
                size = 2.8,
                fill = 'white',
                color = 'black',
                fontface = 'bold',
                alpha = 0.9,
                label.size = 0.2,
                inherit.aes = FALSE,
                show.legend = FALSE
              ) +
              # Add connecting segments
              geom_segment(
                data = highlight_data,
                aes(x = (pseudobulk1 + 1), 
                    y = (pseudobulk2 + 1),
                    xend = (pseudobulk1 + 1) * 1.2, 
                    yend = (pseudobulk2 + 1) * 1.2),
                color = 'gray50',
                size = 0.3,
                alpha = 0.7,
                inherit.aes = FALSE
              )
          } else if (nrow(highlight_data) > 15) {
            message("Too many highlight genes (", nrow(highlight_data), "), skipping labels to avoid crowding")
          } else {
            message("No highlight genes found for labeling")
          }
        } else if (grepl("integrated", plotname)) {
          p <- ggplot(pseudobulk_summary, aes(x = (pseudobulk1 + 1), y = (pseudobulk2 + 1), color = factor(outside_region))) +
            geom_point(size = 0.8, alpha = 0.6) +
            scale_x_continuous(name = paste0(layers_to_compare[i], " Expression"), labels = scales::comma) +
            scale_y_continuous(name = paste0(layers_to_compare[j], " Expression"), labels = scales::comma) +
            labs(title = paste(layers_to_compare[j], "vs", layers_to_compare[i]), 
                 subtitle = paste("Pearson:", round(pearson_value, 2), " | Spearman:", round(spearman_value, 2))) +
            theme_bw() +
            scale_color_manual(values = c("FALSE" = "#7F7F7F",  # Softer gray for unchanged
                                          "TRUE" = "#D62728"),   # Softer red for deviated
                               labels = c(paste0("Unchanged (", round(inside_proportion, 1), "%)"),
                                          paste0("Deviated (", round(outside_proportion, 1), "%)"))) +
            guides(color = guide_legend(title = NULL)) +
            theme(legend.position = "bottom",
                  plot.margin = margin(15, 15, 15, 15, "mm"),
                  text = element_text(size = 10),
                  axis.text = element_text(size = 9),
                  legend.text = element_text(size = 9))
          
          # Add highlighted points
          if (nrow(highlight_data) > 0) {
            p <- p + geom_point(data = highlight_data, 
                                aes(x = (pseudobulk1 + 1), 
                                    y = (pseudobulk2 + 1)), 
                                color = "#1F77B4",
                                size = 2.5, 
                                alpha = 0.9,
                                inherit.aes = FALSE)
          }
          
          # Add labels with simple positioning
          if (nrow(highlight_data) > 0 && nrow(highlight_data) <= 15) {
            p <- p + 
              geom_label(
                data = highlight_data, 
                aes(x = (pseudobulk1 + 1) + 0.2, 
                    y = (pseudobulk2 + 1) + 0.2, 
                    label = Gene), 
                size = 2.8,
                fill = 'white',
                color = 'black',
                fontface = 'bold',
                alpha = 0.9,
                label.size = 0.2,
                inherit.aes = FALSE,
                show.legend = FALSE
              ) +
              geom_segment(
                data = highlight_data,
                aes(x = (pseudobulk1 + 1), 
                    y = (pseudobulk2 + 1),
                    xend = (pseudobulk1 + 1) + 0.2, 
                    yend = (pseudobulk2 + 1) + 0.2),
                color = 'gray50',
                size = 0.3,
                alpha = 0.7,
                inherit.aes = FALSE
              )
          }
        }
        
        # Store the plot in the list
        plots[[paste0("plot_", i, "_", j)]] <- p
        
        # Store the count and proportion in the list
        outside_info[[paste0(layers_to_compare[i], "_vs_", layers_to_compare[j])]] <- list(
          count = outside_count,
          proportion = outside_proportion / 100
        )
      }
    }
  }
  
  # Rest of the code remains the same...
  final_plots <- list()
  
  for (i in 2:length(layers_to_compare)) {
    for (j in 1:(i-1)) {
      plot_key <- paste0("plot_", i, "_", j)
      if (!is.null(plots[[plot_key]])) {
        final_plots <- append(final_plots, list(plots[[plot_key]]))
      }
    }
  }
  
  if (length(final_plots) > 0) {
    n_plots <- length(final_plots)
    n_layers <- length(layers_to_compare)
    
    if (n_layers > 1) {
      ncol_grid <- n_layers - 1
      nrow_grid <- ceiling(n_plots / ncol_grid)
    } else {
      ncol_grid <- 1
      nrow_grid <- 1
    }
    
    message(paste("Creating grid with", n_plots, "plots,", nrow_grid, "rows and", ncol_grid, "columns"))
    
    library(cowplot)
    grid_plot <- plot_grid(plotlist = final_plots, ncol = ncol_grid, nrow = nrow_grid)
    
    final_plot <- plot_grid(grid_plot, 
                            ggdraw() + draw_label(title, size = 14), 
                            ncol = 1, 
                            rel_heights = c(1, 0.05))
    
    ggsave(
      file.path(working_dir, plotname), 
      final_plot, 
      height = 7 * nrow_grid, 
      width = 7 * ncol_grid,
      dpi = 300
    )
    message("Plot saved as PDF: ", file.path(working_dir, plotname))
  } else {
    message("No plots to arrange.")
  }
}

## Individual pseudobulk plots
{
  runid <- "0729_stem"
  query <- "4"
  
  if (query == "3") {
    layers_to_compare <- levels(seurat$combined_sample)
  } else if (query == "4") {
    layers_to_compare <- levels(seurat$merged_sample)
  }
  
  rownames(seurat[['RNA']])[grep("(?i)WNT9A", rownames(seurat[['RNA']]))]
  # 0729_stem
  highlight_genes <- c("ALPL.4", "CD53", "FLT3.3", "DDX56", "Pou-r (POU6F1)", "MEX3C", "Myc (MYCN)", "NOP56", "NOP58", "HSPE1", "CD53", "NFIA", "SPON1", "SHANK2", "NRG1", "Gsc (GSC)", "THBD", "HSPA9", "SAP30L", "SART3", "CD34", "LRP8")
  # 0729_stem_all
  # highlight_genes <- c("ALPL.4", "CD53", "FLT3.3", "DDX56", "Pou-r (POU6F1)", "MEX3C", "Myc (MYCN)", "NOP56", "NOP58", "HSPE1", "CD53", "NFIA", "SPON1", "SHANK2", "NRG1", "THBD", "CD34", "Unc.b (GSC)", "NOTCH1.2", "TP53INP1", "Hes. (HES1)",   "Hes. (HES1).2", "Hes. (HES1).3", "Nr6a (NR6A1)", "TP53I13", "Myc (MYCN)", "Cebpa (CEBPD)", "LRP5", "CD34", "SETD1A", "SOX4", "HES1", "NR6A1", "MYSM1", "ABCG2", "ABCG2.3", "ABCG2.4", "ABCG2.6", "ABCG2.8", "ABCG2.9", "ABCB1", "ABCB1.2", "ABCB1.3", "ABCB1.4", "ABCB1.5", "SPON1", "WNT9A", "CUL4A", "CDK6", "SKI", "METTL3", "MTCH2", "SRF", "PBX1", "ARID4A", "BTG2", "SOX18", "NRG1")
  
  # Debug: Check what genes are in your data
  all_genes <- rownames(seurat[['RNA']])
  message("Total genes in dataset: ", length(all_genes))
  
  # Check which highlight genes exist in the data
  found_genes <- highlight_genes[highlight_genes %in% all_genes]
  missing_genes <- highlight_genes[!highlight_genes %in% all_genes]
  
  message("Highlight genes found in data: ", length(found_genes))
  message("Found genes: ", paste(found_genes, collapse = ", "))
  message("Missing genes: ", paste(missing_genes, collapse = ", "))
  
  k <- 0.85
  b <- 10
  
  # Function to clean filenames
  clean_filename <- function(text) {
    # Replace spaces with underscores
    text <- gsub(" ", "_", text)
    # Replace _. with .
    text <- gsub("_\\.", ".", text)
    # Replace ._ with .
    text <- gsub("\\._", ".", text)
    return(text)
  }
  
  base_plotname <- paste0('pseudobulk.RNA.merged.single.', runid)
  
  for (i in 1:length(layers_to_compare)) {
    for (j in 1:i) {
      if (i != j) {
        if (query == "3") {
          mask1 <- seurat@meta.data$combined_sample == layers_to_compare[i]
          mask2 <- seurat@meta.data$combined_sample == layers_to_compare[j]
        } else if (query == "4") {
          mask1 <- seurat@meta.data$merged_sample == layers_to_compare[i]
          mask2 <- seurat@meta.data$merged_sample == layers_to_compare[j]
        }
        pseudobulk1 <- rowMeans(seurat[['RNA']]$counts[, mask1, drop = FALSE]) * median(rowSums(seurat[['RNA']]$counts > 0))
        pseudobulk2 <- rowMeans(seurat[['RNA']]$counts[, mask2, drop = FALSE]) * median(rowSums(seurat[['RNA']]$counts > 0))
        Genes <- rownames(seurat[['RNA']])
        
        # Calculate Pearson and Spearman correlation values
        pearson_value <- cor(pseudobulk1, pseudobulk2, method = "pearson", use = "complete.obs")
        spearman_value <- cor(pseudobulk1, pseudobulk2, method = "spearman", use = "complete.obs")
        
        # Create the summary data frame and add outside_region directly
        pseudobulk_summary <- data.frame(
          Gene = Genes,
          pseudobulk1 = pseudobulk1,
          pseudobulk2 = pseudobulk2,
          outside_region = !(pseudobulk2 > (k * (pseudobulk1 - b)) & pseudobulk2 < (1 / k * pseudobulk1 + b)) & !(pseudobulk2 > (1 / k * pseudobulk1 - b) & pseudobulk2 < (k * (pseudobulk1 + b)))
        )
        
        # Sort by pseudobulk2 in decreasing order
        pseudobulk_summary <- pseudobulk_summary[order(-pseudobulk_summary$pseudobulk2), ]
        
        # Clean layer names for filename
        layer1_clean <- clean_filename(layers_to_compare[i])
        layer2_clean <- clean_filename(layers_to_compare[j])
        
        # Create individual plot filename
        individual_plotname <- paste0(base_plotname, '.', layer2_clean, '_vs_', layer1_clean, '.pdf')
        
        # Write the sorted data frame to a CSV file
        csv_filename <- paste0(base_plotname, '.', layer2_clean, '_vs_', layer1_clean, '.csv')
        write.csv(pseudobulk_summary, file = file.path(working_dir, "pseudobulk", csv_filename), row.names = FALSE, quote = FALSE)
        
        highlight_data <- pseudobulk_summary[pseudobulk_summary$Gene %in% highlight_genes, ]
        
        # Debug information for each comparison
        message("\n--- Comparison: ", layers_to_compare[j], " vs ", layers_to_compare[i], " ---")
        message("Highlight genes found in this comparison: ", nrow(highlight_data))
        message("Saving as: ", individual_plotname)
        
        # Calculate proportions
        total_genes <- length(pseudobulk1)
        outside_count <- sum(pseudobulk_summary$outside_region, na.rm = TRUE)  # Count TRUE values for outside
        inside_count <- total_genes - outside_count
        outside_proportion <- (outside_count / total_genes) * 100
        inside_proportion <- (inside_count / total_genes) * 100
        
        # Create individual plot with larger size and better formatting
        p <- ggplot(pseudobulk_summary, aes(x = (pseudobulk1 + 1), y = (pseudobulk2 + 1), color = factor(outside_region))) +
          geom_point(size = 1.2, alpha = 0.7) +
          scale_x_log10(name = paste0(layers_to_compare[i], " Expression"), labels = scales::comma) +
          scale_y_log10(name = paste0(layers_to_compare[j], " Expression"), labels = scales::comma) +
          labs(title = paste(layers_to_compare[j], "vs", layers_to_compare[i]), 
               subtitle = paste("Pearson:", round(pearson_value, 3), " | Spearman:", round(spearman_value, 3)),
               caption = paste("Total genes:", total_genes, "| Deviated genes:", outside_count, paste0("(", round(outside_proportion, 1), "%)"))) +
          theme_bw() +
          scale_color_manual(values = c("FALSE" = "#7F7F7F",  # Softer gray for unchanged
                                        "TRUE" = "#D62728"),   # Softer red for deviated
                             labels = c(paste0("Unchanged (", round(inside_proportion, 1), "%)"),
                                        paste0("Deviated (", round(outside_proportion, 1), "%)"))) +
          guides(color = guide_legend(title = NULL)) +
          theme(legend.position = "bottom",
                plot.margin = margin(20, 20, 20, 20, "mm"),
                text = element_text(size = 14),
                axis.text = element_text(size = 12),
                legend.text = element_text(size = 12),
                plot.title = element_text(size = 16, face = "bold"),
                plot.subtitle = element_text(size = 14),
                plot.caption = element_text(size = 12, hjust = 0))
        
        # Add highlighted points with a different color
        if (nrow(highlight_data) > 0) {
          p <- p + geom_point(data = highlight_data, 
                              aes(x = (pseudobulk1 + 1), 
                                  y = (pseudobulk2 + 1)), 
                              color = "#1F77B4",  # Softer blue for highlighted points
                              size = 2, 
                              alpha = 0.9,
                              inherit.aes = FALSE)
        }
        
        # Add labels for individual plots with geom_label_repel
        if (nrow(highlight_data) > 0) {
          message("Adding ", nrow(highlight_data), " gene labels to individual plot")
          
          p <- p + 
            # Add labels with repel and black connecting lines
            geom_label_repel(
              data = highlight_data, 
              aes(x = (pseudobulk1 + 1), 
                  y = (pseudobulk2 + 1), 
                  label = Gene), 
              size = 3.5,
              fill = 'white',
              color = 'black',
              fontface = 'bold',
              alpha = 0.9,
              label.size = 0.3,
              box.padding = 0.5,        # Padding around labels
              point.padding = 0.8,      # Padding between labels and points
              segment.color = 'black',  # Black connecting lines
              segment.size = 0.4,
              segment.alpha = 0.8,
              max.overlaps = Inf,
              min.segment.length = 0,   # Show all segments
              force = 2,
              force_pull = 1,
              nudge_x = 0.1,
              nudge_y = 0.1,
              inherit.aes = FALSE,
              show.legend = FALSE
            )
        } else {
          message("No highlight genes found for labeling")
        }
        
        # Save individual plot
        ggsave(
          file.path(working_dir, individual_plotname), 
          p, 
          height = 10, 
          width = 10,
          dpi = 300
        )
        message("Individual plot saved: ", individual_plotname)
      }
    }
  }
  
  message("\nAll individual plots have been saved!")
}

# Print the counts and proportions of genes outside the region
outside_summary <- do.call(rbind, lapply(names(outside_info), function(comparison) {
  data.frame(
    Comparison = comparison,
    Count = outside_info[[comparison]]$count,
    Proportion = round(outside_info[[comparison]]$proportion * 100, 2)
  )
}))

# Print the table with knitr
knitr::kable(outside_summary,
             col.names = c("Comparison", "Genes Outside the Region", "Proportion (%)"),
             format = "pipe")  

# "pseudobulk.RNA.merged.single.png"
# |Comparison              | Genes Outside the Region| Proportion (%)|
# |:-----------------------|------------------------:|--------------:|
# |STEM_B_vs_ALDH          |                      641|           4.53|
# |CONTROL_vs_ALDH         |                     1154|           8.16|
# |CONTROL_vs_STEM_B       |                     1618|          11.44|
# |STEM_A_vs_ALDH          |                     1248|           8.82|
# |STEM_A_vs_STEM_B        |                     1681|          11.88|
# |STEM_A_vs_CONTROL       |                       41|           0.29|
# |Undetermined_vs_ALDH    |                      402|           2.84|
# |Undetermined_vs_STEM_B  |                      617|           4.36|
# |Undetermined_vs_CONTROL |                      969|           6.85|
# |Undetermined_vs_STEM_A  |                     1064|           7.52|

# "pseudobulk.RNA.merged.genewise.single.png"
# |Comparison              | Genes Outside the Region| Proportion (%)|
# |:-----------------------|------------------------:|--------------:|
# |STEM_B_vs_ALDH          |                     1571|          11.11|
# |CONTROL_vs_ALDH         |                     2503|          17.70|
# |CONTROL_vs_STEM_B       |                     2416|          17.08|
# |STEM_A_vs_ALDH          |                     2632|          18.61|
# |STEM_A_vs_STEM_B        |                     2531|          17.89|
# |STEM_A_vs_CONTROL       |                      612|           4.33|
# |Undetermined_vs_ALDH    |                     1765|          12.48|
# |Undetermined_vs_STEM_B  |                     1663|          11.76|
# |Undetermined_vs_CONTROL |                     1244|           8.80|
# |Undetermined_vs_STEM_A  |                     1453|          10.27|

# "pseudobulk.integrated.merged.single.png"
# |Comparison              | Genes Outside the Region| Proportion (%)|
# |:-----------------------|------------------------:|--------------:|
# |STEM_B_vs_ALDH          |                     2627|          87.57|
# |CONTROL_vs_ALDH         |                     2884|          96.13|
# |CONTROL_vs_STEM_B       |                     2949|          98.30|
# |STEM_A_vs_ALDH          |                     2912|          97.07|
# |STEM_A_vs_STEM_B        |                     2956|          98.53|
# |STEM_A_vs_CONTROL       |                     1251|          41.70|
# |Undetermined_vs_ALDH    |                     2844|          94.80|
# |Undetermined_vs_STEM_B  |                     2828|          94.27|
# |Undetermined_vs_CONTROL |                     2947|          98.23|
# |Undetermined_vs_STEM_A  |                     2956|          98.53|

### GSEA
## pooled
# Extract cluster information
clusters <- seurat$seurat_clusters

# Initialize a data frame to store all top genes info
all_top_genes_info <- data.frame(GeneID = character(), GeneName = character(), Cluster = character(), avg_log2FC = numeric(), p_adj = numeric(), stringsAsFactors = FALSE)

# Loop through each cluster's markers to extract top genes
for (cluster_id in as.numeric(levels(seurat$seurat_clusters))) {
  marker_file <- paste0("/scratch/users/jiamuyu/proj_botryllus/scRNAseq/250125_00_ciona_stemcell/markers/", cluster_id, "_markers.csv")
  
  if (file.exists(marker_file)) {
    markers <- read.csv(marker_file)
    
    # Get top genes based on avg_log2FC
    top_genes <- markers %>%
      arrange(desc(avg_log2FC)) %>%
      head(10)  # Adjust the number of top genes as needed
    
    # Append to data frame
    all_top_genes_info <- rbind(all_top_genes_info, 
                                data.frame(GeneID = top_genes$gene_id, 
                                           GeneName = top_genes$gene_name,
                                           avg_log2FC = top_genes$avg_log2FC,
                                           p_adj = top_genes$p_val_adj,
                                           Cluster = as.character(cluster_id)))
  }
}

geneset <- read.gmt("/scratch/users/jiamuyu/proj_botryllus/scRNAseq/c8.all.v2024.1.Hs.symbols.gmt")

for (cluster_id in unique(all_top_genes_info$Cluster)) {
  
  # Filter marker genes for the current cluster
  cluster_genes <- all_top_genes_info %>%
    filter(Cluster == cluster_id) %>%
    arrange(desc(avg_log2FC))
  
  # Create a named vector for GSEA
  genelist <- cluster_genes$avg_log2FC
  names(genelist) <- cluster_genes$GeneName
  
  gsea_results <- GSEA(genelist, TERM2GENE=geneset, pvalueCutoff = 0.05, minGSSize = 1, maxGSSize = 500000)
  
  # Convert GSEA results to a data frame
  gsea_df <- as.data.frame(gsea_results@result)
  
  # Check if any terms were enriched
  if (nrow(gsea_df) > 0) {
    # Filter significant terms
    significant_terms <- gsea_df[gsea_df$pvalue < 0.05, c("ID", "pvalue", "p.adjust", "core_enrichment")]
    
    # Check if any significant terms remain after filtering
    if (nrow(significant_terms) > 0) {
      # Define the output file name
      output_file <- file.path(working_dir, paste0("markers/gsea_cluster_", cluster_id, ".txt"))
      
      # Write the significant terms to a text file
      write.table(significant_terms, 
                  file = output_file, 
                  sep = "\t", 
                  row.names = FALSE, 
                  quote = FALSE, 
                  col.names = TRUE)
      
      # Print a message indicating completion for the cluster
      cat("Processed cluster:", cluster_id, "with significant terms.\n")
    } else {
      cat("No significant terms for cluster:", cluster_id, "\n")
    }
  } else {
    cat("No terms enriched for cluster:", cluster_id, "\n")
  }
}

all_top_genes_info <- all_top_genes_info[order(all_top_genes_info$avg_log2FC, decreasing = T),]
genelist <- all_top_genes_info$avg_log2FC
names(genelist) <- all_top_genes_info$GeneName



# Run GSEA using the KEGG database as an example
gsea_results <- gseKEGG(geneList = genelist, 
                        organism = 'hsa', 
                        pvalueCutoff = 0.05)

# Extract results for plotting
data <- as.data.frame(gsea_results@result)
data <- data[, c("ID", "NES", "setSize", "pvalue")]

# Create a scatter plot
data <- data[order(data$NES, decreasing = TRUE),]
data$ID <- factor(data$ID, levels = data$ID)

# Convert GSEA results to a data frame
gsea_df <- as.data.frame(gsea_results@result)

# Filter and select using base R
significant_terms <- gsea_df[gsea_df$p.adjust < 0.05, c("ID", "pvalue")]

# View the resulting data frame
print(significant_terms)

# Plot GSEA results
p <- ggplot(data, aes(x = ID, y = NES)) +
  geom_point(aes(size = setSize, alpha = -log10(p.adjust)), 
             shape = 21, stroke = 0.7, fill = "#0000ff", colour = "black") +
  scale_size_continuous(range = c(0.2, 8)) + 
  xlab("Gene Sets") + 
  ylab("Normalized Enrichment Score (NES)") + 
  theme_classic(base_size = 15) +
  theme(axis.line = element_line(color = "black", size = 0.6),
        axis.text = element_text(face = "bold"),
        axis.title = element_text(size = 13))

# Save the plot
ggsave(file.path(working_dir, "GSEA.png"), plot = p, width = 6, height = 8)

### Annotate the clusters
# names(names) <- levels(seurat)
# seurat <- RenameIdents(seurat, names)

# names <- c("0_neural_stem", "1_neural_progenitor", "2_neural_stem", "3_endothelial", 
#            "4_neural_progenitor", "5_endothelial/neural", "6_epithelial_stem", 
#            "7_endocrine", "8_neural_progenitor", "9_mesenchymal", 
#            "10_epithelial", "11_epithelial/endothelial", "12_epithelial", 
#            "13_immune/neural", "14_epithelial", "15_epithelial", 
#            "16_neural-derived_epithelial", "17_hemocyte", "18_epithelial", 
#            "19_epithelial/immune", "20_endocrine/neural", 
#            "21_neural-derived endothelial", "22_immune/neural", 
#            "23_endothelial/neural")

# names <- c("0_multipotent_stem", "1_neural_stem", "2_haematopoietic_stem", "3_neural", 
#            "4_germ", "5_epithelial", "6_primordial_germ", 
#            "7_neural_progenitor", "8_epithelial", "9_neural", 
#            "10_dendritic", "11_glial", "12_neural", 
#            "13_germ", "14_neural_progenitor", "15_epithelial", 
#            "16_?", "17_oligodendrocyte", "18_germ", 
#            "19_myeloid")

names <- c("0_multipotent_stem", "1_neural_stem", "2_lymphoid_progenitor", "3_neural", 
           "4_monocytes", "5_melanocyte+NK", "6_myeloid_progenitor", 
           "7_neural_progenitor", "8_glial", "9_neural+endothelial", 
           "10_endothelial", "11_neural_stem", "12_neural", 
           "13_myeloid", "14_neural_progenitor", "15_monocyte/DC", 
           "16_lymphoid", "17_fibroblast/muscle", "18_myeloid", 
           "19_myeloid")

# Create a named vector to map old cluster IDs to new names
names_mapping <- setNames(names, as.character(0:19))  # Ensure the keys are characters

# Add the annotation to the metadata
seurat@meta.data$annotation <- names_mapping[as.character(seurat$seurat_clusters)]

UMAP_centers_cluster <- tibble(
  UMAP_1 = as.data.frame(seurat@reductions$UMAP@cell.embeddings)$UMAP_1,
  UMAP_2 = as.data.frame(seurat@reductions$UMAP@cell.embeddings)$UMAP_2,
  cluster = seurat@meta.data$annotation
) %>%
  group_by(cluster) %>%
  dplyr::summarize(x = median(UMAP_1), y = median(UMAP_2))

# Create the UMAP plot with repelled labels
plot_umap_by_cluster <- bind_cols(seurat@meta.data, as.data.frame(seurat@reductions$UMAP@cell.embeddings)) %>%
  ggplot(aes(UMAP_1, UMAP_2, color = seurat_clusters)) +
  geom_point(size = 0.2) +
  geom_label_repel(  # Use geom_label_repel instead of geom_label
    data = UMAP_centers_cluster,
    mapping = aes(x, y, label = cluster),
    size = 3,
    fill = 'white',
    color = 'black',
    fontface = 'bold',
    alpha = 0.7,
    label.size = 0,
    show.legend = FALSE,
    box.padding = 1,  # Adjust padding around labels
    point.padding = 0.5,  # Adjust padding between labels and points
    segment.color = NA  # Color of the connecting line
  ) +
  theme_bw() +
  scale_color_manual(
    name = 'Cluster', values = custom_colors$discrete,
    guide = guide_legend(override.aes = list(size = 2))
  ) +
  theme(legend.position = 'right') +
  coord_fixed() +
  annotate(
    geom = 'text', x = Inf, y = -Inf,
    label = paste0('n = ', format(nrow(seurat@meta.data), big.mark = ',', trim = TRUE)),
    vjust = -1.5, hjust = 1.25, color = 'black', size = 2.5
  )

tSNE_centers_cluster <- tibble(
  tSNE_1 = as.data.frame(seurat@reductions$tSNE@cell.embeddings)$tSNE_1,
  tSNE_2 = as.data.frame(seurat@reductions$tSNE@cell.embeddings)$tSNE_2,
  cluster = seurat@meta.data$annotation
) %>%
  group_by(cluster) %>%
  dplyr::summarize(x = median(tSNE_1), y = median(tSNE_2))

plot_tsne_by_cluster <- bind_cols(seurat@meta.data, as.data.frame(seurat@reductions$tSNE@cell.embeddings)) %>%
  ggplot(aes(tSNE_1, tSNE_2, color = seurat_clusters)) +
  geom_point(size = 0.2) +
  geom_label_repel(  # Use geom_label_repel instead of geom_label
    data = tSNE_centers_cluster,
    mapping = aes(x, y, label = cluster),
    size = 3,
    fill = 'white',
    color = 'black',
    fontface = 'bold',
    alpha = 0.7,
    label.size = 0,
    show.legend = FALSE,
    box.padding = 0.5,  # Adjust padding around labels
    point.padding = 0.5,  # Adjust padding between labels and points
    segment.color = NA  # Color of the connecting line
  ) +
  theme_bw() +
  scale_color_manual(
    name = 'Cluster', values = custom_colors$discrete,
    guide = guide_legend(override.aes = list(size = 2))
  ) +
  theme(legend.position = 'right') +
  coord_fixed() +
  annotate(
    geom = 'text', x = Inf, y = -Inf,
    label = paste0('n = ', format(nrow(seurat@meta.data), big.mark = ',', trim = TRUE)),
    vjust = -1.5, hjust = 1.25, color = 'black', size = 2.5
  )

ggsave(
  file.path(working_dir, 'final/UMAP.annotation.pdf'),
  plot_umap_by_cluster + plot_tsne_by_cluster + 
    plot_layout(ncol = 2, widths = c(0.5, 0.5)),
  height = 8,  
  width = 16
)

### UMAP with circles
seurat@meta.data$cell_type <- as.factor(seurat@meta.data$seurat_clusters)
levels(seurat@meta.data$cell_type) <- c(levels(seurat@meta.data$cell_type), "Stem", "Progenitors", "Specialized")
seurat@meta.data$cell_type[seurat@meta.data$cell_type %in% c("0", "1", "2", "5", "11")] <- "Stem"
seurat@meta.data$cell_type[seurat@meta.data$cell_type %in% c("6", "7")] <- "Progenitors"
seurat@meta.data$cell_type[seurat@meta.data$cell_type %in% c("3", "4", "8", "9", "10", "12", "13", "14", "15", "16", "17", "18", "19")] <- "Specialized"
seurat@meta.data$cell_type <- droplevels(seurat@meta.data$cell_type)
desired_order <- c("Stem", "Progenitors", "Specialized")
seurat@meta.data$cell_type <- factor(seurat$cell_type, levels = desired_order)

# UMAP_centers_cluster <- tibble(
#   UMAP_1 = as.data.frame(seurat@reductions$UMAP@cell.embeddings)$UMAP_1,
#   UMAP_2 = as.data.frame(seurat@reductions$UMAP@cell.embeddings)$UMAP_2,
#   cluster = seurat@meta.data$cell_type_singler_blueprintencode_main
# ) %>%
#   group_by(cluster) %>%
#   dplyr::summarize(x = median(UMAP_1), y = median(UMAP_2))

plot_umap_by_celltype <- bind_cols(seurat@meta.data, as.data.frame(seurat@reductions$UMAP@cell.embeddings)) %>%
  ggplot(aes(UMAP_1, UMAP_2, color = cell_type)) +
  geom_point(size = 0.2) +
  theme_bw() +
  scale_color_manual(
    name = 'Cell Stage', values = custom_colors$discrete[c(22, 16, 38)],
    guide = guide_legend(ncol = 1, override.aes = list(size = 2))
  ) +
  theme(legend.position = 'right') +
  coord_fixed() +
  annotate(
    geom = 'text', x = Inf, y = -Inf,
    label = paste0('n = ', format(nrow(seurat@meta.data), big.mark = ',', trim = TRUE)),
    vjust = -1.5, hjust = 1.25, color = 'black', size = 2.5
  )

ggsave(
  file.path(working_dir, 'final/UMAP.single.svg'),
  plot_umap_by_celltype,
  height = 6,
  width = 6
)

# plot_umap_by_sample <- bind_cols(seurat_filtered@meta.data, as.data.frame(seurat_filtered@reductions$UMAP@cell.embeddings)) %>%
#   ggplot(aes(UMAP_1, UMAP_2, color = combined_sample)) +
#   geom_point(size = 0.2) +
#   theme_bw() +
#   scale_color_manual(
#     name = 'Sample', values = custom_colors$discrete[c(1, 12, 16)],
#     guide = guide_legend(ncol = 1, override.aes = list(size = 2))
#   ) +
#   theme(legend.position = 'right') +
#   coord_fixed() +
#   annotate(
#     geom = 'text', x = Inf, y = -Inf,
#     label = paste0('n = ', format(nrow(seurat_filtered@meta.data), big.mark = ',', trim = TRUE)),
#     vjust = -1.5, hjust = 1.25, color = 'black', size = 2.5
#   )
# 
# ggsave(
#   file.path(working_dir, 'final/UMAP_by_3_samples.single.svg'),
#   plot_umap_by_sample,
#   height = 6,
#   width = 6
# )

# UMAP by sample with circles, only3
seurat <- only3_seurat(seurat)

plot_data <- data.frame(
  seurat@meta.data,
  seurat@reductions$UMAP@cell.embeddings
)

maskTable <- generateMask(
  dims = plot_data[, c("UMAP_1", "UMAP_2")], 
  cluster = plot_data$cell_type, 
  minDensity = 0,
  smoothSigma = 0.01
)

plot_umap_by_sample <- ggplot(plot_data, aes(x = UMAP_1, y = UMAP_2)) + 
  geom_point(aes(color = combined_sample), size = 0.2) +
  # Add irregular boundary around Stem cells (black dashed)
  geom_path(
    data = maskTable[maskTable$cluster == "Stem", ], 
    aes(group = group),
    linewidth = 0.8,
    linetype = "dashed",
    color = "black"
  ) +
  # Add irregular boundary around Specialized cells (blue dashed)
  geom_path(
    data = maskTable[maskTable$cluster == "Specialized", ], 
    aes(group = group),
    linewidth = 0.8,
    linetype = "dashed",
    color = custom_colors$discrete[13]
  ) +
  theme_bw() +
  scale_color_manual(
    name = 'Sample', 
    values = custom_colors$discrete[c(1, 12, 16)],
    guide = guide_legend(ncol = 1, override.aes = list(size = 2))
  ) +
  theme(legend.position = 'right') +
  coord_fixed() +
  # Add cell count annotation
  annotate(
    geom = 'text', 
    x = x_range[2], y = y_range[1],
    label = paste0('n = ', format(nrow(seurat@meta.data), big.mark = ',', trim = TRUE)),
    vjust = -1, hjust = 1, color = 'black', size = 2.5
  ) +
  # Add legend for boundary colors
  annotate("segment", 
           x = range(plot_data$UMAP_1)[1], y = range(plot_data$UMAP_2)[1] + 1.5, 
           xend = range(plot_data$UMAP_1)[1] + 1, yend = range(plot_data$UMAP_2)[1] + 1.5,
           color = "black", linetype = "dashed", size = 1) +
  annotate("text", 
           x = range(plot_data$UMAP_1)[1] + 2, y = range(plot_data$UMAP_2)[1] + 1.5, 
           label = "Stem cells", hjust = 0, vjust = 0.5, size = 3) +
  annotate("segment", 
           x = range(plot_data$UMAP_1)[1], y = range(plot_data$UMAP_2)[1], 
           xend = range(plot_data$UMAP_1)[1] + 1, yend = range(plot_data$UMAP_2)[1],
           color = custom_colors$discrete[13], linetype = "dashed", size = 1) +
  annotate("text", 
           x = range(plot_data$UMAP_1)[1] + 2, y = range(plot_data$UMAP_2)[1], 
           label = "Specialized cells", hjust = 0, vjust = 0.5, size = 3)

print(plot_umap_by_sample)

ggsave(
  file.path(working_dir, 'final/UMAP_only3_samples.single.svg'),
  plot_umap_by_sample,
  height = 6,
  width = 6
)

# UMAP by sample with circles, 4 samples
seurat <- purge4_seurat(seurat)

# Generate mask using entire pooled dataset (same boundaries for all plots)
plot_data <- data.frame(
  seurat@meta.data,
  seurat@reductions$UMAP@cell.embeddings
)

maskTable <- generateMask(
  dims = plot_data[, c("UMAP_1", "UMAP_2")], 
  cluster = plot_data$cell_type, 
  minDensity = 0,
  smoothSigma = 0.001
)

# Get coordinate ranges for consistent annotation positioning
x_range <- range(plot_data$UMAP_1)
y_range <- range(plot_data$UMAP_2)

# Create individual plots for each sample
blank_plots <- list()
sample_plots <- list()
samples <- levels(droplevels(plot_data$merged_sample))
i <- 0

for(sample in samples) {
  i <- i + 1
  
  # Filter data for this sample
  sample_data <- plot_data[plot_data$merged_sample == sample, ]
  
  # Create plot for this sample
  p <- ggplot(sample_data, aes(x = UMAP_1, y = UMAP_2)) + 
    geom_point(aes(color = merged_sample), size = 0.2, show.legend = FALSE) +
    # Add the same mask boundaries to each plot
    geom_path(
      data = maskTable[maskTable$cluster == "Stem", ], 
      aes(x = UMAP_1, y = UMAP_2, group = group),
      linewidth = 0.8,
      linetype = "dashed",
      color = "black",
      inherit.aes = FALSE
    ) +
    geom_path(
      data = maskTable[maskTable$cluster == "Progenitors", ], 
      aes(x = UMAP_1, y = UMAP_2, group = group),
      linewidth = 0.8,
      linetype = "dashed",
      color = custom_colors$discrete[16],
      inherit.aes = FALSE
    ) +
    geom_path(
      data = maskTable[maskTable$cluster == "Specialized", ], 
      aes(x = UMAP_1, y = UMAP_2, group = group),
      linewidth = 0.8,
      linetype = "dashed",
      color = custom_colors$discrete[13],
      inherit.aes = FALSE
    ) +
    theme_bw() +
    scale_color_manual(values = custom_colors$discrete[i]) +
    coord_fixed() +
    xlim(x_range) + ylim(y_range) +  # Keep same axis limits for all plots
    ggtitle(sample) +
    theme(
      legend.position = 'none',
      plot.title = element_text(hjust = 0.5, face = "bold")
    ) +
    # Add cell count for this sample
    annotate(
      geom = 'text', 
      x = x_range[2], y = y_range[1],
      label = paste0('n = ', format(nrow(sample_data), big.mark = ',', trim = TRUE)),
      vjust = -1, hjust = 1, color = 'black', size = 2.5
    )
  
  sample_plots[[sample]] <- p
  
  # Create plot for this sample
  p <- ggplot(sample_data, aes(x = UMAP_1, y = UMAP_2)) + 
    geom_point(aes(color = merged_sample), size = 0.2, show.legend = FALSE) +
    theme_bw() +
    scale_color_manual(values = custom_colors$discrete[i]) +
    coord_fixed() +
    xlim(x_range) + ylim(y_range) +  # Keep same axis limits for all plots
    ggtitle(sample) +
    theme(
      legend.position = 'none',
      plot.title = element_text(hjust = 0.5, face = "bold")
    ) +
    # Add cell count for this sample
    annotate(
      geom = 'text', 
      x = x_range[2], y = y_range[1],
      label = paste0('n = ', format(nrow(sample_data), big.mark = ',', trim = TRUE)),
      vjust = -1, hjust = 1, color = 'black', size = 2.5
    )
  
  blank_plots[[sample]] <- p
}

# Create horizontal legend
legend_plot <- ggplot() +
  annotate("segment", 
           x = 1, y = 0.5, xend = 2, yend = 0.5,
           color = "black", linetype = "dashed", size = 1) +
  annotate("text", 
           x = 2.2, y = 0.5, 
           label = "Stem cells", hjust = 0, vjust = 0.5, size = 3) +
  annotate("segment", 
           x = 4, y = 0.5, xend = 5, yend = 0.5,
           color = custom_colors$discrete[16], linetype = "dashed", size = 1) +
  annotate("text", 
           x = 5.2, y = 0.5, 
           label = "Progenitor cells", hjust = 0, vjust = 0.5, size = 3) +
  annotate("segment", 
           x = 7.5, y = 0.5, xend = 8.5, yend = 0.5,
           color = custom_colors$discrete[13], linetype = "dashed", size = 1) +
  annotate("text", 
           x = 8.7, y = 0.5, 
           label = "Specialized cells", hjust = 0, vjust = 0.5, size = 3) +
  xlim(0, 12) + ylim(0, 1) +
  theme_void()

# Combine all plots in a grid
combined_plots_1 <- wrap_plots(blank_plots, ncol = 4)
combined_plots_2 <- wrap_plots(sample_plots, ncol = 4)  # Adjust ncol as needed
final_plot <- combined_plots_1 / combined_plots_2 / legend_plot + plot_layout(heights = c(10, 10, 1))

print(final_plot)

ggsave(
  file.path(working_dir, paste0('final/UMAP_by_sample_circled.', runid, '.pdf')),
  final_plot,
  height = 8,
  width = 14
)

# per sample by 4 samples
plot_data <- data.frame(
  seurat@meta.data,
  seurat@reductions$UMAP@cell.embeddings
)

# Generate mask for cell types
maskTable <- generateMask(
  dims = plot_data[, c("UMAP_1", "UMAP_2")], 
  cluster = plot_data$cell_type, 
  minDensity = 0,
  smoothSigma = 0.01
)

# Create temp_labels that matches facets
temp_labels <- bind_cols(seurat@meta.data, as.data.frame(seurat@reductions$UMAP@cell.embeddings)) %>%
  group_by(merged_sample) %>%
  dplyr::summarise(n = n()) %>%
  ungroup()

# Create the faceted plot with geometric forms
umap_plot <- bind_cols(seurat@meta.data, as.data.frame(seurat@reductions$UMAP@cell.embeddings)) %>%
  ggplot(aes(UMAP_1, UMAP_2, color = merged_sample)) +
  geom_point(size = 0.2, show.legend = FALSE) +
  # Add geometric forms for Stem cells
  geom_path(
    data = maskTable[maskTable$cluster == "Stem", ], 
    aes(x = UMAP_1, y = UMAP_2, group = group),  # Fixed: add x and y aesthetics
    linewidth = 0.8,
    linetype = "dashed",
    color = "black",
    inherit.aes = FALSE
  ) +
  # Add geometric forms for Progenitors
  geom_path(
    data = maskTable[maskTable$cluster == "Progenitors", ], 
    aes(x = UMAP_1, y = UMAP_2, group = group),  # Fixed: add x and y aesthetics
    linewidth = 0.8,
    linetype = "dashed",
    color = custom_colors$discrete[16],
    inherit.aes = FALSE
  ) +
  # Add geometric forms for Specialized cells
  geom_path(
    data = maskTable[maskTable$cluster == "Specialized", ], 
    aes(x = UMAP_1, y = UMAP_2, group = group),  # Fixed: add x and y aesthetics
    linewidth = 0.8,
    linetype = "dashed",
    color = custom_colors$discrete[13],
    inherit.aes = FALSE
  ) +
  geom_text(
    data = temp_labels,
    aes(x = Inf, y = -Inf, label = paste0('n = ', format(n, big.mark = ',', trim = TRUE))),
    vjust = -1.5, hjust = 1.25,
    color = 'black', size = 2.8,
    inherit.aes = FALSE
  ) +
  theme_bw() +
  scale_color_manual(values = custom_colors$discrete) +
  coord_fixed() +
  theme(
    legend.position = 'none'
  ) +
  facet_wrap(~merged_sample, ncol = 4)

# Create horizontal legend at bottom
legend_plot <- ggplot() +
  # Stem cells legend
  annotate("segment", 
           x = 1, y = 0.5, xend = 2, yend = 0.5,
           color = "black", linetype = "dashed", size = 1) +
  annotate("text", 
           x = 2.2, y = 0.5, 
           label = "Stem cells", hjust = 0, vjust = 0.5, size = 3) +
  # Progenitor cells legend
  annotate("segment", 
           x = 4, y = 0.5, xend = 5, yend = 0.5,
           color = custom_colors$discrete[16], linetype = "dashed", size = 1) +
  annotate("text", 
           x = 5.2, y = 0.5, 
           label = "Progenitor cells", hjust = 0, vjust = 0.5, size = 3) +
  # Specialized cells legend
  annotate("segment", 
           x = 7.5, y = 0.5, xend = 8.5, yend = 0.5,
           color = custom_colors$discrete[13], linetype = "dashed", size = 1) +
  annotate("text", 
           x = 8.7, y = 0.5, 
           label = "Specialized cells", hjust = 0, vjust = 0.5, size = 3) +
  xlim(0, 12) + ylim(0, 1) +
  theme_void()

# Combine plots vertically
combined_plot <- umap_plot / legend_plot + plot_layout(heights = c(10, 1))
print(combined_plot)

ggsave(
  file.path(working_dir, 'final/UMAP_by_4_samples.single.svg'),
  plot_umap_by_sample,
  height = 6,
  width = 6
)

## print markers

# # Specify the output file path
# output_file <- file.path(working_dir, "top_10_marker_genes.txt")
# 
# # Open a connection to the output file
# sink(output_file)
# 
# # Write the header
# cat("Cluster_ID\tTop_10_Marker_Genes\n")
# 
# # Iterate over each cluster and write the top markers to the file
# for (cluster_id in names(cluster_markers)) {
#   # Get the top markers for the current cluster
#   top_genes <- cluster_markers[[cluster_id]]
#   
#   # Select only the first 10 genes (if available)
#   top_genes_subset <- head(top_genes, 10)
#   
#   # Create a line with the cluster ID and the top genes
#   line <- paste(cluster_id, paste(top_genes_subset, collapse = ", "), sep = "\t")
#   
#   # Write the line to the file
#   cat(line, "\n")
# }
# 
# # Close the connection
# sink()
# 
# # Notify the user
# cat("Top 10 marker genes for each cluster have been written to", output_file, "\n")

### monocle
seurat.cds <- as.cell_data_set(seurat)
seurat.cds <- cluster_cells(seurat.cds, reduction_method = "UMAP")
seurat.cds <- learn_graph(seurat.cds, use_partition = TRUE)
seurat.cds <- order_cells(seurat.cds, reduction_method = "UMAP")

seurat.cds <- reduce_dimension(seurat.cds, reduction_method = "tSNE")
seurat.cds <- cluster_cells(seurat.cds, reduction_method = "tSNE")
seurat.cds <- learn_graph(seurat.cds, use_partition = TRUE)
seurat.cds <- order_cells(seurat.cds, reduction_method = "tSNE")

# save_monocle_objects(seurat.cds, directory_path = file.path(working_dir, "cytotrace2_result.single.cds"))
save_monocle_objects(seurat.cds, directory_path = file.path(working_dir, "cytotrace2_result.filtered.single.cds"))
save_monocle_objects(seurat.cds, directory_path = file.path(working_dir, "cytotrace2_result.5.singlet.s-g.cds"))

p1 <- plot_cells(seurat.cds, color_cells_by = "partition", reduction_method = "UMAP")

p2 <- plot_cells(seurat.cds, 
                 color_cells_by = "pseudotime", 
                 show_trajectory = TRUE,
                 reduction_method = "UMAP"
)

plot_cells(seurat.cds, color_cells_by = "partition", reduction_method = "tSNE")

plot_cells(seurat.cds,
           color_cells_by = "pseudotime",
           show_trajectory = TRUE,
           reduction_method = "tSNE"
)

ggsave(
  file.path(working_dir, 'monocle.5.singlet.s-g.png'),
  # file.path(working_dir, 'monocle_filtered.png'),
  p1 + p2,
  height = 6,  
  width = 12
)

### CytoTRACE2
## load data
# saveRDS(seurat, file.path(working_dir, "cytotrace2_result.single.rds"))
seurat <- readRDS(file.path(working_dir, "cytotrace2_result.single.rds"))

# saveRDS(seurat, file.path(working_dir, "gfp.single.rds"))
seurat <- readRDS(file.path(working_dir, "gfp.single.rds"))

seurat <- add_sample_metadata(seurat)
seurat <- filtering_seurat(seurat)

# saveRDS(seurat, file.path(working_dir, "cytotrace2_result.filtered.single.rds"))
seurat <- readRDS(file.path(working_dir, "cytotrace2_result.filtered.single.rds"))

# saveRDS(seurat, file.path(working_dir, "gfp.filtered.rds"))
seurat <- readRDS(file.path(working_dir, "gfp.filtered.rds"))

seurat <- subset(seurat, subset = GFP > 0)
seurat <- filtering_seurat(seurat)
seurat <- purge4_seurat(seurat)
seurat <- only3_seurat(seurat)

## not seurat
cytotrace2_result <- cytotrace2(seurat[['RNA']]$counts,   
                                species = "human",
                                is_seurat = FALSE,
                                seed = 42)  

plots <- plotData(cytotrace2_result,
                  annotation = seurat[['seurat_clusters']],
                  expression_data = data,
                  is_seurat = FALSE,
                  pc_dims = pca_cutoff,
                  seed = 14)

## is seurat
cytotrace2_result <- cytotrace2(seurat,   
                                species = "human",
                                is_seurat = TRUE,
                                slot_type = "counts",
                                seed = 42) 

# Create a temporary data frame with sample labels
temp_labels <- seurat@meta.data %>%
  group_by(merged_sample) %>%
  tally()

CytoTRACE2_Score_umap_plot <- bind_cols(seurat@meta.data, as.data.frame(seurat@reductions$UMAP@cell.embeddings)) %>%
  ggplot(aes(UMAP_1, UMAP_2, color = CytoTRACE2_Score)) +
  geom_point(size = 1) +
  geom_text(
    data = temp_labels,
    aes(x = Inf, y = -Inf, label = paste0('n = ', format(n, big.mark = ',', trim = TRUE)), vjust = -1.5, hjust = 1.25),
    color = 'black', size = 2.8
  ) +
  theme_bw() +
  scale_color_viridis(
    name = "CytoTRACE2_Score",
    guide = guide_colorbar(frame.colour = 'black', ticks.colour = 'black'),
    labels = scales::comma,
  ) +
  theme(legend.position = 'right') +
  coord_fixed() +
  facet_wrap(~merged_sample, ncol = 6)

CytoTRACE2_Potency_umap_plot <- bind_cols(seurat@meta.data, as.data.frame(seurat@reductions$UMAP@cell.embeddings)) %>%
  ggplot(aes(UMAP_1, UMAP_2, color = CytoTRACE2_Potency)) +
  geom_text(
    data = temp_labels,
    aes(x = Inf, y = -Inf, label = paste0('n = ', format(n, big.mark = ',', trim = TRUE)), vjust = -1.5, hjust = 1.25),
    color = 'black', size = 2.8
  ) +
  geom_point(size = 1) +
  theme_bw() +
  scale_color_viridis_d(
    name = "CytoTRACE_Potency",
    guide = guide_legend(ncol = 1, override.aes = list(size = 2))
  ) +
  theme(legend.position = 'right') +
  coord_fixed() +
  facet_wrap(~merged_sample, ncol = 6)  

CytoTRACE2_Relative_umap_plot <- bind_cols(seurat@meta.data, as.data.frame(seurat@reductions$UMAP@cell.embeddings)) %>%
  ggplot(aes(UMAP_1, UMAP_2, color = CytoTRACE2_Relative)) +
  geom_text(
    data = temp_labels,
    aes(x = Inf, y = -Inf, label = paste0('n = ', format(n, big.mark = ',', trim = TRUE)), vjust = -1.5, hjust = 1.25),
    color = 'black', size = 2.8
  ) +
  geom_point(size = 1) +
  theme_bw() +
  scale_color_viridis(
    name = "CytoTRACE2_Rank",
    guide = guide_colorbar(frame.colour = 'black', ticks.colour = 'black'),
    labels = scales::comma,
  ) +
  theme(legend.position = 'right') +
  coord_fixed() +
  facet_wrap(~merged_sample, ncol = 6)  

combined_plot <- CytoTRACE2_Score_umap_plot / CytoTRACE2_Relative_umap_plot / CytoTRACE2_Potency_umap_plot 

ggsave(
  file.path(working_dir, 'final/cytotrace.UMAP.single.svg'),
  combined_plot,
  height = 8,  
  width = 10
)

# Create a temporary data frame with sample labels
temp_labels <- seurat@meta.data %>%
  group_by(merged_sample) %>%
  tally()

CytoTRACE2_Score_tSNE_plot <- bind_cols(seurat@meta.data, as.data.frame(seurat@reductions$tSNE@cell.embeddings)) %>%
  ggplot(aes(tSNE_1, tSNE_2, color = CytoTRACE2_Score)) +
  geom_point(size = 1) +
  geom_text(
    data = temp_labels,
    aes(x = Inf, y = -Inf, label = paste0('n = ', format(n, big.mark = ',', trim = TRUE)), vjust = -1.5, hjust = 1.25),
    color = 'black', size = 2.8
  ) +
  theme_bw() +
  scale_color_viridis(
    name = "CytoTRACE2_Score",
    guide = guide_colorbar(frame.colour = 'black', ticks.colour = 'black'),
    labels = scales::comma,
  ) +
  theme(legend.position = 'right') +
  coord_fixed() +
  facet_wrap(~merged_sample, ncol = 6)

CytoTRACE2_Potency_tSNE_plot <- bind_cols(seurat@meta.data, as.data.frame(seurat@reductions$tSNE@cell.embeddings)) %>%
  ggplot(aes(tSNE_1, tSNE_2, color = CytoTRACE2_Potency)) +
  geom_text(
    data = temp_labels,
    aes(x = Inf, y = -Inf, label = paste0('n = ', format(n, big.mark = ',', trim = TRUE)), vjust = -1.5, hjust = 1.25),
    color = 'black', size = 2.8
  ) +
  geom_point(size = 1) +
  theme_bw() +
  scale_color_viridis_d(
    name = "CytoTRACE_Potency",
    guide = guide_legend(ncol = 1, override.aes = list(size = 2))
  ) +
  theme(legend.position = 'right') +
  coord_fixed() +
  facet_wrap(~merged_sample, ncol = 6)  

CytoTRACE2_Relative_tSNE_plot <- bind_cols(seurat@meta.data, as.data.frame(seurat@reductions$tSNE@cell.embeddings)) %>%
  ggplot(aes(tSNE_1, tSNE_2, color = CytoTRACE2_Relative)) +
  geom_text(
    data = temp_labels,
    aes(x = Inf, y = -Inf, label = paste0('n = ', format(n, big.mark = ',', trim = TRUE)), vjust = -1.5, hjust = 1.25),
    color = 'black', size = 2.8
  ) +
  geom_point(size = 1) +
  theme_bw() +
  scale_color_viridis(
    name = "CytoTRACE2_Rank",
    guide = guide_colorbar(frame.colour = 'black', ticks.colour = 'black'),
    labels = scales::comma,
  ) +
  theme(legend.position = 'right') +
  coord_fixed() +
  facet_wrap(~merged_sample, ncol = 6)  

combined_plot <- CytoTRACE2_Score_tSNE_plot / CytoTRACE2_Relative_tSNE_plot / CytoTRACE2_Potency_tSNE_plot 

ggsave(
  file.path(working_dir, 'final/cytotrace.tSNE.single.svg'),
  combined_plot,
  height = 8,  
  width = 10
)

table_samples_by_CytoTRACE2_Potency <- seurat@meta.data %>%
  group_by(merged_sample, CytoTRACE2_Potency) %>%
  dplyr::summarise(count = n()) %>%
  spread(CytoTRACE2_Potency, count, fill = 0) %>%
  ungroup() %>%
  mutate(total_cell_count = rowSums(.[c(2:ncol(.))])) %>%
  dplyr::select(c('merged_sample', 'total_cell_count', everything())) %>%
  arrange(factor(merged_sample, levels = levels(seurat@meta.data$merged_sample)))

knitr::kable(table_samples_by_CytoTRACE2_Potency)
# |sample       | total_cell_count| Differentiated| Unipotent| Oligopotent| Multipotent|
#   |:------------|----------------:|--------------:|---------:|-----------:|-----------:|
#   |ALDH         |              723|            645|        33|          39|           6|
#   |STEM_B_1     |              844|            596|       124|         112|          12|
#   |STEM_B_2     |              953|            711|       106|         121|          15|
#   |CONTROL      |             1144|           1115|        22|           7|           0|
#   |STEM_A       |             1226|           1197|        24|           5|           0|
#   |Undetermined |             3575|           2957|       310|         295|          13|

table_clusters_by_CytoTRACE2_Potency <- seurat@meta.data %>%
  dplyr::rename('cluster' = 'seurat_clusters') %>%
  group_by(cluster, CytoTRACE2_Potency) %>%
  dplyr::summarize(count = n()) %>%
  spread(CytoTRACE2_Potency, count, fill = 0) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(total_cell_count = rowSums(.[c(2:ncol(.))])) %>%
  dplyr::select(c('cluster', 'total_cell_count', everything())) %>%
  arrange(factor(cluster, levels = levels(seurat@meta.data$seurat_clusters)))

knitr::kable(table_clusters_by_CytoTRACE2_Potency)

# |cluster | total_cell_count| Differentiated| Unipotent| Oligopotent| Multipotent|
#   |:-------|----------------:|--------------:|---------:|-----------:|-----------:|
#   |0       |             1595|            632|       397|         525|          41|
#   |1       |             1058|           1051|         7|           0|           0|
#   |2       |              731|            705|        24|           2|           0|
#   |3       |              559|            559|         0|           0|           0|
#   |4       |              544|            519|        19|           6|           0|
#   |5       |              532|            385|       110|          32|           5|
#   |6       |              511|            483|        21|           7|           0|
#   |7       |              387|            386|         1|           0|           0|
#   |8       |              347|            347|         0|           0|           0|
#   |9       |              321|            321|         0|           0|           0|
#   |10      |              320|            310|         9|           1|           0|
#   |11      |              298|            273|        19|           6|           0|
#   |12      |              260|            260|         0|           0|           0|
#   |13      |              188|            188|         0|           0|           0|
#   |14      |              179|            179|         0|           0|           0|
#   |15      |              162|            158|         4|           0|           0|
#   |16      |              149|            143|         6|           0|           0|
#   |17      |              144|            144|         0|           0|           0|
#   |18      |              119|            117|         2|           0|           0|
#   |19      |               61|             61|         0|           0|           0|

# Define your desired order for the samples
sample_order <- c("ALDH", "STEM_B_1", "STEM_B_2", "CONTROL", "STEM_A", "Undetermined") 

# Convert the sample column to a factor with the specified order
seurat@meta.data$sample <- factor(seurat@meta.data$sample, levels = sample_order)

temp_labels <- table_samples_by_CytoTRACE2_Potency %>%
  dplyr::select(merged_sample, total_cell_count) %>%
  mutate(n = total_cell_count)  # Create n for geom_text

p1 <- table_samples_by_CytoTRACE2_Potency %>%
  dplyr::select(-c('total_cell_count')) %>%
  reshape2::melt(id.vars = 'merged_sample') %>%
  dplyr::mutate(sample = factor(sample, levels = sample_order)) %>%
  ggplot(aes(merged_sample, value)) +
  geom_bar(aes(fill = variable), position = 'stack', stat = 'identity') +
  geom_text(
    data = temp_labels,
    aes(x = merged_sample, y = Inf, label = paste0('n = ', format(n, big.mark = ',', trim = TRUE)), vjust = -1),
    color = 'black', size = 2.8
  ) +
  scale_fill_manual(name = 'CytoTRACE2_Potency', values = custom_colors$discrete) +
  scale_y_continuous(name = 'Number of cells', labels = scales::comma, expand = c(0.01, 0)) +
  coord_cartesian(clip = 'off') +
  theme_bw() +
  theme(
    legend.position = 'left',
    plot.title = element_text(hjust = 0.5),
    text = element_text(size = 16),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    plot.margin = margin(t = 20, r = 0, b = 0, l = 0, unit = 'pt')
  )

temp_labels <- seurat@meta.data %>%
  group_by(seurat_clusters) %>%
  tally() %>%
  dplyr::rename('cluster' = seurat_clusters)

p2 <- table_clusters_by_CytoTRACE2_Potency %>%
  dplyr::select(-c('total_cell_count')) %>%
  reshape2::melt(id.vars = 'cluster') %>%
  dplyr::mutate(cluster = factor(cluster, levels = levels(seurat@meta.data$seurat_clusters))) %>%
  ggplot(aes(cluster, value)) +
  geom_bar(aes(fill = variable), position = 'stack', stat = 'identity') +
  geom_text(
    data = temp_labels,
    aes(x = cluster, y = Inf, label = paste0('n = ', format(n, big.mark = ',', trim = TRUE)), vjust = -1),
    color = 'black', size = 2.8
  ) +
  scale_fill_manual(name = 'CytoTRACE2_Potency', values = custom_colors$discrete) +
  scale_y_continuous(labels = scales::comma, expand = c(0.01, 0)) +
  coord_cartesian(clip = 'off') +
  theme_bw() +
  theme(
    legend.position = 'right',
    plot.title = element_text(hjust = 0.5),
    text = element_text(size = 16),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title = element_blank(),
    plot.margin = margin(t = 20, r = 0, b = 0, l = 10, unit = 'pt')
  )

ggsave(
  file.path(working_dir, 'composition_cytotrace.single.png'),
  p1 + p2 +
    plot_layout(ncol = 2, widths = c(
      seurat@meta.data$merged_sample %>% unique() %>% length(),
      seurat@meta.data$seurat_clusters %>% unique() %>% length()
    )),
  width = 18, height = 8
)

## single R
# singler_ref <- BlueprintEncodeData()
# /home/users/jiamuyu/.cache/R/ExperimentHub
singler_ref <- HumanPrimaryCellAtlasData()

singler_results_blueprintencode_main <- SingleR::SingleR(
  test = GetAssayData(seurat, assay = 'integrated', slot = 'scale.data'),
  ref = singler_ref,
  labels = singler_ref@colData@listData$label.main
)

singler_results_blueprintencode_fine <- SingleR::SingleR(
  test = GetAssayData(seurat, assay = 'integrated', slot = 'scale.data'),
  ref = singler_ref,
  labels = singler_ref@colData@listData$label.fine
)

# saveRDS(singler_results_blueprintencode_main, file.path(working_dir, "singler_main.rds"))
# saveRDS(singler_results_blueprintencode_fine, file.path(working_dir, "singler_fine.rds"))

singler_results_blueprintencode_main <- readRDS(file.path(working_dir, "singler_main.rds"))
singler_results_blueprintencode_fine <- readRDS(file.path(working_dir, "singler_fine.rds"))

# saveRDS(singler_results_blueprintencode_main, file.path(working_dir, "singler_main_filtered.rds"))
# saveRDS(singler_results_blueprintencode_fine, file.path(working_dir, "singler_fine_filtered.rds"))

singler_results_blueprintencode_main <- readRDS(file.path(working_dir, "singler_main_filtered.rds"))
singler_results_blueprintencode_fine <- readRDS(file.path(working_dir, "singler_fine_filtered.rds"))

p <- plotScoreHeatmap(
  singler_results_blueprintencode_main,
  show.labels = TRUE,
  annotation_col = data.frame(
    donor = seurat@meta.data$sample,
    row.names = rownames(singler_results_blueprintencode_main)
  )
)

ggsave(
  # file.path(working_dir, 'singler_main.png'), p,
  file.path(working_dir, 'singler_main_filtered.png'), p,
  height = 12, width = 12
)

p <- plotScoreHeatmap(
  singler_results_blueprintencode_fine,
  show.labels = TRUE,
  annotation_col = data.frame(
    donor = seurat@meta.data$sample,
    row.names = rownames(singler_results_blueprintencode_fine)
  )
)

ggsave(
  # file.path(working_dir, 'singler_fine.png'), p,
  file.path(working_dir, 'singler_fine_filtered.png'), p,
  height = 12, width = 12
)

seurat@meta.data$cell_type_singler_blueprintencode_main <- singler_results_blueprintencode_main@listData$labels

singler_scores <- singler_results_blueprintencode_main@listData$scores %>%
  as_tibble() %>%
  dplyr::mutate(assigned_score = NA)

for (i in seq_len(nrow(singler_scores))) {
  singler_scores$assigned_score[i] <- singler_scores[[singler_results_blueprintencode_main@listData$labels[i]]][i]
}

seurat@meta.data$cell_type_singler_blueprintencode_main_score <- singler_scores$assigned_score

seurat@meta.data$cell_type_singler_blueprintencode_fine <- singler_results_blueprintencode_fine@listData$labels

singler_scores <- singler_results_blueprintencode_fine@listData$scores %>%
  as_tibble() %>%
  dplyr::mutate(assigned_score = NA)

for (i in seq_len(nrow(singler_scores))) {
  singler_scores$assigned_score[i] <- singler_scores[[singler_results_blueprintencode_fine@listData$labels[i]]][i]
}

seurat@meta.data$cell_type_singler_blueprintencode_fine_score <- singler_scores$assigned_score

temp_labels <- seurat@meta.data %>%
  group_by(merged_sample) %>%
  tally()

p1 <- ggplot() +
  geom_half_violin(
    data = seurat@meta.data,
    aes(
      x = merged_sample,
      y = cell_type_singler_blueprintencode_main_score,
      fill = sample
    ),
    side = 'l', show.legend = FALSE, trim = FALSE
  ) +
  geom_half_boxplot(
    data = seurat@meta.data,
    aes(
      x = merged_sample,
      y = cell_type_singler_blueprintencode_main_score,
      fill = sample
    ),
    side = 'r', outlier.color = NA, width = 0.4, show.legend = FALSE
  ) +
  geom_text(
    data = temp_labels,
    aes(
      x = merged_sample,
      y = -Inf,
      label = paste0('n = ', format(n, big.mark = ',', trim = TRUE)),
      vjust = -1
    ),
    color = 'black', size = 2.8
  ) +
  scale_color_manual(values = custom_colors$discrete) +
  scale_fill_manual(values = custom_colors$discrete) +
  scale_y_continuous(
    name = 'Assignment score',
    labels = scales::comma,
    limits = c(0,1)
  ) +
  theme_bw() +
  labs(title = 'Samples') +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank()
  )

temp_labels <- seurat@meta.data %>%
  group_by(seurat_clusters) %>%
  tally()

p2 <- ggplot() +
  geom_half_violin(
    data = seurat@meta.data,
    aes(
      x = seurat_clusters,
      y = cell_type_singler_blueprintencode_main_score,
      fill = seurat_clusters
    ),
    side = 'l', show.legend = FALSE, trim = FALSE
  ) +
  geom_half_boxplot(
    data = seurat@meta.data,
    aes(
      x = seurat_clusters,
      y = cell_type_singler_blueprintencode_main_score,
      fill = seurat_clusters
    ),
    side = 'r', outlier.color = NA, width = 0.4, show.legend = FALSE
  ) +
  geom_text(
    data = temp_labels,
    aes(
      x = seurat_clusters,
      y = -Inf,
      label = paste0('n = ', format(n, big.mark = ',', trim = TRUE)),
      vjust = -1
    ),
    color = 'black', size = 2.8
  ) +
  scale_color_manual(values = custom_colors$discrete) +
  scale_fill_manual(values = custom_colors$discrete) +
  scale_y_continuous(
    name = 'Assignment score',
    labels = scales::comma,
    limits = c(0,1)
  ) +
  labs(title = 'Clusters') +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank()
  )

cell_counts <- table(seurat@meta.data$cell_type_singler_blueprintencode_main)

# Identify cell types with more than one observation
valid_cell_types <- names(cell_counts[cell_counts > 1])

# Filter the metadata to include only these cell types
filtered_meta <- seurat@meta.data %>%
  filter(cell_type_singler_blueprintencode_main %in% valid_cell_types)

# Update temp_labels to match filtered data
temp_labels <- filtered_meta %>%
  group_by(cell_type_singler_blueprintencode_main) %>%
  tally()

p3 <- ggplot() +
  gghalves::geom_half_violin(
    data = filtered_meta,
    mapping = aes(
      x = cell_type_singler_blueprintencode_main,
      y = cell_type_singler_blueprintencode_main_score,
      fill = cell_type_singler_blueprintencode_main
    ),
    side = "left",
    show.legend = FALSE,
    trim = FALSE
  ) +
  gghalves::geom_half_boxplot(
    data = filtered_meta,
    mapping = aes(
      x = cell_type_singler_blueprintencode_main,
      y = cell_type_singler_blueprintencode_main_score,
      fill = cell_type_singler_blueprintencode_main
    ),
    side = "right",
    outlier.color = NA, 
    width = 0.4, 
    show.legend = FALSE
  ) +
  geom_text(
    data = temp_labels,
    aes(
      x = cell_type_singler_blueprintencode_main,
      y = 0,  # Using 0 instead of -Inf for better positioning
      label = paste0('n = ', format(n, big.mark = ',', trim = TRUE))
    ),
    vjust = -1,
    color = 'black', 
    size = 2.8
  ) +
  scale_color_manual(values = custom_colors$discrete) +
  scale_fill_manual(values = custom_colors$discrete) +
  scale_y_continuous(
    name = 'Assignment score',
    labels = scales::comma,
    limits = c(0, 1)
  ) +
  labs(title = 'Cell types') +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)  # Angled labels for readability
  )

ggsave(
  file.path(working_dir, 'singler_scores.brain.png'),
  p1 + p2 + p3 + plot_layout(ncol = 1), height = 13, width = 16
)

table_samples_by_cell_type <- seurat@meta.data %>%
  dplyr::group_by(seurat_clusters, cell_type_singler_blueprintencode_main) %>%
  dplyr::summarize(count = n()) %>%
  tidyr::spread(cell_type_singler_blueprintencode_main, count, fill = 0) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(total_cell_count = rowSums(.[c(2:ncol(.))])) %>%
  dplyr::select(c("seurat_clusters", "total_cell_count", dplyr::everything()))

# Prepare the data for the heatmap
heatmap_data <- table_samples_by_cell_type[, -which(names(table_samples_by_cell_type) == "total_cell_count")]

# Reshape to long format
heatmap_data_long <- melt(heatmap_data, id.vars = "seurat_clusters", variable.name = "cell_type", value.name = "count")

# Normalize by row max
heatmap_data_long <- heatmap_data_long %>%
  group_by(seurat_clusters) %>%
  mutate(count = count / sum(count)) %>%
  ungroup()

# Create a matrix for clustering
heatmap_matrix <- acast(heatmap_data_long, seurat_clusters ~ cell_type, value.var = "count")
heatmap_matrix <- sweep(heatmap_matrix, 1, rowSums(heatmap_matrix), FUN = "/")

# Perform hierarchical clustering
dist_matrix <- dist(t(heatmap_matrix))  # Distance on columns
hc <- hclust(dist_matrix)

# Reorder columns based on clustering
ordered_cols <- heatmap_matrix[, hc$order]

# Create a new long format data frame for the reordered heatmap
heatmap_data_ordered <- melt(ordered_cols)
colnames(heatmap_data_ordered) <- c("seurat_clusters", "cell_type", "likelihood")

# Convert seurat_clusters to numeric for proper ordering
heatmap_data_ordered <- heatmap_data_ordered %>%
  mutate(seurat_clusters = as.numeric(as.character(seurat_clusters)))

# Create the heatmap with ordered columns
p <- ggplot(heatmap_data_ordered, aes(x = cell_type, y = factor(seurat_clusters), fill = likelihood)) +  # Use factor for proper ordering
  geom_tile(color = "white") +
  scale_fill_gradient(low = "white", high = "blue") +
  theme_minimal() +
  labs(title = "Per-cluster Cell Types Prediction by SingleR",
       x = "Cell Type",
       y = "Seurat Clusters") +
  scale_y_discrete(expand = c(0, 0)) +  # Ensure y-axis ticks for each row
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.ticks.x = element_line(),
        axis.ticks.y = element_line(),
        panel.background = element_rect(fill = "white"),  
        plot.background = element_rect(fill = "white")
  )  

ggsave(
  # file.path(working_dir, 'singler_per_cluster_heatmap.png'),
  file.path(working_dir, 'singler_per_cluster_heatmap_filtered.png'),
  p,
  height = 6,  
  width = 8
)

heatmap_matrix <- acast(heatmap_data_long, seurat_clusters ~ cell_type, value.var = "count")
heatmap_matrix <- sweep(heatmap_matrix, 2, colSums(heatmap_matrix), FUN = "/")

# Perform hierarchical clustering
dist_matrix <- dist(t(heatmap_matrix))  # Distance on columns
hc <- hclust(dist_matrix)

# Reorder columns based on clustering
ordered_cols <- heatmap_matrix[, hc$order]

# Create a new long format data frame for the reordered heatmap
heatmap_data_ordered <- melt(ordered_cols)
colnames(heatmap_data_ordered) <- c("seurat_clusters", "cell_type", "likelihood")

# Convert seurat_clusters to numeric for proper ordering
heatmap_data_ordered <- heatmap_data_ordered %>%
  mutate(seurat_clusters = as.numeric(as.character(seurat_clusters)))

# Create the heatmap with ordered columns
p <- ggplot(heatmap_data_ordered, aes(x = cell_type, y = factor(seurat_clusters), fill = likelihood)) +  # Use factor for proper ordering
  geom_tile(color = "white") +
  scale_fill_gradient(low = "white", high = "blue") +
  theme_minimal() +
  labs(title = "Per-celltype Cell Types Prediction by SingleR",
       x = "Cell Type",
       y = "Seurat Clusters") +
  scale_y_discrete(expand = c(0, 0)) +  # Ensure y-axis ticks for each row
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.ticks.x = element_line(),
        axis.ticks.y = element_line(),
        panel.background = element_rect(fill = "white"),  
        plot.background = element_rect(fill = "white")
  )  

ggsave(
  # file.path(working_dir, 'singler_per_celltype_heatmap.png'),
  file.path(working_dir, 'singler_per_celltype_heatmap_filtered.png'),
  p,
  height = 6,  
  width = 8
)

table_samples_by_cell_type <- seurat@meta.data %>%
  dplyr::group_by(seurat_clusters, cell_type_singler_blueprintencode_fine) %>%
  dplyr::summarize(count = n()) %>%
  tidyr::spread(cell_type_singler_blueprintencode_fine, count, fill = 0) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(total_cell_count = rowSums(.[c(2:ncol(.))])) %>%
  dplyr::select(c("seurat_clusters", "total_cell_count", dplyr::everything()))

# Prepare the data for the heatmap
heatmap_data <- table_samples_by_cell_type[, -which(names(table_samples_by_cell_type) == "total_cell_count")]

# Reshape to long format
heatmap_data_long <- melt(heatmap_data, id.vars = "seurat_clusters", variable.name = "cell_type", value.name = "count")

# Normalize by row max
heatmap_data_long <- heatmap_data_long %>%
  group_by(seurat_clusters) %>%
  mutate(count = count / sum(count)) %>%
  ungroup()

# Create a matrix for clustering
heatmap_matrix <- acast(heatmap_data_long, seurat_clusters ~ cell_type, value.var = "count")
heatmap_matrix <- sweep(heatmap_matrix, 1, rowSums(heatmap_matrix), FUN = "/")

# attenuate sparsity of columns
heatmap_matrix <- heatmap_matrix[, colSums(heatmap_matrix, na.rm = TRUE) > median(colSums(heatmap_matrix, na.rm = TRUE))]

# secondary normalisation after filtering
# heatmap_matrix <- sweep(heatmap_matrix, 1, rowSums(heatmap_matrix), FUN = "/")

# Perform hierarchical clustering
dist_matrix <- dist(t(heatmap_matrix))  # Distance on columns
hc <- hclust(dist_matrix)

# Reorder columns based on clustering
ordered_cols <- heatmap_matrix[, hc$order]

# Create a new long format data frame for the reordered heatmap
heatmap_data_ordered <- melt(ordered_cols)
colnames(heatmap_data_ordered) <- c("seurat_clusters", "cell_type", "likelihood")

# Convert seurat_clusters to numeric for proper ordering
heatmap_data_ordered <- heatmap_data_ordered %>%
  mutate(seurat_clusters = as.numeric(as.character(seurat_clusters)))

# ## auxiliary visualisation
# # Calculate column sums
# col_sums <- colSums(heatmap_matrix, na.rm = TRUE)
# 
# # Calculate median of column sums
# median_col_sum <- median(col_sums)
# 
# # Create a data frame of column sums for plotting
# col_sums_df <- data.frame(
#   cell_type = names(col_sums),
#   colsum = col_sums
# )
# 
# ggplot(col_sums_df, aes(x = 1, y = colsum)) +
#   geom_violin(fill = "lightblue", alpha = 0.7) +
#   geom_boxplot(width = 0.1, fill = "white", color = "black") +
#   geom_hline(yintercept = median_col_sum, color = "red", linetype = "dashed") +
#   annotate("text", x = 1.2, y = median_col_sum, 
#            label = paste("Median:", round(median_col_sum, 3)), 
#            color = "red", hjust = 0) +
#   theme_minimal() +
#   labs(title = "Distribution of Column Sums",
#        x = "",
#        y = "Column Sum") +
#   theme(axis.text.x = element_blank(),
#         axis.ticks.x = element_blank())

# Create the heatmap with ordered columns
p <- ggplot(heatmap_data_ordered, aes(x = cell_type, y = factor(seurat_clusters), fill = likelihood)) +  # Use factor for proper ordering
  geom_tile(color = "white") +
  scale_fill_gradient(low = "white", high = "blue") +
  theme_minimal() +
  labs(title = "Per-cluster Cell Types Prediction by SingleR",
       x = "Cell Type",
       y = "Seurat Clusters") +
  scale_y_discrete(expand = c(0, 0)) +  # Ensure y-axis ticks for each row
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.ticks.x = element_line(),
        axis.ticks.y = element_line(),
        panel.background = element_rect(fill = "white"),  
        plot.background = element_rect(fill = "white"),
        plot.margin = margin(10, 10, 10, 30, "pt")
  )  

ggsave(
  # file.path(working_dir, 'singler_per_cluster_heatmap_fine.png'),
  file.path(working_dir, 'singler_per_cluster_heatmap_filtered_fine.png'),
  p,
  height = 6,  
  width = 14
)

heatmap_matrix <- acast(heatmap_data_long, seurat_clusters ~ cell_type, value.var = "count")
heatmap_matrix <- sweep(heatmap_matrix, 2, colSums(heatmap_matrix), FUN = "/")
# attenuate sparsity of columns if needed
# heatmap_matrix <- heatmap_matrix[rowSums(heatmap_matrix, na.rm = TRUE) > median(rowSums(heatmap_matrix, na.rm = TRUE)), ]

# ## auxiliary visualisation
# # Calculate column sums
# row_sums <- rowSums(heatmap_matrix, na.rm = TRUE)
# 
# # Calculate median of column sums
# median_row_sum <- median(row_sums)
# 
# # Create a data frame of column sums for plotting
# row_sums_df <- data.frame(
#   cell_type = names(row_sums),
#   rowsum = row_sums
# )
# 
# ggplot(row_sums_df, aes(x = 1, y = rowsum)) +
#   geom_violin(fill = "lightblue", alpha = 0.7) +
#   geom_boxplot(width = 0.1, fill = "white", color = "black") +
#   geom_hline(yintercept = median_row_sum, color = "red", linetype = "dashed") +
#   annotate("text", x = 1.2, y = median_row_sum,
#            label = paste("Median:", round(median_row_sum, 3)),
#            color = "red", hjust = 0) +
#   theme_minimal() +
#   labs(title = "Distribution of Column Sums",
#        x = "",
#        y = "Column Sum") +
#   theme(axis.text.x = element_blank(),
#         axis.ticks.x = element_blank())

# Perform hierarchical clustering
dist_matrix <- dist(t(heatmap_matrix))  # Distance on columns
hc <- hclust(dist_matrix)

# Reorder columns based on clustering
ordered_cols <- heatmap_matrix[, hc$order]

# Create a new long format data frame for the reordered heatmap
heatmap_data_ordered <- melt(ordered_cols)
colnames(heatmap_data_ordered) <- c("seurat_clusters", "cell_type", "likelihood")

# Convert seurat_clusters to numeric for proper ordering
heatmap_data_ordered <- heatmap_data_ordered %>%
  mutate(seurat_clusters = as.numeric(as.character(seurat_clusters)))

# Create the heatmap with ordered columns
p <- ggplot(heatmap_data_ordered, aes(x = cell_type, y = factor(seurat_clusters), fill = likelihood)) +  # Use factor for proper ordering
  geom_tile(color = "white") +
  scale_fill_gradient(low = "white", high = "blue") +
  theme_minimal() +
  labs(title = "Per-celltype Cell Types Prediction by SingleR",
       x = "Cell Type",
       y = "Seurat Clusters") +
  scale_y_discrete(expand = c(0, 0)) +  # Ensure y-axis ticks for each row
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.ticks.x = element_line(),
        axis.ticks.y = element_line(),
        panel.background = element_rect(fill = "white"),  
        plot.background = element_rect(fill = "white"),
        plot.margin = margin(10, 10, 10, 50, "pt")
  )  

ggsave(
  # file.path(working_dir, 'singler_per_celltype_heatmap_fine.png'),
  file.path(working_dir, 'singler_per_celltype_heatmap_filtered_fine.png'),
  p,
  height = 6,  
  width = 24
)

heatmap_matrix <- acast(heatmap_data_long, seurat_clusters ~ cell_type, value.var = "count")

# Perform hierarchical clustering
dist_matrix <- dist(t(heatmap_matrix))  # Distance on columns
hc <- hclust(dist_matrix)

# Reorder columns based on clustering
ordered_cols <- heatmap_matrix[, hc$order]

# Create a new long format data frame for the reordered heatmap
heatmap_data_ordered <- melt(ordered_cols)
colnames(heatmap_data_ordered) <- c("seurat_clusters", "cell_type", "likelihood")

# Convert seurat_clusters to numeric for proper ordering
heatmap_data_ordered <- heatmap_data_ordered %>%
  mutate(seurat_clusters = as.numeric(as.character(seurat_clusters)))

# Create the heatmap with ordered columns
p <- ggplot(heatmap_data_ordered, aes(x = cell_type, y = factor(seurat_clusters), fill = likelihood)) +  # Use factor for proper ordering
  geom_tile(color = "white") +
  scale_fill_gradient(low = "white", high = "blue") +
  theme_minimal() +
  labs(title = "Per-celltype Cell Types Prediction by SingleR",
       x = "Cell Type",
       y = "Seurat Clusters") +
  scale_y_discrete(expand = c(0, 0)) +  # Ensure y-axis ticks for each row
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.ticks.x = element_line(),
        axis.ticks.y = element_line(),
        panel.background = element_rect(fill = "white"),  
        plot.background = element_rect(fill = "white"),
        plot.margin = margin(10, 10, 10, 50, "pt")
  )  

ggsave(
  # file.path(working_dir, 'singler_per_celltype_heatmap_fine.png'),
  file.path(working_dir, 'singler_heatmap_filtered_fine.png'),
  p,
  height = 6,  
  width = 24
)

# Prepare the data for the heatmap
heatmap_data <- table_samples_by_cell_type[, -which(names(table_samples_by_cell_type) == "total_cell_count")]

# Reshape to long format
heatmap_data_long <- melt(heatmap_data, id.vars = "seurat_clusters", variable.name = "cell_type", value.name = "count")

# Normalize by row total
heatmap_data_long <- heatmap_data_long %>%
  group_by(seurat_clusters) %>%
  mutate(count = count / sum(count)) %>%  # Normalize by total for each cluster
  ungroup()

# Create a matrix for clustering
heatmap_matrix <- acast(heatmap_data_long, seurat_clusters ~ cell_type, value.var = "count")

# Perform hierarchical clustering
dist_matrix <- dist(t(heatmap_matrix))  # Distance on columns
hc <- hclust(dist_matrix)

# Reorder columns based on clustering
ordered_cols <- heatmap_matrix[, hc$order]

# Create a new long format data frame for the reordered heatmap
heatmap_data_ordered <- melt(ordered_cols)
colnames(heatmap_data_ordered) <- c("seurat_clusters", "cell_type", "likelihood")

# Convert seurat_clusters to numeric for proper ordering
heatmap_data_ordered <- heatmap_data_ordered %>%
  mutate(seurat_clusters = as.numeric(as.character(seurat_clusters)))

# Create the heatmap with ordered columns
p <- ggplot(heatmap_data_ordered, aes(x = cell_type, y = factor(seurat_clusters), fill = likelihood)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "white", high = "blue") +
  theme_minimal() +
  labs(title = "Per-cluster Cell Types Prediction by SingleR",
       x = "Cell Type",
       y = "Seurat Clusters") +
  scale_y_discrete(expand = c(0, 0)) +  # Ensure y-axis ticks for each row
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.ticks.x = element_line(),
        axis.ticks.y = element_line(),
        panel.background = element_rect(fill = "white"),  
        plot.background = element_rect(fill = "white")
  )  

# Save the heatmap
ggsave(
  file.path(working_dir, 'singler_heatmap.png'),
  p,
  height = 6,  
  width = 8
)

table_samples_by_cell_type %>% knitr::kable()

# |sample | total_cell_count| Astrocyte| B_cell| Chondrocytes| DC| Endothelial_cells| Epithelial_cells| Fibroblasts| Gametocytes| Hepatocytes| HSC_-G-CSF| Keratinocytes| Macrophage| Monocyte| Neuroepithelial_cell| Neurons| Neutrophils| NK_cell| Osteoblasts| Platelets| Smooth_muscle_cells| T_cells| Tissue_stem_cells|
# |:------|----------------:|---------:|------:|------------:|--:|-----------------:|----------------:|-----------:|-----------:|-----------:|----------:|-------------:|----------:|--------:|--------------------:|-------:|-----------:|-------:|-----------:|---------:|-------------------:|-------:|-----------------:|
# |Brain  |              470|        14|      1|           10| 44|                 6|               22|           1|           9|           2|          1|            50|          3|       23|                    2|       6|         155|      25|          76|         1|                   6|       2|                11|


table_clusters_by_cell_type <- seurat@meta.data %>%
  dplyr::group_by(seurat_clusters, cell_type_singler_blueprintencode_main) %>%
  dplyr::rename(cluster = seurat_clusters) %>%
  dplyr::summarize(count = n()) %>%
  tidyr::spread(cell_type_singler_blueprintencode_main, count, fill = 0) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(total_cell_count = rowSums(.[c(2:ncol(.))])) %>%
  dplyr::select(c("cluster", "total_cell_count", dplyr::everything()))

table_clusters_by_cell_type %>% knitr::kable()

# |cluster | total_cell_count| Astrocyte| B_cell| Chondrocytes| DC| Endothelial_cells| Epithelial_cells| Fibroblasts| Gametocytes| Hepatocytes| HSC_-G-CSF| Keratinocytes| Macrophage| Monocyte| Neuroepithelial_cell| Neurons| Neutrophils| NK_cell| Osteoblasts| Platelets| Smooth_muscle_cells| T_cells| Tissue_stem_cells|
# |:-------|----------------:|---------:|------:|------------:|--:|-----------------:|----------------:|-----------:|-----------:|-----------:|----------:|-------------:|----------:|--------:|--------------------:|-------:|-----------:|-------:|-----------:|---------:|-------------------:|-------:|-----------------:|
# |0       |              154|         8|      0|            7|  3|                 1|               18|           1|           6|           1|          0|            12|          1|        4|                    2|       0|          40|       5|          35|         0|                   2|       0|                 8|
# |1       |              119|         0|      1|            1| 22|                 3|                3|           0|           2|           1|          0|            11|          1|       10|                    0|       0|          36|      14|           9|         1|                   1|       2|                 1|
# |2       |               75|         5|      0|            1|  5|                 0|                0|           0|           0|           0|          0|             0|          0|        1|                    0|       6|          25|       0|          28|         0|                   3|       0|                 1|
# |3       |               62|         1|      0|            1|  5|                 1|                0|           0|           0|           0|          0|            19|          0|        4|                    0|       0|          23|       5|           2|         0|                   0|       0|                 1|
# |4       |               36|         0|      0|            0|  5|                 1|                0|           0|           0|           0|          0|             0|          0|        2|                    0|       0|          26|       1|           1|         0|                   0|       0|                 0|
# |5       |               24|         0|      0|            0|  4|                 0|                1|           0|           1|           0|          1|             8|          1|        2|                    0|       0|           5|       0|           1|         0|                   0|       0|                 0|
# 

UMAP_centers_cluster <- tibble(
  UMAP_1 = as.data.frame(seurat@reductions$UMAP@cell.embeddings)$UMAP_1,
  UMAP_2 = as.data.frame(seurat@reductions$UMAP@cell.embeddings)$UMAP_2,
  cluster = seurat@meta.data$cell_type_singler_blueprintencode_main
) %>%
  group_by(cluster) %>%
  dplyr::summarize(x = median(UMAP_1), y = median(UMAP_2))

# Create the UMAP plot with repelled labels
plot_umap_by_cluster <- bind_cols(seurat@meta.data, as.data.frame(seurat@reductions$UMAP@cell.embeddings)) %>%
  ggplot(aes(UMAP_1, UMAP_2, color = seurat_clusters)) +
  geom_point(size = 0.2) +
  geom_label_repel(  # Use geom_label_repel instead of geom_label
    data = UMAP_centers_cluster,
    mapping = aes(x, y, label = cluster),
    size = 4.5,
    fill = 'white',
    color = 'black',
    fontface = 'bold',
    alpha = 0.7,
    label.size = 0,
    show.legend = FALSE,
    box.padding = 0.5,  # Adjust padding around labels
    point.padding = 0.5,  # Adjust padding between labels and points
    segment.color = NA  # Color of the connecting line
  ) +
  theme_bw() +
  scale_color_manual(
    name = 'Cluster', values = custom_colors$discrete,
    guide = guide_legend(override.aes = list(size = 2))
  ) +
  theme(legend.position = 'right') +
  coord_fixed() +
  annotate(
    geom = 'text', x = Inf, y = -Inf,
    label = paste0('n = ', format(nrow(seurat@meta.data), big.mark = ',', trim = TRUE)),
    vjust = -1.5, hjust = 1.25, color = 'black', size = 2.5
  )



tSNE_centers_cluster <- tibble(
  tSNE_1 = as.data.frame(seurat@reductions$tSNE@cell.embeddings)$tSNE_1,
  tSNE_2 = as.data.frame(seurat@reductions$tSNE@cell.embeddings)$tSNE_2,
  cluster = seurat@meta.data$cell_type_singler_blueprintencode_main
) %>%
  group_by(cluster) %>%
  dplyr::summarize(x = median(tSNE_1), y = median(tSNE_2))

plot_tsne_by_cluster <- bind_cols(seurat@meta.data, as.data.frame(seurat@reductions$tSNE@cell.embeddings)) %>%
  ggplot(aes(tSNE_1, tSNE_2, color = seurat_clusters)) +
  geom_point(size = 0.2) +
  geom_label_repel(  # Use geom_label_repel instead of geom_label
    data = tSNE_centers_cluster,
    mapping = aes(x, y, label = cluster),
    size = 4.5,
    fill = 'white',
    color = 'black',
    fontface = 'bold',
    alpha = 0.7,
    label.size = 0,
    show.legend = FALSE,
    box.padding = 0.5,  # Adjust padding around labels
    point.padding = 0.5,  # Adjust padding between labels and points
    segment.color = NA  # Color of the connecting line
  ) +
  theme_bw() +
  scale_color_manual(
    name = 'Cluster', values = custom_colors$discrete,
    guide = guide_legend(override.aes = list(size = 2))
  ) +
  theme(legend.position = 'right') +
  coord_fixed() +
  annotate(
    geom = 'text', x = Inf, y = -Inf,
    label = paste0('n = ', format(nrow(seurat@meta.data), big.mark = ',', trim = TRUE)),
    vjust = -1.5, hjust = 1.25, color = 'black', size = 2.5
  )

ggsave(
  file.path(working_dir, 'singler_annotation_main.brain.png'),
  plot_umap_by_cluster + plot_tsne_by_cluster,
  height = 6,  
  width = 12
)

# ### in silico perturbation
# library(muscat)
# library(dplyr)
# seurat <- readRDS("/home/jyubq/proj_botryllus/expressionCount/230629_00_testis/botryllus_testis_preprocessed.rds")
# expr_matrix <- as.matrix(seurat@assays$RNA$counts)
# metadata <- seurat@meta.data
# gene_of_interest <- c("FUN_009332-T1.exon1", "FUN_009332-T1.exon2", "FUN_009332-T1.exon3", "FUN_009332-T1.exon4", "FUN_009332-T1.exon5", "FUN_009333-T1.exon1")
# expr_matrix[gene_of_interest, ] <- 0
# perturbed_seurat <- CreateSeuratObject(counts = expr_matrix, meta.data = metadata)
# sce <- as.SingleCellExperiment(perturbed_seurat)
# de_results <- muscat::pbDS(sce, design = ~ condition, contrast = "conditionperturbed - conditioncontrol")
# summary(de_results)

# ## Save the updated Seurat object
# saveRDS(seurat, "/home/jyubq/proj_botryllus/expressionCount/230629_00_testis/botryllus_testis_processed.rds")

# ### read loom and check if velocyto works, doesn't work 
# loom <- connect(filename = "/scratch/users/jiamuyu/proj_botryllus/scRNAseq/240804_00_botryllus_ss3_new_velos/multi_input_PAC10_Plate1_A01_ILW92_Testes_A1_Zooid_STDandSZ_1_p01c01r01_7315b_196Days_S5761_dedup_B3ZLZ.loom", mode = "r")

### the 10x tutorial way to convert seurat objects to loom files, functional
seurat.loom <- as.loom(seurat, filename = file.path(working_dir, "cytotrace2_result.single.loom"), verbose = FALSE)
seurat.loom$close_all()

cluster_data <- data.frame(
  Barcode = colnames(seurat),
  Cluster = paste0("Cluster ", seurat@meta.data$seurat_clusters)
)

write.csv(
  cluster_data,
  file = file.path(working_dir, "cluster_identites.single.csv"),
  row.names = FALSE,
  quote = FALSE
)

umap_coords <- seurat@reductions$UMAP@cell.embeddings

umap_data <- data.frame(
  Barcode = rownames(umap_coords),
  `X Coordinate` = umap_coords[, 1],
  `Y Coordinate` = umap_coords[, 2]
)

write.csv(
  umap_data,
  file = file.path(working_dir, "umap_coordinates.single.csv"),
  row.names = FALSE,
  quote = FALSE
)

### extract GFP
# Define the absolute paths to your 10X data directories
base_dir <- "/scratch/groups/ayeletv/scRNA_Ciona_Nov2024/usftp21.novogene.com/03.CountData.250521"

data_dirs <- list(
  "TG_CONTROL" = file.path(base_dir, "TG_CONTROL/TG_CONTROL/outs/filtered_feature_bc_matrix"),
  "TG_STEM_B_1" = file.path(base_dir, "TG_STEM_B_1/TG_STEM_B_1/outs/filtered_feature_bc_matrix"), 
  "TG_STEM_B_2" = file.path(base_dir, "TG_STEM_B_2/TG_STEM_B_2/outs/filtered_feature_bc_matrix"),
  "TG_ALDH" = file.path(base_dir, "TG_ALDH/TG_ALDH/outs/filtered_feature_bc_matrix"),
  "TG_STEM_A" = file.path(base_dir, "TG_STEM_A/TG_STEM_A/outs/filtered_feature_bc_matrix"),
  "Undetermined" = file.path(base_dir, "Undetermined/Undetermined/outs/filtered_feature_bc_matrix")
)

# Sample name mapping (based on your existing cell names)
sample_mapping <- list(
  "TG_CONTROL" = "CONTROL",
  "TG_STEM_B_1" = "STEM_B_1", 
  "TG_STEM_B_2" = "STEM_B_2",
  "TG_ALDH" = "ALDH",
  "TG_STEM_A" = "STEM_A",
  "Undetermined" = "Undetermined"
)

extract_gfp_expression <- function() {
  all_gfp_values <- c()
  all_gfp_names <- c()
  
  for(sample_dir in names(data_dirs)) {
    # Read 10X data
    data_10x <- Read10X(data_dirs[[sample_dir]])
    
    # Check if KY21.GFP.1 exists in the data
    if("KY21.GFP.1" %in% rownames(data_10x)) {
      # Extract GFP expression
      gfp_expr <- data_10x["KY21.GFP.1", ]
      
      # Modify cell barcodes to match your Seurat object format
      sample_name <- sample_mapping[[sample_dir]]
      corrected_names <- paste0(sample_name, "_", names(gfp_expr))
      
      # Combine values and names
      all_gfp_values <- c(all_gfp_values, as.numeric(gfp_expr))
      all_gfp_names <- c(all_gfp_names, corrected_names)
      
      cat("Found KY21.GFP.1 in", sample_dir, "with", length(gfp_expr), "cells\n")
      cat("Sample cell name format:", head(corrected_names, 3), "\n")
    } else {
      cat("KY21.GFP.1 not found in", sample_dir, "\n")
      cat("Available features:\n")
      print(head(rownames(data_10x), 10))
    }
  }
  
  # Create named vector
  names(all_gfp_values) <- all_gfp_names
  return(all_gfp_values)
}

all_gfp <- extract_gfp_expression()

# Check the naming after correction
cat("Sample GFP names after correction:\n")
print(head(names(all_gfp)))

# Get current cell names from your Seurat object
current_cells <- colnames(seurat)

# Match GFP data to your current cells
matched_gfp <- rep(0, length(current_cells))  # Initialize with zeros
names(matched_gfp) <- current_cells

# Find overlapping cells and assign GFP values
overlap_cells <- intersect(names(all_gfp), current_cells)
matched_gfp[overlap_cells] <- all_gfp[overlap_cells]

# Check matching results
cat("Number of overlapping cells:", length(overlap_cells), "\n")
cat("Total cells in Seurat object:", length(current_cells), "\n")
cat("Total cells with GFP data:", length(all_gfp), "\n")

# Add GFP expression to your Seurat object
seurat[["GFP"]] <- matched_gfp

# Check the results
cat("Added GFP expression for", sum(matched_gfp > 0), "out of", length(current_cells), "cells\n")
cat("GFP expression range:", range(matched_gfp), "\n")

# Optional: Create a feature plot to visualize GFP expression
p <- FeaturePlot(seurat, features = "GFP", order = TRUE, reduction = "UMAP", cols = colorRampPalette(c("#FFF2F2", "#990000"))(2), pt.size = 1.5) +
  ggtitle("GFP") +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(
  file.path(working_dir, 'final/featureplot.GFP.pdf'),
  p, 
  height = 6, width = 6
)

### GFP violin
gfp_data <- data.frame(
  combined_sample = seurat@meta.data$merged_sample,
  expr_value = seurat@meta.data$GFP
)

p <- ggstatsplot::ggbetweenstats(
  data = gfp_data,
  x = combined_sample,
  y = expr_value,
  type = "nonparametric",
  pairwise.comparisons = TRUE,
  pairwise.display = "significant",
  centrality.plotting = TRUE,
  centrality.point.args = list(size = 0, color = "darkred"),
  centrality.label.args = list(alpha = 0),
  p.adjust.method = "bonferroni",
  palette = "Set1",
  package = "RColorBrewer",
  violin.args = list(width = 0.8, alpha = 0.7),
  point.args = list(position = ggplot2::position_jitterdodge(jitter.width = 0.8), alpha = 0.4, size = 1.5, stroke = 0, na.rm = TRUE),
  title = "GFP",
  results.subtitle = FALSE,
  xlab = "",
  ylab = "Expression",
  ggplot.component = list(theme(axis.title.y.right = element_blank(), 
                                axis.text.y.right = element_blank(), 
                                axis.ticks.y.right = element_blank()))
)

ggsave(
  file.path(working_dir, 'final/violin.4.GFP.pdf'),
  p, 
  height = 8, width = 4
)


### GFP violin enhanced
gfp_data <- data.frame(
  combined_sample = seurat@meta.data$merged_sample,
  expr_value = seurat@meta.data$GFP
)

# Calculate statistics for circles
circle_stats <- gfp_data %>%
  group_by(combined_sample) %>%
  dplyr::summarise(
    total_cells = n(),
    expressing_cells = sum(expr_value > 0, na.rm = TRUE),
    proportion_expressing = expressing_cells / total_cells,
    mean_expression = mean(expr_value[expr_value > 0], na.rm = TRUE),
    percentage_text = paste0(round(proportion_expressing * 100, 1), "%")
  ) %>%
  mutate(
    mean_expression = ifelse(is.nan(mean_expression), 0, mean_expression)
  )

print(circle_stats)

# Create the main violin plot
p_violin <- ggstatsplot::ggbetweenstats(
  data = gfp_data,
  x = combined_sample,
  y = expr_value,
  type = "nonparametric",
  pairwise.comparisons = TRUE,
  pairwise.display = "significant",
  centrality.plotting = TRUE,
  centrality.point.args = list(size = 0, color = "darkred"),
  centrality.label.args = list(alpha = 0),
  p.adjust.method = "bonferroni",
  palette = "Set1",
  package = "RColorBrewer",
  violin.args = list(width = 0.8, alpha = 0.7),
  point.args = list(position = ggplot2::position_jitterdodge(jitter.width = 0.8), alpha = 0.4, size = 1.5, stroke = 0, na.rm = TRUE),
  title = "GFP",
  results.subtitle = FALSE,
  xlab = "",
  ylab = "Expression",
  ggplot.component = list(theme(axis.title.y.right = element_blank(), 
                                axis.text.y.right = element_blank(), 
                                axis.ticks.y.right = element_blank()))
)

# Create circle plot
create_circle_plot <- function(stats_df) {
  # Calculate color scale limits (non-normalized)
  min_expr <- floor(min(stats_df$mean_expression, na.rm = TRUE))
  max_expr <- ceiling(max(stats_df$mean_expression, na.rm = TRUE))
  
  # Create data for circles
  circle_data <- data.frame()
  text_data <- data.frame()
  
  for(i in 1:nrow(stats_df)) {
    sample_name <- stats_df$combined_sample[i]
    prop <- stats_df$proportion_expressing[i]
    mean_expr <- stats_df$mean_expression[i]  # Use actual expression values
    pct_text <- stats_df$percentage_text[i]
    
    # Outer circle (white, constant size)
    theta <- seq(0, 2*pi, length.out = 100)
    outer_circle <- data.frame(
      x = i + 0.3 * cos(theta),
      y = 0.3 * sin(theta),
      sample = sample_name,
      circle_type = "outer",
      mean_expression = NA
    )
    
    # Inner circle (colored, proportional size)
    inner_radius <- 0.25 * sqrt(prop)  # Scale by sqrt for area proportionality
    inner_circle <- data.frame(
      x = i + inner_radius * cos(theta),
      y = inner_radius * sin(theta),
      sample = sample_name,
      circle_type = "inner",
      mean_expression = mean_expr  # Use actual values, not normalized
    )
    
    # Now both data frames have the same columns
    circle_data <- rbind(circle_data, outer_circle, inner_circle)
    
    # Text data
    text_data <- rbind(text_data, data.frame(
      x = i,
      y = 0,
      label = pct_text,
      sample = sample_name
    ))
  }
  
  # Create the circle plot
  p_circles <- ggplot() +
    # Outer circles (white with black border)
    geom_polygon(data = circle_data[circle_data$circle_type == "outer", ],
                 aes(x = x, y = y, group = sample),
                 fill = "white", color = "#84817a", size = 1) +
    # Inner circles (colored)
    geom_polygon(data = circle_data[circle_data$circle_type == "inner", ],
                 aes(x = x, y = y, group = sample, fill = mean_expression),
                 color = "#0652DD", size = 0.5) +
    geom_text(data = text_data,
              aes(x = x, y = y, label = label),
              size = 3, fontface = "bold", 
              color = "black") +  # Black text layer
    scale_fill_gradient(low = "lightblue", high = "#0652DD", 
                        name = "Mean GFP Expression",
                        na.value = "transparent",
                        limits = c(min_expr, max_expr),  # Non-normalized scale
                        breaks = seq(min_expr, max_expr, by = 1)) +
    scale_x_continuous(breaks = 1:nrow(stats_df), 
                       labels = stats_df$combined_sample) +
    coord_fixed() +
    theme_minimal() +
    theme(
      axis.title = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid = element_blank(),
      # axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.box = "horizontal",
      legend.box.just = "center",
      plot.caption = element_text(hjust = 0.5, size = 10, margin = margin(t = 10))
    ) +
    guides(fill = guide_colorbar(title.position = "top", 
                                 title.hjust = 0.5,
                                 barwidth = 8, 
                                 barheight = 0.8)) +
    labs(caption = "Circle size = % expressing GFP  |  Color = Mean GFP expression")
  
  return(p_circles)
}

p_circles <- create_circle_plot(circle_stats)

ggsave(
  file.path(working_dir, 'final/violin.4.GFP.png'),
  p_violin / p_circles + 
    plot_layout(heights = c(3, 0.5)), 
  height = 10, width = 6
)
