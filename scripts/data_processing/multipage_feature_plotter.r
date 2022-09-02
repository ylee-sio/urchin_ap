library(future)
library(furrr)
library(gridExtra)
library(tictoc)

# read in clustered, saved list of seurat objects
all_stages = readRDS("data/working_data/all_stages.clustered_seurat_obj_list.rds")

# create list of names to be referenced for naming separate files for each stage
stage_names = c(
  "8 Cell",
  "64 Cell",
  "Morula",
  "Early Blastula",
  "Mesenchymal Blastula",
  "Hatched Blastula",
  "Early Gastrula",
  "Late Gastrula"
)

# filter for features of defined genesets which have value > 0 in entire dataset
assay_data_counts = map(
  .x = all_stages,
  .f = function(x)(
    GetAssayData(
      object = x,
      slot = "counts"
      )
  )
)

# find which features have no expression
positive_features_each_stage = map(
  .x = assay_data_counts,
  .f = function(x)(
    which(rowSums(x) > 0) %>% 
      names()
  )
)

# create seurat objects filtered for positive features
stages_filtered_positive_features = map2(
  .x = all_stages,
  .y = positive_features_each_stage,
  .f = function(x,y)(
    subset(
      x = x,
      features = y
    )
  )
)

splg = stages_filtered_positive_features[[8]]
test_smt_geneid = rownames(splg)[rownames(splg) %in% sp_kegg_smts$GeneID %>% which()]
test_smts_df = sp_kegg_smts[sp_kegg_smts$GeneID %in% test_smt_geneid %>% which(),] %>% unique()
test_smts_df_with_plot_names = test_smts_df %>% mutate(plot_name = paste0(Name, "(", GeneID, ")"))

# plots overlays of genes in genesets of interest with a known marker. Takes about one minute
# per 500 genes
test_features = map(
  .x = test_smts_df_with_plot_names$GeneID,
  .f = function(x)(
    FeaturePlot(
      splg,
      features = c(x, "Endo16"),
      blend = TRUE,
      order = TRUE
      ) 
    )
  )

# fixes labels for each feature plots. takes less than one minute

test_features_fixed_labels = map2(
  .x = test_features,
  .y = test_smts_df_with_plot_names$plot_name,
  .f = function(x,y)(
    x = (x[[1]] + ggtitle(y)) + 
      x[[2]] +
      x[[3]] + ggtitle("Overlay") +
      x[[4]]
  )
)

# creates a single pdf file containing all overlays. this will take 6-8 minutes for a 
# target geneset of around 500
pdf(onefile = T, width = 8, height = 8, file = "splg_end16_smt_overlay.pdf")
test_features_fixed_labels
dev.off()


