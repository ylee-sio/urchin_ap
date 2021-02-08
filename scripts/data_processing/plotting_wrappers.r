output_path = "~/Projects/urchin_ap/output/"

# parallel plotting of DimPlots
# requires lists of seurat_objs

plot_multiple_dp = function(cluster_list, file_name){

# creates directory for each output
  output_dir_name = paste0(output_path, file_name, "/")
  dir.create(output_dir_name)
  dimplot_list = future_map(cluster_list, DimPlot)
  plots = CombinePlots(dimplot_list)

# creates pdf of DimPlots
  pdf(file = paste0(output_dir_name, file_name, "_dimplot.pdf", height = 12, width = 12)
  plots
  dev.off()

  return(plots)

}

