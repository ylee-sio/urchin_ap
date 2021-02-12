library(tidyverse)
library(Seurat)
library(ggplotly)

# standardizes output directory
output_path = "~/Projects/urchin_ap/output/"

# standardizes output annotation
annotate_vis = function(output_dir_name){
  annotation = readline(prompt="What is main purpose of this plot? ")
  txtfile = file(paste0(output_dir_name, "annotation.txt"))
  writeLines(annotation, txtfile)
  close(txtfile)
}

# ensures no replacement of directories when name inputs are identical to existing directories
append_vis_num = function(output_dir_name){
  directory_counts = length(list.files(output_dir_name))
  new_directory_label = directory_counts + 1
  new_directory_label = paste0(output_dir_name, "_", as.character(new_directory_label),"/")
  return(new_directory_label)
}

# parallel plotting of DimPlots
# requires lists of seurat_objs
plot_multiple_dimplots = function(cluster_list, file_name){

# creates unique directory for each output  
  unlabeled_output_dir_name = paste0(output_path, file_name)
  output_dir_name = append_vis_num(unlabeled_output_dir_name)
  dir.create(output_dir_name)

  dimplot_list = future_map(cluster_list, DimPlot)
  plots = CombinePlots(dimplot_list)

# non-optional annotation
  annotate_vis(output_dir_name = output_dir_name)
# creates pdf of DimPlots
  pdf(file = paste0(output_dir_name, file_name, "_dimplot.pdf", height = 12, width = 12)
  plots
  dev.off()

  
  return(plots)

}

# generates DotPlots which have labeled features along with option for plotly interactivity and export of interactive plots  
labelled_dotplot = function(scrna_df, feature_df, plot_title, by_stage = F, plot_height, interactive, col.min, col.max, dot.min, cols, file_name){
  
  
  dp = DotPlot(scrna_df, features = feature_df$GeneID, col.min = col.min, col.max = col.max, dot.min = dot.min, cols = cols) + 
    scale_x_discrete(breaks=c(feature_df$GeneID),labels=c(feature_df$Name)) +
    RotatedAxis() +
    ggtitle(as.character(plot_title)) + 
    coord_flip()
  
    if (interactive == T){
      
      dp = ggplotly(dp, height = plot_height)
    }
    if (by_stage == T){
    # in order to create dotplots by stage, scrna_df must be an seurat_obj that has multiple timepoints integrated  
      dp = DotPlot(scrna_df, features = feature_df$GeneID, group.by = "orig.ident") + 
        scale_x_discrete(breaks=c(feature_df$GeneID),
                         labels=c(feature_df$Name)) +
        RotatedAxis() +
        ggtitle(as.character(plot_title)) + 
        coord_flip() 
    }
# creates unique directory for each output
  unlabeled_output_dir_name = paste0(output_path, file_name)
  output_dir_name = append_vis_num(unlabeled_output_dir_name)
  dir.create(output_dir_name)

# non-optional annotation
  annotate_vis(output_dir_name = output_dir_name)

# creates pdf for static DotPlot
  pdf(file = paste0(output_dir_name, file_name, "_dotplot", ".pdf"), height = 24, width=12)
  dp
  dev.off()
  
  htmlwidgets::saveWidget(as_widget(dp), paste0(output_dir_name, file_name,"_index.html"))

  return(dp)

}  
