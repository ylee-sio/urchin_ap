hugo = tibble(
  
  GeneID = c("merged-gene-150", "merged-gene-151", "merged-gene-152", "merged-gene-154",
             "LOC582052", "LOC583249", "LOC592380", "LOC105440521", "LOC590740"),
  Name = c("ABCA1", "ABCA2", "ABCA3", "ABCA5",
           "VLDLR", "VPS26", "PICALM_transcript_variant_1", "PICALM_transcript_variant_2", "CD2AP")
)

t2 = map2(.x = urchin_5_clustered_labelled_seurat_objs,
          .y = stage_names,
         .f = function(x,y)(
           
           lab_DP(x, 
                  hugo,
                  plot_title = y,
                  F, 
                  500,
                  F,
                  -2,
                  3, 
                  0,
                  cols = c(
                    "grey",
                    "blue"
                  )
                  )
           )
         )

pdf("~/sp_transportomics/out.dump/bellen_genes.pdf", onefile = T, width = 8, height = 6)
t2
dev.off()

hugo_without_cd2ap = tibble(
  
  GeneID = c("merged-gene-150", "merged-gene-151", "merged-gene-152", "merged-gene-154",
             "LOC582052", "LOC583249", "LOC592380", "LOC105440521"),
  Name = c("ABCA1", "ABCA2", "ABCA3", "ABCA5",
           "VLDLR", "VPS26", "PICALM_transcript_variant_1", "PICALM_transcript_variant_2")
)

t2 = map2(.x = urchin_5_clustered_labelled_seurat_objs,
          .y = stage_names,
          .f = function(x,y)(
            
            lab_DP(x, 
                   hugo_without_cd2ap,
                   plot_title = y,
                   F, 
                   500,
                   F,
                   -2,
                   3, 
                   0,
                   cols = c(
                     "grey",
                     "blue"
                   )
            )
          )
)

pdf("~/sp_transportomics/out.dump/bellen_genes_without_cd2ap.pdf", onefile = T, width = 8, height = 6)
t2
dev.off()



