#find translated SMTs
#extract SMTs with positive log2FC values for both puromycin detected translation events and selective polysomal translation events

maternal_translational_1 = read_csv("~/sp_transportomics/data_sources/primary/lv_translation_data/chasse_maternal_mrna_translational_stats.csv")

maternal_translational_2 = 
  maternal_translational_1 %>% 
  transmute(transcript_id = `Transcript ID`,
            annotation = Annotation,
            logfc_uf_vs_ufpuro = `logFC (UnF_vs_UnFpuro)`,
            logfc_f_vs_fpuro = `logFC (F_vs_Fpuro)`,
            logfc_f_vs_unf_polysome_recruitment = `logFC (F_vs_UnF)`
            )

maternal_translational_3a = maternal_translational_2[str_which(maternal_translational_2$annotation, "ABC"),]
maternal_translational_3b = maternal_translational_2[str_which(maternal_translational_2$annotation, "solute"),]

maternal_translational_4 = bind_rows(maternal_translational_3a, maternal_translational_3b)
maternal_translational_5 = subset(maternal_translational_4, logfc_f_vs_fpuro > 0 & logfc_f_vs_unf_polysome_recruitment > 0)

#maternal_translational_puro_positive
maternal_translational_puro_positive = read_csv("~/sp_transportomics/data_sources/secondary/maternal_translational_puro_positive.csv")

maternal_translational_puro_positive =
  maternal_translational_puro_positive %>%
  mutate(GeneID = maternal_translational_puro_positive$annotation %>% 
           str_extract_all("(?i)ABC[:alpha:][:digit:]*|(?i)SLC[:digit:]*[:alpha:][:digit:]*|NBC") %>%
           toupper()
)

#number of SMTs translated; confirmed by puromycin assay
maternal_translational_puro_positive_genes = 
  maternal_translational_puro_positive$GeneID %>% 
  unique() 

maternal_translational_puro_positive_genes[c(1,14,26,54,56)] = c("ABCC4", "ABCB1", "ABCB4", "ABCC9", "ABCC1")
puromycin_labelled_smt_assay_transcripts = urchin_5_smts_final_cleaned_with_mito_smts$Name[urchin_5_smts_final_cleaned_with_mito_smts$Name %in% maternal_translational_puro_positive_genes %>% which()]
puromycin_labelled_smt_assay_transcripts_df = urchin_5_smts_final_cleaned_with_mito_smts[which(urchin_5_smts_final_cleaned_with_mito_smts$Name %in% puromycin_labelled_smt_assay_transcripts),]

#maternal_translational_poly_recruitment_positive
maternal_translational_poly_recruitment_positive = read_csv("~/sp_transportomics/data_sources/secondary/maternal_translational_poly_recruitment_positive.csv")

maternal_translational_poly_recruitment_positive = 
  maternal_translational_poly_recruitment_positive %>% 
  mutate(
    GeneID = maternal_translational_poly_recruitment_positive$annotation %>% 
      str_extract_all("(?i)ABC[:alpha:][:digit:]*|(?i)SLC[:digit:]*[:alpha:][:digit:]*|NBC") %>% 
      toupper()
  )

#number of SMTs translated; confirmed by fractionation of polysomes
maternal_translational_poly_recruitment_positive_genes = 
  maternal_translational_poly_recruitment_positive$GeneID %>% 
  unique()

maternal_translational_poly_recruitment_positive_genes[c(1,5,47,52)] = c("ABCC4", "ABCB4", "ABCB1", "ABCC9")
polysomally_recruited_transcripts = urchin_5_smts_final_cleaned_with_mito_smts$Name[urchin_5_smts_final_cleaned_with_mito_smts$Name %in% maternal_translational_poly_recruitment_positive_genes %>% which()]
polysomally_recruited_transcripts_df = urchin_5_smts_final_cleaned_with_mito_smts[which(urchin_5_smts_final_cleaned_with_mito_smts$Name %in% polysomally_recruited_transcripts),]
#Sp1 transcriptional data


#translation factors from kegg
urchin_5_translation_factors = 
  urchin_5_keys[[1]] %>% 
  subset(
    GeneSet == "translation_factors"
  ) %>% 
  select(
    GeneID,
    Name
  )

#smts from kegg that are plasma membrane
urchin_5_smts_final_cleaned = 
  left_join(
    abc_slc_3_geneid_final_cleaned,
    urchin_5_smts,
    by = "Name"
  ) %>% 
  transmute(
    GeneID = GeneID.y,
    Name = Name)

#distirbution of translation factors in the 8-cell stage - pre-cell type merging
urchin_5_sp1_dp = lab_DP(urchin_5_sp1_clustered, 
                         urchin_5_smts_final_cleaned_with_mito_smts,
                         plot_title = "SMT expression profiles of cell types at the 8 cell stage",
                         F, 
                         1000,
                         T,
                         0,
                         3, 
                         0,
                         cols = c("grey",
                                  "blue")
)

#cell type merging and labeling
urchin_5_sp1_cell_type_labelled = RenameIdents(urchin_5_sp1_clustered,
                                               '0' = "MTH1",
                                               '1' = "MTH1",
                                               '3' = "MTH2",
                                               '4' = "MTH2",
                                               '2' = "MTL"
                                               
)

#distribution of smt transcripts at the 8 cell stage; post cell type merging and labeling
urchin_5_sp1_cell_type_labelled_dp = 
  lab_DP(urchin_5_sp1_cell_type_labelled, 
         urchin_5_smts_final_cleaned_with_mito_smts,
         plot_title = "SMT expression profiles of cell types at the 8 cell stage",
         F, 
         4000,
         T,
         0,
         3, 
         0,
         cols = c("grey",
                  "blue")
  )

#RES: Three clusters total; 
#RES: Two types of clusters: high mitochondrial activity, low mitochondrial activity; 
#RES: low mitochondrial activity cluster is greater than 2-fold enriched for 90% of SMTs


#distribution of translation factor transcripts at the 8 cell stage; post cell type merging and labeling
lab_DP(urchin_5_sp1_cell_type_labelled, 
       urchin_5_translation_factors,
       plot_title = "SMT expression profiles of cell types at the 8 cell stage",
       F, 
       1500,
       T,
       0,
       3, 
       0,
       cols = c("grey",
                "blue")
)

#RES: High mitochondrial activity clusters are differentially enriched for translation factors

lab_DP(urchin_5_sp1_cell_type_labelled, 
       puromycin_labelled_smt_assay_transcripts_df,
       plot_title = "SMT expression profiles of cell types at the 8 cell stage",
       F, 
       1500,
       T,
       0,
       3, 
       0,
       cols = c("grey",
                "blue")
)

lab_DP(urchin_5_sp1_cell_type_labelled, 
       polysomally_recruited_transcripts_df,
       plot_title = "SMT expression profiles of cell types at the 8 cell stage",
       F, 
       1500,
       T,
       0,
       3, 
       0,
       cols = c("grey",
                "blue")
)
#FINAL RES: 
