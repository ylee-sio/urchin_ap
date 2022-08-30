setwd("~/spurpuratus_GO/")
library(org.Mm.eg.db)
cluster_pigment <- subset(urchin_5_splg_clustered_unmodified, idents = '14')

expr_pigment <- as.matrix(GetAssayData(cluster_pigment))

# Select genes that are expressed > 0 in at least 75% of cells (somewhat arbitrary definition)
n.gt.pigment <- apply(expr_pigment, 1, function(x)length(which(x > 1)))

expressed.genes_pigment <- rownames(expr_pigment)[which(n.gt.pigment/ncol(expr_pigment) >= 0.75)]

all.genes_pigment <- rownames(expr_pigment)

# define geneList as 1 if gene is in expressed.genes, 0 otherwise
geneList_pigment <- ifelse(all.genes_pigment %in% expressed.genes_pigment, 1, 0)

names(geneList_pigment) <- all.genes_pigment
# Create topGOdata object
GOdata_pigment <- new("topGOdata",
              ontology = "MF", # use biological process ontology
              allGenes = geneList_pigment,
              geneSelectionFun = function(x)(x == 1),
              annot = annFUN.org, mapping = "org.Spurpuratus.eg.db", ID = "symbol")

resultFisher_pigment <- runTest(GOdata_pigment, algorithm = "elim", statistic = "fisher")
GenTable(GOdata_pigment, Fisher = resultFisher_pigment, topNodes = 20, numChar = 60) %>% arrange(desc(Annotated))

#WGCNA
options(stringsAsFactors = F)
datExpr <- t(as.matrix(GetAssayData(urchin_5_splg_clustered_unmodified)))[,VariableFeatures(urchin_5_splg_clustered_unmodified)]  # only use variable genes in analysis

net <- blockwiseModules(datExpr, power = 10,
                        corType = "bicor", # use robust correlation
                        networkType = "signed", minModuleSize = 10,
                        reassignThreshold = 0, mergeCutHeight = 0.15,
                        numericLabels = F, pamRespectsDendro = FALSE,
                        saveTOMs = TRUE,
                        saveTOMFileBase = "TOM",
                        verbose = 3)


tic()
makeOrgPackageFromNCBI(version = "4.2",
                       author = "YL <yoonlee@ucsd.edu>",
                       maintainer = "YL <yoonlee@ucsd.edu>",
                       outputDir = ".",
                       tax_id = "7668",
                       genus = "Strongylocentrotus",
                       species = "purpuratus",
                       verbose = T,
                       rebuildCache = F)

toc()


