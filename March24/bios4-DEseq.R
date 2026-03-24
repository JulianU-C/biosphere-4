## Create an R project (top right corner) in an existing directory and select
## the biosphere_IV-main/March24/ folder you just downloaded and extracted
getwd
## Count Data: Differential expression==========================================
## =============================================================================
library(tximport)
library(readr)
library(DESeq2)
library(tidyverse)
library(EnhancedVolcano)

## load metadata
## -----------------------------------------------------------------------------
md <- read.csv("March24/meta.csv")

# set reference levels
md$tissue <- factor(md$tissue, levels = c("skin", "spleen"))
# Default is alphabetical order, so can just leave alone and skin will be ref.

## load read counts
## -----------------------------------------------------------------------------
countsTableRaw <- read.delim("March24/gene_counts.txt", 
                             comment.char = "#",
                             row.names = 1)

counts <- countsTableRaw %>% dplyr::select(!(Chr:Length))
countsCol <- colnames(counts)
rownames(md) <- countsCol

# check to make sure columns and rows match for metadata and count data
all(rownames(md) %in% colnames(counts))
all(rownames(md) == colnames(counts))
counts <- counts[, rownames(md)]
#write.csv(counts, "counts.csv", row.names = T)

## load emapper annotations
## -----------------------------------------------------------------------------
annoRaw <- read.delim("March24/emapper.emapper/emapper.emapper.annotations",
                      skip = 4,
                      header = TRUE)
# remove last 3 rows
tail(annoRaw)
anno <- annoRaw[-( (nrow(annoRaw) - 2) : nrow(annoRaw) ), ] 
anno$gene_id <- str_extract(anno[[1]], "^MSTRG\\.\\d+")
tail(annoRaw)

# remove duplicates by keeping best hit (highest bit-score)
annoBest <- anno %>%
  group_by(gene_id) %>%
  slice_max(order_by = score, with_ties = FALSE) 

rownames(annoBest) <- annoBest$gene_id

## DESeq 
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = md,
                              design = ~tissue)
# prefilter
smallestGroupsize <- 8
keep <- rowSums(counts(dds) >= 10) >= smallestGroupsize
dds <- dds[keep, ]

# run DESeq
dds <- DESeq(dds)
resultsNames(dds)

# full model
res <- results(dds, name = "tissue_spleen_vs_skin")
res <- res[order(res$padj),]
summary(res)

# extract DEGs
genes <- as.data.frame(res)
genes$gene_id <- rownames(genes)
degs <- genes %>% 
  filter(log2FoldChange > 0.5 | log2FoldChange < -0.5) %>% 
  filter(padj < 0.1)

# add annotations to all genes (not just DEGs)
genes <- as.data.frame(res)
genes$gene_id <- rownames(genes)
genes <- left_join(genes, annoBest, by = "gene_id")
genes <- genes %>% 
  mutate(
    Description = ifelse(is.na(Description) | 
                           Description == "-", "uncharacterized", Description),
    COG_category = ifelse(is.na(COG_category) |
                            COG_category == "-", "uncharacterized", COG_category),
    Preferred_name = ifelse(is.na(Preferred_name) |
                              Preferred_name == "-", "uncharacterized", Preferred_name)
  )

# I want to write out these data to a .csv
write.csv(genes, file = "March24/DEG_data.csv", row.names = F)

# I want to look for my one special gene
gene_of_interest <- genes %>% filter(gene_id == "MSTRG.5605")

## Visualize expression vs. reference
vol1 <- EnhancedVolcano(genes,
                        lab = "",
                        x = "log2FoldChange",
                        y = "padj",
                        pCutoff = 0.001,
                        FCcutoff = 1,
                        pointSize = 1.0,
                        labSize = 2.0,
                        title = "",
                        subtitle = "")

vol1

## Heatmap! _ can't get to run
library(pheatmap) # I like this package and too lazy to do in base/tidyverse
vsd <- vst(dds, blind = FALSE)
topDEGs <- rownames(res)[order(res$padj)][1:100]
dfMF <- as.data.frame(colData(dds)[,c("bd","itra")])
row_anno <- genes %>%
  filter(gene_id %in% topDEGs) %>%
  select(gene_id, COG_category) %>%
  column_to_rownames("gene_id")
(hm1 <- pheatmap(assay(vsd)[topDEGs, ],
                 scale = "row",
                 cluster_rows = TRUE,
                 cluster_cols = FALSE,
                 show_rownames = TRUE,
                 show_colnames = FALSE,
                 cellwidth = 10,
                 cellheight = 10,
                 annotation_row = row_anno,
                 annotation_col = dfMF,
                 annotation_names_col = F,
                 annotation_names_row = F,
                 treeheight_row = 0)
)
