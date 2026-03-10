## Create an R project (top right corner) in an existing directory and select
## the biosphere_IV-main folder you just downloaded and extracted

## 96-well plate data: Bacteria growth curves ==================================
## =============================================================================

# 1/2 plate setup:
  # col1 = Water
  # col2 = E. coli (MClover)
  # col3 = J. liv (dusy basin)
  # col4 = Agrobacterium UM1803
  # col5 = Enterococcus sp.
  # col6 = Media
  # B-D = No STX
  # E-G = STX
  # A and H = Empty wells

# this is the only library we'll need, don't get stuck in library hell
# tidyverse can do all the data manipulation you need
library(tidyverse)

# place all your files into one folder accessible from your working directory
# check your current directory using getwd()
# make life EASY, name your files like "DayX_yourfilename.txt"
files <- list.files(path = "./",
                    pattern = "*.csv",
                    full.names = TRUE)

## Here's a function for converting all my plates to long form
processPlate <- function(file_path) {
  # read the data, in each .txt file the OD readings in 96-well format appear at
  # line 26, and since 96-well plate is only 8 rows long we can set n_max to 8
  raw <- as.data.frame(read_csv(file_path))
  # then we need to cleanup the data
  data1 <- raw %>% select(2:7)
  row.names(data1) <- raw$...1
  data1 <- data1 %>% filter(!(row.names(data1) %in% c('A', "H")))
  # convert your clean plate to long format for analysis
  numbLong <- data1 %>%
    rownames_to_column("row") %>%
    pivot_longer(
      cols = -row,
      names_to = "column",
      values_to = "OD600") %>%
    # this 'mutate' function will create new columns
    mutate(
      column = as.integer(column),
      well = paste0(row, column),
      strain = case_when(column == 1 ~ "Water",
                         column == 2 ~ "E.coli",
                         column == 3 ~ "J.lividum",
                         column == 4 ~ "Agrobacterium UM1803",
                         column == 5 ~ "Enterococcus sp.",
                         column == 6 ~ "Media"),
      toxin = case_when(row %in% c("E", "F", "G") ~ "STX", TRUE ~ "SHAM"),
      file = basename(file_path),
      day = str_extract(file, "Day\\d+"),
      dayNum = str_extract_all(day, "\\d"),
      dose = str_extract(file, "_\\d+.\\d+")
    )
}

# merge everything together, how does it look? 
allData <- map_dfr(files, processPlate)
allData # todo bien?

# if there are any wells that you think are contaminated and need to be removed
# allData <- filter(well %in% c("this well", "this well", "etc."))

################################################################################

# normalize data (set media for each day as baseline)
allData_norm <- allData %>% 
  group_by(day) %>% 
  mutate(
    meanMedia = mean(OD600[strain == "Media"], na.rm = TRUE),
    corrOD600 = OD600 - meanMedia)
allData_norm

# calculate mean and standard deviation
summStats <- allData_norm %>% 
  group_by(strain, toxin, day, dayNum, dose) %>% 
  summarise(
    meanOD = mean(OD600, na.rm = TRUE),
    stdOD = sd(OD600, na.rm = TRUE),
    n = n(),
    seOD = stdOD / sqrt(n),
    .groups = "drop") %>% 
  filter(strain != "Media")
  
summStats # How does it look? Wells are gone! Don't need them anymore

# all I want to do is remove that "_" in front of the dose
summStats$dose <- gsub("_", "", summStats$dose)

# Ah the best part... figure generation, truly a time of self-expression
# figure generation is considered advanced in most coding languages
# good thing for us is R is easy and built for scientific figures
# select a theme for your plots, if unsure just use themebw() or themeclassic()
theme_set(theme_bw())

# Advanced line plot:
AdvLine <- ggplot(summStats, aes(x = as.numeric(dayNum),
                      y = meanOD,
                      color = strain, 
                      group = strain)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = meanOD - stdOD, ymax = meanOD + stdOD),
                width = 0.1, size = 0.75) +
  scale_x_continuous(breaks = c(0, 3, 6, 9)) +
  labs(x = "Incubation time (day)",
       y = "Optical density (600nm)",
       color = NULL, 
       caption = "Error bars: mean OD +/- SD") +
  facet_wrap(~toxin + dose)

AdvLine

# Save your plot as a .png
ggsave("AdvLine.png", plot = AdvLine, 
       width = 7,
       height = 5, 
       units = c("in"),
       dpi = 300)

## Count Data: Differential expression==========================================
## =============================================================================
library(tximport)
library(readr)
library(DESeq2)
library(tidyverse)
library(EnhancedVolcano)

## load metadata
## -----------------------------------------------------------------------------
# NOTE: ok I loaded my metadata as a tibble before (read_csv) and it wouldn't
# subset in dds?? Not sure why but make sure use read.csv for this
md <- read.csv("meta.csv")

# set "no" treatment and "neg" infection as reference levels
md$tissue <- factor(md$tissue, levels = c("skin", "spleen"))
# Default is alphabetical order, so can just leave alone and skin will be ref.

## load read counts
## -----------------------------------------------------------------------------
countsTableRaw <- read.delim("gene_counts.txt", 
                             comment.char = "#",
                             row.names = 1)

counts <- countsTableRaw %>% dplyr::select(!(Chr:Length))
countsCol <- colnames(counts)
rownames(md) <- countsCol

# check to make sure columns and rows match for metadata and count data
all(rownames(md) %in% colnames(counts))
all(rownames(md) == colnames(counts)) # if false run code below
counts <- counts[, rownames(md)]
#write.csv(counts, "counts.csv", row.names = T)

## emapper annotations
#annoRaw <- read.delim("emapper.emapper.annotations",
  #                    skip = 4,
   #                   header = TRUE)
# remove last 3 rows
#tail(annoRaw)
#anno <- annoRaw[-( (nrow(annoRaw) - 2) : nrow(annoRaw) ), ] 
#anno$gene_id <- str_extract(anno[[1]], "^MSTRG\\.\\d+")
#tail(annoRaw)

# remove duplicates by keeping best hit (highest bit-score)
#annoBest <- anno %>%
#  group_by(gene_id) %>%
 # slice_max(order_by = score, with_ties = FALSE) 

#rownames(annoBest) <- annoBest$gene_id

## DESeq on full dataset, then just look at MSTRG.5605
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

sxph <- genes %>% filter(gene_id == "MSTRG.5605")
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
