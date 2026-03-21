## Create an R project (top right corner) in an existing directory and select
## the biosphere_IV-main folder you just downloaded and extracted

## 96-well plate data: Bacteria growth curves ==================================
## =============================================================================

# 1/2 plate setup:
  # col1 = Water
  # col2 = E. coli
  # col3 = J. lividum
  # col4 = Agrobacterium sp.
  # col5 = Enterococcus sp.
  # col6 = Media
  # B-D = No Toxin
  # E-G = Toxin
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
      toxin = case_when(row %in% c("E", "F", "G") ~ "Toxin", TRUE ~ "SHAM"),
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
