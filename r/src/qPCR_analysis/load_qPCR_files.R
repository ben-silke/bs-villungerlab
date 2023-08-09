# Ok it is better to make these reusable, 
# but i think for now its faster for me to just create files where i need them
library(readxl)
library(plotly)
library(htmlwidgets)
library(dplyr)
library(dplyr)

# Import necessary util functions and set the directory correctly.
setwd("~/bs-villungerlab/")
source("~/bs-villungerlab/r/src/pca_utils.R")
source("~/bs-villungerlab/r/src/utils.R")
source("~/bs-villungerlab/r/src/star_utils.R")
source("~/bs-villungerlab/r/src/gene_selection_utils.R")

source("~/bs-villungerlab/r/src/qPCR_analysis/qPCR_utils.R")


# /Users/bsilke/bs-villungerlab/lab_work/qPCR/NOC_qPCR/20230808-p53_time_series-Noc_results_validation.xlsx
EXTENSION = "/Users/bsilke/bs-villungerlab/lab_work/qPCR/NOC_qPCR/"
FILE_NAME = "20230808-p53_time_series-Noc_results_validation.xlsx"

df = read_excel(file.path(EXTENSION, FILE_NAME), sheet='raw_data')

# View(df)
df$Content = NULL
# df$Sample <- gsub("+","p", df$Sample)
df$Sample <- gsub(" 0h","_0h", df$Sample)
df$Sample <- gsub(" 16h","_16h", df$Sample)
df$Sample <- gsub(" 20h","_20h", df$Sample)
df$Sample <- gsub(" 24h","_24h", df$Sample)
df$Sample <- gsub(" 36h","_36h", df$Sample)
df$Sample <- gsub(" 48h","_48h", df$Sample)

df$Cq <- as.numeric(as.character(df$Cq))

control_df <- subset(df, Target=="GAPDH")
# View(control_df)
colnames(control_df)
sample_new_df <- subset(control_df, grepl('new', Sample))
plus_new_df <- subset(control_df, grepl('+', Sample))
p7_new_df <- subset(control_df, grepl('p7', Sample))
# View(control_df)

df
groupd_df <- df %>%
  group_by(Target)
groupd_df

avg_control_df <- control_df %>%
  group_by(Sample) %>%
  summarise(avg_Cq = mean(Cq, na.rm = TRUE))

averaged_df <- df %>%
  group_by(Sample) %>%
  summarise(
    avg_Cq = mean(Cq, na.rm = TRUE),
  )

target_name="NINJ1"


# df_plot <- plot_ly(
#   long_df,
#   x = ~Sample,
#   y = ~avg_ddct2,
#   mode = "markers+lines",
#   hoverinfo = "text",
#   # text = ~paste("T: ", Timepoint, "<br>l2fc: ", log2foldchange, "<br>padj: ", padj, '<br>ID: ', symbol),
#   # split= ~symbol
# ) %>%
#   layout(title="ZM Treatment: Primers")
# df_plot
targets = c('BMF', "FOXM1", "SQSTM1", "NINJ1", "ZMAT3", "PHLDA3", "CCNA2", "CDCA8", "CDC25A", "AURKB", "ARID5B", "ARID5B", "ANKRD1")

full_df <- data.frame()
for (target in targets) {
  long_df <- create_long_df_for_sample(target)
  full_df <- merge(full_df, long_df, all = TRUE)
}

BMF_df <- create_graph_for_sample("new", 'BMF')
SQSTM1_df <- create_graph_for_sample("new", 'SQSTM1_df')
SQSTM1_df
p <- ggplot(long_df, aes(x = timepoint, y = avg_ddct2)) +
  geom_point() +
  geom_line() +
  labs(title = plot_title,
       x = "time(hrs)",
       y = expression(paste(log[2](x), ' fold change'))) +
  theme(plot.title = element_text(hjust = 0.5), # Center the title
        plot.title.position = "plot",
        legend.position = "bottom")
p

