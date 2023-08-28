library(readxl)
library(plotly)
library(htmlwidgets)
library(dplyr)
library("fgsea")

# Import necessary util functions and set the directory correctly.
setwd("~/bs-villungerlab/")
source("~/bs-villungerlab/r/src/pca_utils.R")
source("~/bs-villungerlab/r/src/utils.R")
source("~/bs-villungerlab/r/src/star_utils.R")
source("~/bs-villungerlab/r/src/gene_selection_utils.R")

source("~/bs-villungerlab/r/src/qPCR_utils.R")

# dir.create("lab_work/qPCR/validation")
EXTENSION = "lab_work/qPCR"

noc_df <- read_excel(file.path(EXTENSION, "validation/full_Noc_validation_updated.xlsx"), sheet='full_Noc_validation')
zm_df <- read_excel(file.path(EXTENSION, "validation/full_Noc_validation_updated.xlsx"), sheet='full_Noc_validation')
# View(noc_df)


long_df
df_noc <- noc_df %>%
  group_by(treatment) %>%
  summarise( by=cols_to_keep
  )
# library(dplyr)

filtered_df <- noc_df %>%
  dplyr::select('treatment', 'treatment_avg', 'treatment_sd', 'sample_target')

# View(filtered_df)

distinct_rows <- filtered_df %>%
  distinct()
# View(distinct_rows)

# distinct_rows$cell_type = strsplit(distinct_rows$treatment, split = "_")[[1]]
distinct_rows$cell_type <- sapply(strsplit(distinct_rows$treatment, split = "_"), `[`, 1)
# distinct_rows$ctrl <- sapply(strsplit(distinct_rows$treatment, split = "_"), `[`, 1)

distinct_rows$target_treatment = paste0(distinct_rows$cell_type,'_',distinct_rows$sample_target)

unique(distinct_rows$target_treatment)

sample_dmso <- subset(distinct_rows, grepl('DMSO', treatment))
names = colnames(sample_dmso)
names[2] = 'dmso_treatment_average'
names[3] = 'dmso_treatment_sd'
colnames(sample_dmso) = names
sample_dmso

sample_noc <- subset(distinct_rows, grepl('Noc', treatment))
sample_noc
names = colnames(sample_noc)
names[2] = 'noc_treatment_average'
names[3] = 'noc_treatment_sd'
colnames(sample_noc) = names
sample_noc


df = data.frame(treatments=unique(distinct_rows$target_treatment))
df
merged_df <- df %>%
  left_join(sample_dmso, by = c("treatments" = "target_treatment"))
merged_df$treatment=NULL
merged_df$sample_target=NULL
merged_df$cell_type=NULL

merged_df
merged_df <- merged_df %>%
  left_join(sample_noc, by = c("treatments" = "target_treatment"))


merged_df$fold_change = merged_df$noc_treatment_average / merged_df$dmso_treatment_average
merged_df$fc_sd = (merged_df$noc_treatment_average / merged_df$dmso_treatment_average) * sqrt ((merged_df$noc_treatment_sd/merged_df$noc_treatment_average)^2 + (merged_df$dmso_treatment_sd/merged_df$dmso_treatment_average))
r
view(merged_df)


graphing_df = merged_df
colnames(graphing_df)
graphing_df <- graphing_df %>% dplyr::select('treatments', 'fold_change', 'fc_sd', "cell_type", "sample_target")
graphing_df
colnames(graphing_df)


p <- ggplot(graphing_df, aes(x = sample_target, y = fold_change, fill = cell_type)) + 
  geom_bar(stat="identity", position=position_dodge()) + 
  geom_errorbar(aes(ymin = fold_change - fc_sd, ymax = fold_change + fc_sd), 
                width=.2,                    # Width of the error bars
                position=position_dodge(.9)) + 
  labs(title="Nocodazole Fold Change", 
       y="Fold Change", 
       x="Sample Target") + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # This line slants the x-axis labels


p
print(p)

