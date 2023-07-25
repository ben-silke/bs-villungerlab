from metaprogramming_r import SalmonRFileWriter, StarRFileWriter

import click

@click.command()
@click.option('--treatment', '-t', help='treatment (e.g. ZM)')
@click.option('--output', '-o', help='Where should the files be saved?')
@click.option('--nresults', '-n', help='how many genes do you want graphs for??')
@click.option('--main_time', '-mt', help='What time should serve as the reference to get the most expressed genes?')
# @click.option('--data_location', '-o', help="where should we save data??")

# python3 ../../../metap_script.py -s -o "results/output_encode_1to6" -md -a -r -f "data/encode_output_htseq_counts"

SHORT_TIMES = ["c(0, 8, 12, 16, 24, 48)", [8,12,16,24,48], "c(8, 12, 16, 24, 48)"]
LONG_TIMES = ["c(0, 16, 20, 24, 36, 48)", [16,20,24,36,48], "c(16, 20, 24, 36, 48)"]

TIME_DICT = {
    'ZM': LONG_TIMES,
    'Noc': LONG_TIMES,
    'DHCB': LONG_TIMES,
    'Nutl': SHORT_TIMES,
    'Etop': SHORT_TIMES
}


def create_files(treatment, output, main_time=0, n_results=20):
    pass


def write_output_file(treatment, output, n_results, main_time):
    content = f"""
library(DESeq2)
library(DESeq2)
library("apeglm")
library("ashr")
library(EnhancedVolcano)
library(org.Hs.eg.db)
library("pheatmap")
library("RColorBrewer")
library("PoiClaClu")
library("limma")
library("glue")
library("Rsubread")
library("stringr")
library("DESeq2")
library(dplyr)
library(openxlsx)
library(ggplot2)
library(tidyverse)
library(tidyverse)

setwd("~/bs-villungerlab/")
source("~/bs-villungerlab/r/src/pca_utils.R")
source("~/bs-villungerlab/r/src/utils.R")
source("~/bs-villungerlab/r/src/star_analysis/star_utils.R")
source("~/bs-villungerlab/r/src/star_analysis/gene_selection_utils.R")

load('~/bs-villungerlab/results/{output}/{treatment}_star_data.RData')
dds_{treatment} <- ddseq_{treatment}
"""

    for time in times:
        content = content + f"""
results_{treatment}_t{time}_df <- return_results(dds, "timepoint_t{time}_vs_t0", "_{time}")=
        """
    
    # Append the main
    times = set(times).remove(main_time)
    content = content + f"""
df <- results_{treatment}_t{main_time}_df
"""
    times_list = '''
    dfs <- list('''
    for time in times:
        times_list = times_list + f'results_{treatment}_t{time}_df, '
    times_list = times_list+')'

    content = content + times_list

    content = content + """

#now you can apply filtering
upr_df_sorted <- df[order(-df$log2FoldChange), ]
# set the number of results which you want
upr_top <- head(upr_df_sorted, {n_results})
upr_top_merged_df <- merge_all_data(upr_top, dfs, '{output}/generegulation/{treatment}_gene_regulation_data.csv')
upr_top_long_df <- make_longdf_for_plot(upr_top_merged_df, dfs)
upr_plot <- plot_longdf(upr_top_long_df, "{treatment} upregulated genes")
upr_plot
ggsave(filename = "{output}/generegulation/{treatment}_upregulated_genes.pdf", plot = upr_plot)

downr_df_sorted <- df[order(df$log2FoldChange), ]
# set the number of results which you want
downr_top <- head(upr_df_sorted, {n_results})
downr_top_merged_df <- merge_all_data(downr_top, dfs)
downr_top_long_df <- make_longdf_for_plot(downr_top_merged_df, dfs)
downr_plot <- plot_longdf(upr_top_long_df, "{treatment} downregulated genes")
downr_plot
ggsave(filename = "{output}/generegulation/{treatment}_downregulated_genes.pdf", plot = downr_plot)
"""



if __name__ == '__main__':
    create_files()


# Run to create star scripts:
# python3 metap_script.py -d "python/src" -a -r -s -md -f "data/"