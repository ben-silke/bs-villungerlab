#specify input directories (directories always go on top of a file, that way they can easily be exchanged)
dataset_input_directory <- "C:/Users/vbrennsteiner/OneDrive - CeMM Research Center GmbH/vbrennsteiner_home/gw_filip/data/raw"
plot_output_directory <- "C:/Users/vbrennsteiner/OneDrive - CeMM Research Center GmbH/vbrennsteiner_home/gw_filip/data/processed"
plot_output_directory <- "C:/Users/vbrennsteiner/OneDrive - CeMM Research Center GmbH/vbrennsteiner_home/gw_filip/assets"
dataset_input_name <- "EC_PH_120_M1075-O572-P12664-A01-E04-PepNorm-SumInd-June2022_Proteins.txt"
dataset_output_directory <- "C:/Users/vbrennsteiner/OneDrive - CeMM Research Center GmbH/vbrennsteiner_home/gw_filip/results"

#load necessary packages
library(readr)
library(ggplot2)
library(magrittr)
library(dplyr)

#load data
df <- readr::read_tsv(paste0(dataset_input_directory,"/",dataset_input_name))

#filter data: 
df <- df %>%
  filter(
    `Protein FDR Confidence Combined` == "High",
    `Contaminant` == FALSE,
    `Number of Protein Unique Peptides` >= 1,
    `Number of Peptides` >= 2
    )

#inspect shape of data
dim(df) #EC_PH_120 export has 3526 rows and 221 columns

#extract column names that we need: there should be two comparisons in the data
ratio_colnames <- colnames(df)[grepl("Abundance Ratio log2",colnames(df))]

#schematic of volcanoplotting, for the real plots we use a bigger function with more parameters
for (ratio in ratio_colnames) {
  
  #compute name for pvalue column that belongs to the current ratio
  pvalue <- sub("log2","P-Value",ratio)
  
  #subset dataset: we only need the current columns for plotting
  df_subset <- df %>%
    dplyr::select(Accession, `Gene Symbol`, all_of(ratio), all_of(pvalue))
  
  #rename columns for easier plotting, since dealing with variable names in ggplot is always finicky 
  colnames(df_subset)[3] <- "log2_foldChange"
  colnames(df_subset)[4] <- "p-value"
  
  #before plotting, assign hit status to datapoints
  #case_when is similar to if-statements. The mutate()-function is used to make new columns, the below statement translates to 
  #make a new column called hit_status, which is "upregulated" if log2fc/pval are such and "downregulated" if they are thus, etc.
  #the TRUE ~ statement specifies the default value, it must always be given.
  df_subset %<>%
    mutate(hit_status = case_when(
      log2_foldChange > 1.5 & `p-value` < 0.05 ~ "upregulated",
      log2_foldChange < -1.5 & `p-value` < 0.05 ~ "downregulated",
      TRUE ~ "nonregulated"
    )) %>%
    mutate(label = case_when(
      hit_status == "upregulated" | hit_status == "downregulated" ~ `Gene Symbol`,
      TRUE ~ NA
    ))
  
  #plot volcanoplot
  plt <- ggplot(data = df_subset) + 
    #white background
    theme_bw() +
    #hide legend
    theme(
      legend.position = "none"
    ) +
    #scatterplot: point layer
    geom_point(aes(x = log2_foldChange, y = -log10(`p-value`), color = hit_status)) + 
    #labels; since this layer uses the label column, only hits are labelled
    geom_text(aes(x = log2_foldChange, y = -log10(`p-value`), label = label)) +
    #background lines
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") + 
    geom_vline(xintercept = 1.5, linetype = "dashed") + 
    geom_vline(xintercept = -1.5, linetype = "dashed") + 
    #see geom_point: there, we assign hit_status to color, here we specify what each hit status should look like
    scale_color_manual(
      values = c("nonregulated" = "grey",
                 "upregulated" = "#E95E30",
                 "downregulated" = "#0F5569") 
    ) +
    #slightly adapt the ratio name and use as title
    labs(title = sub("  ","/",ratio)) + 
    #specify limits !!! This can cut off points, use with care !!!
    xlim(c(-3,3)) + 
    ylim(c(0,10))
  
  #show plot
  print(plt)
}

#Actual Plots for FG use pre-written plot function

#source function for volcanoplots
source("C:/Users/vbrennsteiner/OneDrive - CeMM Research Center GmbH/vbrennsteiner_home/bd_r_scripts_curated/msDataVisualizations/scripts/volcanoplotLabelled.R")

#source function for paper-like figure layout
source("C:/Users/vbrennsteiner/OneDrive - CeMM Research Center GmbH/vbrennsteiner_home/bd_r_scripts_curated/figureLayout/scripts/paperLayout.R")

#generate list of plots with selected settings
EC_PH_volcanos <- screeningVolcano(
  dataset_directory = dataset_input_directory,
  dataset_matchstring = dataset_input_name,
  plot_all_given_column = 'all',
  plot_output_directory = plot_output_directory,
  dataset_output_directory = dataset_output_directory,
  print_plot = T,
  return_plot = T,
  return_plot_dataset = F,
  fc_threshold = 1,
  fc_threshold_negative = -1,
  label_hits_sel_both = "both",
  label_hit_down_color = "#0F5569",
  label_hit_up_color = "#E95E30",
  hideDownregulatedLabels = FALSE,
  selected_accessions = c("Q92934", "Q9NQS1", "Q9NZ45", "P07996", "Q99805", "Q9Y4K0", "P10646", "Q9ULB4")
)

#save graphics
for (plot_name in names(EC_PH_volcanos$plots)) {
  plt <- EC_PH_volcanos$plots[plot_name]
  plotResize(
    outdir = plot_output_directory,
    plotname = sprintf("%s_volcanoplot_BindDegsFilter.png",plot_name),
    ggplot_object_list = plt,
    caption = sprintf("**EC_PH_120:** Nostatin Treatment"),
    pbu_width = 2,
    pbu_height = 2,
    showAxisTitle.x = TRUE,
    showAxisTitle.y = TRUE,
    showTitle = TRUE,
    showLegend = FALSE
  )
}

#GSEA functions
dataset_input_directory <- "C:/Users/vbrennsteiner/OneDrive - CeMM Research Center GmbH/vbrennsteiner_home/gw_filip/data/raw"
dataset_input_name <- "EC_PH_120_M1075-O572-P12664-A01-E04-PepNorm-SumInd-June2022_Proteins.txt"
plot_output_directory <- "C:/Users/vbrennsteiner/OneDrive - CeMM Research Center GmbH/vbrennsteiner_home/gw_filip/assets"
dataset_output_directory <- "C:/Users/vbrennsteiner/OneDrive - CeMM Research Center GmbH/vbrennsteiner_home/gw_filip/results"
uniprot_dataset_directory <- 

loadPDdataAndShowRatios <- function(
    dataset_input_directory,
    dataset_input_name
) {
  df <- readr::read_delim(paste0(dataset_input_directory,"/",dataset_input_name))
  colnames(df) <- gsub(" ","_",colnames(df))
  ratio_colnames <- colnames(df)[grepl("Ratio_log2|Ratio_P-Value", colnames(df))]
  print(ratio_colnames)
  
  return(df)
}

basicGSEA <- function(
    PD_protein_dataset,
    ratio_colname,
    FDR_filter = "High",
    Contaminant_filter = FALSE,
    Peptides_filter_gteq = 2,
    Unique_Peptides_filter_gteq = 1,
    ontology = "BP",
    keyType = "ENSEMBL",
    OrgDb = "org.Hs.eg.db",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH"
) {
  
  require(DESeq2)
  require(clusterProfiler)
  require(org.Hs.eg.db)
  require(readr)
  require(ggplot2)
  require(magrittr)
  require(dplyr)
  
  #filter dataset
  PD_protein_dataset %<>%
    filter(
      `Protein_FDR_Confidence_Combined` == FDR_filter,
      `Contaminant` == Contaminant_filter,
      `Number_of_Peptides` >= Peptides_filter_gteq,
      `Number_of_Protein_Unique_Peptides` >= Unique_Peptides_filter_gteq
    )
  
  #drop NA rows
  PD_protein_dataset <- PD_protein_dataset %>%
    dplyr::select(`Ensembl_Gene_ID`, GO_Accessions, Gene_Symbol, all_of(ratio_colname)) %>%
    dplyr::filter(if_all(everything(), function(x)!is.na(x)))
  
  #establish score column 
  PD_protein_dataset$score <- PD_protein_dataset[,ratio_colname] %>% unlist(use.names = FALSE)
  
  #rank decending by score
  PD_protein_dataset %<>%
    arrange(desc(score))
  
  #quick inspection plot to check distribution of score (should be 0-centered and not have too much heteroskedasticity)
  score_hist <- ggplot(data = PD_protein_dataset) +
    theme_bw() +
    geom_histogram(aes(x = score), bins = 100) +
    labs(title = "Distribution of 'Score'")
  print(score_hist)
  
  #another way to show scores
  score_dist <- PD_protein_dataset %>%
    mutate(index = row_number()) %>%
    ggplot() + 
    theme_bw() +
    geom_col(aes(x = index, y = score)) + 
    labs(title = "Ranked Scores") + 
    geom_hline(yintercept = 0) 
  print(score_dist)
    
  #run GSEA
  gene_list <- PD_protein_dataset$score
  names(gene_list) <- PD_protein_dataset$`Ensembl_Gene_ID`
  gse <- gseGO(
    geneList = gene_list,
    ont = ontology,
    keyType = keyType,
    OrgDb = OrgDb,
    pvalueCutoff = pvalueCutoff,
    pAdjustMethod = pAdjustMethod,
    eps = 3e-100
  )
  
  print(gse)
  return(gse)
}

#run analysis
dataset <- loadPDdataAndShowRatios(
  dataset_input_directory,
  dataset_input_name
)

gse <- basicGSEA(
  PD_protein_dataset = dataset,
  ratio_colname = "Abundance_Ratio_log2_NoA-500nM__Control"
)

#extract ensemble IDs and fold changes from gse object
ensemble_fc_df <- data.frame(
  "ensemble_id" = names(gse@geneList),
  "foldchange" = gse@geneList
)

#extract metadata dataframe from gsea object
gse_df <- gse@result

#filter for FDR < 0.01 
gse_df %<>%
  filter(qvalue < 0.01) %>%
  arrange(desc(NES))

#visualize different gsea distributions
available_terms <- gse_df$Description

#function to plot selected gse from available list
gsePlotter <- function(
    gseGoResultObject,
    available_terms,
    current_term
) {
  rdf <- gseGoResultObject@result
  current_term_id <- rdf$ID[rdf$Description == current_term]
  current_term_NES <- rdf$NES[rdf$Description == current_term]

  #Re-derive the number of genes from the original gene set that were found in the data
  expressed_genes <- gseGoResultObject@geneList
  current_set_genes <- gseGoResultObject@geneSets[[current_term_id]]
  found_genes <- expressed_genes[names(expressed_genes) %in% current_set_genes]
  found_genes_df <- data.frame(
    "ID" = names(found_genes),
    "FC" = found_genes
  )
  
  #annotate leading edge genes, which are in "core enrichment"
  leading_edge_genes <- rdf$core_enrichment[rdf$ID == current_term_id] %>% 
    str_split(.,"/") %>% 
    unlist() 
  
  found_genes_df$leading_edge <- FALSE
  found_genes_df$leading_edge[found_genes_df$ID %in% leading_edge_genes] <- TRUE
  
  #some messages to console
  message(sprintf(" --> showing term %s (%s)",current_term_id,current_term))
  message(sprintf(" --> current NES: %f",current_term_NES))
  
  #plot object
  gp <- NA
  if (current_term %in% available_terms) {
    print("term found")
    gp <- gseaplot(gseGoResultObject, geneSetID = current_term_id)
  }
  
  #compile list of plot object and ensembl IDs for output
  output_list <- list(
    "gsea_plot" = list(gp),
    "ensembl_IDs" = list(found_genes_df)
  )
  
  return(output_list)
}

#plot the FDR filtered_terms and their NES
gsePlotAndTableSaver <- function(
    proteomeDiscoverer_result_df,
    gseaGO_result_object,
    available_terms,
    dataset_output_directory,
    plot_output_directory,
    project_savename = "NoA_500nM_1percFDRfiltered"
) {
  #result summary table from gse. This contains information on all enriched terms but 
  #not the individual accessions that are enriched in a dataset. This table is also saved to a file
  gse_df <- gseaGO_result_object@result
  openxlsx::write.xlsx(gse_df,paste0(dataset_output_directory,"/",project_savename,"_GOterms_OverviewTable.xlsx"))
  
  #iterate over provided GO terms, generate GSEA plots and Accesion tables & save
  for (t in available_terms) {
    #gsePlotter returns a plot and a list of accessions.
    gsePlotterOutput <- gsePlotter(gseaGO_result_object,available_terms,current_term = t)
    
    #extract some metrics to show
    NES <- gse_df$NES[gse_df$Description == t] %>% round(digits = 4) 
    qval <- gse_df$qvalue[gse_df$Description == t] %>% signif(digits = 2) %>% format(format = "e") 

    #Extract plot from output
    gseaplot <- gsePlotterOutput[["gsea_plot"]]
    ranked_list_metric_plt <- gseaplot[[1]][[1]] +
      ylab("Ranked \nScore")
    running_enrichment_score_plt <- gseaplot[[1]][[2]] + 
      ylab("Enrichment \nScore")
    #adding plots together with "/" and "+" is part of the patchwork package,
    #which is very convenient to stack and combine different plots into one figure
    plt <- ranked_list_metric_plt / running_enrichment_score_plt +
      patchwork::plot_annotation(
        title = paste0("GO Term: ",t,"\nNES: ",NES,", qvalue: ",qval," (",gse_df$ID[gse_df$Description == t],")")
      )
    
    #save plot to current directory with descriptive name starting with formatted NES score
    ggsave(paste0(plot_output_directory,"/",gsub("\\.","",NES),"_",t,"_GSEA_plot.png"),plt,height = 4,width = 7)
    
    #compile output list of gene names with corresponding foldchange
    found_ids <- gsePlotterOutput[["ensembl_IDs"]][[1]]
    #match uniprot IDs and gene names to ensembl gene IDs
    pd_reference <- proteomeDiscoverer_result_df %>%
      dplyr::select(Accession, `Gene Symbol`, `Ensembl Gene ID`)
    found_ids <- left_join(found_ids,pd_reference,by = c("ID" = "Ensembl Gene ID"))
    openxlsx::write.xlsx(found_ids,paste0(dataset_output_directory,"/",t,project_savename,"_filtered_Accessions.xlsx"))
  }
}

#run functions
undebug(gsePlotAndTableSaver)
gsePlotAndTableSaver(
    proteomeDiscoverer_result_df = df,
    gseaGO_result_object = gse,
    available_terms = available_terms,
    dataset_output_directory = dataset_output_directory,
    plot_output_directory = plot_output_directory,
    project_savename = "NoA_500nM_1percFDRfiltered"
    )


