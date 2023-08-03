###

## This is probably not worth pursing. These values do not behave like this so this won't work.
## Any threshold should be done with a benchmark NES score, such as 3.0=very high

####

# probably better to use up and down because this provides us with a more normal distribution
decrease_file <- "/Users/bsilke/bs-villungerlab/results/output_encode/ZM/iregulon_analysis/Transcription_factors_ZM_3n_decrease.tsv"
increase_file <- "/Users/bsilke/bs-villungerlab/results/output_encode/ZM/iregulon_analysis/Transcription_factors_ZM_3n_increase.tsv"

df_inc <- read.table(increase_file, skip=13, header=FALSE, sep="\t", fill=TRUE)
df_inc$inc = TRUE

df_inc
df_dec <- read.table(decrease_file, skip=13, header=FALSE, sep="\t", fill=TRUE)
df_dec$inc = FALSE



headers = c('Motif', 'id',	'AUC',	'NES',	'ClusterCode',	'Transcription factor',	'Target genes', "increase?")
colnames(df_dec) <- headers
colnames(df_inc) <- headers

df_inc$ClusterCode <- paste0(df_inc$ClusterCode,'_inc')
df_dec$ClusterCode <- paste0(df_dec$ClusterCode,'_dec')

df_dec$NES <- -df_dec$NES

df <- rbind(df_inc, df_dec)

# View(df)

ggplot(df, aes(x=NES)) +
  geom_histogram(binwidth=0.5, fill="blue", color="black") +
  labs(title="Histogram of Values", x="Values", y="Frequency")
mean_nes <- mean(df$NES)
mean_nes
# ahh idk

sd_nes <- sd(df$NES)
sd_nes
df$zscore <- (df$NES - mean_nes) / sd_nes

# This is predicated on the assumption of a normal distribution which may likely not be true,
# However, it could be beneficial to provide a brief insight into splitting these groups
df$pnorm <- pnorm(abs(df$zscore))
df$prob <- 1-df$pnorm

View(df)
transcription_factors <- subset(df, grepl('T', (ClusterCode)))
# View(transcription_factors)

motifs <- subset(df, grepl('M', (ClusterCode)))

# View(motifs)

