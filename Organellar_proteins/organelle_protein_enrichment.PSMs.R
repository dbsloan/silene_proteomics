library (tidyverse)

arabidopsis <- read.csv("../quantification/arabidopsis_psms.csv")
silene <- read.csv("../quantification/silene_psms.csv")

min_total_reads = 5
min_log = 0.1

process_dataframe <- function(df) {
  if (!is.data.frame(df)) {
    stop("Input must be a dataframe.")
  }

  colnames(df)[1] <- "Accession"
  
  df$Accession <- ifelse(
    startsWith(df$Accession, "ATCG"), paste0("plastid_", df$Accession),
    ifelse(startsWith(df$Accession, "ATMG"), paste0("mito_", df$Accession), df$Accession)
  )
  
  #get rid of species name from headers
  colnames(df) <- sub("^[^_]*_", "", colnames(df))
  
  
  df <- df %>%
    filter(startsWith(Accession, "mito") | startsWith(Accession, "plastid")) %>%
    mutate(Genome = ifelse(startsWith(Accession, "mito"), "Mitochondrial", "Chloroplast"))
  
  return(df)
}

arabidopsis_processed <- process_dataframe(arabidopsis)
silene_processed <- process_dataframe(silene)

arabidopsis_processed$species <- "Arabidopsis"
silene_processed$species <- "Silene"

combined_df <- rbind(arabidopsis_processed, silene_processed) 

combined_df$species <- factor (combined_df$species, levels = c("Arabidopsis", "Silene"))

combined_df_filter_low <- combined_df %>% filter (nuclear_1 + nuclear_2 > 0 | chloroplast_1 + chloroplast_2 + mito_1 + mito_2 >= min_total_reads)

combined_df_filter_low %>%
  filter(species != "Agrostemma") %>%
  ggplot(aes(
    x = log10(min_log + (chloroplast_1 + chloroplast_2) / pmax(nuclear_1 + nuclear_2, 1)),
    y = log10(min_log + (mito_1 + mito_2) / pmax(nuclear_1 + nuclear_2, 1)),
    color = Genome,
    shape = Genome,
    size = chloroplast_1 + chloroplast_2 + mito_1 + mito_2 + nuclear_1 + nuclear_2
  )) +
  geom_point(alpha=0.7) +
  facet_wrap(~species) + 
  xlab("log10 (Chloroplast / Total Leaf PSM Ratio)") +
  ylab("log10 (Mitochondrial / Total Leaf PSM Ratio)") +
  scale_color_manual(values=c("chartreuse4", "goldenrod3")) +
  theme_bw() +
  theme(legend.position = "top", 
        legend.margin = margin(t = -5, b = -5),
        axis.title = element_text(size=7, face="bold"), 
        axis.text = element_text(size=6), 
        strip.text = element_text(size=7, face="bold.italic"), 
        legend.title = element_text(size=7, face="bold"), 
        legend.text = element_text(size=6),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
  ) +
  guides(size = "none")

ggsave ("organellar_protein_enrichment.PSMs.pdf", width = 3.25, height = 2.5)

