library (tidyverse)

arabidopsis <- read.csv("../quantification/arabidopsis_abundance.csv")
silene <- read.csv("../quantification/silene_abundance.csv")

process_dataframe <- function(df, floor_percentile = 0.01) {
  if (!is.data.frame(df)) {
    stop("Input must be a dataframe.")
  }

  df_name <- deparse(substitute(df))
  cat ("Dataset:", df_name, "\n")
  
  df$Accession <- ifelse(
    startsWith(df$Accession, "ATCG"), paste0("plastid_", df$Accession),
    ifelse(startsWith(df$Accession, "ATMG"), paste0("mito_", df$Accession), df$Accession)
  )
  
  samples <- c("cp1", "cp2", "mt1", "mt2", "lf1", "lf2")
  median_floor_df <- data.frame(
    median = rep(NA, length(samples)),  # Repeat NA for the length of samples
    floor = rep(NA, length(samples)),   # Repeat NA for the length of samples
    row.names = samples
  )
  
  for (val in samples) {
    abundance_string <- paste0("Abundance_", val)
    found_string <- paste0("Found_", val)
    
    filtered_df <- df %>%
      filter(!is.na(.data[[abundance_string]]) & .data[[found_string]] == "High")
    median_value <- median(filtered_df[[abundance_string]], na.rm = TRUE)
    first_percentile <- quantile(filtered_df[[abundance_string]], probs = floor_percentile, na.rm = TRUE)
    median_floor_df[val, "median"] <- median_value
    median_floor_df[val, "floor"] <- first_percentile
  }
  
  for (val in samples) {
    abundance_string <- paste0("Abundance_", val)
    found_string <- paste0("Found_", val)
    df[[abundance_string]] <- ifelse(
      is.na(df[[abundance_string]]) | df[[abundance_string]] < median_floor_df[val, "floor"],
      median_floor_df[val, "floor"], 
      df[[abundance_string]]
    )
    
    #comment out the following line to include "Peak Found" abundance values in the calculations
    df[[abundance_string]] <- ifelse(df[[found_string]] != "High", median_floor_df[val, "floor"],  df[[abundance_string]])
    
    df[[abundance_string]] <- df[[abundance_string]] / median_floor_df[val, "median"]
  }  
  
  df <- df %>%
    filter(startsWith(Accession, "mito") | startsWith(Accession, "plastid")) %>%
    mutate(Genome = ifelse(startsWith(Accession, "mito"), "Mitochondrial", "Chloroplast"))
  
  print (median_floor_df)
  return(df)
}

arabidopsis_processed <- process_dataframe(arabidopsis)
silene_processed <- process_dataframe(silene)

arabidopsis_processed$species <- "Arabidopsis"
silene_processed$species <- "Silene"

combined_df <- rbind(arabidopsis_processed, silene_processed) 

combined_df$species <- factor (combined_df$species, levels = c("Arabidopsis", "Silene"))

combined_df %>%
  ggplot(aes(x = log10((Abundance_cp1	+ Abundance_cp2) / (Abundance_lf1	+ Abundance_lf2)), 
             y = log10((Abundance_mt1	+ Abundance_mt2) / (Abundance_lf1	+ Abundance_lf2)), 
             color=Genome, shape=Genome,
             size = Abundance_cp1 + Abundance_cp2 + Abundance_mt1 + Abundance_mt2 + Abundance_lf1 + Abundance_lf2)) +
  geom_point(alpha=0.7) +
  facet_wrap(~species) + 
  xlim(c(-3.2,3.2)) +
  ylim(c(-3.2,3.2)) +
  xlab("log10 (Chloroplast / Total Leaf Abundance Ratio)") +
  ylab("log10 (Mitochondrial / Total Leaf Abundance Ratio)") +
  scale_color_manual(values=c("chartreuse4", "goldenrod3")) +
  theme_bw() +
  theme(legend.position = "top", 
        axis.title = element_text(size=7, face="bold"), 
        axis.text = element_text(size=6), 
        strip.text = element_text(size=7, face="bold.italic"), 
        legend.title = element_text(size=7, face="bold"), 
        legend.text = element_text(size=6),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
  ) +
  guides(size = "none")

ggsave ("organelle_genome_enrichment.abundance.pdf", width = 4.5, height = 3.5)

