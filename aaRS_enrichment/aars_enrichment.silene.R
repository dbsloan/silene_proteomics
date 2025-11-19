library(tidyverse)
library(ggforce)

log_floor = 0.1

aars_df = read.csv("silene_aaRS_list.csv")
psm_df = read.csv("../quantification/silene_psms.csv")

joined_df <- aars_df %>%
  left_join(psm_df, by = c("Accession" = "Master.Protein.Accessions"))

ggplot(data=joined_df, aes(
    x=log10(log_floor + (conica_chloroplast_1 + conica_chloroplast_2)/(conica_nuclear_1 + conica_nuclear_2)), 
    y=log10(log_floor + (conica_mito_1 + conica_mito_2)/(conica_nuclear_1 + conica_nuclear_2)), 
    color=Targeting)) +
  geom_point(alpha=0.5) +
  geom_text(aes(
    label=Type, 
    x=log10(log_floor + (conica_chloroplast_1 + conica_chloroplast_2)/(conica_nuclear_1 + conica_nuclear_2)), 
    y=log10(log_floor + (conica_mito_1 + conica_mito_2)/(conica_nuclear_1 + conica_nuclear_2))), 
    nudge_y = 0.05, size=2) +
  xlab("log10 (Chloroplast / Total Leaf PSM Ratio)") +
  ylab("log10 (Mitochondrial / Total Leaf PSM Ratio)") +
  theme_bw()


joined_df <- joined_df %>%
  mutate(
    Mitochondrial = ifelse(grepl("M", Target), 1, 0),
    Cytosolic = ifelse(grepl("N", Target), 1, 0),
    Chloroplast = ifelse(grepl("P", Target), 1, 0)
  )

pivot_longer(joined_df,cols=c(Mitochondrial,Cytosolic,Chloroplast),names_to = "Subcellular",values_to = "Presence") %>%
  mutate(Subcellular = factor(Subcellular, levels = c("Cytosolic", "Chloroplast", "Mitochondrial"))) %>%
  ggplot() +
  geom_arc_bar(aes(
    x0=log10(log_floor + (conica_chloroplast_1 + conica_chloroplast_2)/pmax(conica_nuclear_1 + conica_nuclear_2, 1)), 
    y0=log10(log_floor + (conica_mito_1 + conica_mito_2)/pmax(conica_nuclear_1 + conica_nuclear_2, 1)), 
    r0=0,r=0.025,amount=Presence,fill=Subcellular),
    stat="pie", alpha=0.7, linetype="blank") +
  coord_equal() +
  #comment out the following line to exclude aaRS text labels
  geom_text(aes(
    label=Type, 
    x=log10(log_floor + (conica_chloroplast_1 + conica_chloroplast_2)/pmax(conica_nuclear_1 + conica_nuclear_2,1)), 
    y=log10(log_floor + (conica_mito_1 + conica_mito_2)/pmax(conica_nuclear_1 + conica_nuclear_2,1))), 
    nudge_y = 0.03, size=2) +
  scale_fill_manual(values = c("Mitochondrial"="goldenrod3","Cytosolic"="dodgerblue3","Chloroplast"="chartreuse4")) +
  xlab("log10 (Chloroplast / Total Leaf PSM Ratio)") +
  ylab("log10 (Mitochondrial / Total Leaf PSM Ratio)") +
  theme_bw() +
  theme(legend.position = "right", 
        axis.title = element_text(size=7, face="bold"), 
        axis.text = element_text(size=6), 
        legend.title = element_blank(), 
        legend.text = element_text(size=6),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
  ) 

ggsave("aars_enrichment.silene.pdf", width=6, height=3.5)

write.csv(joined_df, "aars_enrichment.silene.csv", row.names = FALSE)
