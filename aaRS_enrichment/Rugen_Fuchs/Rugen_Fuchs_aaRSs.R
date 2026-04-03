library(tidyverse)
library(ggforce)

aars_df = read.csv("aars_enrichment.arabidopsis.Rugen_Fuchs.csv")


pivot_longer(aars_df,cols=c(Mitochondrial,Cytosolic,Chloroplast),names_to = "Subcellular",values_to = "Presence") %>%
  mutate(Subcellular = factor(Subcellular, levels = c("Cytosolic", "Chloroplast", "Mitochondrial"))) %>%
  ggplot() +
  geom_arc_bar(aes(
    y0=log10(Fuchs2020 + 1), 
    x0=log10(Rugen2024 + 1e-7), 
    r0=0,r=0.04,amount=Presence,fill=Subcellular),
    stat="pie", alpha=0.7, linetype="blank") +
  coord_equal() +
  #comment out the following line to exclude aaRS text labels
  geom_text(aes(
    label=Type, 
    y=log10(Fuchs2020 + 1), 
    x=log10(Rugen2024 + 1e-7)), 
    nudge_y = 0.03, size=2) +
  scale_fill_manual(values = c("Mitochondrial"="goldenrod3","Cytosolic"="dodgerblue3","Chloroplast"="chartreuse4")) +
  ylab("Copies per Mitochondrion (Fuchs et al. 2020)") +
  xlab("Proportion of Mitochondrial Proteome (Rugen et al. 2024)") +
  theme_bw() +
  theme(legend.position = "right", 
        axis.title = element_text(size=8, face="bold"), 
        axis.text = element_text(size=7), 
        legend.title = element_blank(), 
        legend.text = element_text(size=8),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
  ) +
  xlim(-7,-3) +
  ylim(0,3)

ggsave("Rugen_Fuchs_aaRS.pdf", width=5, height=6)

