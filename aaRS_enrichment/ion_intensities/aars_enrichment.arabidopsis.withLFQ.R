library(tidyverse)
library(ggforce)

log_floor = 0.01

joined_df = read.csv("aars_enrichment.arabidopsis.withLFQ.noPF.csv")

pivot_longer(joined_df,cols=c(Mitochondrial,Cytosolic,Chloroplast),names_to = "Subcellular",values_to = "Presence") %>%
  mutate(Subcellular = factor(Subcellular, levels = c("Cytosolic", "Chloroplast", "Mitochondrial"))) %>%
  ggplot() +
  geom_arc_bar(aes(
    x0=log10(log_floor + (arabidopsis_chloroplast_1_LFQ + arabidopsis_chloroplast_2_LFQ)/(arabidopsis_nuclear_1_LFQ + arabidopsis_nuclear_2_LFQ)), 
    y0=log10(log_floor + (arabidopsis_mito_1_LFQ + arabidopsis_mito_2_LFQ)/(arabidopsis_nuclear_1_LFQ + arabidopsis_nuclear_2_LFQ)), 
    r0=0,r=0.04,amount=Presence,fill=Subcellular),
    stat="pie", alpha=0.7, linetype="blank") +
  coord_equal() +
  #comment out the following line to exclude aaRS text labels
  geom_text(aes(
    label=Type, 
    x=log10(log_floor + (arabidopsis_chloroplast_1_LFQ + arabidopsis_chloroplast_2_LFQ)/(arabidopsis_nuclear_1_LFQ + arabidopsis_nuclear_2_LFQ)), 
    y=log10(log_floor + (arabidopsis_mito_1_LFQ + arabidopsis_mito_2_LFQ)/(arabidopsis_nuclear_1_LFQ + arabidopsis_nuclear_2_LFQ))), 
    nudge_y = 0.03, size=2) +
  scale_fill_manual(values = c("Mitochondrial"="goldenrod3","Cytosolic"="dodgerblue3","Chloroplast"="chartreuse4")) +
  xlab("log10 (Chloroplast / Total Leaf Ion Intensity Ratio)") +
  ylab("log10 (Mitochondrial / Total Leaf Ion Intensity Ratio)") +
  theme_bw() +
  theme(legend.position = "right", 
        axis.title = element_text(size=7, face="bold"), 
        axis.text = element_text(size=6), 
        legend.title = element_blank(), 
        legend.text = element_text(size=6),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
  ) 

ggsave("aars_enrichment.arabidopsis.withLFQ.noPF.pdf", width=6, height=3.5)


