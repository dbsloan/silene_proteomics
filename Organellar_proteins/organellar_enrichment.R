library(tidyverse)

enrichment_LFQ = read.csv("organellar_enrichment_LFQ_total.csv")
enrichment_LFQ$Type = factor(enrichment_LFQ$Type, levels=c("Chloroplast", "Mitochondrial", "Leaf"))


ggplot (data = enrichment_LFQ, aes(x=Cp_Abundance_Scaled, y=Mt_Abundance_Scaled, color=Type, shape=Species))+
  geom_point(alpha=0.5, size=3)+
  scale_x_log10(limits=c(0.01,0.5)) +
  scale_y_log10(limits=c(5e-5,1e-1)) +
  theme_bw()+
  theme(
    axis.title = element_text(size=8, face="bold"), 
    axis.text = element_text(size=7), 
    legend.title = element_text(size=8, face="bold"), 
    legend.text = element_text(size=7),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )+
  scale_color_manual(values=c("chartreuse4", "goldenrod3", "black"))+
  xlab("Plastid-Encoded Peptide Abundance")+
  ylab("Mitochondrial-Encoded Peptide Abundance")
  
ggsave("organellar_enrichment_LFQ_total.pdf", width=4.25, height=3)


enrichment_PSM = read.csv("organellar_enrichment_PSM_total.csv")
enrichment_PSM$Type = factor(enrichment_LFQ$Type, levels=c("Chloroplast", "Mitochondrial", "Leaf"))


ggplot (data = enrichment_PSM, aes(x=Cp_Abundance_Scaled, y=Mt_Abundance_Scaled, color=Type, shape=Species))+
  geom_point(alpha=0.5, size=3)+
  scale_x_log10(limits=c(0.01,0.2)) +
  scale_y_log10(limits=c(1e-3,1e-1)) +
  theme_bw()+
  theme(
    axis.title = element_text(size=8, face="bold"), 
    axis.text = element_text(size=7), 
    legend.title = element_text(size=8, face="bold"), 
    legend.text = element_text(size=7),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )+
  scale_color_manual(values=c("chartreuse4", "goldenrod3", "black"))+
  xlab("Proportion Plastid-Encoded PSMs")+
  ylab("Proportion Mitochondrial-Encoded PSMs")

ggsave("organellar_enrichment_PSM_total.pdf", width=4.25, height=3)


enrichment_LFQ_by_gene = read.csv("organellar_enrichment_LFQ.csv")

mt_max = 2
mt_min = -2
cp_max = 2
cp_min = -2

enrichment_LFQ_by_gene %>%
  filter(UniquePeptides >= 5) %>%
#  ggplot (aes(x=log10(cp/leaf), y=log10(mito/leaf), color=Genome, shape=Species, size=mito+leaf+cp))+
  ggplot (aes(pmin(cp_max, pmax(cp_min,x=log10(cp_noPF/leaf_noPF))), y=pmin(mt_max,pmax(mt_min,log10(mito_noPF/leaf_noPF))), color=Genome, shape=Species, size=mito_noPF+leaf_noPF+cp_noPF))+
  geom_point(alpha=0.5)+
#  facet_wrap(~Species)+
  theme_bw()+
  theme(
    axis.title = element_text(size=8, face="bold"), 
    axis.text = element_text(size=7), 
    legend.title = element_text(size=8, face="bold"), 
    legend.text = element_text(size=7),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )+
  scale_color_manual(values=c("goldenrod3", "chartreuse4"))+
  xlab("Chloroplast Enrichment Ratio")+
  ylab("Mitochondrial Enrichment Ratio")

ggsave("organellar_enrichment_LFQ.pdf", width=5, height=3)


enrichment_PSM_by_gene = read.csv("organellar_enrichment_PSMs.csv")

mt_max = 1
mt_min = -1
cp_max = 1
cp_min = -1


enrichment_PSM_by_gene %>%
  filter(UniquePeptides >= 5) %>%
  ggplot (aes(pmin(cp_max, pmax(cp_min,x=log10(cp/leaf))), y=pmin(mt_max,pmax(mt_min,log10(mito/leaf))), color=Genome, shape=Species, size=mito+leaf+cp))+
  geom_point(alpha=0.5)+
#  facet_wrap(~Species)+
  theme_bw()+
  theme(
    axis.title = element_text(size=8, face="bold"), 
    axis.text = element_text(size=7), 
    legend.title = element_text(size=8, face="bold"), 
    legend.text = element_text(size=7),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )+
  scale_color_manual(values=c("goldenrod3", "chartreuse4"))+
  xlab("Chloroplast Enrichment Ratio")+
  ylab("Mitochondrial Enrichment Ratio")

ggsave("organellar_enrichment_PSM.pdf", width=4.25, height=3)

  