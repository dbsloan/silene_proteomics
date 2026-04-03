library(tidyverse)

mito_anticorrelation = read.csv("Cyto_Organellar_Anticorrelation.csv")

ggplot (data=mito_anticorrelation, aes(y=log10(FuchsCytosolic+1), x=log10(FuchsOrganellar+1)))+
  geom_smooth(method = "lm", fill = "gray80", color="gray50", alpha = 0.5) +
  geom_point(alpha=0.7) +
  geom_text(aes(
    label=Type, 
    nudge_y = 0.1)) +
  ylab('Cytosolic aaRS Copies per Mitochondrion (Fuchs et al. 2020)') +
  xlab('Chloroplast aaRS Copies per Mitochondrion (Fuchs et al. 2020)') +
  theme_bw() +
  theme(legend.position = "right", 
        axis.title = element_text(size=8, face="bold"), 
        axis.text = element_text(size=7), 
        legend.title = element_blank(), 
        legend.text = element_text(size=8),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
  )


cor.test(
  x = log10(mito_anticorrelation$FuchsCytosolic + 1),
  y = log10(mito_anticorrelation$FuchsOrganellar + 1)
)