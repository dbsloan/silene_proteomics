library(tidyverse)

phers = read.csv ("PheRS.csv")

phers$Sample = factor (phers$Sample, levels = c("Total Leaf", "Isolated Chloroplasts", "Isolated Mitochondria"))
phers$Targeting = factor (phers$Targeting, levels = c("Cytosolic", "Chloroplast", "Mitochondrial"))

  
ggplot(data=phers, aes(x=Rep, y=PSMs, fill=Targeting)) +
  geom_bar(stat = "identity") +
  facet_wrap(~Sample) +
  theme_classic() +
  scale_fill_manual(values = c("Mitochondrial"="goldenrod3","Cytosolic"="dodgerblue3","Chloroplast"="chartreuse4")) +
  xlab ("Sample") +
  ylab ("# PSMs") +
  theme(legend.position = "right", 
        axis.title = element_text(size=7, face="bold"), 
        axis.text = element_text(size=6), 
        strip.text = element_text(size=6), 
        legend.title = element_blank(), 
        legend.text = element_text(size=6),
  ) 

ggsave("PheRS.pdf", width=4, height=2.75)