library(tidyverse)
library(GGally)

#df <- read.csv("LFQ_conica.proteins_cleaned.noPF.csv")
df <- read.csv("LFQ_arabidopsis.proteins_cleaned.noPF.csv")

if (dim(df)[2]==7){
  df <- df[,2:7]
}
df <- df %>% mutate(across(everything(), ~log10(.x + 1e6)))

colnames(df) <- c("Cp 1", "Cp 2", "Mito 1", "Mito 2", "Leaf 1", "Leaf 2")

axis_min = 6
axis_max = 13

cor_panel <- function(data, mapping, ...){
  
  x <- eval_data_col(data, mapping$x)
  y <- eval_data_col(data, mapping$y)
  
  pearson  <- cor(x, y, method = "pearson",  use="complete.obs")
  spearman <- cor(x, y, method = "spearman", use="complete.obs")
  
  label <- paste0(
    "Pearson: ",  sprintf("%.3f", pearson),
    "\nSpearman: ", sprintf("%.3f", spearman)
  )
  
  ggplot(data = data, mapping = mapping) +
    annotate("text", x = mean(x,na.rm=TRUE), y = mean(y,na.rm=TRUE), label = label, size=2) +
    theme_void()
}

my_points <- function(data, mapping, ...) {
  
  xvar <- as_label(mapping$x)
  yvar <- as_label(mapping$y)
  
  color <- "gray70"
  
  if ((xvar == "Cp 1" & yvar == "Cp 2")) {
    color <- "chartreuse4"
  }else if ((xvar == "Mito 1" & yvar == "Mito 2")) {
    color <- "goldenrod3"
  }else if ((xvar == "Leaf 1" & yvar == "Leaf 2")) {
    color <- "black"
  }
  
  
  ggplot(data = data, mapping = mapping) +
    geom_point(color = color, alpha = 0.15, size=0.5) +
    theme_bw()+
    xlim(axis_min,axis_max) +
    ylim(axis_min,axis_max)
}

my_density <- function(data, mapping, ...){
  
  x <- eval_data_col(data, mapping$x)
  
  ggplot(data = data, mapping = mapping) +
    geom_density(fill="gray80") +
    xlim(axis_min,axis_max) +
    ylim(0, 0.8) +
    theme_bw()
}

ggpairs(
  df,
  lower = list(continuous = my_points),
  upper = list(continuous = cor_panel),
  diag  = list(continuous = my_density)
) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  xlab ("log10 Abundance") +
  ylab ("log10 Abundance")

ggsave("Arabidopsis.CorrMatrix.LFQ.pdf", width=6.7, height=6.7)


