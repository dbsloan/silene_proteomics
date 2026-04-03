library(tidyverse)

arabidopsis = read.csv("aars_enrichment.arabidopsis.csv")
silene = read.csv("aars_enrichment.silene.csv")

#treat ArgRS 3-way enzyme as organellar 
arabidopsis <- arabidopsis %>%
  mutate(Targeting = if_else(Targeting == "Cyto+Organellar" & Type == "ArgRS", "Organellar", Targeting))

arabidopsis$Targeting <- ifelse(startsWith(arabidopsis$Targeting, "Cyto"), "Cyto", "Organellar")
silene$Targeting <- ifelse(startsWith(silene$Targeting, "Cyto"), "Cyto", "Organellar")


arabidopsis_sum <- arabidopsis %>%
  group_by(Type, Targeting) %>%
  summarize(
    arabidopsis_mito_1 = sum(arabidopsis_mito_1, na.rm = TRUE),
    arabidopsis_mito_2 = sum(arabidopsis_mito_2, na.rm = TRUE),
    arabidopsis_nuclear_1 = sum(arabidopsis_nuclear_1, na.rm = TRUE),
    arabidopsis_nuclear_2 = sum(arabidopsis_nuclear_2, na.rm = TRUE),
    .groups = "drop"
  )

silene_sum <- silene %>%
  group_by(Type, Targeting) %>%
  summarize(
    conica_mito_1 = sum(conica_mito_1, na.rm = TRUE),
    conica_mito_2 = sum(conica_mito_2, na.rm = TRUE),
    conica_nuclear_1 = sum(conica_nuclear_1, na.rm = TRUE),
    conica_nuclear_2 = sum(conica_nuclear_2, na.rm = TRUE),
    .groups = "drop"
  )

arabidopsis_sum <- arabidopsis_sum %>%
  mutate(mito_sum = arabidopsis_mito_1 + arabidopsis_mito_2)
arabidopsis_sum <- arabidopsis_sum %>%
  mutate(leaf_sum = pmax(arabidopsis_nuclear_1 + arabidopsis_nuclear_2, 1))

silene_sum <- silene_sum %>%
  mutate(mito_sum = conica_mito_1 + conica_mito_2)
silene_sum <- silene_sum %>%
  mutate(leaf_sum = pmax(conica_nuclear_1 + conica_nuclear_2, 1))

combined_df <- bind_rows(
  arabidopsis_sum %>% mutate(species = "Arabidopsis"),
  silene_sum %>% mutate(species = "Silene")
)

combined_df <- combined_df %>% select(-arabidopsis_mito_1, -arabidopsis_mito_2, -arabidopsis_nuclear_1, -arabidopsis_nuclear_2, -conica_mito_1, -conica_mito_2, -conica_nuclear_1, -conica_nuclear_2)

combined_df <- combined_df %>%
  mutate(mito_enrich = mito_sum / leaf_sum)

write.csv(combined_df, "differential_enrichment/combined_df_psm_counts.csv", row.names = FALSE)

combined_df <- combined_df %>% select(-mito_sum, -leaf_sum)


wide_df <- combined_df %>%
  pivot_wider(
    names_from = Targeting,
    values_from = mito_enrich,
    values_fill = 0
  )

wide_df <- wide_df %>%
  mutate(cyto_weight = replace(
    Cyto / (Cyto + Organellar),
    is.nan(Cyto / (Cyto + Organellar)),
    0
  ))

custom_order <- c(
  "AlaRS", "LeuRS", "ThrRS", "ValRS", 
  "ArgRS", "IleRS","AsnRS", "AspRS",
  "CysRS", "GluRS", "HisRS", "PheRS",
  "SerRS", "GlnRS", "GlyRS", "LysRS",
  "MetRS", "ProRS", "TrpRS", "TyrRS"
)

wide_df <- wide_df %>% mutate(Type = factor(Type, levels = rev(custom_order)))

ggplot(wide_df, aes(x = species, y = Type, fill = cyto_weight)) +
  geom_tile(color = "white") +
  scale_fill_gradient(
    low = "goldenrod3",
    high = "dodgerblue3",
    name = "Enrichment Bias"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 7),
    axis.text.y = element_text(size = 7),
    axis.title = element_blank(),
    panel.grid = element_blank()
  )

ggsave("aaRS_differential_enrichment.pdf", width=3.25, height=4)
