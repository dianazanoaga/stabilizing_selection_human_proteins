# ============================================================
# Pie chart of the 48 phenotypes 
# one slice per trait, colored by trait category
# ============================================================

library(dplyr)
library(tibble)
library(ggplot2)

rm(list = ls())

trait_categories <- tribble(
  ~trait,               ~category,
  "albumin",            "Liver Function",
  "ALP",                "Liver Function",
  "ALT",                "Liver Function",
  "AST",                "Liver Function",
  "bilirubin",          "Liver Function",
  "GGT",                "Liver Function",
  "total_protein",      "Liver Function",
  "ApoA",               "Lipid Metabolism",
  "ApoB",               "Lipid Metabolism",
  "cholesterol",        "Lipid Metabolism",
  "HDL",                "Lipid Metabolism",
  "LDL",                "Lipid Metabolism",
  "TG",                 "Lipid Metabolism",
  "LPA",                "Lipid Metabolism",
  "glucose",            "Glycaemic Control",
  "HbA1c",              "Glycaemic Control",
  "Cr",                 "Renal Function",
  "cystatinC",          "Renal Function",
  "urea",               "Renal Function",
  "urate",              "Renal Function",
  "systolic_BP",        "Cardiovascular Vitals",
  "diastolic_BP",       "Cardiovascular Vitals",
  "heart_rate",         "Cardiovascular Vitals",
  "BMI",                "Anthropometry",
  "BMR",                "Anthropometry",
  "body_fat_mass",      "Anthropometry",
  "height",             "Anthropometry",
  "weight",             "Anthropometry",
  "waist_circumference","Anthropometry",
  "hip_circumference",  "Anthropometry",
  "Ca",                 "Mineral Metabolism",
  "phosphate",          "Mineral Metabolism",
  "vitamin_D",          "Mineral Metabolism",
  "IGF1",               "Endocrine Function",
  "SHBG",               "Endocrine Function",
  "CRP",                "Inflammation",
  "RBC_count",          "Haematology",
  "WBC_count",          "Haematology",
  "platelet_count",     "Haematology",
  "lymphocyte_count",   "Haematology",
  "monocyte_count",     "Haematology",
  "neutrophil_count",   "Haematology",
  "reticulocyte_count", "Haematology",
  "MCH",                "Haematology",
  "FVC",                "Respiratory Function",
  "grip_strength_right","Muscle Strength",
  "fluid_intelligence", "Cognitive Function",
  "neuroticism",        "Cognitive Function"
)

## clockwise order of categories around the pie + colors
cat_colors <- c(
  "Anthropometry"         = "#E0352B",
  "Cardiovascular Vitals" = "#F08A0C",
  "Cognitive Function"    = "#BE8A0A",
  "Endocrine Function"    = "#86A10A",
  "Glycaemic Control"     = "#159A0A",
  "Haematology"           = "#1DB33A",
  "Lipid Metabolism"      = "#139A6A",
  "Liver Function"        = "#16A09A",
  "Mineral Metabolism"    = "#1588CC",
  "Inflammation"          = "#0F5BAA",
  "Respiratory Function"  = "#BA0BCF",
  "Renal Function"        = "#DA0A8C",
  "Muscle Strength"       = "#F5076E"
)

pie_df <- trait_categories %>%
  mutate(category = factor(category, levels = names(cat_colors)),
         label = gsub("_", " ", trait),
         label = recode(label, "cystatinC" = "cystatin C")) %>%
  arrange(category, label) %>%
  mutate(pos = row_number(),
         theta = 360 * (pos - 0.5) / n(),          # clockwise from 12 o'clock
         left = theta > 180,
         angle = ifelse(left, 270 - theta, 90 - theta),
         hjust = ifelse(left, 0, 1),
         pos = factor(pos))

p <- ggplot(pie_df, aes(x = pos, y = 1, fill = category)) +
  geom_col(width = 1, color = "white", linewidth = 0.4) +
  geom_text(aes(y = 0.96, label = label, angle = angle, hjust = hjust),
            size = 3.4, fontface = "bold") +
  coord_polar(theta = "x", start = 0, clip = "off") +
  scale_y_continuous(limits = c(-0.03, 1), expand = c(0, 0)) +
  scale_fill_manual(values = cat_colors, breaks = sort(names(cat_colors))) +
  guides(fill = guide_legend(ncol = 2)) +
  labs(fill = "Trait category") +
  theme_void() +
  theme(legend.title = element_text(face = "bold", size = 14),
        legend.text = element_text(size = 11),
        legend.key.size = unit(0.6, "cm"),
        plot.background = element_rect(fill = "white", color = NA))

p

ggsave("~/selection_proteins/scripts_figure_may_26/figure_may_26/trait_categories_pie.pdf",
       p, width = 12, height = 7)
ggsave("~/selection_proteins/scripts_figure_may_26/figure_may_26/trait_categories_pie.png",
       p, width = 12, height = 7, dpi = 300)
