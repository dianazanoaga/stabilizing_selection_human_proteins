library(readxl)
library(ggplot2)
library(patchwork)
library(ggpubr)
library(ggplot2)
library(patchwork)
library(ggpubr)

# load ldpred2 prot no cis 

ldpred2_prot_no_cis <- read_excel(
  "~/selection_proteins/data_last/process/ldpred2-auto_proteins_coefficients_no_cis.xlsx"
)

ldpred2_prot_all <-  read_excel("~/selection_proteins/data_last/process/ldpred2-auto_proteins_coefficients.xlsx")

ld_a <- merge(ldpred2_prot_all, ldpred2_prot_no_cis, by.y = "pheno", by.x = "Gene")
ld_a <- ld_a[!duplicated(ld_a$Gene),] #2915

cor.test(ld_a$alpha_est_1.x, ld_a$alpha_est_1.y, method = "spearman")

cor.test(ld_a$h2_est_1.x, ld_a$h2_est_1.y, method = "spearman")

cor.test(ld_a$p_est_1.x, ld_a$p_est_1.y, method = "spearman")

theme_manuscript <- function() {
  theme_classic(base_size = 14) +
    theme(
      axis.title = element_text(face = "bold", size = 16),
      axis.text = element_text(face = "bold", size = 12, colour = "black"),
      plot.title = element_text(face = "bold", hjust = 0.5, size = 20),
      strip.text = element_text(face = "bold"),
      plot.tag = element_text(face = "bold", size = 18),
      axis.line = element_line(linewidth = 0.8),
      axis.ticks = element_line(linewidth = 0.8)
    )
}

my_col <- "#2C7FB8"

p1 <- ggplot(ld_a, aes(alpha_est_1.x, alpha_est_1.y)) +
  geom_point(
    colour = my_col,
    alpha = 0.35,
    size = 1.5
  ) +
  geom_smooth(
    method = "lm",
    colour = "#08306B",
    fill = "#9ECAE1",
    linewidth = 1
  ) +
  stat_cor(
    method = "spearman",
    aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")),
    label.x.npc = 0.35,
    label.y.npc = 0.05,
    size = 6
  ) +
  labs(
    x = expression(bold(alpha~"(cis+trans variants)")),
    y = expression(bold(alpha~"(only trans variants)"))
  ) +
  scale_x_log10() +
  scale_y_log10() +
  theme_manuscript() +
  labs(tag = "A")

p2 <- ggplot(ld_a, aes(h2_est_1.x, h2_est_1.y)) +
  geom_point(
    colour = my_col,
    alpha = 0.35,
    size = 1.5
  ) +
  geom_smooth(
    method = "lm",
    colour = "#08306B",
    fill = "#9ECAE1",
    linewidth = 1
  ) +
  stat_cor(
    method = "spearman",
    aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")),
    label.x.npc = 0.4,
    label.y.npc = 0.05,
    size = 6
  ) +
  labs(
    x = expression(bold(h^2~"(cis+trans variants)")),
    y = expression(bold(h^2~"(only trans variants)"))
  ) +
  theme_manuscript() +
  labs(tag = "B")

p3 <- ggplot(ld_a, aes(p_est_1.x, p_est_1.y)) +
  geom_point(
    colour = my_col,
    alpha = 0.35,
    size = 1.5
  ) +
  geom_smooth(
    method = "lm",
    colour = "#08306B",
    fill = "#9ECAE1",
    linewidth = 1
  ) +
  stat_cor(
    method = "spearman",
    aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")),
    label.x.npc = 0.35,
    label.y.npc = 0.05,
    size = 6
  ) +
  labs(
    x = expression(bold(p~"(cis+trans variants)")),
    y = expression(bold(p~"(only trans variants)"))
  ) +
  scale_x_log10() +
  scale_y_log10() +
  theme_manuscript() +
  labs(tag = "C")

final_plot <-
  p1 + p2 + p3 +
  plot_layout(ncol = 3) +
  plot_annotation(
    title = "LDpred2-auto parameter estimates with and without cis variants",
    theme = theme(
      plot.title = element_text(
        face = "bold",
        size = 22,
        hjust = 0.5,
        margin = margin(b = 15)
      )
    )
  )
final_plot

