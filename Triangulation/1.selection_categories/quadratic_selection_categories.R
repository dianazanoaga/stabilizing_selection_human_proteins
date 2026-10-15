# Load libraries
library(data.table)
library(ggplot2)

# =========================
# Load OBS protein-level data only
# =========================

prot_lev_nchild <- fread("~/selection_proteins/data_last/process/p_nchild_edu.tsv")
prot_lev_long <- fread("~/selection_proteins/data_last/process/p_longevity.tsv")
prot_lev_compscore <- fread("~/selection_proteins/data_last/process/p_health_score.tsv")

# =========================
# Rename relevant columns
# =========================

colnames(prot_lev_nchild)[11] <- "alpha_nchild"
colnames(prot_lev_long)[11] <- "alpha_long"
colnames(prot_lev_compscore)[11] <- "alpha_health_score"

colnames(prot_lev_nchild)[14] <- "pval_nchild"
colnames(prot_lev_long)[14] <- "pval_long"
colnames(prot_lev_compscore)[14] <- "pval_health_score"

colnames(prot_lev_nchild)[20] <- "dist"
colnames(prot_lev_long)[20] <- "dist"
colnames(prot_lev_compscore)[20] <- "dist"


colnames(prot_lev_nchild)[10] <-  "PValues_PRS"
colnames(prot_lev_long)[10] <-  "PValues_PRS"
colnames(prot_lev_compscore)[10] <-  "PValues_PRS"


colnames(prot_lev_nchild)[15] <-  "pval_nchild_adj"
colnames(prot_lev_long)[15] <-  "pval_long_adj"
colnames(prot_lev_compscore)[15] <-  "pval_health_score_adj"

colnames(prot_lev_nchild)[16] <-  "PValues_PRS_adj"
colnames(prot_lev_long)[16] <-  "PValues_PRS_adj"
colnames(prot_lev_compscore)[16] <-  "PValues_PRS_adj"

# =========================
# Helper functions
# =========================

getcol <- function(dt, col) {
  if (col %in% names(dt)) dt[[col]] else rep(NA, nrow(dt))
}

classify_effect_obs <- function(dt,
                                alpha_col,
                                pval_col,
                                prs_pval_col = "PValues_PRS",
                                dist_col = "dist",
                                pthr = 0.05) {
  
  dt <- as.data.table(dt)
  
  alpha <- getcol(dt, alpha_col)
  pval  <- getcol(dt, pval_col)
  dist  <- getcol(dt, dist_col)
  prs_p <- getcol(dt, prs_pval_col)
  
  category <- rep("Not significant", nrow(dt))
  
  # =========================
  # Directional linear only
  # linear significant
  # quadratic non-significant
  # =========================
  
  category[
    !is.na(prs_p) &
      prs_p < pthr &
      !is.na(pval) &
      pval >= pthr
  ] <- "Directional - linear only"
  
  # =========================
  # Significant quadratic patterns
  # =========================
  
  sig_quad <- !is.na(pval) & pval < pthr & !is.na(alpha)
  
  category[sig_quad & dist == 1 & alpha < 0] <-
    "Stabilizing"
  
  category[sig_quad & dist == 1 & alpha > 0] <-
    "Disruptive"
  
  category[sig_quad & dist == 0 & alpha < 0] <-
    "Directional - non linear (stabilizing)"
  
  category[sig_quad & dist == 0 & alpha > 0] <-
    "Directional - non linear (disruptive)"
  
  category
}
make_counts_obs <- function(dt,
                            trait,
                            alpha_col,
                            pval_col,
                            prs_pval_col = "PValues_PRS",
                            dist_col = "dist",
                            pthr = 0.05) {
  
  cat <- classify_effect_obs(
    dt,
    alpha_col = alpha_col,
    pval_col = pval_col,
    prs_pval_col = prs_pval_col,
    dist_col = dist_col,
    pthr = pthr
  )
  
  out <- data.table(
    trait = trait,
    category = cat
  )
  
  out <- out[, .N, by = .(trait, category)]
  out
}
# =========================
# Build OBS counts only
# =========================
counts_list <- list(
  
  make_counts_obs(
    prot_lev_compscore,
    trait = "Health Score",
    alpha_col = "alpha_health_score",
    pval_col = "pval_health_score",
    prs_pval_col = "PValues_PRS"
  ),
  
  make_counts_obs(
    prot_lev_long,
    trait = "Longevity",
    alpha_col = "alpha_long",
    pval_col = "pval_long",
    prs_pval_col = "PValues_PRS"
  ),
  
  make_counts_obs(
    prot_lev_nchild,
    trait = "N.children",
    alpha_col = "alpha_nchild",
    pval_col = "pval_nchild",
    prs_pval_col = "PValues_PRS"
  )
)

counts_dt <- rbindlist(counts_list, fill = TRUE)

# =========================
# Category order
# Not significant last = grey on top
# =========================
counts_dt[, category := factor(
  category,
  levels = c(
    "Not significant",
    "Directional - linear only", 
    "Directional - non linear (disruptive)",
    "Disruptive",
    "Directional - non linear (stabilizing)",
    "Stabilizing"
  )
)]

counts_dt[, trait := factor(
  trait,
  levels = c("Health Score", "Longevity", "N.children")
)]

# =========================
# Colors
# using your previous green/red palette
# =========================

selection_cols <- c(
  "Stabilizing" = "aquamarine4",
  "Directional - non linear (stabilizing)" = "darkseagreen3",
  "Disruptive" = "tomato3",
  "Directional - non linear (disruptive)" = "salmon1",
  "Directional - linear only" = "khaki2",
  "Not significant" = "gray95"
)


# =========================
# Cleaner manuscript colors
# =========================

# shorter legend labels
legend_labs <- c(
  "Not significant" = "Not significant",
  "Directional - linear only" = "Directional linear only",
  "Directional - non linear (disruptive)" = "Directional non-linear (disruptive)",
  "Disruptive" = "Disruptive",
  "Directional - non linear (stabilizing)" = "Directional non-linear (stabilizing)",
  "Stabilizing" = "Stabilizing"
)

category_levels <- levels(counts_dt$category)
p_obs <- ggplot(counts_dt, aes(x = trait, y = N, fill = category)) +
  geom_col(
    width = 0.6,
    color = "white",
    linewidth = 0.35
  )  +
  scale_fill_manual(
    values = selection_cols,
    breaks = category_levels,
    labels = legend_labs,
    drop = FALSE
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.05))
  ) +
  labs(
    x = NULL,
    y = "Number of proteins",
    fill = "Selection pattern",
    title = "Selection patterns inferred from regression analyses of protein abundance on fitness proxies",
    subtitle = "Proteins were classified using a nominal significance threshold of p < 0.05."
  ) +
  theme_classic(base_size = 12, base_family = "Arial") +
  theme(
    axis.line = element_line(linewidth = 0.7, color = "black"),
    axis.ticks = element_line(linewidth = 0.7, color = "black"),
    plot.title = element_text(face = "bold", size = 18, hjust = 0),
    plot.subtitle = element_text(size = 16, hjust = 0),
    axis.text.x = element_text(
      size = 18,
      color = "black"
    ),
    axis.text.y = element_text(
      size = 16,
      color = "black"
    ),
    axis.title.y = element_text(
      size = 18,
      margin = margin(r = 8)
    ),
    
    legend.position = "right",
    legend.title = element_text(
      size = 16,
      face = "bold"
    ),
    legend.text = element_text(size = 14.5),
    legend.key.size = unit(0.5, "cm"),
    legend.spacing.y = unit(0.2, "cm"),
    
    plot.margin = margin(8, 8, 8, 8)
  )

p_obs


## adj 

counts_list <- list(
  
  make_counts_obs(
    prot_lev_compscore,
    trait = "Health Score",
    alpha_col = "alpha_health_score",
    pval_col = "pval_health_score_adj",
    prs_pval_col = "PValues_PRS_adj"
  ),
  
  make_counts_obs(
    prot_lev_long,
    trait = "Longevity",
    alpha_col = "alpha_long",
    pval_col = "pval_long_adj",
    prs_pval_col = "PValues_PRS_adj"
  ),
  
  make_counts_obs(
    prot_lev_nchild,
    trait = "N.children",
    alpha_col = "alpha_nchild",
    pval_col = "pval_nchild_adj",
    prs_pval_col = "PValues_PRS_adj"
  )
)

counts_dt <- rbindlist(counts_list, fill = TRUE)

# =========================
# Category order
# Not significant last = grey on top
# =========================
counts_dt[, category := factor(
  category,
  levels = c(
    "Not significant",
    "Directional - linear only", 
    "Directional - non linear (disruptive)",
    "Disruptive",
    "Directional - non linear (stabilizing)",
    "Stabilizing"
  )
)]

counts_dt[, trait := factor(
  trait,
  levels = c("Health Score", "Longevity", "N.children")
)]

# =========================
# Colors
# using your previous green/red palette
# =========================

selection_cols <- c(
  "Stabilizing" = "aquamarine4",
  "Directional - non linear (stabilizing)" = "darkseagreen3",
  "Disruptive" = "tomato3",
  "Directional - non linear (disruptive)" = "salmon1",
  "Directional - linear only" = "khaki2",
  "Not significant" = "gray95"
)


# =========================
# Cleaner manuscript colors
# =========================

# shorter legend labels
legend_labs <- c(
  "Not significant" = "Not significant",
  "Directional - linear only" = "Directional linear only",
  "Directional - non linear (disruptive)" = "Directional non-linear (disruptive)",
  "Disruptive" = "Disruptive",
  "Directional - non linear (stabilizing)" = "Directional non-linear (stabilizing)",
  "Stabilizing" = "Stabilizing"
)

p_obs <- ggplot(counts_dt, aes(x = trait, y = N, fill = category)) +
  geom_col(
    width = 0.6,
    color = "white",
    linewidth = 0.35
  )  +
  scale_fill_manual(
    values = selection_cols,
    breaks = category_levels,
    labels = legend_labs,
    drop = FALSE
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.05))
  ) +
  labs(
    x = NULL,
    y = "Number of proteins",
    fill = "Selection pattern",
    title = "Selection patterns inferred from regression analyses of protein abundance on fitness proxies",
    subtitle = "Proteins were classified using a nominal significance threshold of FDR p< 0.05."
  ) +
  theme_classic(base_size = 12, base_family = "Arial") +
  theme(
    axis.line = element_line(linewidth = 0.7, color = "black"),
    axis.ticks = element_line(linewidth = 0.7, color = "black"),
    plot.title = element_text(face = "bold", size = 13.5, hjust = 0),
    axis.text.x = element_text(
      size = 11.5,
      color = "black"
    ),
    axis.text.y = element_text(
      size = 10,
      color = "black"
    ),
    axis.title.y = element_text(
      size = 12.5,
      margin = margin(r = 8)
    ),
    
    legend.position = "right",
    legend.title = element_text(
      size = 13,
      face = "bold"
    ),
    legend.text = element_text(size = 10),
    legend.key.size = unit(0.5, "cm"),
    legend.spacing.y = unit(0.15, "cm"),
    
    plot.margin = margin(8, 8, 8, 8)
  )

p_obs

