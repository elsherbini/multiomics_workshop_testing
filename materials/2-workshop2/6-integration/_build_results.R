# Precompute the Module 6 integration results.
#
# Run by hand from the website root, after modules 4 and 5:
#   RENV_CONFIG_AUTOLOADER_ENABLED=FALSE \
#     Rscript materials/2-workshop2/6-integration/_build_results.R
#
# One number per person from each table, then the three models the Friday
# discussion walks through. Nothing exotic: this is the honest baseline every
# fancier integration method has to beat.

suppressPackageStartupMessages({library(tidyverse); library(msigdbr); library(pROC)})
here <- "materials/2-workshop2/6-integration"
out <- file.path(here, "results"); dir.create(out, showWarnings = FALSE, recursive = TRUE)
source("materials/2-workshop2/R/count_tables.R")

d <- read_counts("instructional")
groups <- read_csv("materials/2-workshop2/4-differential/results/sample_groups.csv",
                   show_col_types = FALSE)

# Host score: mean z-scored log-CPM of two inflammatory Hallmark sets.
hallmark <- msigdbr(species = "Homo sapiens", collection = "H")
inflam <- hallmark |>
  filter(gs_name %in% c("HALLMARK_INFLAMMATORY_RESPONSE", "HALLMARK_TNFA_SIGNALING_VIA_NFKB")) |>
  pull(ensembl_gene) |> unique()

cpm <- sweep(d$genes, 2, colSums(d$genes), "/") * 1e6
expressed <- rownames(cpm)[rowMeans(cpm) >= 5]
lcpm <- log2(cpm[expressed, ] + 1)
set_genes <- intersect(inflam, rownames(lcpm))
cat("inflammatory genes used for the score:", length(set_genes), "of", length(inflam), "\n")
z <- t(scale(t(lcpm[set_genes, ])))
host_score <- colMeans(z, na.rm = TRUE)

per_sample <- groups |>
  select(sample_id, participant_id, time_point, outcome, community, lacto_fraction) |>
  mutate(anaerobe_fraction = 1 - lacto_fraction,
         host_score = host_score[sample_id])
write_csv(per_sample, file.path(out, "scores_per_sample.csv"))

per_person <- per_sample |>
  group_by(participant_id, outcome) |>
  summarise(micro_score = mean(anaerobe_fraction),
            host_score = mean(host_score), .groups = "drop") |>
  mutate(preterm = as.integer(outcome == "preterm"),
         micro_z = as.numeric(scale(micro_score)),
         host_z = as.numeric(scale(host_score)))
write_csv(per_person, file.path(out, "scores_per_person.csv"))

fits <- list(
  `microbiome only` = glm(preterm ~ micro_z, per_person, family = binomial),
  `host only`       = glm(preterm ~ host_z, per_person, family = binomial),
  `both`            = glm(preterm ~ micro_z + host_z, per_person, family = binomial)
)

model_tbl <- imap_dfr(fits, function(f, nm) {
  r <- roc(per_person$preterm, predict(f, type = "link"), quiet = TRUE)
  ci <- as.numeric(ci.auc(r))
  tibble(model = nm, auc = as.numeric(auc(r)), auc_lo = ci[1], auc_hi = ci[3],
         aic = AIC(f))
})
lrt <- anova(fits$`microbiome only`, fits$both, test = "Chisq")
model_tbl$lrt_p_vs_micro <- c(NA, NA, lrt$`Pr(>Chi)`[2])
write_csv(model_tbl, file.path(out, "model_comparison.csv"))
print(model_tbl)

coef_tbl <- imap_dfr(fits, function(f, nm) {
  broom_like <- summary(f)$coefficients
  as_tibble(broom_like, rownames = "term") |>
    mutate(model = nm, odds_ratio = exp(Estimate)) |>
    rename(estimate = Estimate, std_error = `Std. Error`, p_value = `Pr(>|z|)`)
})
write_csv(coef_tbl, file.path(out, "model_coefficients.csv"))

roc_tbl <- imap_dfr(fits, function(f, nm) {
  r <- roc(per_person$preterm, predict(f, type = "link"), quiet = TRUE)
  tibble(model = nm, fpr = 1 - r$specificities, tpr = r$sensitivities)
})
write_csv(roc_tbl, file.path(out, "roc_curves.csv"))

cat("correlation between the two scores:",
    sprintf("%.3f", cor(per_person$micro_z, per_person$host_z)), "\n")
cat("module 6 results written to", out, "\n")
