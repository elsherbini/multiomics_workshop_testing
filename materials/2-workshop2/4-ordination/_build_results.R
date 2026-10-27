# Precompute the Module 4 ordination and PERMANOVA results.
#
# Run by hand from the website root:
#   RENV_CONFIG_AUTOLOADER_ENABLED=FALSE \
#     Rscript materials/2-workshop2/4-ordination/_build_results.R
#
# The coordinates are taken from ord_plot()$data, so the figures on the slides
# are exactly what the microViz pipeline shown on those slides produces.
#
# One pipeline, three knobs, both tables:
#   tax_transform(trans) |> dist_calc(dist) |> ord_calc("PCoA") |> ord_plot()
# Bray-Curtis on proportions, Jaccard on presence/absence, Euclidean on CLR
# (which is Aitchison distance). The gene table takes the third one.

suppressPackageStartupMessages({library(tidyverse); library(microViz); library(phyloseq)})
here <- "materials/2-workshop2/4-ordination"
out <- file.path(here, "results"); dir.create(out, showWarnings = FALSE, recursive = TRUE)
source("materials/2-workshop2/R/count_tables.R")

d <- read_counts("instructional")
groups <- read_csv("materials/2-workshop2/3-differential/results/sample_groups.csv",
                   show_col_types = FALSE)

add_groups <- function(ps) {
  sd <- as(sample_data(ps), "data.frame") |>
    rownames_to_column("sid") |>
    left_join(groups |> select(sample_id, community, lacto_fraction, bv), by = c("sid" = "sample_id")) |>
    column_to_rownames("sid")
  sample_data(ps) <- sample_data(sd)
  ps
}

ps_sp <- add_groups(as_counts(table = "microbiome"))
ps_gn <- add_groups(as_counts(table = "rnaseq", min_count = 10, min_samples = 100))
cat("species:", ntaxa(ps_sp), " genes kept:", ntaxa(ps_gn), "\n")

# --- ordinations, three ways on the species table and one on the genes -------
# ord_plot's data frame keeps the sample id in a column named after the
# phyloseq sample names; pull it out robustly instead of guessing.
ord_tbl <- function(ps, rank, trans, dist, method, label) {
  pl <- ps |>
    tax_transform(trans, rank = rank) |>
    dist_calc(dist) |>
    ord_calc(method = method) |>
    ord_plot(colour = "community", shape = "outcome")
  df <- pl$data
  ax <- names(df)[1:2]
  tibble(sample_id = rownames(df),
         axis1 = df[[ax[1]]], axis2 = df[[ax[2]]],
         axis1_name = ax[1], axis2_name = ax[2],
         ordination = label)
}

ords <- bind_rows(
  ord_tbl(ps_sp, "species", "identity", "bray",      "PCoA", "species / Bray-Curtis"),
  ord_tbl(ps_sp, "species", "binary",   "jaccard",   "PCoA", "species / Jaccard"),
  ord_tbl(ps_sp, "species", "clr",      "euclidean", "PCoA", "species / Aitchison"),
  ord_tbl(ps_gn, "symbol",  "clr",      "euclidean", "PCoA", "genes / Aitchison")
) |>
  left_join(groups |> select(sample_id, time_point, outcome, community, lacto_fraction, bv, nugent),
            by = "sample_id") |>
  left_join(d$metadata |> select(sample_id, rnaseq_batch), by = "sample_id")
write_csv(ords, file.path(out, "ordinations.csv"))

# --- PERMANOVA on both tables, same variables --------------------------------
perm <- function(ps, rank, trans, dist, label) {
  res <- ps |>
    tax_transform(trans, rank = rank) |>
    dist_calc(dist) |>
    dist_permanova(variables = c("community", "outcome", "time_point", "rnaseq_batch"),
                   n_perms = 199, seed = 1, by = "margin") |>
    perm_get()
  as_tibble(as.data.frame(res), rownames = "term") |>
    mutate(table = label) |>
    rename(p_value = `Pr(>F)`)
}
permanova <- bind_rows(
  perm(ps_sp, "species", "identity", "bray", "species / Bray-Curtis"),
  perm(ps_gn, "symbol", "clr", "euclidean", "genes / Aitchison")
)
write_csv(permanova, file.path(out, "permanova.csv"))
print(permanova |> select(table, term, R2, p_value))

# --- how much of each table's variance the first axes hold -------------------
eig <- function(ps, rank, trans, dist, method, label) {
  o <- ps |> tax_transform(trans, rank = rank) |> dist_calc(dist) |> ord_calc(method = method)
  ev <- microViz::ord_get(o)
  v <- if (!is.null(ev$values$Relative_eig)) ev$values$Relative_eig else {
    e <- if (!is.null(ev$CA$eig)) ev$CA$eig else ev$eig
    e / sum(e)
  }
  tibble(table = label, axis = seq_along(v), variance_share = as.numeric(v)) |> head(10)
}
write_csv(bind_rows(
  eig(ps_sp, "species", "identity", "bray", "PCoA", "species / Bray-Curtis"),
  eig(ps_gn, "symbol", "clr", "euclidean", "PCoA", "genes / Aitchison")
), file.path(out, "scree.csv"))

cat("module 4 results written to", out, "\n")
