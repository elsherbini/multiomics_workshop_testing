# Precompute the summaries the Module 2 slides plot.
#
# Run by hand when the data or this script changes, from the website root:
#   RENV_CONFIG_AUTOLOADER_ENABLED=FALSE \
#     Rscript materials/2-workshop2/2-count-table-anatomy/_build_results.R
#
# Writes results/*.csv next to the slides. The slides read those with readr
# only, so rendering the site never needs Bioconductor and renv.lock stays
# as it is. Because the slides do not compute anything, there is no cache to
# invalidate: edit a csv, re-render, and the new numbers are on the slide.

suppressPackageStartupMessages({library(tidyverse); library(vegan)})
here <- "materials/2-workshop2/2-count-table-anatomy"
out <- file.path(here, "results"); dir.create(out, showWarnings = FALSE, recursive = TRUE)
source("materials/2-workshop2/R/count_tables.R")

d <- read_counts("instructional")
md <- d$metadata

# --- per-sample summaries, the same five numbers on both tables ---------------
summarise_table <- function(counts, label) {
  lib <- colSums(counts)
  prop <- sweep(counts, 2, lib, "/")
  tibble(
    sample_id        = colnames(counts),
    table            = label,
    library_size     = lib,
    features_total   = nrow(counts),
    features_detected = colSums(counts > 0),
    zero_fraction    = colMeans(counts == 0),
    shannon          = vegan::diversity(t(prop), index = "shannon"),
    top_feature_share = apply(prop, 2, max)
  )
}

sample_qc <- bind_rows(
  summarise_table(d$species, "species"),
  summarise_table(d$genes, "genes")
) |> left_join(md |> select(sample_id, time_point, outcome, rnaseq_batch, nugent), by = "sample_id")
write_csv(sample_qc, file.path(out, "sample_qc.csv"))

# --- per-feature summaries ----------------------------------------------------
feature_summary <- function(counts, label, info, id_col, name_col) {
  prop <- sweep(counts, 2, colSums(counts), "/")
  tibble(
    feature    = rownames(counts),
    table      = label,
    prevalence = rowMeans(counts > 0),
    mean_share = rowMeans(prop),
    max_share  = apply(prop, 1, max)
  ) |>
    left_join(info |> select(feature = all_of(id_col), name = all_of(name_col)), by = "feature")
}

features <- bind_rows(
  feature_summary(d$species, "species", d$taxonomy, "taxon", "Genus"),
  feature_summary(d$genes, "genes", d$gene_info, "ensembl_gene_id", "symbol")
)
write_csv(features, file.path(out, "feature_summary.csv.gz"))

# --- a corner of each table, for the "what does it look like" slides ----------
corner <- function(counts, n_f = 6, n_s = 5) {
  as_tibble(counts[seq_len(n_f), seq_len(n_s)], rownames = "feature")
}
write_csv(corner(d$species), file.path(out, "corner_species.csv"))
write_csv(corner(d$genes), file.path(out, "corner_genes.csv"))

# --- filtering: how many features survive each rule --------------------------
filter_curve <- function(counts, label) {
  prev <- rowMeans(counts > 0)
  tibble(table = label,
         min_prevalence = seq(0, 0.9, by = 0.05)) |>
    rowwise() |>
    mutate(features_kept = sum(prev >= min_prevalence),
           share_of_counts_kept = sum(counts[prev >= min_prevalence, ]) / sum(counts)) |>
    ungroup()
}
write_csv(bind_rows(filter_curve(d$species, "species"), filter_curve(d$genes, "genes")),
          file.path(out, "filter_curve.csv"))

cat("module 2 results written to", out, "\n")
cat("samples:", nrow(md), " species:", nrow(d$species), " genes:", nrow(d$genes), "\n")
