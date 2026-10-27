# Precompute the Module 5 enrichment results.
#
# Run by hand from the website root, AFTER module 3's build script, whose
# ranked results this reads:
#   RENV_CONFIG_AUTOLOADER_ENABLED=FALSE \
#     Rscript materials/2-workshop2/5-enrichment/_build_results.R
#
# The same ranking-and-sets machinery twice: Hallmark sets over 13,738 ranked
# genes, then taxonomic ranks over 53 ranked species.

suppressPackageStartupMessages({library(tidyverse); library(fgsea); library(msigdbr)})
here <- "materials/2-workshop2/5-enrichment"
out <- file.path(here, "results"); dir.create(out, showWarnings = FALSE, recursive = TRUE)
source("materials/2-workshop2/R/count_tables.R")
m3 <- "materials/2-workshop2/3-differential/results"
set.seed(1)

hallmark <- msigdbr(species = "Homo sapiens", collection = "H")
hallmark_sets <- split(hallmark$ensembl_gene, hallmark$gs_name)
cat("Hallmark sets:", length(hallmark_sets), "\n")

ranking <- function(file) {
  res <- read_csv(file.path(m3, file), show_col_types = FALSE) |> filter(!is.na(stat))
  setNames(res$stat, res$feature)
}

run_fgsea <- function(stat, sets, label, min_size = 15, max_size = 500) {
  fgsea(sets, stat, minSize = min_size, maxSize = max_size, nproc = 1) |>
    as_tibble() |>
    mutate(ranking = label,
           leadingEdge = map_chr(leadingEdge, ~ paste(head(.x, 12), collapse = ";"))) |>
    arrange(padj)
}

# --- gene sets over the gene ranking ----------------------------------------
stat_out  <- ranking("genes_by_outcome.csv.gz")
stat_comm <- ranking("genes_by_community.csv.gz")

gsea <- bind_rows(
  run_fgsea(stat_out,  hallmark_sets, "preterm vs term"),
  run_fgsea(stat_comm, hallmark_sets, "anaerobe vs Lactobacillus")
)
write_csv(gsea, file.path(out, "gsea_hallmark.csv"))
print(gsea |> filter(ranking == "preterm vs term") |> select(pathway, NES, padj, size) |> head(8))

# --- the ranked list itself, so the slides can draw a running score ---------
top_set <- gsea |> filter(ranking == "preterm vs term") |> slice_min(padj, n = 1) |> pull(pathway)
cat("top set for the running-score figure:", top_set, "\n")
sym <- read_counts("instructional")$gene_info |> select(feature = ensembl_gene_id, symbol)
ranked <- tibble(feature = names(stat_out), stat = stat_out) |>
  arrange(desc(stat)) |>
  mutate(rank = row_number(),
         in_top_set = feature %in% hallmark_sets[[top_set]],
         set = top_set) |>
  left_join(sym, by = "feature")
write_csv(ranked, file.path(out, "gene_ranking_outcome.csv.gz"))

# Running enrichment score for that one set, precomputed so the slide is a plot.
running_score <- function(stat, set) {
  s <- sort(stat, decreasing = TRUE)
  hit <- names(s) %in% set
  nr <- sum(abs(s[hit]))
  inc <- ifelse(hit, abs(s) / nr, -1 / sum(!hit))
  tibble(rank = seq_along(s), hit = hit, running = cumsum(inc))
}
write_csv(running_score(stat_out, hallmark_sets[[top_set]]) |> mutate(set = top_set),
          file.path(out, "running_score.csv.gz"))

# --- how lopsided is each ranking? the check module 5 teaches ---------------
#
# The naive check, "is the distribution of statistics centred", passes both of
# our rankings. It is the wrong check. GSEA weights genes by the magnitude of
# their statistic, so what matters is whether the STRONGEST effects point one
# way. We therefore report the share positive as a function of how large an
# effect you insist on.
balance_inputs <- list(`preterm vs term` = stat_out,
                       `anaerobe vs Lactobacillus` = stat_comm)

tail_balance <- imap_dfr(balance_inputs, function(st, nm) {
  tibble(min_abs_stat = seq(0, 8, by = 0.25)) |>
    rowwise() |>
    mutate(n_genes = sum(abs(st) >= min_abs_stat),
           share_positive = mean(st[abs(st) >= min_abs_stat] > 0)) |>
    ungroup() |>
    filter(n_genes >= 50) |>
    mutate(ranking = nm)
})
write_csv(tail_balance, file.path(out, "ranking_balance.csv"))

balance_summary <- imap_dfr(balance_inputs, function(st, nm) {
  res <- read_csv(file.path(m3, if (nm == "preterm vs term") "genes_by_outcome.csv.gz"
                            else "genes_by_community.csv.gz"), show_col_types = FALSE) |>
    filter(!is.na(stat))
  top500 <- res |> slice_max(abs(stat), n = 500)
  tibble(ranking = nm,
         n_genes = length(st),
         share_positive_all = mean(st > 0),
         n_significant = sum(res$padj < 0.05, na.rm = TRUE),
         share_positive_significant = mean(res$stat[which(res$padj < 0.05)] > 0),
         share_positive_top500 = mean(top500$stat > 0))
})
write_csv(balance_summary, file.path(out, "ranking_balance_summary.csv"))
print(balance_summary)

# --- taxonomic ranks as sets over the species ranking -----------------------
sp <- read_csv(file.path(m3, "species_by_outcome.csv"), show_col_types = FALSE) |> filter(!is.na(stat))
stat_sp <- setNames(sp$stat, sp$feature)
genus_sets  <- split(sp$feature, sp$Genus)  |> keep(~ length(.x) >= 2)
family_sets <- split(sp$feature, sp$Family) |> keep(~ length(.x) >= 2)
cat("genus sets with 2+ species:", length(genus_sets),
    " family sets:", length(family_sets), "\n")

tax_sets <- bind_rows(
  run_fgsea(stat_sp, genus_sets,  "genus",  min_size = 2, max_size = 50),
  run_fgsea(stat_sp, family_sets, "family", min_size = 2, max_size = 50)
)
write_csv(tax_sets, file.path(out, "tax_set_enrichment.csv"))
print(tax_sets |> select(ranking, pathway, size, NES, padj) |> head(10))

write_csv(tibble(feature = names(stat_sp), stat = stat_sp) |>
            arrange(desc(stat)) |> mutate(rank = row_number()) |>
            left_join(sp |> select(feature, Genus, Family, species), by = "feature"),
          file.path(out, "species_ranking_outcome.csv"))

cat("module 5 results written to", out, "\n")
