# Precompute the Module 3 differential-testing results.
#
# Run by hand from the website root:
#   RENV_CONFIG_AUTOLOADER_ENABLED=FALSE \
#     Rscript materials/2-workshop2/3-differential/_build_results.R
#
# Four contrasts, all at week_18 so there is one sample per person and no
# repeated-measures question to answer in a 90-minute block.

suppressPackageStartupMessages({library(tidyverse); library(DESeq2)})
here <- "materials/2-workshop2/3-differential"
out <- file.path(here, "results"); dir.create(out, showWarnings = FALSE, recursive = TRUE)
source("materials/2-workshop2/R/count_tables.R")

d <- read_counts("instructional")

# The grouping participants build in Module 2: is this community mostly
# Lactobacillus or mostly anaerobes? Defined from the species table only.
lacto <- intersect(d$taxonomy$taxon[!is.na(d$taxonomy$Genus) & d$taxonomy$Genus == "Lactobacillus"],
                   rownames(d$species))
lacto_fraction <- colSums(d$species[lacto, ]) / colSums(d$species)

md <- d$metadata |>
  mutate(lacto_fraction = lacto_fraction[sample_id],
         community = if_else(lacto_fraction >= 0.5, "lactobacillus", "anaerobe"),
         bv = if_else(nugent >= 7, "bv", "no_bv"))
write_csv(md |> select(sample_id, participant_id, time_point, outcome, nugent,
                       lacto_fraction, community, bv),
          file.path(out, "sample_groups.csv"))

w18 <- md |> filter(time_point == "week_18")

deseq_contrast <- function(counts, meta, group_col, hi, lo, min_count, min_frac) {
  meta <- meta |> filter(.data[[group_col]] %in% c(hi, lo))
  m <- counts[, meta$sample_id, drop = FALSE]
  keep <- rowSums(m >= min_count) >= ceiling(min_frac * ncol(m))
  design <- as.formula(paste("~", group_col))
  dds <- DESeqDataSetFromMatrix(m[keep, ], as.data.frame(meta), design)
  dds <- DESeq(dds, quiet = TRUE)
  res <- results(dds, contrast = c(group_col, hi, lo))
  list(n_in = sum(keep), n_tested = ncol(m),
       res = as_tibble(as.data.frame(res), rownames = "feature"))
}

label_genes <- function(res) res |>
  left_join(d$gene_info |> select(feature = ensembl_gene_id, symbol, gene_name), by = "feature")
label_species <- function(res) res |>
  left_join(d$taxonomy |> select(feature = taxon, Genus, Family, species), by = "feature")

runs <- list()

# Taught contrast: group from the species table, test the gene table.
g_comm <- deseq_contrast(d$genes, w18, "community", "anaerobe", "lactobacillus", 10, 0.25)
write_csv(label_genes(g_comm$res), file.path(out, "genes_by_community.csv.gz"))
runs$genes_by_community <- g_comm

# Exercise contrast: clinical outcome, tested on the species table.
s_out <- deseq_contrast(d$species, w18, "outcome", "preterm", "term", 1, 0.10)
write_csv(label_species(s_out$res), file.path(out, "species_by_outcome.csv"))
runs$species_by_outcome <- s_out

# Shown for comparison: the same question asked of the gene table.
g_out <- deseq_contrast(d$genes, w18, "outcome", "preterm", "term", 10, 0.25)
write_csv(label_genes(g_out$res), file.path(out, "genes_by_outcome.csv.gz"))
runs$genes_by_outcome <- g_out

# Shown for comparison: a crude clinical proxy instead of the 16S-derived group.
g_bv <- deseq_contrast(d$genes, w18 |> filter(!is.na(nugent)), "bv", "bv", "no_bv", 10, 0.25)
write_csv(label_genes(g_bv$res), file.path(out, "genes_by_bv.csv.gz"))
runs$genes_by_bv <- g_bv

summary_tbl <- imap_dfr(runs, function(r, nm) {
  tibble(contrast = nm,
         features_in = r$n_in,
         samples = r$n_tested,
         n_sig = sum(r$res$padj < 0.05, na.rm = TRUE),
         n_sig_lfc1 = sum(r$res$padj < 0.05 & abs(r$res$log2FoldChange) > 1, na.rm = TRUE))
})
write_csv(summary_tbl, file.path(out, "contrast_summary.csv"))
print(summary_tbl)
cat("module 3 results written to", out, "\n")
