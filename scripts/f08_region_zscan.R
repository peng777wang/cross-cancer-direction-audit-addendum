# Regional z scan at the reversed loci.
#
# Fine-mapping the reversed arm is uninformative by construction (no own-GWAS
# signal -> no credible set), so we test the localization-artifact hypothesis
# directly instead: does the reversed arm's own data contain ANY variant with
# a sign opposing the MTAG lead and |z| > 2 anywhere in the region?
#
# No LD matrix and no susie_rss is used, so this analysis cannot be affected
# by reference-panel LD mismatch.
#
# Usage: Rscript f08_region_zscan.R

ROOT <- normalizePath(Sys.getenv("CGS_ROOT", unset = getwd()), mustWork = FALSE)
OUT  <- file.path(ROOT, "f08_region_zscan.csv")
Z_CUT <- 2.0

loci <- trimws(readLines(file.path(ROOT, "loci273.txt")))
loci <- loci[nzchar(loci)]

res <- read.csv(file.path(ROOT, "pleiotropy_resource.csv"))
res$cls <- sub(" .*", "", as.character(res$class_mtag))
res$mtag_significant <- as.logical(res$mtag_significant)

f07 <- read.csv(file.path(ROOT, "f07_p0_credsets.csv"))
rev_tags <- unique(f07$locus)

out <- list()
for (tag in rev_tags) {
  i   <- as.integer(sub("L", "", tag))
  snp <- loci[i]
  fz  <- file.path(ROOT, sprintf("reg_%s_z.tsv", tag))
  if (!file.exists(fz)) next
  z <- read.table(fz, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)

  # the reversed ARM is the one whose own GWAS is not significant; the other
  # arm of the same locus is the concordant (own-significant) cancer
  res$own_gwas_significant <- as.logical(res$own_gwas_significant)
  arms <- res[res$snp == snp & res$cls == "reversed" & res$mtag_significant &
              !res$own_gwas_significant, ]
  for (a in seq_len(nrow(arms))) {
    tr  <- arms$trait[a]
    col <- paste0("z_", tr)
    if (!col %in% names(z)) next
    v <- z[[col]]
    v <- v[!is.na(v)]
    if (!length(v)) next

    # MTAG direction for this cancer at the lead variant
    lead_mtag <- arms$z_mtag[a]
    lead_own  <- arms$z_raw[a]
    sgn <- sign(lead_mtag)

    max_abs <- max(abs(v))
    max_signed <- v[which.max(abs(v))]
    opp <- v[sgn * v < 0]
    n_opp <- sum(abs(opp) > Z_CUT)

    out[[length(out) + 1]] <- data.frame(
      locus = tag, lead_snp = snp, arm = tr,
      n_variants = length(v),
      lead_own_z = round(lead_own, 3),
      lead_mtag_z = round(lead_mtag, 3),
      own_z_at_lead = round(lead_own, 3),
      max_abs_own_z = round(max_abs, 3),
      max_own_z_signed = round(max_signed, 3),
      opposes_lead = sgn * max_signed < 0,
      n_opposing_gt2 = n_opp,
      max_abs_opposing = if (length(opp)) round(max(abs(opp)), 3) else 0)
  }
}

out <- do.call(rbind, out)
write.csv(out, OUT, row.names = FALSE)

cat("loci:", nrow(out), "\n")
cat("median |z_own| at lead:", round(median(abs(out$own_z_at_lead)), 2), "\n")
cat("median max|z_own| in region:", round(median(out$max_abs_own_z), 2),
    " (range ", round(min(out$max_abs_own_z), 2), "-", round(max(out$max_abs_own_z), 2), ")\n", sep = "")
cat("arms with any opposing variant |z|>2:", sum(out$n_opposing_gt2 > 0), "of", nrow(out), "\n")
cat("median count of opposing |z|>2 variants:", median(out$n_opposing_gt2), "\n")
cat("max opposing |z| across all arms:", round(max(out$max_abs_opposing), 2), "\n")
cat("F08_DONE\n")
