# Fine-mapping the reversed loci (P0 tier) with crash-safe incremental output.
#
# * 42 reversed loci, only their MTAG-significant arms (~84 fits)
# * two-stage signal number: L = 3, escalate to L = 10 only if L = 3 saturates
# * every finished locus is appended to the CSV immediately, and loci already
#   in that CSV are skipped on restart, so a crash costs at most one locus
#
# Usage: Rscript f07_susie_p0.R [n_loci]

EXTRA_R_LIB <- Sys.getenv("CGS_R_LIB", unset = "")
if (nzchar(EXTRA_R_LIB)) .libPaths(c(EXTRA_R_LIB, .libPaths()))
suppressMessages(library(susieR))

ROOT <- normalizePath(Sys.getenv("CGS_ROOT", unset = getwd()), mustWork = FALSE)
OUTCSV <- file.path(ROOT, "f07_p0_credsets.csv")
L_SMALL <- 3
L_LARGE <- 10

cat("Rfast loaded:", requireNamespace("Rfast", quietly = TRUE), "\n")

args <- commandArgs(trailingOnly = TRUE)
NMAX <- if (length(args)) as.integer(args[1]) else 42L

loci <- trimws(readLines(file.path(ROOT, "loci273.txt")))
loci <- loci[nzchar(loci)]

res_tab <- read.csv(file.path(ROOT, "pleiotropy_resource.csv"))
res_tab$cls <- sub(" .*", "", as.character(res_tab$class_mtag))
res_tab$mtag_significant <- as.logical(res_tab$mtag_significant)

rev_snps <- unique(res_tab$snp[res_tab$cls == "reversed"])
idx_rev  <- head(which(loci %in% rev_snps), NMAX)
cat("P0 loci to consider:", length(idx_rev), "\n")

done <- character(0)
if (file.exists(OUTCSV)) {
  prev <- read.csv(OUTCSV)
  done <- unique(prev$locus)
  cat("resume: ", nrow(prev), " rows on disk from ", length(done), " loci\n", sep = "")
} else {
  prev <- NULL
}

COLS <- c("locus","index","lead_snp","trait","n_snp","n_cs","cs_sizes",
          "top_pip","top_snp","min_r2_in_cs","escalated","converged",
          "n_iter","saturated")

fit_one <- function(tag, i, tr) {
  fR <- file.path(ROOT, sprintf("ld_%s_R.tsv", tag))
  fm <- file.path(ROOT, sprintf("reg_%s_meta.tsv", tag))
  fz <- file.path(ROOT, sprintf("reg_%s_z.tsv", tag))
  fn <- file.path(ROOT, sprintf("reg_%s_n.txt", tag))
  if (!all(file.exists(fR, fm, fz, fn))) return(NULL)

  R    <- as.matrix(read.table(fR, header = TRUE, row.names = 1, check.names = FALSE))
  meta <- read.table(fm, header = TRUE, stringsAsFactors = FALSE)
  z    <- read.table(fz, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  nn   <- scan(fn, quiet = TRUE, what = "")
  nvec <- as.numeric(nn[seq(2, length(nn), 2)]); names(nvec) <- nn[seq(1, length(nn), 2)]

  col <- paste0("z_", tr)
  if (!col %in% names(z)) return(NULL)
  d <- data.frame(snp = z$snp, z = z[[col]])
  d <- d[!is.na(d$z), ]
  if (nrow(d) < 50) return(NULL)

  idx <- match(d$snp, meta$snp); ok <- !is.na(idx)
  d <- d[ok, ]; idx <- idx[ok]
  Rs <- (R[idx, idx] + t(R[idx, idx])) / 2
  Rs <- Rs + diag(max(0, -min(eigen(Rs, symmetric = TRUE,
                                    only.values = TRUE)$values) + 1e-4), nrow(Rs))

  s <- try(susie_rss(z = d$z, R = Rs, n = nvec[[tr]], L = L_SMALL,
                     max_iter = 1000), silent = TRUE)
  if (inherits(s, "try-error")) return(NULL)
  cs <- susie_get_cs(s, Xcorr = Rs, min_abs_corr = 0.5, coverage = 0.95)
  escalated <- FALSE
  if (length(cs$cs) >= L_SMALL) {
    escalated <- TRUE
    s2 <- try(susie_rss(z = d$z, R = Rs, n = nvec[[tr]], L = L_LARGE,
                        max_iter = 1000), silent = TRUE)
    if (!inherits(s2, "try-error")) {
      s <- s2
      cs <- susie_get_cs(s, Xcorr = Rs, min_abs_corr = 0.5, coverage = 0.95)
    }
  }

  pip <- s$pip
  min_r2 <- NA_real_
  if (length(cs$cs)) {
    first <- cs$cs[[which.max(vapply(cs$cs, function(x) max(pip[x]), 0))]]
    if (length(first) > 1) {
      sub <- Rs[first, first]; diag(sub) <- NA
      min_r2 <- suppressWarnings(min(sub, na.rm = TRUE))
    }
  }
  data.frame(locus = tag, index = i, lead_snp = loci[i], trait = tr,
             n_snp = nrow(d), n_cs = length(cs$cs),
             cs_sizes = paste(vapply(cs$cs, length, 0L), collapse = ";"),
             top_pip = round(max(pip, na.rm = TRUE), 4),
             top_snp = names(pip)[which.max(pip)],
             min_r2_in_cs = round(min_r2, 3),
             escalated = escalated, converged = s$converged, n_iter = s$niter,
             saturated = length(cs$cs) >= L_LARGE)
}

for (i in idx_rev) {
  tag <- sprintf("L%03d", i)
  if (tag %in% done) { cat("skip (done)", tag, "\n", flush = TRUE); next }

  arms <- res_tab[res_tab$snp == loci[i] & res_tab$cls == "reversed" &
                  res_tab$mtag_significant, ]
  rows <- list()
  for (a in seq_len(nrow(arms))) {
    r <- try(fit_one(tag, i, arms$trait[a]), silent = TRUE)
    if (!inherits(r, "try-error") && !is.null(r)) rows[[length(rows) + 1]] <- r
  }
  if (!length(rows)) { cat(tag, "no rows\n", flush = TRUE); next }

  chunk <- do.call(rbind, rows)[, COLS]
  write.table(chunk, OUTCSV, sep = ",", row.names = FALSE, col.names = !file.exists(OUTCSV),
              append = file.exists(OUTCSV), quote = TRUE)
  cat(sprintf("%s %s  +%d rows  (file now %d loci)\n", tag, loci[i], nrow(chunk),
              length(unique(read.csv(OUTCSV)$locus))), flush = TRUE)
}

if (file.exists(OUTCSV)) {
  out <- read.csv(OUTCSV)
  conc <- do.call(rbind, lapply(split(out, out$locus), function(g) if (nrow(g) == 2)
    data.frame(locus = g$locus[1], arm1 = g$trait[1], arm2 = g$trait[2],
               same_top_snp = g$top_snp[1] == g$top_snp[2]) else NULL))
  if (!is.null(conc)) {
    write.csv(conc, file.path(ROOT, "f07_lead_concordance.csv"), row.names = FALSE)
    cat("same top SNP across arms:", sum(conc$same_top_snp), "of", nrow(conc), "\n")
  }
  cat("F07_DONE rows =", nrow(out), "\n")
}
