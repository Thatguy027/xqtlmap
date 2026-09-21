## xqtl_pipeline.R -----------------------------------------------------------
## Reusable functions for xQTL allele-frequency mapping from GATK
## ASEReadCounter tables. Sourced by run_xqtl.R; nothing here runs on load.
##
## Input contract
##   <experiment>/aser/<sample_name>.table[.gz]   one per sample
##   <experiment>/sample_sheet.tsv                columns:
##       sample_name  parent1  parent2  condition  contrast_to
##     plus any number of optional extra columns (e.g. timepoint). Extra
##     columns become part of the sample label, in sheet order, inserted
##     between parent2 and condition.
##
##   contrast_to is a comma/semicolon-separated list of OTHER sample_names.
##   Each entry produces one contrast with this row's sample on the left (L)
##   and the named sample on the right (R). Pairs are taken as written --
##   reciprocals are not added and not deduplicated.
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
  library(vcfR)
  library(AlphaSimR)
  library(Rfast)
  library(ggpubr)
  library(xQTLStats)
})

SHEET_REQUIRED <- c("sample_name", "parent1", "parent2", "condition", "contrast_to")

## Columns that are per-sample *measurements*, not design factors. They are read
## off the sheet and used, but kept out of the label -- and therefore out of every
## cache key and output filename -- so that recording a pool size does not rename
## a contrast. Anything not listed here and not required still becomes a label
## component, which is what makes `timepoint` and `stage` show up in filenames.
SHEET_NONLABEL <- c("n_worms")

## -- small utilities --------------------------------------------------------

msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")

`%||%` <- function(a, b) if (is.null(a)) b else a

## Resolve a path that may be absolute or relative to `base`.
resolve_path <- function(p, base) {
  if (grepl("^(/|~)", p)) path.expand(p) else file.path(base, p)
}

## Read a .qs or .qs2 serialized object, whichever the extension says.
read_serialized <- function(path) {
  if (!file.exists(path)) stop("serialized object not found: ", path, call. = FALSE)
  if (grepl("\\.qs2$", path)) {
    if (!requireNamespace("qs2", quietly = TRUE))
      stop("package 'qs2' is required to read ", path, call. = FALSE)
    qs2::qs_read(path)
  } else if (grepl("\\.qs$", path)) {
    if (!requireNamespace("qs", quietly = TRUE))
      stop("package 'qs' is required to read ", path, call. = FALSE)
    qs::qread(path)
  } else if (grepl("\\.rds$", path, ignore.case = TRUE)) {
    readRDS(path)
  } else {
    stop("don't know how to read ", path, " (expected .qs, .qs2 or .rds)", call. = FALSE)
  }
}

## Source the xQTLSims helper functions this pipeline depends on.
## NOTE: helperFxs.R defines createFounderPop twice; sourcing the whole file
## leaves the later (gmap-aware) definition in scope, which is the one we want.
load_xqtlsims <- function(xQTLSims.dir) {
  source.dir <- file.path(xQTLSims.dir, "R")
  for (f in c("simWormCrosses.R", "helperFxs.R")) {
    p <- file.path(source.dir, f)
    if (!file.exists(p)) stop("missing xQTLSims source: ", p, call. = FALSE)
    source(p, local = FALSE)
  }
  invisible(TRUE)
}

## -- sample sheet -----------------------------------------------------------

## Read + validate the sample sheet. Returns a tibble with the required
## columns plus `label` (used for filenames and column suffixes) and
## `extra_cols` recorded as an attribute.
## Accepts .tsv or .csv, sniffed from the header line rather than the extension:
## a table exported from a spreadsheet is routinely named .tsv and comma separated.
read_sample_sheet <- function(path) {
  if (!file.exists(path)) stop("metadata table not found: ", path, call. = FALSE)
  hdr <- readLines(path, n = 1L, warn = FALSE)
  delim <- if (lengths(regmatches(hdr, gregexpr("\\t", hdr))) >=
               lengths(regmatches(hdr, gregexpr(",", hdr)))) "\t" else ","
  sheet <- readr::read_delim(path, delim = delim, progress = FALSE,
                             col_types = readr::cols(.default = readr::col_character()))

  missing <- setdiff(SHEET_REQUIRED, names(sheet))
  if (length(missing))
    stop("sample sheet is missing required column(s): ", paste(missing, collapse = ", "),
         call. = FALSE)

  ## trim whitespace everywhere; turn "" into NA
  sheet <- sheet %>%
    dplyr::mutate(dplyr::across(dplyr::everything(), ~ {
      x <- stringr::str_trim(dplyr::coalesce(.x, ""))
      dplyr::na_if(x, "")
    }))

  ## required fields other than contrast_to must be present
  for (col in setdiff(SHEET_REQUIRED, "contrast_to")) {
    bad <- which(is.na(sheet[[col]]))
    if (length(bad))
      stop("sample sheet row(s) ", paste(bad, collapse = ", "),
           " have an empty '", col, "'", call. = FALSE)
  }

  dup <- sheet$sample_name[duplicated(sheet$sample_name)]
  if (length(dup))
    stop("duplicated sample_name in sample sheet: ",
         paste(unique(dup), collapse = ", "), call. = FALSE)

  ## label = sample_name _ parent1 _ parent2 _ <extra cols, in sheet order> _ condition
  extra_cols <- setdiff(names(sheet), c(SHEET_REQUIRED, SHEET_NONLABEL))
  label_cols <- c("sample_name", "parent1", "parent2", extra_cols, "condition")
  sheet$label <- apply(sheet[, label_cols, drop = FALSE], 1,
                       function(r) paste(r[!is.na(r)], collapse = "_"))

  if (anyDuplicated(sheet$label))
    stop("sample labels are not unique: ",
         paste(unique(sheet$label[duplicated(sheet$label)]), collapse = ", "), call. = FALSE)

  attr(sheet, "extra_cols") <- extra_cols
  sheet
}

## Expand contrast_to into a table of L/R sample pairs, in sheet order.
build_contrast_table <- function(sheet) {
  rows <- list()
  for (i in seq_len(nrow(sheet))) {
    raw <- sheet$contrast_to[i]
    if (is.na(raw)) next
    targets <- stringr::str_trim(stringr::str_split(raw, "[,;]")[[1]])
    targets <- targets[nzchar(targets)]
    for (tgt in targets) {
      if (identical(tgt, sheet$sample_name[i]))
        stop("sample '", tgt, "' lists itself in contrast_to", call. = FALSE)
      j <- match(tgt, sheet$sample_name)
      if (is.na(j))
        stop("contrast_to for sample '", sheet$sample_name[i],
             "' references unknown sample '", tgt, "'", call. = FALSE)
      rows[[length(rows) + 1L]] <- tibble::tibble(
        L_sample = sheet$sample_name[i], R_sample = tgt,
        L_label  = sheet$label[i],       R_label  = sheet$label[j],
        L_row = i, R_row = j
      )
    }
  }
  if (!length(rows))
    return(tibble::tibble(L_sample = character(), R_sample = character(),
                          L_label = character(), R_label = character(),
                          L_row = integer(), R_row = integer()))
  ct <- dplyr::bind_rows(rows)

  dup <- duplicated(ct[, c("L_sample", "R_sample")])
  if (any(dup)) {
    warning("dropping ", sum(dup), " exactly duplicated contrast(s) in the sample sheet",
            call. = FALSE)
    ct <- ct[!dup, , drop = FALSE]
  }
  ct
}

## Locate the ASEReadCounter table for each sample. Accepts .table / .table.gz
## (and .tsv variants), which is why the old hardcoded ".table" suffix broke
## once the files were gzipped.
find_aser_files <- function(aser.dir, sample_names) {
  if (!dir.exists(aser.dir)) stop("aser directory not found: ", aser.dir, call. = FALSE)
  exts <- c(".table", ".table.gz", ".tsv", ".tsv.gz", ".txt", ".txt.gz")
  paths <- vapply(sample_names, function(sn) {
    cand <- file.path(aser.dir, paste0(sn, exts))
    hit <- cand[file.exists(cand)]
    if (!length(hit))
      stop("no ASEReadCounter table for sample '", sn, "' in ", aser.dir,
           " (looked for ", paste0(sn, exts, collapse = ", "), ")", call. = FALSE)
    if (length(hit) > 1L)
      warning("multiple files for sample '", sn, "'; using ", basename(hit[1]), call. = FALSE)
    hit[1]
  }, character(1))
  stats::setNames(paths, sample_names)
}

## -- count tables -----------------------------------------------------------

ASER_COLS <- c("contig", "position", "refCount", "altCount")

## Build a phased p1/p2 count table for one sample.
## Reimplements xQTLSims::makeCountTables for a single row so that file paths
## are explicit (no sample.dir / sample.suffix globals) and gzipped input works.
make_count_table_one <- function(aser.path, p1, p2, vcf, gt, gmap, X.drop = FALSE) {
  founderPop <- createFounderPop(vcf, gt, c(p1, p2), gmap, X.drop = X.drop)
  genMap <- getGenMap(founderPop)

  scounts <- readr::read_tsv(aser.path, progress = FALSE, show_col_types = FALSE)
  missing <- setdiff(ASER_COLS, names(scounts))
  if (length(missing))
    stop(basename(aser.path), " is missing column(s): ", paste(missing, collapse = ", "),
         call. = FALSE)

  scounts <- tibble::tibble(
    id  = paste0(scounts$contig, "_", scounts$position),
    ref = scounts$refCount,
    alt = scounts$altCount
  )
  ## defensive: a duplicated site would fan out in the left_join below
  if (anyDuplicated(scounts$id)) {
    warning(basename(aser.path), ": dropping ", sum(duplicated(scounts$id)),
            " duplicated site(s)", call. = FALSE)
    scounts <- dplyr::distinct(scounts, id, .keep_all = TRUE)
  }

  ## Keep every genetic-map marker; sites absent from the ASE table get 0/0.
  ## (GATK ASEReadCounter emits different site sets for different BAMs.)
  scounts <- dplyr::left_join(genMap, scounts, by = "id")
  n.missing <- sum(is.na(scounts$ref))
  scounts$ref[is.na(scounts$ref)] <- 0
  scounts$alt[is.na(scounts$alt)] <- 0
  names(scounts)[1] <- "ID"

  countdf <- phaseBiparental(scounts, p1, founderPop, genMap)
  attr(countdf, "p1") <- p1
  attr(countdf, "p2") <- p2
  attr(countdf, "n.markers") <- nrow(countdf)
  attr(countdf, "n.missing") <- n.missing
  countdf
}

## Per-sample count tables + AFDs, cached to <cache.dir>/<label>.rds so that
## re-running after a sample-sheet edit only recomputes what changed.
compute_sample_afds <- function(sheet, aser.paths, vcf, gt, gmap, cfg, cache.dir,
                                force = FALSE) {
  dir.create(cache.dir, showWarnings = FALSE, recursive = TRUE)
  countdfs <- list(); afds <- list()

  for (i in seq_len(nrow(sheet))) {
    lab <- sheet$label[i]; sn <- sheet$sample_name[i]
    cache.file <- file.path(cache.dir, paste0(lab, ".rds"))

    ## calcAFD's sample.size is the number of individuals in the pool; it caps the
    ## effective n behind afd.se via n = min(sample.size * sel.strength, depth-based
    ## n). A single config-wide value is wrong whenever pools differ in size, so an
    ## n_worms column on the sheet overrides it per sample.
    afd.params <- cfg$afd.params
    if ("n_worms" %in% names(sheet) && !is.na(sheet$n_worms[i]) &&
        nzchar(sheet$n_worms[i])) {
      afd.params$sample.size <- as.numeric(sheet$n_worms[i])
    } else if ("n_worms" %in% names(sheet)) {
      msg("  ", lab, ": n_worms not recorded, falling back to sample.size=",
          cfg$sample.size)
    }

    if (!force && file.exists(cache.file)) {
      cached <- readRDS(cache.file)
      if (identical(cached$params, afd.params) &&
          identical(cached$aser.mtime, as.numeric(file.mtime(aser.paths[[sn]])))) {
        msg("  ", lab, ": cached")
        countdfs[[lab]] <- cached$countdf; afds[[lab]] <- cached$afd
        next
      }
      msg("  ", lab, ": cache stale, recomputing")
    }

    msg("  ", lab, ": building count table")
    countdf <- make_count_table_one(aser.paths[[sn]], sheet$parent1[i], sheet$parent2[i],
                                    vcf, gt, gmap, X.drop = cfg$X.drop)
    msg("  ", lab, ": ", attr(countdf, "n.markers"), " markers (",
        attr(countdf, "n.missing"), " with no ASE coverage); calcAFD")

    afd <- calcAFD(countdf, experiment.name = lab,
                   sample.size  = afd.params$sample.size,
                   sel.strength = cfg$sel.strength,
                   bin.width    = cfg$bin.width,
                   eff.length   = cfg$eff.length,
                   uchr         = cfg$uchr)

    saveRDS(list(countdf = countdf, afd = afd, params = afd.params,
                 aser.mtime = as.numeric(file.mtime(aser.paths[[sn]]))),
            cache.file)
    countdfs[[lab]] <- countdf; afds[[lab]] <- afd
  }
  list(countdfs = countdfs, afds = afds)
}

## -- naming -----------------------------------------------------------------

## Contrast basename, built from the sample sheet columns directly rather than
## by re-parsing the label with positional str_split (which the old script did
## and which silently breaks whenever the column set changes).
contrast_basename <- function(sheet, L_row, R_row, bin.width) {
  L <- sheet[L_row, ]; R <- sheet[R_row, ]
  parents <- paste0(L$parent1, "_", L$parent2)
  if (!(identical(L$parent1, R$parent1) && identical(L$parent2, R$parent2))) {
    warning("contrasting samples from different crosses: ", L$label, " vs ", R$label,
            call. = FALSE)
    parents <- paste0(L$parent1, "_", L$parent2, "_vs_", R$parent1, "_", R$parent2)
  }
  ## keep the historical F<tp_L>-<tp_R> element when a timepoint column exists
  tp <- ""
  if ("timepoint" %in% names(sheet) && !is.na(L$timepoint) && !is.na(R$timepoint))
    tp <- paste0("_F", L$timepoint, "-", R$timepoint)

  paste0(parents, tp, "_contrast_", L$condition, "-", R$condition, "_", bin.width)
}

## -- contrasts --------------------------------------------------------------

## LOD-drop support intervals, one row per chromosome, sorted by LOD.
lod_intervals <- function(results, lod.drop = 2) {
  results %>%
    dplyr::group_by(chrom) %>%
    dplyr::filter(LOD > max(LOD) - lod.drop) %>%
    dplyr::mutate(lcon = min(physical.position), rcon = max(physical.position)) %>%
    dplyr::arrange(dplyr::desc(LOD)) %>%
    dplyr::distinct(chrom, lcon, rcon, .keep_all = TRUE) %>%
    dplyr::mutate(marker = paste0(chrom, ":", lcon, "-", rcon)) %>%
    dplyr::select(marker, physical.position, LOD) %>%
    dplyr::ungroup()
}

run_contrasts <- function(contrasts, sheet, afds, plots, cfg, out.dir) {
  dir.create(out.dir, showWarnings = FALSE, recursive = TRUE)
  out <- list()

  for (k in seq_len(nrow(contrasts))) {
    Ll <- contrasts$L_label[k]; Rl <- contrasts$R_label[k]
    base <- contrast_basename(sheet, contrasts$L_row[k], contrasts$R_row[k], cfg$bin.width)
    msg("  ", Ll, " vs ", Rl, " -> ", base)

    results <- calcContrastStats(results = list(afds[[Ll]], afds[[Rl]]),
                                 L = paste0("_", Ll), R = paste0("_", Rl))
    ## replace the underflow-prone p/LOD columns with log-space versions
    results <- recompute_lod(results, cfg$lod.convention)

    interval_df <- lod_intervals(results, cfg$LOD.drop)

    utils::write.table(interval_df, file.path(out.dir, paste0(base, ".tsv")),
                       col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
    utils::write.table(results, file.path(out.dir, paste0(base, "_plot_DF.tsv")),
                       col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

    sc <- plotContrast(results, suffix1 = Ll, suffix2 = Rl)
    s  <- plotSummary(results, effective.n.tests = cfg$effective.n.tests)

    panel <- ggpubr::ggarrange(plots[[Ll]], plots[[Rl]], s, nrow = 3)
    ggplot2::ggsave(file.path(out.dir, paste0(base, ".png")), plot = panel,
                    height = cfg$contrast.plot.height, width = cfg$contrast.plot.width)

    pub <- write_pub_outputs(results, base, out.dir, cfg,
                             title = paste0(sheet$condition[contrasts$L_row[k]], " vs ",
                                            sheet$condition[contrasts$R_row[k]]))

    out[[base]] <- list(results = results, intervals = interval_df,
                        peak.intervals = pub$intervals,
                        contrast.plot = sc, summary.plot = s, pub.plot = pub$plot)
  }
  out
}

## -- LOD in log space -------------------------------------------------------
##
## xQTLStats::calcContrastStats computes
##     p   = 2 * pnorm(abs(z), lower.tail = FALSE)
##     LOD = PvalToLOD(p) = qchisq(2*p, df = 1, lower.tail = FALSE) / (2*log(10))
##
## `p` underflows to exactly 0 once |z| >~ 38.5 (the double-precision floor is
## ~5e-324), at which point qchisq(0, lower.tail = FALSE) returns Inf. So LOD
## does not "cap" at 300 -- it silently becomes Inf, and -log10(p) saturates at
## ~308. That is not just cosmetic: any `LOD > max(LOD) - drop` filter evaluates
## to `LOD > Inf`, which is FALSE everywhere, so the entire chromosome carrying
## the strongest signal drops out of the interval table.
##
## Fix: never form `p` on the linear scale. pnorm(log.p = TRUE) and
## qchisq(log.p = TRUE) work far into the tail, so both LOD and -log10(p) can be
## evaluated exactly for any |z|.
##
## convention:
##   "package" -- bit-for-bit what PvalToLOD() returns, just without the
##                underflow. Use this to stay comparable with earlier results.
##   "chisq"   -- the textbook 1-df LOD, z^2 / (2*ln10).
##
## The two differ by exactly log10(2) = 0.30103 asymptotically, because
## PvalToLOD() expects a ONE-tailed p (it doubles its argument to get the
## chi-square tail probability) but calcContrastStats hands it a TWO-tailed p.
## The package's LOD is therefore deflated by ~0.3 relative to z^2/(2*ln10).
recompute_lod <- function(results, convention = c("package", "chisq")) {
  convention <- match.arg(convention)
  z <- results$z

  ## log of the two-tailed normal p-value; exact for any |z|
  log_p <- log(2) + stats::pnorm(abs(z), lower.tail = FALSE, log.p = TRUE)

  results$log10p    <- log_p / log(10)
  results$neglog10p <- -log_p / log(10)

  results$LOD <- if (convention == "chisq") {
    z^2 / (2 * log(10))
  } else {
    ## PvalToLOD(p) with its argument doubled, all in log space
    ## qchisq's argument 2*p exceeds 1 for |z| <= qnorm(.25, lower = FALSE)
    ## = 0.6745, giving NaN. PvalToLOD() maps that to 0; do the same.
    l <- suppressWarnings(
      stats::qchisq(log(2) + log_p, df = 1, lower.tail = FALSE, log.p = TRUE)) /
      (2 * log(10))
    l[is.nan(l)] <- 0
    l
  }
  results$LOD[is.na(z)] <- NA_real_
  attr(results, "lod.convention") <- convention
  results
}

## Genome-wide significance threshold on the LOD scale, matching `convention`.
lod_threshold <- function(alpha = 0.05, effective.n.tests = 2000,
                          convention = c("package", "chisq")) {
  convention <- match.arg(convention)
  p <- alpha / effective.n.tests
  if (convention == "chisq") stats::qchisq(p,     df = 1, lower.tail = FALSE) / (2 * log(10))
  else                       stats::qchisq(2 * p, df = 1, lower.tail = FALSE) / (2 * log(10))
}

## -- support intervals ------------------------------------------------------

## Contiguous LOD-drop support interval around each chromosome's peak.
##
## Differs deliberately from lod_intervals(): that one (inherited from the
## original script) takes min/max of ALL positions above the threshold on a
## chromosome, so two separated peaks get merged into one interval spanning the
## gap between them. This walks outward from the peak and stops at the first
## marker that falls below peak - lod.drop, which is what a support interval
## normally means.
## `lod.drop.frac`, when set, makes the drop a fraction of each chromosome's own
## peak instead of an absolute number of LOD units. Calibrated on a cross whose
## causal gene was known, where the causal allele sat 18.9 LOD below a 334.7
## peak (5.6%). An absolute drop of that size is destructive on weak peaks
## -- 20 LOD below a LOD-31 peak covers an entire chromosome -- whereas the
## proportional form transfers across peak strengths.
peak_intervals <- function(results, lod.drop = 1.5, threshold = NULL,
                           lod.drop.frac = NULL) {
  res <- results[is.finite(results$LOD), c("chrom", "physical.position", "LOD")]
  if (!nrow(res)) return(NULL)

  out <- lapply(split(res, res$chrom), function(d) {
    d <- d[order(d$physical.position), ]
    i <- which.max(d$LOD)
    drop <- if (!is.null(lod.drop.frac)) d$LOD[i] * lod.drop.frac else lod.drop
    above <- d$LOD > (d$LOD[i] - drop)
    ## keep only the run of TRUEs containing the peak
    run <- cumsum(c(TRUE, diff(above) != 0))
    sel <- which(run == run[i])
    tibble::tibble(
      chrom         = as.character(d$chrom[1]),
      peak.position = d$physical.position[i],
      peak.LOD      = d$LOD[i],
      lod.drop.used = drop,
      lcon          = d$physical.position[min(sel)],
      rcon          = d$physical.position[max(sel)],
      marker        = paste0(d$chrom[1], ":", d$physical.position[min(sel)],
                             "-", d$physical.position[max(sel)])
    )
  })
  out <- dplyr::bind_rows(out)
  out$width.kb   <- (out$rcon - out$lcon) / 1e3
  out$significant <- if (is.null(threshold)) TRUE else out$peak.LOD > threshold
  out %>% dplyr::arrange(dplyr::desc(peak.LOD))
}

## -- publication figure -----------------------------------------------------

## Clean single-panel LOD trace: one facet per chromosome, the QTL support
## interval shaded under the curve, a dashed genome-wide threshold.
plot_lod_trace <- function(results, intervals, threshold = NULL, title = NULL,
                           lod.drop = 1.5, lod.drop.frac = NULL,
                           shade.only.significant = TRUE,
                           y.trans = "identity", base_size = 11) {

  d <- results[is.finite(results$LOD), c("chrom", "physical.position", "LOD")]
  d$chrom <- factor(d$chrom, levels = intersect(
    c(as.character(utils::as.roman(1:5)), "X"), unique(as.character(d$chrom))))
  d <- d[!is.na(d$chrom), ]
  d$pos.mb <- d$physical.position / 1e6

  shade <- intervals
  if (shade.only.significant && "significant" %in% names(shade))
    shade <- shade[shade$significant, , drop = FALSE]
  if (!is.null(shade) && nrow(shade))
    shade$chrom <- factor(shade$chrom, levels = levels(d$chrom))

  ## the trace restricted to each support interval, for the ribbon
  band <- NULL
  if (!is.null(shade) && nrow(shade)) {
    band <- dplyr::bind_rows(lapply(seq_len(nrow(shade)), function(k) {
      s <- shade[k, ]
      d[as.character(d$chrom) == as.character(s$chrom) &
          d$physical.position >= s$lcon & d$physical.position <= s$rcon, ]
    }))
    if (nrow(band)) band$chrom <- factor(band$chrom, levels = levels(d$chrom))
  }

  p <- ggplot2::ggplot(d, ggplot2::aes(x = pos.mb, y = LOD))

  ## Support intervals are shaded at their TRUE width, under the trace only.
  ## A full-height band was tried and rejected: it gives a wide, weak peak more
  ## visual weight than a narrow, very strong one, which inverts the reading.
  if (!is.null(band) && nrow(band)) {
    p <- p + ggplot2::geom_ribbon(data = band,
                                  ggplot2::aes(ymin = 0, ymax = LOD),
                                  fill = "#C4302B", alpha = 0.55)
  }

  p <- p + ggplot2::geom_line(linewidth = 0.45, colour = "grey12")

  if (!is.null(threshold))
    p <- p + ggplot2::geom_hline(yintercept = threshold, linetype = "dashed",
                                 linewidth = 0.35, colour = "grey40")

  if (!is.null(shade) && nrow(shade))
    p <- p + ggplot2::geom_point(data = shade, inherit.aes = FALSE,
                                 ggplot2::aes(x = peak.position / 1e6, y = peak.LOD),
                                 colour = "#C4302B", size = 1.1)

  ylab <- if (identical(y.trans, "sqrt")) "LOD (sqrt scale)" else "LOD"

  p +
    ggplot2::facet_grid(. ~ chrom, scales = "free_x", space = "free_x", switch = "x") +
    ggplot2::scale_y_continuous(trans = y.trans,
                                expand = ggplot2::expansion(mult = c(0, 0.08))) +
    ggplot2::scale_x_continuous(breaks = seq(0, 25, by = 5),
                                expand = ggplot2::expansion(mult = 0.03),
                                guide = ggplot2::guide_axis(check.overlap = TRUE)) +
    ggplot2::labs(x = "Position (Mb)", y = ylab, title = title) +
    ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      panel.spacing.x  = grid::unit(5, "pt"),
      strip.placement  = "outside",
      strip.background = ggplot2::element_blank(),
      strip.text       = ggplot2::element_text(face = "bold", size = base_size),
      axis.line.x      = ggplot2::element_line(linewidth = 0.3),
      axis.ticks.x     = ggplot2::element_line(linewidth = 0.3),
      axis.text.x      = ggplot2::element_text(size = base_size - 3.5),
      axis.title.x     = ggplot2::element_text(margin = ggplot2::margin(t = 4)),
      plot.title       = ggplot2::element_text(face = "bold", size = base_size, hjust = 0),
      plot.title.position = "plot",
      plot.margin      = ggplot2::margin(6, 10, 4, 6)
    )
}

## Write a publication figure + its interval table for one contrast.
write_pub_outputs <- function(results, base, out.dir, cfg, title = NULL) {
  thr <- lod_threshold(0.05, cfg$effective.n.tests, cfg$lod.convention)
  iv  <- peak_intervals(results, lod.drop = cfg$pub.lod.drop, threshold = thr,
                        lod.drop.frac = cfg$pub.lod.drop.frac)

  utils::write.table(iv, file.path(out.dir, paste0(base, "_intervals.tsv")),
                     col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

  p <- plot_lod_trace(results, iv, threshold = thr, title = title,
                      lod.drop = cfg$pub.lod.drop,
                      lod.drop.frac = cfg$pub.lod.drop.frac,
                      shade.only.significant = isTRUE(cfg$pub.shade.only.significant),
                      y.trans = cfg$pub.y.trans %||% "identity")

  for (fmt in cfg$pub.formats)
    ggplot2::ggsave(file.path(out.dir, paste0(base, "_LOD.", fmt)), p,
                    width = cfg$pub.plot.width, height = cfg$pub.plot.height,
                    dpi = 300)
  list(plot = p, intervals = iv, threshold = thr)
}
