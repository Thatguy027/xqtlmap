## run.R ---------------------------------------------------------------------
## The body of the `xqtlmap` command: validate inputs, build per-sample allele
## frequency deviations, plot them, and run every contrast the metadata asks for.
##
## Sourced by bin/xqtlmap, which has already set:
##   PKG_DIR   the repository root
##   opts      parsed --key=value flags
##   ASER_DIR  METADATA  OUT_DIR
## ---------------------------------------------------------------------------

## -- config ------------------------------------------------------------------

config.file <- opts$config %||% file.path(PKG_DIR, "config", "default.yml")
if (!file.exists(config.file)) stop("config not found: ", config.file, call. = FALSE)
cfg <- yaml::yaml.load_file(config.file)

## A config beside the metadata table overrides the defaults for this dataset
## only, so a cross at a different AIL generation needs no flag and no edit to
## the shared file. Precedence: command line > dataset > defaults.
ds.config <- file.path(dirname(METADATA), "xqtlmap.yml")
if (is.null(opts$config) && file.exists(ds.config)) {
  over <- yaml::yaml.load_file(ds.config)
  unknown <- setdiff(names(over), names(cfg))
  if (length(unknown))
    warning("keys in ", ds.config, " are not in the default config: ",
            paste(unknown, collapse = ", "), call. = FALSE)
  for (k in names(over)) cfg[[k]] <- over[[k]]
  message("[config] overlaid ", ds.config, ": ", paste(names(over), collapse = ", "))
}

coerce_like <- function(value, template) {
  if (is.logical(template)) return(as.logical(value))
  if (is.numeric(template)) return(as.numeric(value))
  if (is.list(template) || length(template) > 1) return(trimws(strsplit(value, ",")[[1]]))
  value
}
reserved <- c("aser", "metadata", "out", "config", "force", "check", "help")
for (k in setdiff(names(opts), reserved))
  cfg[[k]] <- if (k %in% names(cfg)) coerce_like(opts[[k]], cfg[[k]]) else opts[[k]]

cfg$uchr <- as.character(unlist(cfg$uchr))
## the parameters that invalidate a cached per-sample AFD when they change
cfg$afd.params <- cfg[c("sample.size", "sel.strength", "bin.width", "eff.length",
                        "uchr", "X.drop", "map.expansion.factor", "expandX")]

## -- reference data ----------------------------------------------------------

ref.dir <- path.expand(cfg$reference.dir %||% "")
if (!nzchar(ref.dir) || !dir.exists(ref.dir))
  stop("reference.dir not found: '", ref.dir, "'\n",
       "  Set it in config/default.yml, in an xqtlmap.yml beside your metadata,\n",
       "  or with --reference.dir=/path/to/xQTLSims. See README, 'Reference data'.",
       call. = FALSE)

gmap.file <- resolve_path(cfg$gmap.file, ref.dir)
vcf.file  <- resolve_path(cfg$vcf.file,  ref.dir)
gt.file   <- resolve_path(cfg$gt.file,   ref.dir)
for (f in c(gmap.file, vcf.file, gt.file))
  if (!file.exists(f)) stop("required reference file not found: ", f, call. = FALSE)

out.dir   <- OUT_DIR
plot.dir  <- file.path(out.dir, "plots")
cache.dir <- file.path(out.dir, "cache")

## -- validate ----------------------------------------------------------------

msg("aser     : ", ASER_DIR)
msg("metadata : ", METADATA)
msg("output   : ", out.dir)

sheet      <- read_sample_sheet(METADATA)
aser.paths <- find_aser_files(ASER_DIR, sheet$sample_name)
contrasts  <- build_contrast_table(sheet)

msg(nrow(sheet), " sample(s), ", nrow(contrasts), " contrast(s)")
cat("\n")
print(as.data.frame(sheet[, c(setdiff(names(sheet), "label"), "label")]), row.names = FALSE)

if (nrow(contrasts)) {
  plan <- data.frame(
    L = contrasts$L_label, R = contrasts$R_label,
    output = vapply(seq_len(nrow(contrasts)),
                    function(k) contrast_basename(sheet, contrasts$L_row[k],
                                                  contrasts$R_row[k], cfg$bin.width),
                    character(1)))
  cat("\ncontrasts:\n"); print(plan, row.names = FALSE)
} else {
  cat("\ncontrasts: none (contrast_to is empty for every sample)\n")
}
cat("\naser files:\n")
for (sn in names(aser.paths)) cat("  ", sn, " -> ", basename(aser.paths[[sn]]), "\n", sep = "")
cat("\n")

if (isTRUE(opts$check)) {
  msg("--check: inputs validated, stopping before any computation.")
  quit(status = 0)
}

dir.create(plot.dir,  showWarnings = FALSE, recursive = TRUE)
dir.create(cache.dir, showWarnings = FALSE, recursive = TRUE)
readr::write_tsv(sheet, file.path(out.dir, "metadata_resolved.tsv"))

## -- load reference ----------------------------------------------------------

msg("sourcing xQTLSims helpers from ", file.path(ref.dir, "R"))
load_xqtlsims(ref.dir)

msg("loading genetic map (expansion factor ", cfg$map.expansion.factor, ")")
gmap <- restructureGeneticMap(gmap.file, expansion.factor = cfg$map.expansion.factor,
                              expandX = cfg$expandX)

msg("loading VCF (", basename(vcf.file), ") -- this is the memory-hungry step")
vcf <- read_serialized(vcf.file)
msg("loading GT matrix (", basename(gt.file), ")")
gt  <- read_serialized(gt.file)

## -- per-sample allele frequency deviations ----------------------------------

msg("computing per-sample count tables and allele-frequency deviations")
sa       <- compute_sample_afds(sheet, aser.paths, vcf, gt, gmap, cfg, cache.dir,
                                force = isTRUE(opts$force))
countdfs <- sa$countdfs
afds     <- sa$afds

rm(vcf, gt); invisible(gc())

if (isTRUE(cfg$save.rda)) {
  tosave <- list(afds, countdfs)
  save(tosave, file = file.path(out.dir, "af_counts.rda"))
  msg("wrote ", file.path(out.dir, "af_counts.rda"))
}

## -- allele frequency plots --------------------------------------------------

msg("plotting allele frequencies, one panel per sample")
plots <- lapply(names(afds), function(lab) plotIndividualExperiment(afds[[lab]], lab))
names(plots) <- names(afds)
for (lab in names(plots)) {
  f <- file.path(plot.dir, paste0(lab, "_", cfg$bin.width, "_af.png"))
  ggplot2::ggsave(f, plots[[lab]], width = cfg$individual.plot.width)
  msg("  wrote plots/", basename(f))
}

## -- contrasts ---------------------------------------------------------------

if (nrow(contrasts)) {
  msg("running contrasts")
  invisible(run_contrasts(contrasts, sheet, afds, plots, cfg, plot.dir))
} else {
  msg("no contrasts requested; contrast_to is empty for every sample")
}

msg("done. outputs in ", out.dir)
