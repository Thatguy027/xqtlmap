#!/usr/bin/env Rscript
## Tests for metadata parsing, labelling and contrast expansion.
##
##   Rscript tests/test_metadata.R
##
## These exercise everything that can go wrong before the VCF is loaded, which
## is where nearly every real failure has been. They need no reference data.
## ---------------------------------------------------------------------------

PKG_DIR <- dirname(dirname(normalizePath(
  sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))][1]))))
`%||%` <- function(a, b) if (is.null(a)) b else a
suppressWarnings(suppressMessages(source(file.path(PKG_DIR, "R", "pipeline.R"))))

PASS <- 0L; FAIL <- 0L
ok <- function(what, expr) {
  got <- tryCatch({ force(expr); TRUE }, error = function(e) conditionMessage(e))
  if (isTRUE(got)) { PASS <<- PASS + 1L; cat("  ok   ", what, "\n") }
  else { FAIL <<- FAIL + 1L; cat("  FAIL ", what, "\n         ", got, "\n") }
}
errors <- function(what, pattern, expr) {
  m <- tryCatch({ force(expr); NA_character_ }, error = function(e) conditionMessage(e))
  if (!is.na(m) && grepl(pattern, m)) { PASS <<- PASS + 1L; cat("  ok   ", what, "\n") }
  else { FAIL <<- FAIL + 1L
         cat("  FAIL ", what, "\n          expected /", pattern, "/, got: ",
             if (is.na(m)) "no error" else m, "\n", sep = "") }
}
tmp <- function(txt) { f <- tempfile(fileext = ".tsv"); writeLines(txt, f); f }

HDR <- "sample_name\tparent1\tparent2\tcondition\tcontrast_to"

cat("\nmetadata parsing\n")
ok("reads the shipped example", {
  s <- read_sample_sheet(file.path(PKG_DIR, "examples", "metadata.tsv"))
  stopifnot(nrow(s) == 7, "n_worms" %in% names(s))
})
ok("reads a comma-separated table named .tsv", {
  f <- tmp(c("sample_name,parent1,parent2,condition,contrast_to",
             "A,P1,P2,ctl,", "B,P1,P2,trt,A"))
  s <- read_sample_sheet(f); stopifnot(nrow(s) == 2, s$condition[2] == "trt")
})
errors("missing required column", "missing required column",
       read_sample_sheet(tmp(c("sample_name\tparent1\tcondition", "A\tP1\tctl"))))
errors("blank required field", "have an empty 'parent2'",
       read_sample_sheet(tmp(c(HDR, "A\tP1\t\tctl\t"))))
errors("duplicate sample_name", "duplicated sample_name",
       read_sample_sheet(tmp(c(HDR, "A\tP1\tP2\tctl\t", "A\tP1\tP2\ttrt\t"))))

cat("\nlabels\n")
ok("extra columns join the label, in sheet order", {
  f <- tmp(c("sample_name\tparent1\tparent2\ttimepoint\tcondition\tcontrast_to",
             "A\tP1\tP2\t1\tctl\t"))
  stopifnot(read_sample_sheet(f)$label == "A_P1_P2_1_ctl")
})
ok("n_worms is used but stays out of the label", {
  f <- tmp(c("sample_name\tparent1\tparent2\tcondition\tcontrast_to\tn_worms",
             "A\tP1\tP2\tctl\t\t5000"))
  s <- read_sample_sheet(f)
  stopifnot(s$label == "A_P1_P2_ctl", s$n_worms == "5000",
            !"n_worms" %in% attr(s, "extra_cols"))
})

cat("\ncontrasts\n")
ok("contrast_to expands to one row per target", {
  f <- tmp(c(HDR, "A\tP1\tP2\tctl\t", "B\tP1\tP2\tpos\t", "C\tP1\tP2\ttrt\tA,B"))
  ct <- build_contrast_table(read_sample_sheet(f))
  stopifnot(nrow(ct) == 2, all(ct$L_sample == "C"), setequal(ct$R_sample, c("A", "B")))
})
ok("semicolons and stray spaces separate targets too", {
  f <- tmp(c(HDR, "A\tP1\tP2\tctl\t", "B\tP1\tP2\tpos\t", "C\tP1\tP2\ttrt\tA; B"))
  stopifnot(nrow(build_contrast_table(read_sample_sheet(f))) == 2)
})
ok("empty contrast_to yields no contrasts", {
  f <- tmp(c(HDR, "A\tP1\tP2\tctl\t", "B\tP1\tP2\ttrt\t"))
  stopifnot(nrow(build_contrast_table(read_sample_sheet(f))) == 0)
})
errors("unknown contrast target", "references unknown sample",
       build_contrast_table(read_sample_sheet(tmp(c(HDR, "A\tP1\tP2\ttrt\tZ")))))
errors("self-contrast", "lists itself",
       build_contrast_table(read_sample_sheet(tmp(c(HDR, "A\tP1\tP2\ttrt\tA")))))

cat("\naser discovery\n")
ok("finds tables by sample name across extensions", {
  d <- file.path(tempdir(), paste0("aser", as.integer(runif(1, 1e5, 1e6))))
  dir.create(d); file.create(file.path(d, c("A.table", "B.tsv.gz")))
  p <- find_aser_files(d, c("A", "B"))
  stopifnot(basename(p[["A"]]) == "A.table", basename(p[["B"]]) == "B.tsv.gz")
})
errors("missing aser table names the sample", "no ASEReadCounter table for sample 'B'", {
  d <- file.path(tempdir(), paste0("aser", as.integer(runif(1, 1e5, 1e6))))
  dir.create(d); file.create(file.path(d, "A.table"))
  find_aser_files(d, c("A", "B"))
})

cat("\nlod\n")
ok("lod_threshold rises with the number of tests", {
  stopifnot(lod_threshold(0.05, 20000) > lod_threshold(0.05, 2000))
})
ok("recompute_lod survives |z| past the underflow point", {
  r <- data.frame(z = c(1, 40, 80))
  l <- recompute_lod(r, "package")$LOD
  stopifnot(all(is.finite(l)), l[3] > l[2], l[2] > l[1])
})

cat(sprintf("\n%d passed, %d failed\n", PASS, FAIL))
quit(status = if (FAIL) 1L else 0L)
