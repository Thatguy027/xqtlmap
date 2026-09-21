# xqtlmap

Allele-frequency mapping for bulk-segregant xQTL experiments in *C. elegans*.

Give it a directory of GATK **ASEReadCounter** tables and a **metadata table**. It phases
every site to the two parents of each cross, estimates the allele-frequency deviation
along each chromosome, and writes a plot per sample plus every contrast the metadata asks
for — LOD trace, support intervals, effect sizes.

```
xqtlmap --aser aser/ --metadata metadata.tsv --out results/
```

---

## Install

R ≥ 4.2, and:

```r
install.packages(c("tidyverse", "vcfR", "AlphaSimR", "Rfast", "ggpubr", "yaml", "qs2"))
remotes::install_github("joshsbloom/xQTLStats")
```

Then clone and put the command on your PATH:

```bash
git clone https://github.com/Thatguy027/xqtlmap.git
ln -s "$PWD/xqtlmap/bin/xqtlmap" ~/bin/xqtlmap
```

`xqtlmap` checks for every package at startup and names the ones you are missing rather
than failing at the first `library()` call.

### Reference data

Three files are needed and are **not** in this repository — they are large and licensed
separately. They ship with [xQTLSims](https://github.com/joshsbloom/xQTLSims):

| file | what it is |
| --- | --- |
| `data/geneticMapXQTLsnplist.rds` | genetic map, measured on an F10 AIL |
| `data/WI.20220216.vcf.qs2` | filtered CaeNDR isotype VCF |
| `data/WI.20220216.vcf.GT.qs2` | numeric genotype matrix for the same sites |

Point `reference.dir` at your clone, in `config/default.yml` or with
`--reference.dir=/path/to/xQTLSims`. `xqtlmap` also sources a handful of helper functions
from that clone's `R/` directory, so it needs the repository, not just the data files.

---

## Input 1 — the ASER directory

One table per sample, from `gatk ASEReadCounter`. Files are matched **by `sample_name`**,
with or without an extension, gzipped or not: `S1.table`, `S1.table.gz`, `S1.tsv`,
`S1.txt.gz` all resolve for a sample called `S1`. Each needs the columns
`contig`, `position`, `refCount`, `altCount`.

## Input 2 — the metadata table

TSV or CSV — the delimiter is sniffed from the header line, so a comma-separated file
named `.tsv` works. Five columns are required:

| column | meaning |
| --- | --- |
| `sample_name` | matches the ASER file, and must be unique |
| `parent1` | strain the allele frequency is phased **to** |
| `parent2` | the other parent of this cross |
| `condition` | treatment, genotype, timepoint label — whatever this pool is |
| `contrast_to` | comma- or semicolon-separated `sample_name`s to contrast this sample **against**. Leave blank for none. |

Any other column you add becomes part of the sample's label, in the order the columns
appear, so `timepoint` and `stage` below produce `S3_XZ1516_ECA191_1_L1_let-363`. That
label names the cache entry and every output file.

```tsv
sample_name	parent1	parent2	timepoint	stage	condition	contrast_to	n_worms
S1	XZ1516	ECA191	1	L1	HT115		100000
S2	XZ1516	ECA191	1	L1	pos-1		15000
S3	XZ1516	ECA191	1	L1	let-363	S1,S2	25000
S4	XZ1516	ECA191	1	L1	sbp-1	S1,S2	31000
```

**`contrast_to` is directional.** `contrast_to = S1` on row `S3` produces
`afd(S3) − afd(S1)`, so a positive effect means the `parent1` allele is at higher
frequency in **S3**. Contrast each knockdown against its control, not the reverse.

**`n_worms` is special.** It is the one extra column that does *not* enter the label. It
sets `sample.size` per sample, which caps the effective *n* behind every standard error:

```
n      = min(n_worms × sel.strength, total reads / eff.length)
afd.se = sqrt(afd × (1 − afd) / n)
```

Leave it out entirely and every pool is assumed to be `sample.size` from the config. Leave
a single cell blank and that row falls back to the config value with a message. Use it
whenever pools differ in size by more than about twofold — without it a small pool gets a
standard error far too tight and dominates anything it is pooled with.

See `examples/metadata.tsv`.

---

## Output

```
results/
├── metadata_resolved.tsv           the table as parsed, with labels
├── cache/<label>.rds               per-sample counts + AFD, reused across runs
└── plots/
    ├── <label>_10000_af.png        allele frequency along each chromosome
    ├── <cross>_contrast_<A>-<B>_10000.png        LOD trace + effect size
    ├── <cross>_contrast_<A>-<B>_10000_intervals.tsv
    └── <cross>_contrast_<A>-<B>_10000_plot_DF.tsv
```

The interval table carries one row per chromosome clearing the threshold:
`peak.position`, `peak.LOD`, `lcon`/`rcon` (support interval), `width.kb`.

The cache is keyed on the label **and** on the parameters that would change the answer
(`bin.width`, `sample.size`, `sel.strength`, `eff.length`, `map.expansion.factor`,
`uchr`, `X.drop`, `expandX`) plus the mtime of the ASER file. Edit a condition name and
only that sample recomputes; change `bin.width` and everything does. `--force` ignores it.

---

## Configuration

Precedence: **command line > `xqtlmap.yml` beside your metadata > `config/default.yml`**.

```bash
xqtlmap --aser aser/ --metadata metadata.tsv --out results/ --bin.width=5000
```

An `xqtlmap.yml` next to the metadata table is picked up automatically, which is the right
place for anything true of one dataset rather than of the tool:

```yaml
map.expansion.factor: 0.186   # these pools are AIL generation 9
```

### `map.expansion.factor` is about the map, not your pools

`restructureGeneticMap` multiplies the shipped map by this factor to normalise it to F2.
That map was measured on an **F10 AIL** and runs ~5.3× the F2 map, which is where the
default 0.2 comes from. So it encodes the *map's* generation.

A population at a different AIL generation accumulated a different amount of recombination,
and the same physical interval spans **fewer** cM at a lower generation — the factor
contracts, it does not expand. Scaling by the ratio of AIL expansions over F2
(Darvasi & Soller 1995, expansion = (t+4)/6):

| pools | factor |
| --- | --- |
| F8 | 0.171 |
| F9 | 0.186 |
| F10 | 0.200 (default) |
| F12 | 0.229 |

Use the generation of the **pools you sequenced**, which is not always the number on the
tube: if cross construction stopped at generation 8 and treatment ran two more generations,
the pools are generation 9 and 10.

### LOD convention

`xQTLStats` forms `p = 2*pnorm(|z|)` on the linear scale, which underflows to zero around
|z| ≈ 38.5 and silently turns LOD into `Inf` — taking the top chromosome out of every
interval table. `xqtlmap` recomputes LOD in log space. `lod.convention: package` is
bit-identical to `xQTLStats::PvalToLOD` minus the underflow; `chisq` is the textbook
1-df LOD, `z²/(2 ln 10)`, which is the other exactly `log10(2)` higher.

---

## Checking before you run

```bash
xqtlmap --aser aser/ --metadata metadata.tsv --out results/ --check
```

Parses the metadata, resolves every ASER file, expands every contrast and prints the plan
and the output names — then stops, before the VCF is loaded. Nearly every real failure is
a metadata problem, and this catches them in a second rather than twenty minutes in.

## Tests

```bash
Rscript tests/test_metadata.R
```

16 tests over parsing, labelling, contrast expansion, file discovery and the LOD
recomputation. They need no reference data.

---

## Notes and limits

- **Both parents must be in the VCF**, by the names in `parent1` and `parent2`. They are
  CaeNDR isotype names, so check spelling against the release you are using.
- **A contrast between samples from different crosses** is allowed but warned about; the
  allele frequencies are phased to different parents and the difference is rarely what you
  want.
- **Support intervals are LOD-drop intervals**, not confidence intervals. Treat them as
  regions of interest.
- **`afd.se` carries no pool-to-pool variance term** — it is binomial in the effective *n*
  above. Pooling several samples per arm with a fixed-effects meta-analysis therefore
  inflates LOD substantially. If you need to pool, use random effects.

## Credits

Built on [xQTLStats](https://github.com/joshsbloom/xQTLStats) and
[xQTLSims](https://github.com/joshsbloom/xQTLSims) by Josh Bloom. Genetic map and variant
data from [CaeNDR](https://caendr.org).

MIT licensed.
