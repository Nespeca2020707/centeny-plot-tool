#!/usr/bin/env Rscript
# =============================================================================
#  motif_report.R  —  Motif Model Composition Report
#  Version: 2.1
# =============================================================================
#
#  PURPOSE
#  -------
#  From --sample_sheet's centeny/model1/model2/model3 BED data, produces
#  TWO outputs in the output directory:
#
#    1. whole_genome_motif_abundance.html
#       Per-model composition: one stacked bar PER SAMPLE (haplotypes
#       shown separately, not merged), each summing to 100%, broken down
#       by that model's own categories (strand for centeny; distance/
#       monomer-count bins for model1/model2; top-N k-patterns for
#       model3). Hover a segment for its exact percentage; a click-to-
#       reveal table underneath each chart gives the exact counts too.
#
#    2. <output>_karyotypes.html
#       Per model: a karyotype-style grid (rows = genome — haplotypes of
#       the same individual/assembly are shown TOGETHER here, split at
#       the first "_" in the sample label, as up to N side-by-side bars
#       per chromosome column; columns = chromosome). Hover a segment for
#       its exact percentage; a click-to-reveal table gives exact counts.
#       Also includes, per model, a collapsible Euclidean/average-linkage
#       dendrogram+heatmap clustering individual SAMPLES (haplotypes kept
#       separate here too) by how similar their category percentages are.
#
#  NOT IN THIS RELEASE: CENP-B_conservation_profiles.html — the FASTA-
#  scanning-based per-reference top-motif tables, nucleotide-conservation
#  heatmap+dendrogram, and cross-genome sequence-similarity clustering
#  (which used find_motif_profile.py) are deferred to a later release.
#  This script currently only needs --sample_sheet; no FASTA/genome_sheet
#  input is read.
#
#  USAGE
#  -----
#    Rscript motif_report.R \
#      --sample_sheet samples.tsv \
#      --output       motif_report
#
#  SAMPLE SHEET
#  ------------
#  Same format the main pipeline uses — tab-separated, no header, 2 or 3
#  columns: <sample_label> <bed_path> [<model>]. <model> is optional:
#  centeny | model1 | model2 | model3 | auto (default — detected from
#  column 4). Genome/haplotype grouping for the karyotype grid is derived
#  automatically from the label (split at the first "_", e.g. "RPE1_hap1"
#  -> genome "RPE1", haplotype "hap1").
# =============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
})

# =============================================================================
#  1.  COMMAND-LINE ARGUMENT PARSING
# =============================================================================

option_list <- list(
  make_option(c("--workdir"), type = "character", default = getwd(),
    help = "Working directory [default: current directory]"),

  make_option(c("--sample_sheet"), type = "character", default = NULL,
    help = "REQUIRED. Same sample sheet format as the main pipeline."),

  make_option(c("--output"), type = "character", default = "motif_report",
    help = "Output base name [default: motif_report]"),

  make_option(c("--composition_threshold"), type = "double", default = 0.5,
    help = paste0("For the model1/model2 composition charts, categories ",
                  "below this percentage (of that model's total features, ",
                  "globally) are folded into an 'Other' bucket, matching ",
                  "typical practice for readability [default: 0.5]")),

  make_option(c("--kpattern_top_n"), type = "integer", default = 10,
    help = paste0("For model3's composition chart/karyotype grid, the N ",
                  "most abundant k-patterns (globally) get their own ",
                  "category; everything else is folded into 'Other' ",
                  "[default: 10]"))
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$sample_sheet)) stop("--sample_sheet is required.", call. = FALSE)

setwd(opt$workdir)

output_base <- tools::file_path_sans_ext(opt$output)
out_dir <- dirname(output_base)
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
#  2.  SMALL SHARED HELPERS (self-contained — this script doesn't source
#      centeny_plot_tool.R, to keep the two tools independently runnable)
# =============================================================================

clean_chr <- function(x) {
  x <- gsub("^chr", "", x, ignore.case = TRUE)
  x <- gsub("_hap[0-9]+$|_(pat|mat)(ernal)?$", "", x, ignore.case = TRUE)
  x
}

natural_sort_chroms <- function(chr_vec) {
  nm <- unique(chr_vec)
  is_num <- grepl("^[0-9]+$", nm)
  c(nm[is_num][order(as.integer(nm[is_num]))], sort(nm[!is_num]))
}

# Mirrors centeny_plot_tool.R's model auto-detection (column 4 prefix),
# duplicated here in miniature since only the model LABEL is needed, not
# any colour/legend logic.
detect_model_label <- function(v4) {
  v4 <- as.character(v4)
  if (any(grepl("^#", v4)))                                return("model3")
  if (any(grepl("^D[0-9]+$", v4)) || any(v4 == "DNA"))      return("model1")
  if (any(grepl("^M[0-9]+$", v4)))                          return("model2")
  "centeny"
}

resolve_model_tag_light <- function(tag) {
  tag <- tolower(trimws(tag))
  if (tag %in% c("", "auto"))            return("auto")
  if (tag %in% c("centeny", "plain"))    return("centeny")
  if (tag %in% c("model1", "1"))         return("model1")
  if (tag %in% c("model2", "2"))         return("model2")
  if (tag %in% c("model3", "3"))         return("model3")
  "auto"
}

html_escape <- function(x) {
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;",  x, fixed = TRUE)
  x <- gsub(">", "&gt;",  x, fixed = TRUE)
  x
}

# JS string-literal escaper (for embedding R character vectors as JS
# array-of-strings literals in generated <script> blocks).
js_str <- function(x) paste0('"', gsub('(["\\\\])', '\\\\\\1', x), '"')

MODEL_ORDER   <- c("centeny", "model1", "model2", "model3")
MODEL_TITLES  <- c(centeny = "CenTeny (strand)", model1 = "Model 1 (distance, bp)",
                    model2 = "Model 2 (distance, monomers)", model3 = "Model 3 (k-pattern)")

# "RPE1_hap1" -> genome group "RPE1", haplotype "1" — used for the
# composition-percentage chart's sample-level rows and the karyotype
# grid's genome/haplotype grouping.
extract_genome_group <- function(label) sub("_.*", "", label)

extract_hap <- function(label) {
  hap <- rep("1", length(label))
  has_hap <- grepl("_hap[0-9]+$", label, ignore.case = TRUE)
  hap[has_hap] <- sub(".*_hap([0-9]+)$", "\\1", label[has_hap], ignore.case = TRUE)
  has_pm  <- grepl("_(pat|mat)(ernal)?$", label, ignore.case = TRUE)
  hap[has_pm] <- tolower(sub(".*_(pat|mat)(ernal)?$", "\\1", label[has_pm], ignore.case = TRUE))
  hap
}

# BED col9 "204,0,0" -> "#CC0000" (uppercase hex, stable dictionary key).
# Malformed/missing values fall back to a neutral grey rather than erroring.
rgb_string_to_hex <- function(v9) {
  parts <- strsplit(as.character(v9), ",", fixed = TRUE)
  vapply(parts, function(p) {
    if (length(p) != 3) return("#808080")
    p <- suppressWarnings(as.integer(trimws(p)))
    p <- suppressWarnings(as.integer(trimws(p)))
    if (any(is.na(p)) || any(p < 0) || any(p > 255)) return("#808080")
    sprintf("#%02X%02X%02X", p[1], p[2], p[3])
  }, character(1))
}

parse_model_value <- function(v4, model) {
  v4 <- as.character(v4)
  if (model == "model1") {
    out <- suppressWarnings(as.numeric(sub("^D", "", v4)))
    out[!grepl("^D[0-9]+$", v4)] <- NA_real_
    return(out)
  }
  if (model == "model2") {
    # M0 in the file means 1 monomer between motifs, M1 means 2, etc.
    out <- suppressWarnings(as.numeric(sub("^M", "", v4)) + 1)
    out[!grepl("^M[0-9]+$", v4)] <- NA_real_
    return(out)
  }
  rep(NA_real_, length(v4))
}

format_value_label <- function(lo, hi, model) {
  if (model == "model2") {
    lo_i <- round(lo); hi_i <- round(hi)
    if (isTRUE(all.equal(lo, hi)))
      return(sprintf("every-%d-monomer%s", lo_i, if (lo_i == 1) "" else "s"))
    return(sprintf("every-%d..%d-monomers", lo_i, hi_i))
  }
  fmt_one <- function(x) {
    if (x >= 1e6) return(sprintf("%.2f Mb", x / 1e6))
    if (x >= 1e3) return(sprintf("%.1f kb", x / 1e3))
    sprintf("%d bp", as.integer(round(x)))
  }
  if (isTRUE(all.equal(lo, hi))) return(fmt_one(lo))
  paste0(fmt_one(lo), "-", fmt_one(hi))
}

# Builds, ONCE PER MODEL and GLOBALLY across every sample, a lookup from
# each row's colour (model1/model2/centeny) or literal k-pattern
# (model3) to a final display category — so the same category always
# gets the same label/colour everywhere in the report (composition chart
# AND karyotype grid), and low-abundance categories are consistently
# folded into "Other" rather than each bar defining its own ad hoc set.
build_category_lookup <- function(rows, model, threshold_pct = 0.5, kpattern_top_n = 10) {
  other_colour <- "#B0B0B0"
  if (nrow(rows) == 0) return(NULL)

  if (model == "centeny") {
    key_to_label <- c("#CC0000" = "Reverse strand", "#0000CC" = "Forward strand")
    present <- unique(ifelse(rows$HEXKEY %in% names(key_to_label),
                              key_to_label[rows$HEXKEY], "Other"))
    categories <- intersect(c("Forward strand", "Reverse strand", "Other"), present)
    colours <- c("Forward strand" = "#0000CC", "Reverse strand" = "#CC0000",
                 "Other" = other_colour)[categories]
    return(list(by = "hex", key_to_label = key_to_label, default_label = "Other",
                categories = categories, colours = colours))
  }

  if (model %in% c("model1", "model2")) {
    vals <- parse_model_value(rows$V4, model)
    dt <- data.table(HEXKEY = rows$HEXKEY, VAL = vals)
    grp <- dt[, .(n = .N,
                   lo = suppressWarnings(min(VAL, na.rm = TRUE)),
                   hi = suppressWarnings(max(VAL, na.rm = TRUE))), by = HEXKEY]
    grp$is_na <- is.infinite(grp$lo)
    grp$label <- ifelse(grp$is_na, "No value",
                         mapply(function(lo, hi) format_value_label(lo, hi, model), grp$lo, grp$hi))
    grp$pct <- grp$n / sum(grp$n) * 100
    other_label <- sprintf("Other (<%.1f%%)", threshold_pct)
    grp$final_label <- ifelse(!grp$is_na & grp$pct < threshold_pct, other_label, grp$label)
    key_to_label <- setNames(grp$final_label, grp$HEXKEY)

    cat_agg <- grp[, .(n = sum(n), rep_val = mean(lo[!is_na], na.rm = TRUE)), by = final_label]
    is_other <- cat_agg$final_label %in% c(other_label, "No value")
    ord <- order(is_other, cat_agg$rep_val)
    categories <- cat_agg$final_label[ord]
    is_other_ord <- is_other[ord]

    n_colour_cats <- sum(!is_other_ord)
    pal <- if (n_colour_cats > 0) grDevices::hcl.colors(n_colour_cats, "Dynamic") else character(0)
    colours <- setNames(rep(other_colour, length(categories)), categories)
    colours[categories[!is_other_ord]] <- pal

    return(list(by = "hex", key_to_label = key_to_label, default_label = other_label,
                categories = categories, colours = colours))
  }

  if (model == "model3") {
    freq <- sort(table(rows$V4), decreasing = TRUE)
    top_pat <- names(freq)[seq_len(min(kpattern_top_n, length(freq)))]
    pal <- grDevices::hcl.colors(max(1, length(top_pat)), "Dynamic")[seq_along(top_pat)]
    key_to_label <- setNames(top_pat, top_pat)
    has_other <- length(freq) > length(top_pat)
    categories <- if (has_other) c(top_pat, "Other") else top_pat
    colours <- setNames(c(pal, if (has_other) other_colour), categories)
    return(list(by = "v4", key_to_label = key_to_label, default_label = "Other",
                categories = categories, colours = colours))
  }

  NULL
}

apply_category_lookup <- function(dt, lookup) {
  key <- if (lookup$by == "hex") dt$HEXKEY else as.character(dt$V4)
  lab <- unname(lookup$key_to_label[key])
  lab[is.na(lab)] <- lookup$default_label
  lab
}


# =============================================================================
#  3.  PARSE SAMPLE SHEET (2-3 columns: label, bed, [model])
# =============================================================================

cat("Reading sample sheet:", opt$sample_sheet, "\n")
raw_sample <- readLines(opt$sample_sheet)
raw_sample <- raw_sample[!grepl("^\\s*#", raw_sample)]
raw_sample <- raw_sample[nchar(trimws(raw_sample)) > 0]

sample_df <- do.call(rbind, lapply(raw_sample, function(line) {
  parts <- strsplit(line, "\t")[[1]]
  if (!length(parts) %in% c(2, 3))
    stop(paste0("Malformed sample sheet line (expected 2-3 tab-separated ",
                "columns): ", line), call. = FALSE)
  label <- trimws(parts[1])
  data.frame(
    label     = label,
    bed       = trimws(parts[2]),
    model_tag = resolve_model_tag_light(if (length(parts) >= 3) parts[3] else "auto"),
    genome_group = extract_genome_group(label),
    hap          = extract_hap(label),
    stringsAsFactors = FALSE
  )
}))

missing_bed <- sample_df$bed[!file.exists(sample_df$bed)]
if (length(missing_bed) > 0)
  stop(paste("Missing BED file(s):\n", paste(missing_bed, collapse = "\n")), call. = FALSE)

cat(sprintf("  %d sample(s) across %d genome(s)\n", nrow(sample_df),
            length(unique(sample_df$genome_group))))

#  4.  LOAD ALL SAMPLE BED DATA ONCE (shared by every section below)
# =============================================================================

cat("\nLoading sample BED data for composition/karyotype analysis...\n")
sample_beds <- vector("list", nrow(sample_df))
for (i in seq_len(nrow(sample_df))) {
  dt <- fread(sample_df$bed[i], header = FALSE)
  dt$V1 <- clean_chr(dt$V1)
  dt$HEXKEY <- rgb_string_to_hex(dt$V9)
  model_label <- if (sample_df$model_tag[i] == "auto") {
    detect_model_label(dt$V4)
  } else {
    sample_df$model_tag[i]
  }
  dt$MODEL        <- model_label
  dt$SAMPLE_LABEL <- sample_df$label[i]
  dt$GENOME_GROUP <- sample_df$genome_group[i]
  dt$HAP          <- sample_df$hap[i]
  sample_beds[[i]] <- dt
}
all_rows <- rbindlist(sample_beds, use.names = TRUE, fill = TRUE)

# =============================================================================
#  5.  PER-MODEL GLOBAL CATEGORY LOOKUP, COMPOSITION %, AND GRID CELLS
# =============================================================================

cat("Classifying features into categories per model...\n")
category_lookups   <- list()
composition_pct    <- list()   # model -> matrix[sample_label, category]  (haplotypes NOT merged)
grid_cells         <- list()   # model -> data.frame(genome_group,hap,chrom,category,n,n_total,pct)
sample_cluster_mat <- list()   # model -> matrix[sample_label, category]  (same as composition_pct; kept
                                # as its own object for clarity at the clustering call site below)

for (model in MODEL_ORDER) {
  rows <- all_rows[all_rows$MODEL == model, ]
  if (nrow(rows) == 0) { cat(sprintf("  [%s] no features found - skipping\n", model)); next }

  lookup <- build_category_lookup(rows, model,
                                   threshold_pct  = opt$composition_threshold,
                                   kpattern_top_n = opt$kpattern_top_n)
  if (is.null(lookup)) next
  category_lookups[[model]] <- lookup
  rows$CATEGORY <- apply_category_lookup(rows, lookup)
  cat(sprintf("  [%s] %d categories: %s\n", model, length(lookup$categories),
              paste(utils::head(lookup$categories, 6), collapse = ", ")))

  # ---- composition % per SAMPLE (haplotypes shown separately, per request) ---
  tab <- table(rows$SAMPLE_LABEL, factor(rows$CATEGORY, levels = lookup$categories))
  pct <- prop.table(tab, margin = 1) * 100
  composition_pct[[model]]    <- pct
  sample_cluster_mat[[model]] <- pct

  # ---- per (genome_group, hap, chromosome) cells (karyotype grid) --------
  cell_n <- rows[, .N, by = .(GENOME_GROUP, HAP, V1, CATEGORY)]
  setnames(cell_n, c("genome_group", "hap", "chrom", "category", "n"))
  cell_tot <- cell_n[, .(n_total = sum(n)), by = .(genome_group, hap, chrom)]
  cell_n <- merge(cell_n, cell_tot, by = c("genome_group", "hap", "chrom"))
  cell_n$pct <- cell_n$n / cell_n$n_total * 100
  grid_cells[[model]] <- as.data.frame(cell_n)
}

# =============================================================================
#  6.  SHARED TOOLTIP (one <div>+CSS+JS block, reused by every page)
# =============================================================================

TOOLTIP_DIV_CSS <- paste0(
  '<div id="globalTooltip" style="position:fixed;pointer-events:none;',
  'background:rgba(20,20,20,0.92);color:#fff;font-size:11px;padding:4px 7px;',
  'border-radius:4px;display:none;z-index:50;white-space:nowrap;"></div>'
)
TOOLTIP_JS <- paste0(
  '<script>\n',
  '(function(){\n',
  '  var tip = document.getElementById("globalTooltip");\n',
  '  document.addEventListener("mousemove", function(e){\n',
  '    var el = e.target.closest ? e.target.closest("[data-tip]") : null;\n',
  '    if (!el) { tip.style.display = "none"; return; }\n',
  '    tip.style.display = "block";\n',
  '    tip.style.left = (e.clientX + 12) + "px";\n',
  '    tip.style.top  = (e.clientY + 12) + "px";\n',
  '    tip.innerHTML = el.getAttribute("data-tip");\n',
  '  });\n',
  '})();\n',
  '</script>'
)

# =============================================================================
#  7.  SVG BUILDERS  (hover = data-tip attribute, read by the shared JS above)
# =============================================================================

# One bar per SAMPLE (haplotypes shown separately — NOT merged), stacked by
# that model's own categories, each bar summing to 100%.
build_composition_pct_svg <- function(pct_matrix, categories, colours,
                                       bar_w = 60, gap = 26, left_pad = 46, plot_h = 380) {
  groups <- rownames(pct_matrix)
  n <- length(groups)
  if (n == 0) return("<p><em>No data.</em></p>")
  top_pad  <- 10
  height   <- top_pad + plot_h + 90
  width    <- left_pad + n * (bar_w + gap) + 40

  axis <- character(0)
  for (yv in seq(0, 100, by = 20)) {
    yy <- top_pad + plot_h - yv / 100 * plot_h
    axis <- c(axis, sprintf('<line x1="%d" y1="%.2f" x2="%.2f" y2="%.2f" stroke="#eee"/>',
                             left_pad - 4, yy, left_pad + n * (bar_w + gap), yy))
    axis <- c(axis, sprintf('<text x="%d" y="%.2f" font-size="10" text-anchor="end">%d</text>',
                             left_pad - 8, yy + 3, yv))
  }

  bars <- character(0)
  for (i in seq_len(n)) {
    g <- groups[i]
    x <- left_pad + (i - 1) * (bar_w + gap)
    y_cursor <- top_pad + plot_h
    for (cat in categories) {
      v <- pct_matrix[g, cat]
      if (is.na(v) || v <= 0) next
      h <- v / 100 * plot_h
      y_cursor <- y_cursor - h
      tip <- html_escape(sprintf("%s<br>%s: %.2f%%", g, cat, v))
      bars <- c(bars, sprintf(
        '<rect x="%.2f" y="%.2f" width="%d" height="%.2f" fill="%s" data-tip="%s"></rect>',
        x, y_cursor, bar_w, h, colours[[cat]], tip))
    }
    bars <- c(bars, sprintf(
      '<text x="%.2f" y="%d" font-size="11" text-anchor="end" transform="rotate(-35 %.2f %d)">%s</text>',
      x + bar_w / 2, top_pad + plot_h + 18, x + bar_w / 2, top_pad + plot_h + 18, html_escape(g)))
  }

  # Explicit pixel width/height (not scaled to fit the page) — with many
  # samples this needs to be genuinely wide; the caller wraps it in a
  # horizontally-scrollable container so bar/text size stays constant and
  # readable no matter how many samples there are.
  sprintf('<svg width="%.0f" height="%.0f" viewBox="0 0 %.0f %.0f" xmlns="http://www.w3.org/2000/svg" font-family="-apple-system,Helvetica,Arial,sans-serif">%s%s</svg>',
          width, height, width, height, paste(axis, collapse = ""), paste(bars, collapse = ""))
}

# Grid: rows = genome (haplotypes combined under one row, split at the first
# "_" in the sample label), columns = chromosome (natural sort, union across
# genomes), up to N thin side-by-side bars per cell (one per haplotype),
# each itself stacked to 100% by category. Hover via data-tip.
build_karyotype_grid_svg <- function(cell_df, categories, colours,
                                      row_h = 92, col_w = 42, bar_w = 15, bar_gap = 3,
                                      left_pad = 120, top_pad = 46) {
  if (nrow(cell_df) == 0) return("<p><em>No data.</em></p>")
  genome_groups <- unique(cell_df$genome_group)
  chrom_order   <- natural_sort_chroms(cell_df$chrom)
  n_rows <- length(genome_groups); n_cols <- length(chrom_order)
  width  <- left_pad + n_cols * col_w + 30
  height <- top_pad + n_rows * row_h + 20

  elems <- character(0)
  for (j in seq_along(chrom_order)) {
    cx <- left_pad + (j - 1) * col_w + col_w / 2
    elems <- c(elems, sprintf('<text x="%.1f" y="%d" font-size="10" text-anchor="middle">chr%s</text>',
                               cx, top_pad - 10, html_escape(chrom_order[j])))
  }

  for (i in seq_along(genome_groups)) {
    g  <- genome_groups[i]
    ry <- top_pad + (i - 1) * row_h
    elems <- c(elems, sprintf('<text x="%d" y="%.1f" font-size="11" font-weight="600" text-anchor="end">%s</text>',
                               left_pad - 10, ry + row_h / 2, html_escape(g)))
    for (j in seq_along(chrom_order)) {
      chrom <- chrom_order[j]
      sub <- cell_df[cell_df$genome_group == g & cell_df$chrom == chrom, ]
      haps <- sort(unique(sub$hap))
      if (length(haps) == 0) next
      grp_w <- length(haps) * bar_w + (length(haps) - 1) * bar_gap
      cx0 <- left_pad + (j - 1) * col_w + (col_w - grp_w) / 2
      for (k in seq_along(haps)) {
        h <- haps[k]
        bx <- cx0 + (k - 1) * (bar_w + bar_gap)
        hap_rows <- sub[sub$hap == h, ]
        y_cursor <- ry + row_h - 6
        for (cat in categories) {
          v <- hap_rows$pct[hap_rows$category == cat]
          v <- if (length(v) == 0) 0 else v
          if (v <= 0) next
          seg_h <- v / 100 * (row_h - 12)
          y_cursor <- y_cursor - seg_h
          tip <- html_escape(sprintf("%s chr%s hap%s<br>%s: %.2f%%", g, chrom, h, cat, v))
          elems <- c(elems, sprintf(
            '<rect x="%.2f" y="%.2f" width="%d" height="%.2f" fill="%s" data-tip="%s"></rect>',
            bx, y_cursor, bar_w, seg_h, colours[[cat]], tip))
        }
      }
    }
  }
  sprintf('<svg width="%.0f" height="%.0f" viewBox="0 0 %.0f %.0f" xmlns="http://www.w3.org/2000/svg" font-family="-apple-system,Helvetica,Arial,sans-serif">%s</svg>',
          width, height, width, height, paste(elems, collapse = ""))
}

build_legend_html <- function(categories, colours) {
  items <- vapply(categories, function(cat) sprintf(
    '<span style="display:inline-flex;align-items:center;gap:4px;margin-right:14px;font-size:11px;">
       <span style="width:10px;height:10px;background:%s;display:inline-block;border:1px solid #0002;"></span>%s</span>',
    colours[[cat]], html_escape(cat)), character(1))
  paste0('<div style="margin:4px 0 10px 0;">', paste(items, collapse = ""), "</div>")
}

# Exact counts/percentages behind a composition/grid figure, as a plain
# collapsible table (click-to-reveal, matching the plot it documents).
build_counts_table_details <- function(df, cols, col_labels, summary_text) {
  rows <- apply(df[, cols, drop = FALSE], 1, function(r) {
    paste0("<tr>", paste0("<td>", html_escape(as.character(r)), "</td>", collapse = ""), "</tr>")
  })
  paste0(
    '<details><summary>', summary_text, '</summary>',
    '<div style="max-height:420px;overflow:auto;border:1px solid #eee;margin-top:6px;">',
    '<table class="motiftbl"><thead><tr>',
    paste0("<th>", html_escape(col_labels), "</th>", collapse = ""),
    "</tr></thead><tbody>", paste(rows, collapse = ""), "</tbody></table></div></details>"
  )
}

# =============================================================================
#  8.  DENDROGRAM + HEATMAP  (Euclidean distance, average-linkage hclust)
#       Used for the per-model sample-composition clustering shown on the
#       karyotypes page. White = 0; colour ramps to dark navy at the row's
#       own max the same way as the main centeny_plot_tool.R clustering
#       feature, for visual consistency across the whole toolchain.
# =============================================================================

compute_dendrogram_segments <- function(hc) {
  n <- length(hc$order)
  leaf_pos <- setNames(seq_len(n), as.character(hc$order))
  node_pos <- numeric(nrow(hc$merge))
  get_pos <- function(idx) if (idx < 0) unname(leaf_pos[as.character(-idx)]) else node_pos[idx]
  get_ht  <- function(idx) if (idx < 0) 0 else hc$height[idx]
  segs <- vector("list", nrow(hc$merge) * 3)
  k <- 0L
  for (i in seq_len(nrow(hc$merge))) {
    a <- hc$merge[i, 1]; b <- hc$merge[i, 2]
    pa <- get_pos(a); pb <- get_pos(b)
    ha <- get_ht(a);  hb <- get_ht(b)
    hm <- hc$height[i]
    pm <- (pa + pb) / 2
    node_pos[i] <- pm
    k <- k + 1L; segs[[k]] <- c(pa, ha, pa, hm)
    k <- k + 1L; segs[[k]] <- c(pb, hb, pb, hm)
    k <- k + 1L; segs[[k]] <- c(pa, hm, pb, hm)
  }
  list(segments = segs, order = hc$order, max_height = max(hc$height), n = n)
}

# Returns list(html=<div+canvas+script>, hc=<hclust object>) or NULL if
# there isn't enough data to cluster (fewer than 3 rows).
build_dendro_heatmap_widget <- function(mat, id_prefix, value_fmt = "%.2f",
                                         dist_method = "euclidean", hclust_method = "average") {
  if (is.null(mat) || nrow(mat) < 3) return(NULL)
  mat[is.na(mat)] <- 0
  d  <- dist(mat, method = dist_method)
  hc <- hclust(d, method = hclust_method)
  dendro  <- compute_dendrogram_segments(hc)
  ordered <- mat[hc$order, , drop = FALSE]
  n_rows <- nrow(ordered); n_cols <- ncol(ordered)
  mx <- max(ordered, na.rm = TRUE); if (!is.finite(mx) || mx <= 0) mx <- 1

  row_labels_js <- paste0("[", paste(js_str(rownames(ordered)), collapse = ","), "]")
  col_labels_js <- paste0("[", paste(js_str(colnames(ordered)), collapse = ","), "]")
  values_js <- paste0("[", paste(sprintf("%.5f", as.vector(t(ordered))), collapse = ","), "]")
  seg_js <- paste0("[", paste(vapply(dendro$segments, function(s)
                    sprintf("[%.4f,%.4f,%.4f,%.4f]", s[1], s[2], s[3], s[4]),
                    character(1)), collapse = ","), "]")

  row_h  <- max(15, min(24, round(3000 / n_rows)))
  cell_w <- max(10, min(40, round(900 / max(1, n_cols))))
  dendro_w <- 150

  html <- paste0(
    '<div id="', id_prefix, '_wrap" style="overflow:auto;border:1px solid #eee;">',
    '<canvas id="', id_prefix, '_cv"></canvas></div>',
    '<script>\n(function(){\n',
    '"use strict";\n',
    'var ROWS=', row_labels_js, ';\n',
    'var COLS=', col_labels_js, ';\n',
    'var VALS=', values_js, ';\n',
    'var SEGS=', seg_js, ';\n',
    'var MAXH=', dendro$max_height, ';\n',
    'var MAXV=', mx, ';\n',
    'var ROW_H=', row_h, ', CELL_W=', cell_w, ', DENDRO_W=', dendro_w, ';\n',
    'var LABEL_W=200, HEADER_H=170;\n',
    'var nRows=ROWS.length, nCols=COLS.length;\n',
    'var canvas=document.getElementById("', id_prefix, '_cv");\n',
    'var ctx=canvas.getContext("2d");\n',
    'var totalW=LABEL_W+DENDRO_W+nCols*CELL_W+20;\n',
    'var totalH=HEADER_H+nRows*ROW_H+10;\n',
    'var dpr=window.devicePixelRatio||1;\n',
    'canvas.style.width=totalW+"px"; canvas.style.height=totalH+"px";\n',
    'canvas.width=Math.round(totalW*dpr); canvas.height=Math.round(totalH*dpr);\n',
    'ctx.setTransform(dpr,0,0,dpr,0,0);\n',
    'function val(r,c){ return VALS[r*nCols+c]; }\n',
    'function colour(v){\n',
    '  if (v<=0) return "#FFFFFF";\n',
    '  var t=Math.max(0,Math.min(1, v/MAXV));\n',
    '  var stops=[[214,234,248],[76,114,176],[8,48,107]];\n',
    '  var seg=t<0.5?0:1; var lt=t<0.5?t/0.5:(t-0.5)/0.5;\n',
    '  var a=stops[seg], b=stops[seg+1];\n',
    '  var r=Math.round(a[0]+(b[0]-a[0])*lt), g=Math.round(a[1]+(b[1]-a[1])*lt), bch=Math.round(a[2]+(b[2]-a[2])*lt);\n',
    '  return "rgb("+r+","+g+","+bch+")";\n',
    '}\n',
    'function heightToX(h){ return DENDRO_W-(h/(MAXH||1))*(DENDRO_W-8); }\n',
    'function posToY(p){ return HEADER_H+(p-0.5)*ROW_H; }\n',
    'function draw(){\n',
    '  ctx.clearRect(0,0,totalW,totalH);\n',
    '  ctx.strokeStyle="#666"; ctx.lineWidth=1;\n',
    '  SEGS.forEach(function(s){\n',
    '    ctx.beginPath();\n',
    '    ctx.moveTo(LABEL_W+heightToX(s[1]), posToY(s[0]));\n',
    '    ctx.lineTo(LABEL_W+heightToX(s[3]), posToY(s[2]));\n',
    '    ctx.stroke();\n',
    '  });\n',
    '  ctx.font="9px sans-serif"; ctx.textAlign="left"; ctx.textBaseline="middle"; ctx.fillStyle="#222";\n',
    '  for (var r=0;r<nRows;r++){\n',
    '    var label=ROWS[r];\n',
    '    if (ctx.measureText(label).width>LABEL_W-8){\n',
    '      while (label.length>1 && ctx.measureText(label+"\\u2026").width>LABEL_W-8) label=label.slice(0,-1);\n',
    '      label+="\\u2026";\n',
    '    }\n',
    '    ctx.fillText(label, 2, HEADER_H+r*ROW_H+ROW_H/2);\n',
    '  }\n',
    '  ctx.save();\n',
    '  var gx0=LABEL_W+DENDRO_W;\n',
    '  for (var c=0;c<nCols;c++){\n',
    '    var cx=gx0+c*CELL_W+CELL_W/2;\n',
    '    ctx.save(); ctx.translate(cx,HEADER_H-6); ctx.rotate(-Math.PI/3);\n',
    '    ctx.textAlign="left"; ctx.textBaseline="middle";\n',
    '    var label=COLS[c]; if (label.length>26) label=label.slice(0,25)+"\\u2026";\n',
    '    ctx.fillText(label,0,0);\n',
    '    ctx.restore();\n',
    '  }\n',
    '  ctx.restore();\n',
    '  for (var r=0;r<nRows;r++){\n',
    '    for (var c=0;c<nCols;c++){\n',
    '      ctx.fillStyle=colour(val(r,c));\n',
    '      ctx.fillRect(gx0+c*CELL_W, HEADER_H+r*ROW_H, CELL_W-1, ROW_H-2);\n',
    '    }\n',
    '  }\n',
    '}\n',
    'draw();\n',
    'var tip = document.getElementById("globalTooltip");\n',
    'canvas.addEventListener("mousemove", function(e){\n',
    '  var rect=canvas.getBoundingClientRect();\n',
    '  var mx=e.clientX-rect.left, my=e.clientY-rect.top;\n',
    '  var gx0=LABEL_W+DENDRO_W;\n',
    '  if (mx<gx0||my<HEADER_H){ if (tip) tip.style.display="none"; return; }\n',
    '  var c=Math.floor((mx-gx0)/CELL_W), r=Math.floor((my-HEADER_H)/ROW_H);\n',
    '  if (c<0||c>=nCols||r<0||r>=nRows){ if (tip) tip.style.display="none"; return; }\n',
    '  if (tip){\n',
    '    tip.style.display="block";\n',
    '    tip.style.left=(e.clientX+12)+"px"; tip.style.top=(e.clientY+12)+"px";\n',
    '    tip.innerHTML = ROWS[r] + "<br>" + COLS[c] + "<br>" + val(r,c).toFixed(2);\n',
    '  }\n',
    '});\n',
    'canvas.addEventListener("mouseleave", function(){ if (tip) tip.style.display="none"; });\n',
    '})();\n</script>'
  )
  list(html = html, hc = hc)
}

# =============================================================================
#  9.  ASSEMBLE THE TWO OUTPUT HTML FILES
# =============================================================================

cat("\nBuilding HTML reports...\n")

REPORT_CSS <- paste0(
  'html,body{margin:0;padding:0;font-family:-apple-system,Helvetica,Arial,sans-serif;color:#222;}',
  '#page{max-width:1100px;margin:0 auto;padding:16px 20px 60px 20px;}',
  'h1{font-size:20px;} h2{font-size:17px;margin-top:28px;border-bottom:2px solid #eee;padding-bottom:4px;}',
  'h3{font-size:13px;color:#555;margin:18px 0 6px 0;}',
  '#toc{position:sticky;top:0;background:#fff;padding:8px 0;border-bottom:1px solid #ddd;font-size:13px;z-index:5;}',
  '#toc a{margin-right:4px;}',
  'table.motiftbl{border-collapse:collapse;width:100%;font-size:12px;}',
  'table.motiftbl th,table.motiftbl td{border-bottom:1px solid #eee;padding:4px 8px;text-align:left;}',
  'table.motiftbl thead th{position:sticky;top:0;background:#fafafa;}',
  'details summary{cursor:pointer;font-size:13px;color:#345;padding:6px 0;}',
  'code{font-size:11px;background:#f6f6f6;padding:1px 4px;border-radius:3px;}',
  'hr{border:none;border-top:1px solid #eee;margin:24px 0;}'
)

# ---- File 1: whole_genome_motif_abundance.html ---------------------------
# Per-model composition %, one bar per SAMPLE (haplotypes shown separately).

abundance_sections <- vapply(MODEL_ORDER, function(model) {
  if (is.null(composition_pct[[model]])) return("")
  lookup <- category_lookups[[model]]
  svg <- build_composition_pct_svg(composition_pct[[model]], lookup$categories, lookup$colours)

  rows <- all_rows[all_rows$MODEL == model, ]
  rows$CATEGORY <- apply_category_lookup(rows, lookup)
  n_tab <- as.data.frame(rows[, .N, by = .(SAMPLE_LABEL, CATEGORY)])
  n_tab$pct <- composition_pct[[model]][cbind(n_tab$SAMPLE_LABEL, n_tab$CATEGORY)]
  n_tab <- n_tab[order(n_tab$SAMPLE_LABEL, -n_tab$pct), ]
  n_tab$pct <- sprintf("%.2f%%", n_tab$pct)
  names(n_tab) <- c("Sample", "Category", "Count", "Percent")
  table_html <- build_counts_table_details(
    n_tab, c("Sample", "Category", "Count", "Percent"),
    c("Sample", "Category", "Count", "Percent"),
    "Show feature counts &amp; percentages table")

  paste0(
    '<h3>', MODEL_TITLES[[model]], '</h3>',
    build_legend_html(lookup$categories, lookup$colours),
    '<div style="overflow:auto;border:1px solid #eee;">', svg, '</div>',
    table_html
  )
}, character(1))

abundance_html <- paste0(
  '<!DOCTYPE html><html><head><meta charset="utf-8">',
  '<title>Whole-genome motif abundance</title>',
  '<style>', REPORT_CSS, '</style></head><body><div id="page">',
  '<h1>Whole-genome motif abundance</h1>',
  '<p style="font-size:12px;color:#666;">Percentage of feature occurrence per model. ',
  "Hover a segment for its exact percentage.</p>",
  paste(abundance_sections, collapse = ""),
  TOOLTIP_DIV_CSS,
  '</div>', TOOLTIP_JS, '</body></html>'
)
abundance_path <- file.path(out_dir, "whole_genome_motif_abundance.html")
writeLines(abundance_html, abundance_path)
cat(sprintf("Report written: %s (%.1f MB)\n", abundance_path, file.info(abundance_path)$size / 1e6))

# ---- File 2: <output>_karyotypes.html --------------------------------
# One section per model: karyotype grid (hover + counts table) plus a
# collapsible sample-composition clustering dendrogram+heatmap.

karyo_section_html <- vapply(MODEL_ORDER, function(model) {
  if (is.null(grid_cells[[model]])) return("")
  lookup <- category_lookups[[model]]
  svg <- build_karyotype_grid_svg(grid_cells[[model]], lookup$categories, lookup$colours)

  gc <- grid_cells[[model]]
  chrom_order_factor <- factor(gc$chrom, levels = natural_sort_chroms(gc$chrom))
  tbl <- gc[order(gc$genome_group, gc$hap, chrom_order_factor, -gc$pct),
            c("genome_group", "hap", "chrom", "category", "n", "pct")]
  tbl$pct <- sprintf("%.2f%%", tbl$pct)
  names(tbl) <- c("Genome", "Haplotype", "Chromosome", "Category", "Count", "Percent")
  table_html <- build_counts_table_details(
    tbl, names(tbl), names(tbl), "Show feature counts &amp; percentages table")

  cluster_block <- {
    widget <- build_dendro_heatmap_widget(sample_cluster_mat[[model]], paste0("clust_", model))
    if (is.null(widget)) {
      '<p><em>Not enough samples to cluster (need at least 3).</em></p>'
    } else {
      paste0(
        '<details><summary>Show sample-similarity clustering</summary>',
        '<p style="font-size:12px;color:#666;margin-top:6px;">Euclidean distance on each ',
        'sample\u2019s category percentages (all chromosomes combined), average-linkage (UPGMA) ',
        'hierarchical clustering.</p>',
        '<div style="margin-top:4px;">', widget$html, '</div></details>'
      )
    }
  }

  paste0(
    '<section id="karyo_', model, '"><h2>', MODEL_TITLES[[model]], '</h2>',
    build_legend_html(lookup$categories, lookup$colours),
    '<div style="overflow:auto;border:1px solid #eee;">', svg, '</div>',
    table_html, cluster_block,
    '</section><hr>'
  )
}, character(1))

karyo_toc <- paste0(
  '<nav id="toc"><strong>Models:</strong> ',
  paste(sprintf('<a href="#karyo_%s">%s</a>', MODEL_ORDER, MODEL_TITLES[MODEL_ORDER]),
        collapse = " &middot; "),
  '</nav>'
)

karyo_html <- paste0(
  '<!DOCTYPE html><html><head><meta charset="utf-8">',
  '<title>CENP-B Karyotype Composition</title>',
  '<style>', REPORT_CSS, '</style></head><body><div id="page">',
  '<h1>CENP-B Karyotype Composition</h1>',
  '<p style="font-size:12px;color:#666;">Each bar is a chromosome/',
  'haplotype\u2019s own feature composition, stacked to 100%. Hover a segment ',
  'for its exact percentage.</p>',
  karyo_toc,
  karyo_section_html,
  TOOLTIP_DIV_CSS,
  '</div>', TOOLTIP_JS, '</body></html>'
)
karyo_path <- paste0(output_base, "_karyotypes.html")
writeLines(karyo_html, karyo_path)
cat(sprintf("Karyotype composition report written: %s (%.1f MB)\n", karyo_path,
            file.info(karyo_path)$size / 1e6))

cat("\nDone.\n")
