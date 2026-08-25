#!/usr/bin/env Rscript
# =============================================================================
#  centeny_plot_tool.R  —  Centeny Karyotype Plotter
#  Version: 1.0
# =============================================================================
#
#  USAGE
#  -----
#    Rscript centeny_plot_tool.R \
#      --cytoband      /path/to/chm13v2.0_cytobands_allchrs.bed \
#      --sample_sheet  samples.tsv \
#      --output        centeny_plot
#
#  OUTPUTS 
#  -------
#  With the default --format html:
#  1. <output>_<group>.html                     — interactive canvas viewer
#  2. <output>_<group>_metrics.tsv              — per-sample x per-chromosome metrics
#  3. <output>_<group>_cross_sample_summary.tsv — cross-assembly comparison
#
#  With a static --format (tiff/png/svg/pdf):
#  1. <output>_<group>.<fmt>                    — karyotype figure
#  2. <output>_<group>.html                     — also written if --html true
#  3. <output>_<group>_metrics.tsv              — per-sample x per-chromosome metrics
#  4. <output>_<group>_cross_sample_summary.tsv — cross-assembly comparison
#  where <group> is one of: centeny, model1, model2, model3
#
#  SAMPLE SHEET FORMAT
#  -------------------
#  Tab-separated, no header, 2 or 3 columns:
#    <sample_label>  <path_to_bed>  [<model>]
#  <model> is optional: centeny | model1 | model2 | model3 | auto (default).
#  Rows of different models can be freely mixed — see grouping note above.
#
#  GENOME SIZES FILE (--genome_sizes, required for --coord_mode percent)
#  -----------------------------------------------------------------------
#  Tab-separated, no header, 2 columns:
#    <sample_label>  <path_to_fai>
# =============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(karyoploteR)
})

# =============================================================================
#  1.  COMMAND-LINE ARGUMENT PARSING
# =============================================================================

option_list <- list(

  make_option(c("--workdir"), type = "character", default = getwd(),
    help = "Working directory [default: current directory]"),

  make_option(c("--cytoband"), type = "character", default = NULL,
    metavar = "FILE",
    help = "REQUIRED. Path to CHM13v2.0 cytoband BED file."),

  make_option(c("--sample_sheet"), type = "character", default = NULL,
    metavar = "FILE",
    help = "REQUIRED. TSV with: <label> <centeny_bed_path>"),

  make_option(c("--output"), type = "character", default = "centeny_plot.tiff",
    metavar = "FILE",
    help = "Output file name [default: centeny_plot.tiff]"),

  make_option(c("--chromosome"), type = "character", default = NULL,
    help = "Optional: plot only one chromosome (e.g. 7, X, Y). Default: all."),

  make_option(c("--jpeg_res"), type = "integer", default = 1200,
    help = "Resolution in DPI — PNG only [default: 1200]"),

  make_option(c("--jpeg_width"), type = "integer", default = 9000,
    help = "Width in pixels — PNG only [default: 9000]"),

  make_option(c("--jpeg_height"), type = "integer", default = 21600,
    help = "Height in pixels — PNG only [default: 21600]"),

  make_option(c("--format"), type = "character", default = "html",
    help = paste0("Output format: tiff, png, svg, pdf, html [default: html]. ",
                  "'html' produces ONLY the interactive canvas viewer (no ",
                  "static image) — use --html true alongside a different ",
                  "--format to get both.")),

  make_option(c("--plot_width"), type = "double", default = 7.5,
    help = "Plot width in inches — svg/pdf only [default: 7.5]"),

  make_option(c("--plot_height"), type = "double", default = 18.0,
    help = "Plot height in inches — svg/pdf only [default: 18.0]"),

  make_option(c("--chr_cex"), type = "double", default = 0.8,
    help = "Font size for chromosome labels [default: 0.8]"),

  make_option(c("--label_cex"), type = "double", default = 0.65,
    help = "Font size for sample track labels [default: 0.65]"),

  make_option(c("--tick_cex"), type = "double", default = 0.5,
    help = paste0("Font size for the genomic coordinate numbers under ",
                  "each ideogram (the 0 10 20 30... Mb axis) [default: 0.5]")),

  make_option(c("--track_height"), type = "double", default = 1.5,
    help = "Height of each sample track [default: 1.5]"),

  make_option(c("--track_gap"), type = "double", default = 5.0,
    help = "Vertical gap between sample tracks [default: 5.0]"),

  make_option(c("--major_tick"), type = "integer", default = 10000000,
    help = "Major tick label interval [default: 10 Mb]"),

  make_option(c("--minor_tick"), type = "integer", default = 5000000,
    help = "Minor tick interval [default: 5 Mb]"),

  make_option(c("--model"), type = "character", default = "auto",
    help = paste0("CENP-B model BED type: 'auto' (default; detected ",
                  "per-file from column 4), '1' (basepair distances), ",
                  "'2' (monomer-unit distances), or '3' (k-pattern ",
                  "clusters). 'auto' is recommended unless column 4 ",
                  "does not follow the D../M../#.. convention.")),

  make_option(c("--top_n_patterns"), type = "integer", default = 10,
    help = paste0("Model 3 only. When a k-pattern BED has a single ",
                  "colour throughout (i.e. not yet colour-coded), the N ",
                  "most abundant k-patterns are auto-coloured and every ",
                  "other pattern is grouped as 'Other' [default: 10]")),

  make_option(c("--legend_cex"), type = "double", default = 0.6,
    help = "Font size for the colour legend [default: 0.6]"),

  make_option(c("--legend_max_entries"), type = "integer", default = 30,
    help = paste0("Cap the colour legend at this many entries; the least ",
                  "abundant remaining colours are collapsed into a single ",
                  "'+ N more' entry so the legend never overflows the page ",
                  "[default: 30]")),

  make_option(c("--ideogram_gap"), type = "double", default = 15,
    help = paste0("Vertical gap between each chromosome ideogram and its ",
                  "feature track(s), in karyoploteR plot-param units. ",
                  "Increase for more breathing room [default: 15]")),

  make_option(c("--coord_mode"), type = "character", default = "bp",
    help = paste0("X-axis coordinate mode: 'bp' (real genomic positions, ",
                  "default) or 'percent' (each chromosome rescaled to 0-",
                  "100%% of its own length). Percent mode is useful for ",
                  "cross-species/cross-assembly comparison where chromosome ",
                  "identity/length correspondence isn't assumed — it lets ",
                  "you ask 'where along the chromosome' without claiming ",
                  "homology to a specific position. Requires --genome_sizes ",
                  "(each sample's OWN chromosome lengths — percent-of-what ",
                  "is meaningless without them). The cytoband ideogram ",
                  "itself stays on the shared --cytoband file's lengths ",
                  "either way, since it's one shared reference banding row, ",
                  "not per-sample.")),

  make_option(c("--genome_sizes"), type = "character", default = NULL,
    help = paste0("REQUIRED when --coord_mode percent. TSV, no header, 2 ",
                  "columns: <sample_label> <path_to_fai>, one row per ",
                  "sample in --sample_sheet, giving that sample's OWN ",
                  "chromosome lengths (a samtools faidx .fai index of that ",
                  "sample's FASTA). Needed because different samples can ",
                  "be different assemblies/individuals with genuinely ",
                  "different chromosome lengths — percent mode divides by ",
                  "each sample's own length, not the shared --cytoband ",
                  "file's length, or percentages would be wrong for any ",
                  "sample whose true chromosome size differs from the ",
                  "cytoband reference's.")),

  make_option(c("--cluster"), type = "character", default = "false",
    help = paste0("Also run an optional clustering analysis per model ",
                  "group: for every sample x chromosome, tabulate the ",
                  "proportion of features in each colour/category, then ",
                  "hierarchically cluster to see which chromosomes (within ",
                  "and across samples) share similar feature composition. ",
                  "Writes a heatmap+dendrogram as both PDF and interactive ",
                  "HTML. 'true'/'false' [default: false]")),

  make_option(c("--html"), type = "character", default = "false",
    help = paste0("Also write the interactive HTML alongside a non-html ",
                  "--format (e.g. --format pdf --html true gets you both). ",
                  "Not needed when --format html (the default) is used — ",
                  "that already produces the HTML on its own. ",
                  "'true'/'false' [default: false]")),

  make_option(c("--html_max_points"), type = "integer", default = 150000,
    help = paste0("Hard cap on total features rendered in the interactive ",
                  "HTML (summed across all chromosomes and sample tracks ",
                  "in a group); above this, points are spatially-binned ",
                  "downsampled (dense regions thinned, sparse regions kept ",
                  "intact) to keep the file small and the browser ",
                  "responsive. Set to 0 to disable the cap entirely and ",
                  "render every single point — file size and browser load ",
                  "time scale directly with your data size in that case. ",
                  "[default: 150000]"))
)

parser <- OptionParser(option_list = option_list)
opt    <- parse_args(parser)

# =============================================================================
#  2.  VALIDATE
# =============================================================================

if (is.null(opt$cytoband))          stop("--cytoband is required.",     call. = FALSE)
if (is.null(opt$sample_sheet))      stop("--sample_sheet is required.", call. = FALSE)
if (!file.exists(opt$cytoband))     stop("Cytoband file not found.",    call. = FALSE)
if (!file.exists(opt$sample_sheet)) stop("Sample sheet not found.",     call. = FALSE)

fmt <- tolower(trimws(opt$format))
if (!fmt %in% c("tiff", "png", "svg", "pdf", "html"))
  stop("--format must be one of: tiff, png, svg, pdf, html", call. = FALSE)

opt$model <- tolower(trimws(opt$model))
if (!opt$model %in% c("auto", "1", "2", "3"))
  stop("--model must be one of: auto, 1, 2, 3", call. = FALSE)

opt$html <- tolower(trimws(opt$html)) %in% c("true", "t", "yes", "y", "1")

opt$coord_mode <- tolower(trimws(opt$coord_mode))
if (!opt$coord_mode %in% c("bp", "percent"))
  stop("--coord_mode must be one of: bp, percent", call. = FALSE)
if (opt$coord_mode == "percent" && (is.null(opt$genome_sizes) || !nzchar(opt$genome_sizes)))
  stop(paste0("--coord_mode percent requires --genome_sizes (each sample's own ",
              "chromosome lengths) — see --help for the file format."), call. = FALSE)
if (!is.null(opt$genome_sizes) && nzchar(opt$genome_sizes) && !file.exists(opt$genome_sizes))
  stop(paste0("--genome_sizes file not found: ", opt$genome_sizes), call. = FALSE)

opt$cluster <- tolower(trimws(opt$cluster)) %in% c("true", "t", "yes", "y", "1")

# =============================================================================
#  3.  PARSE SAMPLE SHEET
#  ------------------------------------------------------------------------
#  Two OR three tab-separated columns, no header:
#    <label>  <bed_path>  [<model>]
#  <model> is optional. Accepted values (case-insensitive):
#    centeny | plain          -> plain centeny BED, colours trusted as-is
#    model1  | 1               -> distance-in-bp   BED (col4 = "D<bp>")
#    model2  | 2               -> distance-in-monomers BED (col4 = "M<n>")
#    model3  | 3               -> k-pattern cluster BED (col4 = "#d1_d2_..")
#    auto | (blank)            -> detect automatically from column 4
# =============================================================================

cat("Reading sample sheet:", opt$sample_sheet, "\n")

raw_sheet <- readLines(opt$sample_sheet)
raw_sheet <- raw_sheet[!grepl("^\\s*#", raw_sheet)]
raw_sheet <- raw_sheet[nchar(trimws(raw_sheet)) > 0]

resolve_model_tag <- function(tag) {
  tag <- tolower(trimws(tag))
  if (tag %in% c("", "auto"))            return("auto")
  if (tag %in% c("centeny", "plain"))    return("generic")
  if (tag %in% c("model1", "1"))         return("distance_bp")
  if (tag %in% c("model2", "2"))         return("monomer")
  if (tag %in% c("model3", "3"))         return("kpattern")
  stop(sprintf(paste0("Unrecognised model tag '%s' in sample sheet. Use ",
                       "one of: centeny, model1, model2, model3, auto ",
                       "(or leave the 3rd column blank)."), tag),
       call. = FALSE)
}

sample_df <- do.call(rbind, lapply(raw_sheet, function(line) {
  parts <- strsplit(line, "\t")[[1]]
  if (!length(parts) %in% c(2, 3))
    stop(paste0("Malformed sample sheet line (expected 2 or 3 tab-",
                "separated columns): ", line), call. = FALSE)
  data.frame(
    label     = trimws(parts[1]),
    bed       = trimws(parts[2]),
    model_tag = resolve_model_tag(if (length(parts) == 3) parts[3] else "auto"),
    stringsAsFactors = FALSE
  )
}))

missing <- sample_df$bed[!file.exists(sample_df$bed)]
if (length(missing) > 0)
  stop(paste("Missing BED files:\n", paste(missing, collapse = "\n")), call. = FALSE)

# =============================================================================
#  4.  HELPERS
# =============================================================================

rgb_to_hex <- function(rgb_string) {
  v <- suppressWarnings(as.numeric(unlist(strsplit(rgb_string, ","))))
  if (any(is.na(v)) || length(v) != 3) return("#808080")
  grDevices::rgb(v[1], v[2], v[3], maxColorValue = 255)
}

load_centeny_bed <- function(path, label) {
  dt <- fread(path, header = FALSE)
  if (ncol(dt) < 9) stop("BED must be BED9.", call. = FALSE)
  dt$HEX <- sapply(dt$V9, rgb_to_hex)
  dt$V1  <- gsub("^chr", "", dt$V1, ignore.case = TRUE)
  dt$V1  <- gsub("_hap[0-9]+$|_(pat|mat)(ernal)?$", "", dt$V1,
                  ignore.case = TRUE)
  dt
}

compute_track_positions <- function(n, base_start = 12,
                                    track_height = 1, gap = 1.15) {
  data.frame(
    r0 = base_start + (seq_len(n) - 1) * (track_height + gap),
    r1 = base_start + (seq_len(n) - 1) * (track_height + gap) + track_height
  )
}

pick_legend_row_window <- function(chr_lengths, chromosomes, span = 5) {
  rel <- chr_lengths[chromosomes] / max(chr_lengths[chromosomes])
  n <- length(chromosomes)
  span <- max(1, min(span, n))
  best_start <- 1; best_avg <- Inf
  for (s in seq_len(n - span + 1)) {
    avg <- mean(rel[s:(s + span - 1)])
    if (avg < best_avg) { best_avg <- avg; best_start <- s }
  }
  c(best_start, best_start + span - 1)
}

# -----------------------------------------------------------------------
#  Model detection & colour-legend helpers
# -----------------------------------------------------------------------

MODEL_LOOKUP <- c("1" = "distance_bp", "2" = "monomer", "3" = "kpattern")

detect_model_type <- function(v4) {
  v4 <- as.character(v4)
  if (any(grepl("^#", v4)))                                return("kpattern")
  if (any(grepl("^D[0-9]+$", v4)) || any(v4 == "DNA"))      return("distance_bp")
  if (any(grepl("^M[0-9]+$", v4)))                          return("monomer")
  "generic"
}

parse_model_value <- function(v4, model) {
  v4 <- as.character(v4)
  if (model == "distance_bp") {
    out <- suppressWarnings(as.numeric(sub("^D", "", v4)))
    out[!grepl("^D[0-9]+$", v4)] <- NA_real_
    return(out)
  }
  if (model == "monomer") {
    out <- suppressWarnings(as.numeric(sub("^M", "", v4)) + 1)
    out[!grepl("^M[0-9]+$", v4)] <- NA_real_
    return(out)
  }
  rep(NA_real_, length(v4))
}

trunc_label <- function(x, n = 28) {
  ifelse(nchar(x) > n, paste0(substr(x, 1, n - 1), "..."), x)
}

format_value_label <- function(lo, hi, unit) {
  fmt_one <- function(x) {
    if (unit == "bp") {
      if (x >= 1e6) return(sprintf("%.2f Mb", x / 1e6))
      if (x >= 1e3) return(sprintf("%.1f kb", x / 1e3))
      return(sprintf("%d bp", as.integer(round(x))))
    }
    sprintf("%d monomers", as.integer(round(x)))
  }
  if (isTRUE(all.equal(lo, hi))) return(fmt_one(lo))
  paste0(fmt_one(lo), " - ", fmt_one(hi))
}

legend_from_value_bins <- function(dt, model) {
  vals <- parse_model_value(dt$V4, model)
  unit <- if (model == "distance_bp") "bp" else "monomer"

  d <- data.table(HEX = dt$HEX, VAL = vals)
  agg <- d[, .(
              n      = .N,
              lo     = suppressWarnings(min(VAL, na.rm = TRUE)),
              hi     = suppressWarnings(max(VAL, na.rm = TRUE)),
              mean_v = suppressWarnings(mean(VAL, na.rm = TRUE))
            ), by = HEX]

  agg$is_na_only <- is.infinite(agg$lo) | is.nan(agg$mean_v) | is.na(agg$mean_v)

  span_threshold <- if (unit == "bp") 2e6 else 1000
  is_wide <- !agg$is_na_only & (agg$hi - agg$lo) > span_threshold
  if (sum(is_wide) > 0) {
    wide      <- agg[is_wide]
    canon_hex <- wide$HEX[which.max(wide$n)]
    merged <- data.table(
      HEX = canon_hex, n = sum(wide$n),
      lo = min(wide$lo), hi = max(wide$hi),
      mean_v = mean(c(min(wide$lo), max(wide$hi))), is_na_only = FALSE
    )
    agg <- rbind(agg[!is_wide], merged)
  }

  agg <- agg[order(is_na_only, mean_v)]

  labels <- vapply(seq_len(nrow(agg)), function(i) {
    if (agg$is_na_only[i]) {
      if (model == "distance_bp") "No distance (first motif / DNA)"
      else                        "No monomer distance (first motif)"
    } else if ((agg$hi[i] - agg$lo[i]) > span_threshold) {
      "Background (wide range of values)"
    } else {
      format_value_label(agg$lo[i], agg$hi[i], unit)
    }
  }, character(1))

  data.frame(label = labels, colour = agg$HEX, n = agg$n, stringsAsFactors = FALSE)
}

legend_from_existing_colours <- function(dt) {
  agg <- dt[, .(n = .N,
                top_pattern = names(sort(table(V4), decreasing = TRUE))[1]),
            by = HEX]
  agg <- agg[order(-n)]
  data.frame(
    label  = trunc_label(agg$top_pattern),
    colour = agg$HEX,
    n      = agg$n,
    stringsAsFactors = FALSE
  )
}

recolor_kpattern_topn <- function(dt, top_n = 10, other_colour = "#B0B0B0") {
  freq    <- sort(table(dt$V4), decreasing = TRUE)
  keep_n  <- min(top_n, length(freq))
  top_pat <- names(freq)[seq_len(keep_n)]
  pal     <- grDevices::hcl.colors(keep_n, palette = "Dynamic")
  colour_map <- setNames(pal, top_pat)

  dt$HEX <- other_colour
  for (p in top_pat) dt$HEX[dt$V4 == p] <- colour_map[[p]]

  other_patterns    <- setdiff(names(freq), top_pat)
  n_other_patterns  <- length(other_patterns)
  n_other_features  <- sum(dt$V4 %in% other_patterns)

  legend_df <- data.frame(
    label  = trunc_label(top_pat),
    colour = unname(colour_map),
    n      = as.integer(freq[top_pat]),
    stringsAsFactors = FALSE
  )
  if (n_other_patterns > 0) {
    legend_df <- rbind(legend_df, data.frame(
      label  = sprintf("Other (%d patterns)", n_other_patterns),
      colour = other_colour,
      n      = n_other_features,
      stringsAsFactors = FALSE
    ))
  }
  list(dt = dt, legend = legend_df)
}

apply_group_kpattern_colouring <- function(grp_bed_list, row_ids, top_n = 10,
                                            other_colour = "#B0B0B0") {
  if (length(row_ids) == 0) return(list(bed_list = grp_bed_list, legend = NULL))

  combined_v4 <- unlist(lapply(row_ids, function(rid) grp_bed_list[[rid]]$V4), use.names = FALSE)
  freq    <- sort(table(combined_v4), decreasing = TRUE)
  keep_n  <- min(top_n, length(freq))
  top_pat <- names(freq)[seq_len(keep_n)]
  pal     <- grDevices::hcl.colors(keep_n, palette = "Dynamic")
  colour_map <- setNames(pal, top_pat)

  for (rid in row_ids) {
    dt  <- grp_bed_list[[rid]]
    hex <- rep(other_colour, nrow(dt))
    m   <- match(dt$V4, top_pat)
    hit <- !is.na(m)
    hex[hit] <- colour_map[m[hit]]
    dt$HEX <- hex
    grp_bed_list[[rid]] <- dt
  }

  other_patterns   <- setdiff(names(freq), top_pat)
  n_other_patterns <- length(other_patterns)
  n_other_features <- sum(combined_v4 %in% other_patterns)

  legend_df <- data.frame(
    label  = trunc_label(top_pat),
    colour = unname(colour_map),
    n      = as.integer(freq[top_pat]),
    stringsAsFactors = FALSE
  )
  if (n_other_patterns > 0) {
    legend_df <- rbind(legend_df, data.frame(
      label  = sprintf("Other (%d patterns)", n_other_patterns),
      colour = other_colour,
      n      = n_other_features,
      stringsAsFactors = FALSE
    ))
  }
  list(bed_list = grp_bed_list, legend = legend_df)
}

legend_from_generic_colours <- function(dt) {
  if (length(unique(dt$HEX)) <= 1) return(NULL)
  strand_map <- c("#CC0000" = "Reverse strand", "#0000CC" = "Forward strand")
  agg <- dt[, .N, by = HEX][order(-N)]
  labels <- ifelse(agg$HEX %in% names(strand_map), strand_map[agg$HEX], "Other")
  data.frame(label = labels, colour = agg$HEX, n = agg$N, stringsAsFactors = FALSE)
}

cap_legend_entries <- function(legend_df, max_entries = 30) {
  no_op <- list(legend_df = legend_df, dropped_colours = character(0),
                placeholder_colour = "#D9D9D9")
  if (is.null(legend_df) || nrow(legend_df) <= max_entries) return(no_op)
  ord   <- order(-legend_df$n)
  keep  <- legend_df[ord[seq_len(max_entries - 1)], , drop = FALSE]
  dropped <- legend_df[ord[-seq_len(max_entries - 1)], , drop = FALSE]
  dropped_n <- sum(dropped$n)
  dropped_k <- nrow(dropped)
  placeholder <- "#D9D9D9"
  extra <- data.frame(
    label  = sprintf("+ %d more colour categories", dropped_k),
    colour = placeholder,
    n      = dropped_n,
    stringsAsFactors = FALSE
  )
  list(legend_df = rbind(keep, extra),
       dropped_colours = dropped$colour,
       placeholder_colour = placeholder)
}

# -----------------------------------------------------------------------
#  Interactive HTML builder
#  ------------------------------------------------------------------------

base64_encode_raw <- function(raw_bytes) {
  b64_chars <- c(LETTERS, letters, 0:9, "+", "/")
  n   <- length(raw_bytes)
  pad <- (3 - n %% 3) %% 3
  if (pad > 0) raw_bytes <- c(raw_bytes, as.raw(rep(0, pad)))
  m <- matrix(as.integer(raw_bytes), nrow = 3)
  idx <- rbind(
    bitwShiftR(m[1, ], 2),
    bitwOr(bitwShiftL(bitwAnd(m[1, ], 3L),  4L), bitwShiftR(m[2, ], 4)),
    bitwOr(bitwShiftL(bitwAnd(m[2, ], 15L), 2L), bitwShiftR(m[3, ], 6)),
    bitwAnd(m[3, ], 63L)
  )
  out <- paste(b64_chars[as.vector(idx) + 1L], collapse = "")
  if (pad == 1) out <- paste0(substr(out, 1, nchar(out) - 1), "=")
  if (pad == 2) out <- paste0(substr(out, 1, nchar(out) - 2), "==")
  out
}

pack_points <- function(chr_idx, track_idx, colour_idx, pos) {
  n <- length(pos)
  if (n == 0) return(raw(0))
  pos_raw <- writeBin(as.integer(pos), raw(), size = 4, endian = "little")
  pos_mat <- matrix(pos_raw, nrow = 4)
  track_raw <- writeBin(as.integer(track_idx), raw(), size = 2, endian = "little")
  track_mat <- matrix(track_raw, nrow = 2)
  out <- raw(n * 8L)
  base <- (seq_len(n) - 1L) * 8L
  out[base + 1L] <- as.raw(chr_idx)
  out[base + 2L] <- track_mat[1, ]
  out[base + 3L] <- track_mat[2, ]
  out[base + 4L] <- as.raw(colour_idx)
  out[base + 5L] <- pos_mat[1, ]
  out[base + 6L] <- pos_mat[2, ]
  out[base + 7L] <- pos_mat[3, ]
  out[base + 8L] <- pos_mat[4, ]
  out
}

giemsa_colour <- function(stain) {
  switch(as.character(stain),
    "gneg"    = "#FFFFFF",
    "gpos25"  = "#C8C8C8",
    "gpos50"  = "#969696",
    "gpos75"  = "#5A5A5A",
    "gpos100" = "#000000",
    "gpos"    = "#000000",
    "acen"    = "#D92523",
    "gvar"    = "#CFF0F0",
    "stalk"   = "#647FA4",
    "#BEBEBE"
  )
}

js_str <- function(x) paste0('"', gsub('(["\\\\])', '\\\\\\1', x), '"')

js_cytoband_array <- function(cytobands, chromosomes, chr_lengths) {
  rows <- character(0)
  for (ci in seq_along(chromosomes)) {
    this_chr <- chromosomes[ci]
    cb <- cytobands[cytobands$chr == this_chr, ]
    if (nrow(cb) == 0) {
      rows <- c(rows, sprintf("[%d,0,%d,%s]", ci - 1L,
                               as.integer(chr_lengths[this_chr]), js_str("#BEBEBE")))
    } else {
      seg <- sprintf("[%d,%d,%d,%s]", ci - 1L,
                      as.integer(round(cb$start)), as.integer(round(cb$end)),
                      js_str(vapply(cb$gieStain, giemsa_colour, character(1))))
      rows <- c(rows, seg)
    }
  }
  paste0("[", paste(rows, collapse = ","), "]")
}

spatially_downsample <- function(sub, chr_len, target_n, n_bins = 1000) {
  n <- nrow(sub)
  if (target_n >= n) return(sub)

  n_bins <- max(1L, min(n_bins, n))
  bin_w  <- max(1, chr_len / n_bins)
  mids   <- (sub$V2 + sub$V3) / 2
  bins   <- pmin(n_bins - 1L, floor(mids / bin_w))

  idx_by_bin <- split(seq_len(n), bins)
  counts_by_bin <- lengths(idx_by_bin)

  lo <- 0; hi <- max(counts_by_bin)
  for (iter in seq_len(40)) {
    mid <- (lo + hi) / 2
    tot <- sum(pmin(counts_by_bin, mid))
    if (tot > target_n) hi <- mid else lo <- mid
  }
  cap <- max(1, ceiling(lo))

  keep_idx <- unlist(lapply(idx_by_bin, function(idxs) {
    if (length(idxs) <= cap) return(idxs)
    sample(idxs, cap)
  }), use.names = FALSE)
  sub[keep_idx, ]
}

build_interactive_html <- function(grp_sample_df, grp_bed_list, chromosomes,
                                    chr_lengths_all, cytobands, legend_df,
                                    group_friendly, html_path, max_points = 150000,
                                    coord_mode = "bp", sample_chr_lengths = list()) {

  n_tracks <- nrow(grp_sample_df)
  if (n_tracks > 65535) {
    stop(sprintf(paste0("This group has %d sample tracks, which exceeds the interactive ",
                         "HTML's hard limit of 65535 (track index is packed as a 2-byte ",
                         "integer). Split the sample sheet into smaller groups, or use a ",
                         "static --format instead for this run."), n_tracks), call. = FALSE)
  }
  percent_mode <- identical(coord_mode, "percent")

  if (is.null(legend_df) || nrow(legend_df) == 0) {
    legend_df <- data.frame(label = "features", colour = "#333333",
                             n = 1, stringsAsFactors = FALSE)
  }
  legend_lookup <- setNames(seq_len(nrow(legend_df)) - 1L, legend_df$colour)
  uncategorised_idx <- nrow(legend_df)  # reserved extra slot, added below if needed

  chr_lengths_bp <- as.numeric(chr_lengths_all[chromosomes])
  names(chr_lengths_bp) <- chromosomes
  chr_lengths_js <- if (percent_mode) setNames(rep(100000L, length(chromosomes)), chromosomes) else chr_lengths_bp

  total_rows <- sum(vapply(grp_bed_list, function(d) sum(d$V1 %in% chromosomes), integer(1)))
  unlimited <- max_points <= 0
  sample_frac <- if (!unlimited && total_rows > max_points) max_points / total_rows else 1
  if (unlimited) {
    cat(sprintf("  HTML point cap disabled (--html_max_points 0): rendering all %s points\n",
                format(total_rows, big.mark = ",")))
  } else if (sample_frac < 1) {
    cat(sprintf("  Downsampling for HTML: %s points -> ~%s (cap=%s)\n",
                format(total_rows, big.mark = ","),
                format(round(total_rows * sample_frac), big.mark = ","),
                format(max_points, big.mark = ",")))
  }

  chr_idx_l   <- vector("list", n_tracks * length(chromosomes))
  track_idx_l <- vector("list", n_tracks * length(chromosomes))
  colour_idx_l <- vector("list", n_tracks * length(chromosomes))
  pos_l       <- vector("list", n_tracks * length(chromosomes))
  slot <- 0L
  any_uncategorised <- FALSE

  set.seed(42)
  for (ti in seq_len(n_tracks)) {
    dat <- grp_bed_list[[grp_sample_df$row_id[ti]]]
    sample_lens <- if (percent_mode) sample_chr_lengths[[grp_sample_df$label[ti]]] else NULL
    for (ci in seq_along(chromosomes)) {
      chr <- chromosomes[ci]
      sub <- dat[dat$V1 == chr, ]
      if (nrow(sub) == 0) next
      if (sample_frac < 1) {
        sub <- spatially_downsample(sub, chr_lengths_bp[chr],
                                     target_n = max(1, round(nrow(sub) * sample_frac)))
      }
      cidx <- legend_lookup[sub$HEX]
      if (any(is.na(cidx))) { any_uncategorised <- TRUE; cidx[is.na(cidx)] <- uncategorised_idx }
      slot <- slot + 1L
      chr_idx_l[[slot]]    <- rep.int(ci - 1L, nrow(sub))
      track_idx_l[[slot]]  <- rep.int(ti - 1L, nrow(sub))
      colour_idx_l[[slot]] <- as.integer(cidx)
      mids <- (sub$V2 + sub$V3) / 2
      pos_l[[slot]] <- if (percent_mode) {
        as.integer(round(mids / sample_lens[chr] * 100 * 1000))
      } else {
        as.integer(round(mids))
      }
    }
  }
  keep <- seq_len(slot)
  chr_idx    <- unlist(chr_idx_l[keep],    use.names = FALSE)
  track_idx  <- unlist(track_idx_l[keep],  use.names = FALSE)
  colour_idx <- unlist(colour_idx_l[keep], use.names = FALSE)
  pos        <- unlist(pos_l[keep],        use.names = FALSE)

  if (any_uncategorised) {
    legend_df <- rbind(legend_df, data.frame(
      label = "(uncategorised)", colour = "#999999", n = NA_integer_,
      stringsAsFactors = FALSE))
  }

  points_b64 <- base64_encode_raw(pack_points(chr_idx, track_idx, colour_idx, pos))

  cytobands_js_input <- cytobands
  if (percent_mode) {
    cytobands_js_input <- copy(cytobands)
    cytobands_js_input$start <- as.integer(round(
      cytobands_js_input$start / chr_lengths_bp[cytobands_js_input$chr] * 100 * 1000))
    cytobands_js_input$end   <- as.integer(round(
      cytobands_js_input$end   / chr_lengths_bp[cytobands_js_input$chr] * 100 * 1000))
  }
  cyto_js    <- js_cytoband_array(cytobands_js_input, chromosomes, chr_lengths_js)
  names_js   <- paste0("[", paste(js_str(paste0("chr", chromosomes)), collapse = ","), "]")
  lens_js    <- paste0("[", paste(as.integer(chr_lengths_js), collapse = ","), "]")
  tracks_js  <- paste0("[", paste(js_str(grp_sample_df$label), collapse = ","), "]")
  legend_js  <- paste0("[", paste(sprintf('{"label":%s,"colour":%s}',
                                           js_str(legend_df$label), js_str(legend_df$colour)),
                                   collapse = ","), "]")
  coord_mode_js <- js_str(coord_mode)

  html <- paste0(
'<!DOCTYPE html><html><head><meta charset="utf-8"><title>', group_friendly, '</title>
<style>
  html,body{margin:0;padding:0;background:#fff;font-family:-apple-system,Helvetica,Arial,sans-serif;color:#222;}
  #page{max-width:1500px;margin:0 auto;padding:12px 12px 220px 12px;}
  h1{font-size:16px;margin:6px 0 12px 2px;}
  .chrRow{margin-bottom:4px;}
  .chrCanvas{width:100%;height:', max(90, round(28 + n_tracks * 16)), 'px;
             display:block;cursor:grab;border-bottom:1px solid #eee;}
  #zoomBar{position:fixed;top:10px;right:14px;z-index:11;background:rgba(255,255,255,0.95);
           border:1px solid #ccc;border-radius:8px;padding:4px 6px;box-shadow:0 1px 6px rgba(0,0,0,0.15);}
  #zoomBar button{width:26px;height:26px;margin:0 2px;cursor:pointer;font-size:14px;border-radius:4px;
                   border:1px solid #ccc;background:#fafafa;}
  #zoomBar button:hover{background:#eee;}
  #legendPanel{position:fixed;right:14px;bottom:14px;max-width:340px;max-height:70vh;
               overflow:auto;background:rgba(255,255,255,0.97);border:1px solid #ccc;
               border-radius:8px;padding:8px 10px;font-size:11px;box-shadow:0 1px 6px rgba(0,0,0,0.18);}
  #legendPanel.minimized{max-height:none;overflow:visible;padding:6px 10px;}
  #legendPanel.minimized #legendItems{display:none;}
  #legendHead{display:flex;align-items:center;justify-content:space-between;gap:10px;cursor:pointer;}
  #legendHead h2{font-size:12px;margin:2px 0;}
  #legendToggle{font-size:13px;border:1px solid #ccc;border-radius:4px;background:#fafafa;
                width:20px;height:20px;line-height:18px;text-align:center;cursor:pointer;flex:0 0 auto;}
  .legendItem{display:flex;align-items:center;gap:6px;padding:2px 0;cursor:pointer;user-select:none;}
  .legendItem.off{opacity:0.35;}
  .swatch{width:11px;height:11px;flex:0 0 auto;border:1px solid #0002;}
  #tooltip{position:fixed;pointer-events:none;background:rgba(20,20,20,0.92);color:#fff;
           font-size:11px;padding:4px 7px;border-radius:4px;display:none;z-index:20;white-space:nowrap;}
  #hint{position:fixed;left:14px;bottom:14px;font-size:11px;color:#888;
        background:rgba(255,255,255,0.85);padding:3px 8px;border-radius:4px;}
</style></head><body>
<div id="zoomBar">
  <button id="zoomInAll" title="Zoom in all">+</button>
  <button id="zoomOutAll" title="Zoom out all">\u2013</button>
  <button id="zoomResetAll" title="Reset all">\u21ba</button>
</div>
<div id="page"><h1>', group_friendly, '</h1><div id="rows"></div></div>
<div id="legendPanel">
  <div id="legendHead"><h2>Colour key</h2><span id="legendToggle" title="Minimize">\u2212</span></div>
  <div id="legendItems"></div>
</div>
<div id="tooltip"></div>
<div id="hint">scroll or +/\u2013 = zoom &middot; drag = pan &middot; click legend to toggle</div>
<script>
(function(){
"use strict";
var CHR_NAMES = ', names_js, ';
var CHR_LENS  = ', lens_js, ';
var TRACKS    = ', tracks_js, ';
var LEGEND    = ', legend_js, ';
var CYTO      = ', cyto_js, ';
var N_TRACKS  = ', n_tracks, ';
var COORD_MODE = ', coord_mode_js, ';
var B64       = "', points_b64, '";

function b64ToBytes(b64) {
  var bin = atob(b64), len = bin.length, bytes = new Uint8Array(len);
  for (var i = 0; i < len; i++) bytes[i] = bin.charCodeAt(i);
  return bytes;
}
var raw = b64ToBytes(B64);
var nPts = Math.floor(raw.length / 8);
var dv = new DataView(raw.buffer);
var chrIdx = new Uint8Array(nPts), trackIdx = new Uint16Array(nPts),
    colIdx = new Uint8Array(nPts), pos = new Int32Array(nPts);
for (var i = 0; i < nPts; i++) {
  var o = i * 8;
  chrIdx[i]   = raw[o];
  trackIdx[i] = dv.getUint16(o+1, true);
  colIdx[i]   = raw[o+3];
  pos[i]      = dv.getInt32(o+4, true);
}

var byChr = [];
for (var c = 0; c < CHR_NAMES.length; c++) byChr.push([]);
for (var i = 0; i < nPts; i++) byChr[chrIdx[i]].push(i);

var hidden = new Set();
var legendItemsEl = document.getElementById("legendItems");
LEGEND.forEach(function(l, idx) {
  var div = document.createElement("div");
  div.className = "legendItem";
  div.innerHTML = \'<span class="swatch" style="background:\' + l.colour + \'"></span><span>\' + l.label + "</span>";
  div.onclick = function() {
    if (hidden.has(idx)) hidden.delete(idx); else hidden.add(idx);
    div.classList.toggle("off");
    states.forEach(redraw);
  };
  legendItemsEl.appendChild(div);
});

var legendPanel = document.getElementById("legendPanel");
var legendToggle = document.getElementById("legendToggle");
document.getElementById("legendHead").addEventListener("click", function() {
  var min = legendPanel.classList.toggle("minimized");
  legendToggle.textContent = min ? "+" : "\u2212";
  legendToggle.title = min ? "Expand" : "Minimize";
});

// Left-hand gutter reserved for track-name labels, kept fully separate
// from the plotted data area so features never render on top of names.
var GUTTER = 92;

function niceStep(bpRange, pxWidth) {
  var target = bpRange / (pxWidth / 70);
  var mag = Math.pow(10, Math.floor(Math.log10(target)));
  var norm = target / mag;
  var step = norm < 1.5 ? 1 : norm < 3.5 ? 2 : norm < 7.5 ? 5 : 10;
  return step * mag;
}
function fmtBp(v) {
  v = Math.round(v);
  if (COORD_MODE === "percent") return (v / 1000).toFixed(2) + "%";
  if (Math.abs(v) >= 1e6) return (v/1e6).toFixed(2) + " Mb";
  if (Math.abs(v) >= 1e3) return (v/1e3).toFixed(1) + " kb";
  return v + " bp";
}

var rowsEl = document.getElementById("rows");
var tooltip = document.getElementById("tooltip");
var states = [];

CHR_NAMES.forEach(function(name, ci) {
  var row = document.createElement("div");
  row.className = "chrRow";
  var label = document.createElement("div");
  label.style.fontSize = "11px"; label.style.fontWeight = "600"; label.style.margin = "6px 0 2px 2px";
  label.textContent = name;
  var canvas = document.createElement("canvas");
  canvas.className = "chrCanvas";
  row.appendChild(label); row.appendChild(canvas);
  rowsEl.appendChild(row);

  var state = {
    ci: ci, len: CHR_LENS[ci], canvas: canvas, ctx: canvas.getContext("2d"),
    viewStart: 0, viewEnd: CHR_LENS[ci], idx: byChr[ci],
    panning: false, panStartX: 0, panStartView: [0,0]
  };
  states.push(state);

  function resize() {
    var cssW = canvas.clientWidth || 800;
    var cssH = parseInt(getComputedStyle(canvas).height, 10) || 120;
    var dpr = window.devicePixelRatio || 1;
    canvas.width = Math.round(cssW * dpr);
    canvas.height = Math.round(cssH * dpr);
    state.ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    redraw(state);
  }
  window.addEventListener("resize", resize);
  setTimeout(resize, 0);

  canvas.addEventListener("wheel", function(e) {
    e.preventDefault();
    var rect = canvas.getBoundingClientRect();
    var plotW = Math.max(1, rect.width - GUTTER);
    var mx = e.clientX - rect.left - GUTTER;
    var frac = Math.max(0, Math.min(1, mx / plotW));
    var bpAt = state.viewStart + frac * (state.viewEnd - state.viewStart);
    var factor = e.deltaY < 0 ? 0.82 : 1/0.82;
    var newLen = Math.max(30, Math.min(state.len, (state.viewEnd - state.viewStart) * factor));
    state.viewStart = bpAt - frac * newLen;
    state.viewEnd   = state.viewStart + newLen;
    clampView(state);
    redraw(state);
  }, { passive: false });

  canvas.addEventListener("mousedown", function(e) {
    state.panning = true; state.panStartX = e.clientX;
    state.panStartView = [state.viewStart, state.viewEnd];
    canvas.style.cursor = "grabbing";
  });
  window.addEventListener("mousemove", function(e) {
    if (!state.panning) { hoverPick(state, e); return; }
    var rect = canvas.getBoundingClientRect();
    var plotW = Math.max(1, rect.width - GUTTER);
    var dxBp = (e.clientX - state.panStartX) / plotW * (state.panStartView[1] - state.panStartView[0]);
    state.viewStart = state.panStartView[0] - dxBp;
    state.viewEnd   = state.panStartView[1] - dxBp;
    clampView(state);
    redraw(state);
  });
  window.addEventListener("mouseup", function() {
    state.panning = false; canvas.style.cursor = "grab";
  });
  canvas.addEventListener("dblclick", function() {
    state.viewStart = 0; state.viewEnd = state.len; redraw(state);
  });
  canvas.addEventListener("mouseleave", function() { tooltip.style.display = "none"; });
});

function zoomAllBy(factor) {
  states.forEach(function(s) {
    var mid = (s.viewStart + s.viewEnd) / 2;
    var newLen = Math.max(30, Math.min(s.len, (s.viewEnd - s.viewStart) * factor));
    s.viewStart = mid - newLen / 2;
    s.viewEnd   = mid + newLen / 2;
    clampView(s);
    redraw(s);
  });
}
document.getElementById("zoomInAll").onclick    = function() { zoomAllBy(0.7); };
document.getElementById("zoomOutAll").onclick   = function() { zoomAllBy(1/0.7); };
document.getElementById("zoomResetAll").onclick = function() {
  states.forEach(function(s) { s.viewStart = 0; s.viewEnd = s.len; redraw(s); });
};

function clampView(s) {
  var len = s.viewEnd - s.viewStart;
  if (s.viewStart < 0) { s.viewStart = 0; s.viewEnd = len; }
  if (s.viewEnd > s.len) { s.viewEnd = s.len; s.viewStart = s.len - len; }
  if (s.viewStart < 0) s.viewStart = 0;
}

function hoverPick(state, e) {
  var rect = state.canvas.getBoundingClientRect();
  var mx = e.clientX - rect.left, my = e.clientY - rect.top;
  if (mx < GUTTER || my < 0 || mx > rect.width || my > rect.height) { tooltip.style.display = "none"; return; }
  var ideoH = 20, trackH = (rect.height - ideoH - 6) / N_TRACKS;
  if (my < ideoH) { tooltip.style.display = "none"; return; }
  var trackAt = Math.floor((my - ideoH) / trackH);
  var plotW = Math.max(1, rect.width - GUTTER);
  var bpAt = state.viewStart + ((mx - GUTTER) / plotW) * (state.viewEnd - state.viewStart);
  var bpTol = (state.viewEnd - state.viewStart) / plotW * 5;
  var best = -1, bestD = 1e18;
  var idxArr = state.idx;
  for (var k = 0; k < idxArr.length; k++) {
    var i = idxArr[k];
    if (trackIdx[i] !== trackAt) continue;
    if (pos[i] < state.viewStart - bpTol || pos[i] > state.viewEnd + bpTol) continue;
    var d = Math.abs(pos[i] - bpAt);
    if (d < bestD) { bestD = d; best = i; }
  }
  if (best === -1 || bestD > bpTol * 3) { tooltip.style.display = "none"; return; }
  tooltip.style.display = "block";
  tooltip.style.left = (e.clientX + 12) + "px";
  tooltip.style.top  = (e.clientY + 12) + "px";
  tooltip.innerHTML = TRACKS[trackAt] + "<br>" + CHR_NAMES[state.ci] + ":" +
                       fmtBp(pos[best]) + "<br>" + LEGEND[colIdx[best]].label;
}

function redraw(state) {
  var ctx = state.ctx, cssW = state.canvas.clientWidth, cssH = state.canvas.clientHeight;
  ctx.clearRect(0, 0, cssW, cssH);
  var vs = state.viewStart, ve = state.viewEnd, span = ve - vs;
  var plotW = Math.max(1, cssW - GUTTER);
  function xOf(p) { return GUTTER + (p - vs) / span * plotW; }

  var ideoY0 = 2, ideoY1 = 16;
  CYTO.forEach(function(seg) {
    if (seg[0] !== state.ci) return;
    if (seg[2] < vs || seg[1] > ve) return;
    var x0 = xOf(Math.max(seg[1], vs)), x1 = xOf(Math.min(seg[2], ve));
    ctx.fillStyle = seg[3];
    ctx.fillRect(x0, ideoY0, Math.max(0.6, x1 - x0), ideoY1 - ideoY0);
  });
  ctx.strokeStyle = "#333"; ctx.lineWidth = 1;
  ctx.strokeRect(GUTTER + 0.5, ideoY0 + 0.5, plotW - 1, ideoY1 - ideoY0 - 1);

  var step = niceStep(span, plotW);
  ctx.fillStyle = "#555"; ctx.font = "9px sans-serif"; ctx.textAlign = "center";
  var t0 = Math.ceil(vs / step) * step;
  for (var t = t0; t < ve; t += step) {
    var x = xOf(t);
    ctx.fillText(fmtBp(t), x, ideoY1 + 10);
    ctx.strokeStyle = "#eee"; ctx.beginPath(); ctx.moveTo(x, ideoY1); ctx.lineTo(x, cssH); ctx.stroke();
  }

  var trackTop = ideoY1 + 14, trackH = (cssH - trackTop) / N_TRACKS;

  // Track-name gutter: drawn on an opaque background of its own so
  // feature marks scrolled/zoomed under it never show through, and a
  // vertical divider makes the boundary with the plot area clear.
  ctx.fillStyle = "#fff";
  ctx.fillRect(0, trackTop, GUTTER - 4, cssH - trackTop);
  ctx.font = "9px sans-serif"; ctx.textAlign = "left"; ctx.fillStyle = "#333";
  for (var ti = 0; ti < N_TRACKS; ti++) {
    var label = TRACKS[ti];
    if (ctx.measureText(label).width > GUTTER - 10) {
      while (label.length > 1 && ctx.measureText(label + "\u2026").width > GUTTER - 10) {
        label = label.slice(0, -1);
      }
      label += "\u2026";
    }
    ctx.fillText(label, 2, trackTop + ti * trackH + trackH * 0.72);
  }
  ctx.strokeStyle = "#ccc"; ctx.beginPath();
  ctx.moveTo(GUTTER - 3, trackTop); ctx.lineTo(GUTTER - 3, cssH); ctx.stroke();

  var idxArr = state.idx;
  ctx.save();
  ctx.beginPath(); ctx.rect(GUTTER, 0, plotW, cssH); ctx.clip();
  ctx.lineWidth = 1.3;
  for (var k = 0; k < idxArr.length; k++) {
    var i = idxArr[k];
    if (hidden.has(colIdx[i])) continue;
    var p = pos[i];
    if (p < vs || p > ve) continue;
    var x = xOf(p);
    var yMid = trackTop + trackIdx[i] * trackH + trackH * 0.55;
    ctx.strokeStyle = LEGEND[colIdx[i]].colour;
    ctx.beginPath();
    ctx.moveTo(x, yMid - trackH * 0.32);
    ctx.lineTo(x, yMid + trackH * 0.32);
    ctx.stroke();
  }
  ctx.restore();
}
})();
</script>
</body></html>
')

  writeLines(html, html_path)
  invisible(html_path)
}

prepare_model_and_legend <- function(dt, forced_model = "auto",
                                      top_n_patterns = 10) {
  model <- if (forced_model == "auto") detect_model_type(dt$V4) else forced_model
  dt$MODEL <- model

  legend <- NULL
  needs_group_recolor <- FALSE
  if (model == "kpattern") {
    needs_group_recolor <- TRUE
  } else if (model %in% c("distance_bp", "monomer")) {
    legend <- legend_from_value_bins(dt, model)
  } else {
    legend <- legend_from_generic_colours(dt)
  }

  list(dt = dt, model = model, legend = legend,
       needs_group_recolor = needs_group_recolor)
}

kpAddBaseNumbers_manual <- function(kp, chr, chr_length,
                                    tick.dist    = 5e6,
                                    label.dist   = 10e6,
                                    cex          = 0.5,
                                    tick.y       = 0.02,
                                    label.y      = 0.05,
                                    unit_divisor = 1e6,
                                    unit_suffix  = "") {
  tick.pos  <- seq(0, chr_length, by = tick.dist)
  label.pos <- seq(0, chr_length, by = label.dist)
  kpSegments(kp, chr = chr, x0 = tick.pos, x1 = tick.pos,
             y0 = 0, y1 = tick.y, data.panel = "ideogram")
  kpText(kp, chr = chr, x = label.pos, y = label.y,
         labels = paste0(round(label.pos / unit_divisor, 1), unit_suffix),
         data.panel = "ideogram", cex = cex, pos = 3)
}

# -----------------------------------------------------------------------
#  Optional clustering analysis: heatmap + dendrogram
#  ------------------------------------------------------------------------
#  For every (sample, chromosome) pair in a group, tabulates what
#  proportion of its features fall into each legend colour/category, then
#  hierarchically clusters those rows to see which chromosomes — within a
#  sample and across different samples — share similar feature
#  composition (same frequency and type of feature). Opt-in via
#  --cluster true. Zero extra dependencies: hclust/heatmap() are base R.
# -----------------------------------------------------------------------

# Computes (leaf-position, height) line segments for a dendrogram from an
# hclust object — the same layout base R's plot.hclust uses internally —
# so both the PDF (via heatmap()) and the hand-built HTML viewer draw the
# exact same tree without re-implementing clustering logic in JS.
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

# Self-contained heatmap + dendrogram HTML viewer: dendrogram on the left
# (pre-computed segments, drawn on canvas), heatmap grid on the right
# (cell colour = proportion), hover tooltip with the exact value. Same
# zero-dependency canvas approach as the main interactive plot.
build_clustering_html <- function(mat, hc, group_friendly, html_path, row_colours = NULL) {
  dendro   <- compute_dendrogram_segments(hc)
  ordered  <- mat[hc$order, , drop = FALSE]
  n_rows   <- nrow(ordered)
  n_cols   <- ncol(ordered)
  ordered_row_colours <- if (!is.null(row_colours)) row_colours[hc$order] else rep("#888888", n_rows)

  row_labels_js  <- paste0("[", paste(js_str(rownames(ordered)), collapse = ","), "]")
  row_colours_js <- paste0("[", paste(js_str(ordered_row_colours), collapse = ","), "]")
  col_labels_js <- paste0("[", paste(js_str(colnames(ordered)), collapse = ","), "]")
  values_js     <- paste0("[", paste(sprintf("%.4f", as.vector(t(ordered))), collapse = ","), "]")
  seg_js        <- paste0("[", paste(vapply(dendro$segments, function(s)
                      sprintf("[%.4f,%.4f,%.4f,%.4f]", s[1], s[2], s[3], s[4]),
                      character(1)), collapse = ","), "]")

  # Fixed, generous per-row height (doesn't shrink to squeeze into a fixed
  # budget any more) so rows stay readable regardless of how many there
  # are — the container scrolls instead.
  row_h  <- max(15, min(24, round(3000 / n_rows)))
  cell_w <- max(14, min(46, round(900 / n_cols)))
  dendro_w <- 160
  chr_strip_w <- 10

  html <- paste0(
'<!DOCTYPE html><html><head><meta charset="utf-8"><title>', group_friendly, ' - clustering</title>
<style>
  html,body{margin:0;padding:0;background:#fff;font-family:-apple-system,Helvetica,Arial,sans-serif;color:#222;}
  #page{padding:14px;}
  h1{font-size:15px;margin:2px 0 10px 2px;}
  #wrap{overflow:auto;border:1px solid #eee;}
  canvas{display:block;}
  #tooltip{position:fixed;pointer-events:none;background:rgba(20,20,20,0.92);color:#fff;
           font-size:11px;padding:4px 7px;border-radius:4px;display:none;z-index:20;white-space:nowrap;}
  #hint{color:#888;font-size:11px;margin-top:6px;}
</style></head><body>
<div id="page">
<h1>', group_friendly, ' - chromosome pattern clustering</h1>
<div id="wrap"><canvas id="cv"></canvas></div>
<div id="hint">rows: sample | chromosome (coloured strip = chromosome identity) &middot; columns: feature category &middot; hover a cell for its exact proportion</div>
</div>
<div id="tooltip"></div>
<script>
(function(){
"use strict";
var ROWS = ', row_labels_js, ';
var ROW_COLOURS = ', row_colours_js, ';
var COLS = ', col_labels_js, ';
var VALS = ', values_js, ';   // row-major, length ROWS.length * COLS.length
var SEGS = ', seg_js, ';      // [pos0,height0,pos1,height1] per dendrogram line
var MAXH = ', dendro$max_height, ';
var ROW_H = ', row_h, ', CELL_W = ', cell_w, ', DENDRO_W = ', dendro_w, ', CHR_STRIP_W = ', chr_strip_w, ';
var LABEL_W = 210, HEADER_H = 90;

var nRows = ROWS.length, nCols = COLS.length;
var canvas = document.getElementById("cv");
var ctx = canvas.getContext("2d");
var totalW = LABEL_W + CHR_STRIP_W + DENDRO_W + nCols * CELL_W + 20;
var totalH = HEADER_H + nRows * ROW_H + 10;
var dpr = window.devicePixelRatio || 1;
canvas.style.width = totalW + "px";
canvas.style.height = totalH + "px";
canvas.width = Math.round(totalW * dpr);
canvas.height = Math.round(totalH * dpr);
ctx.setTransform(dpr, 0, 0, dpr, 0, 0);

function val(r, c) { return VALS[r * nCols + c]; }
function colour(v) {
  // White is reserved for EXACT zero ("no features of this category on
  // this chromosome"). Anything above zero starts at a clearly visible
  // light blue instead of fading toward white, so low-but-nonzero values
  // stay legible; ramps on to dark navy at the max value.
  if (v <= 0) return "#FFFFFF";
  var t = Math.max(0, Math.min(1, v));
  var stops = [[214,234,248],[76,114,176],[8,48,107]];
  var seg = t < 0.5 ? 0 : 1;
  var lt = t < 0.5 ? t / 0.5 : (t - 0.5) / 0.5;
  var a = stops[seg], b = stops[seg+1];
  var r = Math.round(a[0] + (b[0]-a[0])*lt);
  var g = Math.round(a[1] + (b[1]-a[1])*lt);
  var bch = Math.round(a[2] + (b[2]-a[2])*lt);
  return "rgb(" + r + "," + g + "," + bch + ")";
}

function heightToX(h) { return DENDRO_W - (h / (MAXH || 1)) * (DENDRO_W - 8); }
function posToY(p) { return HEADER_H + (p - 0.5) * ROW_H; }

function draw() {
  ctx.clearRect(0, 0, totalW, totalH);

  // dendrogram
  ctx.strokeStyle = "#666"; ctx.lineWidth = 1;
  SEGS.forEach(function(s) {
    ctx.beginPath();
    ctx.moveTo(LABEL_W + CHR_STRIP_W + heightToX(s[1]), posToY(s[0]));
    ctx.lineTo(LABEL_W + CHR_STRIP_W + heightToX(s[3]), posToY(s[2]));
    ctx.stroke();
  });

  // per-row chromosome-identity colour strip + matching label text colour
  for (var r = 0; r < nRows; r++) {
    ctx.fillStyle = ROW_COLOURS[r];
    ctx.fillRect(LABEL_W, HEADER_H + r * ROW_H, CHR_STRIP_W - 2, ROW_H - 2);
  }

  // row labels
  ctx.font = "9px sans-serif"; ctx.textAlign = "left";
  ctx.textBaseline = "middle";
  for (var r = 0; r < nRows; r++) {
    var label = ROWS[r];
    if (ctx.measureText(label).width > LABEL_W - 8) {
      while (label.length > 1 && ctx.measureText(label + "\u2026").width > LABEL_W - 8) {
        label = label.slice(0, -1);
      }
      label += "\u2026";
    }
    ctx.fillStyle = ROW_COLOURS[r];
    ctx.fillText(label, 2, HEADER_H + r * ROW_H + ROW_H / 2);
  }

  // column headers (rotated)
  ctx.save();
  ctx.font = "9px sans-serif"; ctx.fillStyle = "#222";
  var gx0 = LABEL_W + CHR_STRIP_W + DENDRO_W;
  for (var c = 0; c < nCols; c++) {
    var cx = gx0 + c * CELL_W + CELL_W / 2;
    ctx.save();
    ctx.translate(cx, HEADER_H - 6);
    ctx.rotate(-Math.PI / 3);
    ctx.textAlign = "left"; ctx.textBaseline = "middle";
    var label = COLS[c];
    if (label.length > 26) label = label.slice(0, 25) + "\u2026";
    ctx.fillText(label, 0, 0);
    ctx.restore();
  }
  ctx.restore();

  // heatmap grid — a slightly larger inset than before so the extra row
  // height reads as clear spacing between rows, not just bigger cells
  for (var r = 0; r < nRows; r++) {
    for (var c = 0; c < nCols; c++) {
      ctx.fillStyle = colour(val(r, c));
      ctx.fillRect(gx0 + c * CELL_W, HEADER_H + r * ROW_H, CELL_W - 1, ROW_H - 2);
    }
  }
}
draw();

var tooltip = document.getElementById("tooltip");
canvas.addEventListener("mousemove", function(e) {
  var rect = canvas.getBoundingClientRect();
  var mx = e.clientX - rect.left, my = e.clientY - rect.top;
  var gx0 = LABEL_W + CHR_STRIP_W + DENDRO_W;
  if (mx < gx0 || my < HEADER_H) { tooltip.style.display = "none"; return; }
  var c = Math.floor((mx - gx0) / CELL_W);
  var r = Math.floor((my - HEADER_H) / ROW_H);
  if (c < 0 || c >= nCols || r < 0 || r >= nRows) { tooltip.style.display = "none"; return; }
  tooltip.style.display = "block";
  tooltip.style.left = (e.clientX + 12) + "px";
  tooltip.style.top  = (e.clientY + 12) + "px";
  tooltip.innerHTML = ROWS[r] + "<br>" + COLS[c] + "<br>" + (val(r,c)*100).toFixed(1) + "% of features";
});
canvas.addEventListener("mouseleave", function() { tooltip.style.display = "none"; });
})();
</script>
</body></html>
')

  writeLines(html, html_path)
  invisible(html_path)
}

run_clustering_analysis <- function(grp_sample_df, grp_bed_list, chromosomes,
                                     legend_df, group_friendly, out_base_grp) {
  if (is.null(legend_df) || nrow(legend_df) < 2) {
    cat("  Skipping clustering: fewer than 2 colour categories in this group.\n")
    return(invisible(NULL))
  }
  cat("  Building feature-pattern matrix for clustering...\n")

  cats <- sub("  \\(n=.*\\)$", "", legend_df$label)
  hex_to_cat <- setNames(cats, legend_df$colour)

  noise_regex <- paste0(
    "^(Background \\(wide range of values\\)",
    "|No distance \\(first motif / DNA\\)",
    "|No monomer distance \\(first motif\\)",
    "|Other \\([0-9]+ patterns\\)",
    "|\\+ [0-9]+ more colour categories",
    "|\\(uncategorised\\))$"
  )
  is_noise  <- grepl(noise_regex, cats)
  cats_keep <- unique(cats[!is_noise])
  if (length(cats_keep) < 2) {
    cat("  Skipping clustering: fewer than 2 informative (non-background) categories.\n")
    return(invisible(NULL))
  }
  noise_cats <- unique(cats[is_noise])

  n_chr <- length(chromosomes)
  hues  <- ((seq_len(n_chr) - 1) * 137.508) %% 360
  chr_palette <- setNames(grDevices::hcl(h = hues, c = 80, l = 55), chromosomes)

  rows <- list()
  row_names <- character(0)
  row_chr   <- character(0)
  for (i in seq_len(nrow(grp_sample_df))) {
    lbl <- grp_sample_df$label[i]
    dat <- grp_bed_list[[grp_sample_df$row_id[i]]]
    for (chr in chromosomes) {
      sub <- dat[dat$V1 == chr, ]
      if (nrow(sub) == 0) next
      cat_lab <- unname(hex_to_cat[sub$HEX])
      cat_lab[is.na(cat_lab)] <- "(uncategorised)"
      cat_lab <- cat_lab[!(cat_lab %in% noise_cats)]
      if (length(cat_lab) == 0) next  # nothing informative left on this chromosome
      tab  <- table(factor(cat_lab, levels = cats_keep))
      prop <- as.numeric(tab) / sum(tab)
      rows[[length(rows) + 1]] <- prop
      row_names <- c(row_names, paste0(lbl, " | chr", chr))
      row_chr   <- c(row_chr, chr)
    }
  }
  if (length(rows) < 3) {
    cat("  Skipping clustering: fewer than 3 sample x chromosome rows with data.\n")
    return(invisible(NULL))
  }

  mat <- do.call(rbind, rows)
  colnames(mat) <- cats_keep
  mat <- mat[, colSums(mat, na.rm = TRUE) > 0, drop = FALSE]
  rownames(mat) <- row_names

  if (ncol(mat) < 2) {
    cat("  Skipping clustering: fewer than 2 non-empty categories after filtering.\n")
    return(invisible(NULL))
  }

  d  <- dist(mat, method = "euclidean")
  hc <- hclust(d, method = "average")

  row_colours <- unname(chr_palette[row_chr])

  # ---- PDF: base-R heatmap(), zero extra dependency -----------------------
  # White is reserved for exact zero ("no features of this category on
  # this chromosome") — everything above zero starts at a clearly visible
  # light blue
  n_ramp   <- 63
  col_vec  <- c("#FFFFFF", colorRampPalette(c("#D6EAF8", "#4C72B0", "#08306B"))(n_ramp))
  mx       <- max(mat, na.rm = TRUE)
  if (!is.finite(mx) || mx <= 0) mx <- 1
  breaks   <- c(-1, 0, seq(mx / n_ramp, mx, length.out = n_ramp))

  pdf_path <- paste0(out_base_grp, "_clustering.pdf")
  pdf_h <- max(6, 0.26 * nrow(mat) + 3)
  pdf_w <- max(8, 0.30 * ncol(mat) + 9)
  pdf(pdf_path, width = pdf_w, height = pdf_h)
  tryCatch({
    heatmap(mat, Rowv = as.dendrogram(hc), Colv = NA, scale = "none",
            col = col_vec, breaks = breaks,
            RowSideColors = row_colours,
            margins = c(10, 14),
            cexRow = max(0.25, min(0.9, 30 / nrow(mat))),
            cexCol = max(0.4, min(0.8, 25 / ncol(mat))),
            main = paste0(group_friendly, "\nChromosome pattern clustering"))
  }, error = function(e) {
    cat("  NOTE: clustering heatmap PDF failed (", conditionMessage(e), ")\n", sep = "")
  })
  dev.off()
  cat(sprintf("  Clustering heatmap (PDF)  -> %s\n", pdf_path))

  html_path <- paste0(out_base_grp, "_clustering.html")
  build_clustering_html(mat, hc, group_friendly, html_path, row_colours = row_colours)
  cat(sprintf("  Clustering heatmap (HTML) -> %s\n", html_path))

  invisible(list(pdf = pdf_path, html = html_path))
}

# =============================================================================
#  5.  LOAD DATA & RESOLVE MODEL PER ROW
# =============================================================================

cat("Loading cytoband/genome definition:", opt$cytoband, "\n")

load_genome_from_cytobands <- function(path) {
  cb <- fread(path, header = FALSE)
  if (ncol(cb) < 5) stop("Cytoband file must have >= 5 columns (chr,start,end,name,gieStain).", call. = FALSE)
  colnames(cb)[1:5] <- c("chr", "start", "end", "name", "gieStain")
  cb$chr <- gsub("^chr", "", cb$chr, ignore.case = TRUE)
  cb$chr <- gsub("_hap[0-9]+$|_(pat|mat)(ernal)?$", "", cb$chr, ignore.case = TRUE)

  lens <- cb[, .(len = max(end)), by = chr]
  chr_lengths <- setNames(as.numeric(lens$len), lens$chr)

  # Natural sort: numeric chromosome names ascending first (1,2,...,22,23)
  nm     <- names(chr_lengths)
  is_num <- grepl("^[0-9]+$", nm)
  ord    <- c(nm[is_num][order(as.integer(nm[is_num]))], sort(nm[!is_num]))
  chr_lengths <- chr_lengths[ord]

  list(cytobands = cb, chr_lengths = chr_lengths, chromosomes = names(chr_lengths))
}

genome_info  <- load_genome_from_cytobands(opt$cytoband)
cytobands    <- genome_info$cytobands
all_chromosomes <- genome_info$chromosomes
cat(sprintf("  Detected %d chromosomes: %s\n", length(all_chromosomes),
            paste(head(all_chromosomes, 6), collapse = ", ")),
    if (length(all_chromosomes) > 6) sprintf("... (+%d more)\n", length(all_chromosomes) - 6) else "\n")

if (!is.null(opt$chromosome) && !(opt$chromosome %in% all_chromosomes)) {
  stop(sprintf("--chromosome '%s' not found in the cytoband file. Available: %s",
               opt$chromosome, paste(all_chromosomes, collapse = ", ")), call. = FALSE)
}

cat("Loading CenTeny BED files...\n")
sample_df$row_id <- paste0("row", seq_len(nrow(sample_df)))

bed_list <- lapply(seq_len(nrow(sample_df)), function(i) {
  cat("  Loaded:", sample_df$label[i], "(", sample_df$model_tag[i], ")\n")
  load_centeny_bed(sample_df$bed[i], sample_df$label[i])
})
names(bed_list) <- sample_df$row_id

# Guard against chromosomes present in a BED file but absent from the
# cytoband-derived genome
cat("Validating chromosome names against the cytoband-derived genome...\n")
any_dropped <- FALSE
for (idx in seq_len(nrow(sample_df))) {
  rid <- sample_df$row_id[idx]
  lbl <- sample_df$label[idx]
  dt  <- bed_list[[rid]]
  bad <- !(dt$V1 %in% all_chromosomes)
  if (any(bad)) {
    any_dropped <- TRUE
    bad_chrs <- sort(unique(dt$V1[bad]))
    cat(sprintf(paste0("  WARNING: '%s' has %d feature(s) on chromosome(s) not ",
                        "present in the cytoband file (%s) - excluded from ",
                        "plotting/metrics.\n"),
                lbl, sum(bad), paste(bad_chrs, collapse = ", ")))
    bed_list[[rid]] <- dt[!bad, ]
  }
}
if (any_dropped) {
  cat(paste0("  NOTE: if these are DIFFERENT SPECIES sharing one sample sheet, a ",
             "single cytoband file cannot correctly represent all of them (chr23/",
             "24 in apes has no human equivalent, and same-numbered autosomes ",
             "differ in true length between species, which also silently skews ",
             "--coord_mode percent). Run each species as its own --cytoband + ",
             "--sample_sheet pass instead — the pipeline supports any species' ",
             "cytoband already, see the reply above this run.\n"))
}

# ---- per-sample chromosome lengths (percent mode only) --------------------
sample_chr_lengths <- list()
if (opt$coord_mode == "percent") {
  cat("Loading per-sample genome sizes:", opt$genome_sizes, "\n")
  raw_sizes <- readLines(opt$genome_sizes)
  raw_sizes <- raw_sizes[!grepl("^\\s*#", raw_sizes)]
  raw_sizes <- raw_sizes[nchar(trimws(raw_sizes)) > 0]

  genome_sizes_df <- do.call(rbind, lapply(raw_sizes, function(line) {
    parts <- strsplit(line, "\t")[[1]]
    if (length(parts) != 2)
      stop(paste0("Malformed --genome_sizes line (expected 2 tab-separated ",
                  "columns: <sample_label> <path_to_fai>): ", line), call. = FALSE)
    data.frame(label = trimws(parts[1]), fai = trimws(parts[2]), stringsAsFactors = FALSE)
  }))

  missing_fai <- genome_sizes_df$fai[!file.exists(genome_sizes_df$fai)]
  if (length(missing_fai) > 0)
    stop(paste("Missing .fai file(s) referenced in --genome_sizes:\n",
               paste(missing_fai, collapse = "\n")), call. = FALSE)

  needed_labels  <- unique(sample_df$label)
  missing_labels <- setdiff(needed_labels, genome_sizes_df$label)
  if (length(missing_labels) > 0)
    stop(paste0("--coord_mode percent: the following sample label(s) from ",
                "--sample_sheet have no matching row in --genome_sizes:\n  ",
                paste(missing_labels, collapse = "\n  ")), call. = FALSE)

  for (i in seq_len(nrow(genome_sizes_df))) {
    lbl <- genome_sizes_df$label[i]
    fai_path <- genome_sizes_df$fai[i]
    fai <- fread(fai_path, header = FALSE)

    if (ncol(fai) < 2 || nrow(fai) == 0) {
      stop(sprintf(paste0(
        "'%s' (for sample '%s') doesn't look like a valid .fai index — expected ",
        ">=2 tab-separated columns (name, length, ...) but got %d column(s) and ",
        "%s row(s). The most common cause: this path points at the FASTA file ",
        "itself, not its .fai index. Generate one with:\n  samtools faidx %s\n",
        "...then point --genome_sizes at the resulting '%s.fai' file instead."),
        fai_path, lbl, ncol(fai), format(nrow(fai), big.mark = ","),
        fai_path, fai_path), call. = FALSE)
    }

    lens_num <- suppressWarnings(as.numeric(fai$V2))
    if (any(is.na(lens_num))) {
      stop(sprintf(paste0(
        "'%s' (for sample '%s'): column 2 isn't all numeric lengths as expected ",
        "in a .fai index. Check this is really a `samtools faidx` output, not a ",
        "different file."), fai_path, lbl), call. = FALSE)
    }

    # same chromosome-name cleanup as load_centeny_bed(), so names match
    chrom_clean <- gsub("^chr", "", as.character(fai$V1), ignore.case = TRUE)
    chrom_clean <- gsub("_hap[0-9]+$|_(pat|mat)(ernal)?$", "", chrom_clean, ignore.case = TRUE)
    sample_chr_lengths[[lbl]] <- setNames(lens_num, chrom_clean)
  }
  cat(sprintf("  Loaded chromosome lengths for %d sample(s)\n", length(sample_chr_lengths)))

  # A sample's own .fai might legitimately not cover every chromosome the
  # BED references (e.g. no assembled Y) — drop those rows with a clear
  # warning rather
  cat("Validating chromosome names against each sample's own --genome_sizes...\n")
  any_dropped_gs <- FALSE
  for (idx in seq_len(nrow(sample_df))) {
    rid  <- sample_df$row_id[idx]
    lbl  <- sample_df$label[idx]
    dt   <- bed_list[[rid]]
    lens <- sample_chr_lengths[[lbl]]
    bad  <- !(dt$V1 %in% names(lens))
    if (any(bad)) {
      any_dropped_gs <- TRUE
      bad_chrs <- sort(unique(dt$V1[bad]))
      cat(sprintf(paste0("  WARNING: '%s' has %d feature(s) on chromosome(s) not ",
                          "present in its own --genome_sizes .fai (%s) - excluded ",
                          "from plotting/metrics.\n"),
                  lbl, sum(bad), paste(bad_chrs, collapse = ", ")))
      bed_list[[rid]] <- dt[!bad, ]
    }
  }
  if (any_dropped_gs) {
    cat(paste0("  NOTE: this usually means that sample's assembly genuinely lacks ",
               "that chromosome (e.g. no assembled Y), or its .fai uses a different ",
               "chromosome-naming convention than its BED file.\n"))
  }
}

cat("Detecting/applying CENP-B model & building colour legend...\n")
legend_entries  <- list()
resolved_models <- character(nrow(sample_df))
pending_kpattern_recolor <- character(0)
for (idx in seq_len(nrow(sample_df))) {
  rid        <- sample_df$row_id[idx]
  lbl        <- sample_df$label[idx]
  row_forced <- sample_df$model_tag[idx]
  if (row_forced == "auto" && opt$model != "auto") {
    row_forced <- MODEL_LOOKUP[[opt$model]]
  }
  prep <- prepare_model_and_legend(bed_list[[rid]],
                                    forced_model   = row_forced,
                                    top_n_patterns = opt$top_n_patterns)
  bed_list[[rid]]       <- prep$dt
  legend_entries[[rid]] <- prep$legend
  resolved_models[idx]  <- prep$model
  if (isTRUE(prep$needs_group_recolor)) {
    pending_kpattern_recolor <- c(pending_kpattern_recolor, rid)
  }
  cat(sprintf("  %-20s model=%-12s colours=%d\n",
              lbl, prep$model, length(unique(prep$dt$HEX))))
}
sample_df$resolved_model <- resolved_models

# =============================================================================
#  6.  GROUP BY MODEL
# =============================================================================

group_map <- c(generic     = "centeny",
                distance_bp = "model1",
                monomer     = "model2",
                kpattern    = "model3")

group_friendly_map <- c(generic     = "CenTeny map",
                         distance_bp = "Model 1 - distance (bp)",
                         monomer     = "Model 2 - distance (monomers)",
                         kpattern    = "Model 3 - k-pattern clusters")

groups <- split(seq_len(nrow(sample_df)), sample_df$resolved_model)

cat(sprintf("\nDetected %d group(s): %s\n\n", length(groups),
            paste(group_map[names(groups)], collapse = ", ")))

output_base <- if (startsWith(opt$output, "/")) {
  tools::file_path_sans_ext(opt$output)
} else {
  file.path(opt$workdir, tools::file_path_sans_ext(opt$output))
}

out_dir <- dirname(output_base)
if (!dir.exists(out_dir)) {
  cat("Creating output directory:", out_dir, "\n")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
}
if (!dir.exists(out_dir))
  stop(paste0("Could not create output directory: ", out_dir,
              "\nCheck the path and your write permissions."), call. = FALSE)

chr_lengths_all <- genome_info$chr_lengths

# =============================================================================
#  7.  RENDER ONE GROUP
# =============================================================================

render_group <- function(grp_sample_df, grp_bed_list, grp_legend_entries,
                          group_tag, group_friendly) {

  n_tracks <- nrow(grp_sample_df)
  cat(sprintf("=== Rendering %s (%d sample track%s) ===\n",
              group_friendly, n_tracks, if (n_tracks == 1) "" else "s"))

  # ---- combine + aggregate legend across all samples in this group -------
  legend_df <- do.call(rbind, grp_legend_entries)
  if (!is.null(legend_df) && nrow(legend_df) > 0) {
    legend_dt <- as.data.table(legend_df)
    legend_dt <- legend_dt[order(-n)]
    legend_df <- as.data.frame(legend_dt[, .(colour = colour[1], n = sum(n)), by = label])
    legend_df <- legend_df[order(-legend_df$n), ]
  }
  legend_df <- cap_legend_entries(legend_df, max_entries = opt$legend_max_entries)
  cap_dropped_colours     <- legend_df$dropped_colours
  cap_placeholder_colour  <- legend_df$placeholder_colour
  legend_df <- legend_df$legend_df

  # cap_legend_entries() only shrinks the DISPLAYED legend — the features
  # themselves still carry their real original colours. Without this,
  # every capped-out category's points keep rendering in their own
  # distinct colour that no longer appears anywhere in the legend.
  if (length(cap_dropped_colours) > 0) {
    for (rid in names(grp_bed_list)) {
      dt  <- grp_bed_list[[rid]]
      hit <- dt$HEX %in% cap_dropped_colours
      if (any(hit)) {
        dt$HEX[hit] <- cap_placeholder_colour
        grp_bed_list[[rid]] <- dt
      }
    }
  }

  if (!is.null(legend_df) && nrow(legend_df) > 0 && group_tag == "centeny") {
    legend_df$label <- paste0(legend_df$label, "  (n=", format(legend_df$n, big.mark = ","), ")")
  }

  legend_ncol <- 1
  legend_rows <- 0
  if (!is.null(legend_df) && nrow(legend_df) > 0) {
    rows_per_col <- 14
    legend_ncol  <- max(1, ceiling(nrow(legend_df) / rows_per_col))
    legend_rows  <- ceiling(nrow(legend_df) / legend_ncol)
  }

  # ---- coordinate mode (bp vs percent) -------------------------------------
  # Percent mode rescales every chromosome to 0-100% of its own length
  percent_mode <- identical(opt$coord_mode, "percent")

  rescale_to_percent <- function(start, end, chr, chr_lengths) {
    len <- chr_lengths[chr]
    if (any(is.na(len))) {
      bad <- sort(unique(chr[is.na(len)]))
      stop(paste0("Internal error: percent-mode rescale hit chromosome(s) with ",
                  "no known length (", paste(bad, collapse = ", "), "). This ",
                  "should have been filtered out already — please report this."),
           call. = FALSE)
    }
    list(start = start / len * 100, end = end / len * 100)
  }

  cytobands_plot <- cytobands
  if (percent_mode) {
    cytobands_plot <- copy(cytobands)
    r <- rescale_to_percent(cytobands_plot$start, cytobands_plot$end,
                             cytobands_plot$chr, chr_lengths_all)
    cytobands_plot$start <- r$start
    cytobands_plot$end   <- r$end
  }

  grp_bed_list_plot <- grp_bed_list
  if (percent_mode) {
    grp_bed_list_plot <- lapply(grp_sample_df$row_id, function(rid) {
      dt  <- grp_bed_list[[rid]]
      lbl <- grp_sample_df$label[grp_sample_df$row_id == rid][1]
      lens <- sample_chr_lengths[[lbl]]
      dt2 <- copy(dt)
      r <- rescale_to_percent(dt2$V2, dt2$V3, dt2$V1, lens)
      dt2$V2 <- r$start
      dt2$V3 <- r$end
      dt2
    })
    names(grp_bed_list_plot) <- grp_sample_df$row_id
  }

  # ---- layout --------------------------------------------------------------
  pp <- getDefaultPlotParams(plot.type = 1)
  pp$ideogramheight <- 30
  pp$topmargin      <- 20
  pp$bottommargin   <- 20
  pp$leftmargin     <- 0.25
  pp$rightmargin    <- 0.01
  # dedicated gap between each ideogram and its feature track(s)
  pp$data1inmargin  <- opt$ideogram_gap
  # Track (sample) label font auto-shrinks as more tracks are stacked into
  # the same row, so names like "RPE1_hap1" stop overlapping each other
  # once there are a dozen-plus samples.
  label_scale <- min(1, 9 / n_tracks)

  if (!is.null(opt$chromosome)) {

    chromosomes      <- opt$chromosome
    cat("  Plotting only chromosome:", chromosomes, "\n")
    pp$ideogramheight <- 12
    pp$topmargin      <- 10
    legend_bottommargin_extra <- if (legend_rows > 0) 14 + legend_rows * 5.5 else 0
    pp$bottommargin   <- 10 + legend_bottommargin_extra * 0.5
    pp$leftmargin     <- 0.25
    pp$rightmargin    <- 0.01
    pp$data1inmargin  <- opt$ideogram_gap * 0.5

    track_height_use <- 1
    track_gap_use    <- 2
    scale_factor     <- 1.8
    chr_cex_use      <- opt$chr_cex   * scale_factor
    label_cex_use    <- opt$label_cex * 2.5 * label_scale
    tick_cex         <- opt$tick_cex * scale_factor
    tick_y_offset    <- 0.10
    label_y_offset   <- 0.60

    track_pos <- compute_track_positions(n_tracks,
                   track_height = track_height_use,
                   gap          = track_gap_use,
                   base_start   = 2)

    pp$data1max <- max(track_pos$r1) + opt$ideogram_gap * 0.5
    pp$data2max <- 0

  } else {

    chromosomes    <- all_chromosomes
    tick_cex       <- opt$tick_cex
    tick_y_offset  <- 0.06
    label_y_offset <- 0.15
    chr_cex_use    <- opt$chr_cex
    label_cex_use  <- opt$label_cex * label_scale

    track_pos <- compute_track_positions(n_tracks,
                   track_height = opt$track_height,
                   gap          = opt$track_gap,
                   base_start   = 12)

    pp$data1max <- max(track_pos$r1) + opt$ideogram_gap
    pp$data2max <- 0
    if (percent_mode) {
      pp$bottommargin <- pp$bottommargin + (if (legend_rows > 0) 14 + legend_rows * 5.5 else 0)
    }
  }

   draw_once <- function() {
    genome_lengths_plot <- if (percent_mode) {
      setNames(rep(100, length(chromosomes)), chromosomes)
    } else {
      chr_lengths_all[chromosomes]
    }
    genome_df <- data.frame(chr = chromosomes, start = 0,
                             end = as.numeric(genome_lengths_plot))

    kp <- plotKaryotype(
      plot.type   = 1,
      genome      = genome_df,
      plot.params = pp,
      cex         = chr_cex_use,
      cytobands   = cytobands_plot,
      chromosomes = chromosomes
    )

    for (chr in chromosomes) {
      if (percent_mode) {
        kpAddBaseNumbers_manual(kp, chr, 100,
                                tick.dist    = 10, label.dist = 20,
                                cex          = tick_cex,
                                tick.y       = tick_y_offset,
                                label.y      = label_y_offset,
                                unit_divisor = 1, unit_suffix = "%")
      } else {
        kpAddBaseNumbers_manual(kp, chr, chr_lengths_all[chr],
                                tick.dist  = opt$minor_tick,
                                label.dist = opt$major_tick,
                                cex        = tick_cex,
                                tick.y     = tick_y_offset,
                                label.y    = label_y_offset)
      }
    }

    for (i in seq_len(n_tracks)) {
      lbl <- grp_sample_df$label[i]
      dat <- grp_bed_list_plot[[grp_sample_df$row_id[i]]]
      r0  <- track_pos$r0[i]
      r1  <- track_pos$r1[i]
      kpPlotRegions(kp, data.panel = 1, r0 = r0, r1 = r1,
                    data = dat, col = dat$HEX)
      kpAddLabels(kp, labels = lbl, r0 = r0, r1 = r1,
                  data.panel = 1, cex = label_cex_use,
                  side = "left", label.margin = 0.06)
    }

    if (!is.null(legend_df) && nrow(legend_df) > 0) {
      if (is.null(opt$chromosome) && !percent_mode) {
        n_chr_total <- length(chromosomes)
        win         <- pick_legend_row_window(chr_lengths_all, chromosomes, span = 5)
        row_top     <- win[1]
        row_bottom  <- win[2]
        y_top    <- 1 - (row_top - 1) / n_chr_total
        y_bottom <- 1 - row_bottom / n_chr_total
        span     <- y_top - y_bottom
        lx <- grconvertX(0.60, from = "ndc", to = "user")
        ly <- grconvertY(y_top - 0.12 * span, from = "ndc", to = "user")
        xjust <- 0; yjust <- 1
      } else {
        lx <- grconvertX(0.985, from = "ndc", to = "user")
        ly <- grconvertY(0.015, from = "ndc", to = "user")
        xjust <- 1; yjust <- 0
      }
      legend(x = lx, y = ly, xjust = xjust, yjust = yjust,
             legend = legend_df$label,
             fill   = legend_df$colour,
             cex    = opt$legend_cex,
             bty    = "n",
             ncol   = legend_ncol,
             xpd    = NA,
             title  = paste0("CENP-B motif colour key - ", group_friendly))
    }
    invisible(kp)
  }

  open_device <- function(path, this_fmt) {
    if (this_fmt == "tiff") {
      tiff_w <- 15
      tiff_h <- if (!is.null(opt$chromosome)) 4 else 24
      tiff(path, res = 1200, width = tiff_w, height = tiff_h, units = "in",
           compression = "lzw", antialias = "subpixel")
    } else if (this_fmt == "png") {
      png(path, res = opt$jpeg_res, width = opt$jpeg_width,
          height = opt$jpeg_height, units = "px", type = "cairo")
    } else if (this_fmt == "svg") {
      svg(path, width = opt$plot_width, height = opt$plot_height)
    } else if (this_fmt == "pdf") {
      pdf(path, width = opt$plot_width, height = opt$plot_height,
          useDingbats = FALSE)
    }
  }

  out_base_grp <- paste0(output_base, "_", group_tag)
  want_html <- (fmt == "html") || isTRUE(opt$html)

  # ---- primary output (--format), unless format is html -------------------
  primary_path <- NA_character_
  if (fmt != "html") {
    primary_path <- paste0(out_base_grp, ".", fmt)
    cat(sprintf("  Opening %s device -> %s\n", toupper(fmt), primary_path))
    open_device(primary_path, fmt)
    draw_once()
    dev.off()

    if (fmt == "svg") {
      svgz_path <- paste0(out_base_grp, ".svgz")
      svg_bytes <- readBin(primary_path, what = "raw",
                            n = file.info(primary_path)$size)
      con <- gzfile(svgz_path, open = "wb", compression = 9)
      writeBin(svg_bytes, con)
      close(con)
      file.remove(primary_path)
      primary_path <- svgz_path
      cat(sprintf("  SVG compressed: %s (%.1f MB)\n", basename(primary_path),
                  file.info(primary_path)$size / 1e6))
    } else {
      cat(sprintf("  Plot saved: %s\n", primary_path))
    }
  }

  # ---- interactive canvas HTML
  html_path <- NA_character_
  if (want_html) {
    html_path <- paste0(out_base_grp, ".html")
    cat(sprintf("  Building interactive HTML -> %s\n", html_path))
    build_interactive_html(grp_sample_df, grp_bed_list, chromosomes,
                            chr_lengths_all, cytobands, legend_df,
                            group_friendly, html_path,
                            max_points = opt$html_max_points,
                            coord_mode = opt$coord_mode,
                            sample_chr_lengths = sample_chr_lengths)
    cat(sprintf("  HTML saved: %s (%.1f MB)\n", html_path,
                file.info(html_path)$size / 1e6))
  }
  if (fmt == "html") primary_path <- html_path

  # ---- quantitative metrics ------------------------------------------------
  chr_len_lookup <- setNames(as.numeric(chr_lengths_all[chromosomes]), chromosomes)

  metrics_list <- lapply(seq_len(n_tracks), function(i) {
    lbl <- grp_sample_df$label[i]
    dat <- grp_bed_list[[grp_sample_df$row_id[i]]]

    do.call(rbind, lapply(chromosomes, function(chr) {

      sub     <- dat[dat$V1 == chr, ]
      chr_len <- chr_len_lookup[chr]
      model_lbl <- if ("MODEL" %in% names(dat) && nrow(dat) > 0) dat$MODEL[1] else NA_character_

      if (nrow(sub) == 0) {
        return(data.frame(
          sample            = lbl,
          model             = model_lbl,
          chromosome        = chr,
          n_blocks          = 0L,
          total_length_bp   = 0L,
          mean_block_size   = NA_real_,
          median_block_size = NA_real_,
          span_bp           = NA_real_,
          centromere_start  = NA_real_,
          centromere_end    = NA_real_,
          centromere_mid    = NA_real_,
          chr_length_bp     = chr_len,
          rel_position      = NA_real_,
          dominant_colour   = NA_character_,
          n_colours         = 0L,
          stringsAsFactors  = FALSE
        ))
      }

      block_sizes   <- sub$V3 - sub$V2
      cen_start     <- min(sub$V2)
      cen_end       <- max(sub$V3)
      cen_mid       <- (cen_start + cen_end) / 2
      colour_counts <- sort(table(sub$HEX), decreasing = TRUE)

      data.frame(
        sample            = lbl,
        model             = model_lbl,
        chromosome        = chr,
        n_blocks          = nrow(sub),
        total_length_bp   = sum(block_sizes),
        mean_block_size   = round(mean(block_sizes), 1),
        median_block_size = round(median(block_sizes), 1),
        span_bp           = cen_end - cen_start,
        centromere_start  = cen_start,
        centromere_end    = cen_end,
        centromere_mid    = round(cen_mid, 0),
        chr_length_bp     = chr_len,
        rel_position      = round(cen_mid / chr_len, 4),
        dominant_colour   = names(colour_counts)[1],
        n_colours         = length(colour_counts),
        stringsAsFactors  = FALSE
      )
    }))
  })

  metrics_df   <- do.call(rbind, metrics_list)
  metrics_path <- paste0(out_base_grp, "_metrics.tsv")
  fwrite(metrics_df, metrics_path, sep = "\t", quote = FALSE, na = "NA")
  cat(sprintf("  Per-sample metrics       -> %s\n", metrics_path))

  summary_df <- do.call(rbind, lapply(chromosomes, function(chr) {
    sub <- metrics_df[metrics_df$chromosome == chr & !is.na(metrics_df$centromere_mid), ]
    if (nrow(sub) < 2) return(NULL)
    mids <- sub$centromere_mid / 1e6
    data.frame(
      chromosome        = chr,
      chr_length_Mb     = round(chr_len_lookup[chr] / 1e6, 2),
      n_samples         = nrow(sub),
      mean_mid_Mb       = round(mean(mids),  3),
      sd_mid_Mb         = round(sd(mids),    3),
      cv_mid            = round(sd(mids) / mean(mids), 4),
      range_mid_Mb      = round(diff(range(mids)), 3),
      mean_total_len_Mb = round(mean(sub$total_length_bp / 1e6), 3),
      sd_total_len_Mb   = round(sd(sub$total_length_bp   / 1e6), 3),
      mean_n_blocks     = round(mean(sub$n_blocks), 1),
      stringsAsFactors  = FALSE
    )
  }))

  summary_path <- paste0(out_base_grp, "_cross_sample_summary.tsv")
  fwrite(summary_df, summary_path, sep = "\t", quote = FALSE, na = "NA")
  cat(sprintf("  Cross-sample summary     -> %s\n\n", summary_path))

  # ---- optional clustering analysis (heatmap + dendrogram) ----------------
  cluster_out <- NULL
  if (isTRUE(opt$cluster)) {
    cluster_out <- run_clustering_analysis(grp_sample_df, grp_bed_list, chromosomes,
                                            legend_df, group_friendly, out_base_grp)
  }

  list(plot = primary_path, html = html_path,
       metrics = metrics_path, summary = summary_path,
       cluster_pdf  = if (!is.null(cluster_out)) cluster_out$pdf  else NA_character_,
       cluster_html = if (!is.null(cluster_out)) cluster_out$html else NA_character_)
}

# =============================================================================
#  8.  RUN ONE render_group() PER MODEL GROUP
# =============================================================================

all_outputs <- list()
for (gm in names(groups)) {
  idxs <- groups[[gm]]
  grp_sample_df     <- sample_df[idxs, , drop = FALSE]
  grp_bed_list      <- bed_list[grp_sample_df$row_id]
  grp_legend_entries <- legend_entries[grp_sample_df$row_id]
  names(grp_bed_list) <- names(grp_legend_entries) <- grp_sample_df$row_id

  # Model 3: recolour every not-yet-colour-coded file in this group
  # TOGETHER (not one file at a time), so the same #hashID pattern gets
  # the same colour everywhere
  recolor_rows <- intersect(grp_sample_df$row_id, pending_kpattern_recolor)
  if (length(recolor_rows) > 0) {
    rc <- apply_group_kpattern_colouring(grp_bed_list, recolor_rows,
                                          top_n = opt$top_n_patterns)
    grp_bed_list <- rc$bed_list
    for (rid in recolor_rows) grp_legend_entries[[rid]] <- NULL
    if (!is.null(rc$legend)) grp_legend_entries[[recolor_rows[1]]] <- rc$legend
  }

  group_tag      <- unname(group_map[gm])
  group_friendly <- unname(group_friendly_map[gm])

  all_outputs[[group_tag]] <- render_group(grp_sample_df, grp_bed_list,
                                            grp_legend_entries,
                                            group_tag, group_friendly)
}

# =============================================================================
#  9.  SUMMARY
# =============================================================================

cat("================ ALL OUTPUTS ================\n")
for (gt in names(all_outputs)) {
  o <- all_outputs[[gt]]
  if (!is.na(o$html) && identical(o$plot, o$html)) {
    cat(sprintf("[%s]\n  html:    %s\n  metrics: %s\n  summary: %s\n",
                gt, o$html, o$metrics, o$summary))
  } else {
    cat(sprintf("[%s]\n  plot:    %s\n  html:    %s\n  metrics: %s\n  summary: %s\n",
                gt, o$plot,
                if (is.na(o$html)) "(disabled, --html false)" else o$html,
                o$metrics, o$summary))
  }
  if (!is.na(o$cluster_pdf)) {
    cat(sprintf("  cluster: %s\n  cluster: %s\n", o$cluster_pdf, o$cluster_html))
  }
  cat("\n")
}

cat(sprintf("Done! All outputs written to:\n  %s\n\n", dirname(output_base)))
