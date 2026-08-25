# =============================================================================
#  CenTeny Karyotype Plot — Snakemake Pipeline
# =============================================================================
#
#  HOW TO RUN
#  ----------
#  Dry run:
#    snakemake --cores 1 --dry-run
#
#  Full run:
#    snakemake --cores 1
#
# =============================================================================

import os

configfile: "config.yaml"

def parse_sample_sheet(tsv_path):
    samples = []
    with open(tsv_path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) not in (2, 3):
                raise ValueError(
                    f"Malformed sample sheet line (expected 2 or 3 tab-separated columns):\n"
                    f"  '{line}'\nEach row must be: <label><TAB><bed_path>[<TAB><model>]\n"
                    f"where <model> is one of: centeny, model1, model2, model3, auto"
                )
            label, bed = parts[0].strip(), parts[1].strip()
            model_tag  = parts[2].strip() if len(parts) == 3 else "auto"
            samples.append((label, bed, model_tag))
    if not samples:
        raise ValueError(f"No samples found in '{tsv_path}'.")
    return samples

SAMPLE_SHEET = os.path.join(config["workdir"], config["sample_sheet"])
SAMPLES      = parse_sample_sheet(SAMPLE_SHEET)

OUTPUT_BASE = os.path.splitext(config.get("output_file", "centeny_plot"))[0]
OUTPUT_DIR  = os.path.join(config["workdir"], "centeny_outputs")

# ---------------------------------------------------------------------------
rule all:
    input:
        [OUTPUT_DIR] + ([os.path.join(config["workdir"], "whole_genome_motif_abundance.html"),
                         os.path.join(config["workdir"],
                                      config["motif_report_output"] + "_karyotypes.html")]
                         if config.get("motif_report_output") else [])

# ---------------------------------------------------------------------------
rule validate_inputs:
    input:
        cytoband     = config["cytoband"],
        bed_files    = [s[1] for s in SAMPLES],
        sample_sheet = SAMPLE_SHEET
    output:
        touch(os.path.join(config["workdir"], ".inputs_validated"))
    run:
        missing = []
        if not os.path.isfile(input.cytoband):
            missing.append(f"  Cytoband file : {input.cytoband}")
        for label, bed, model_tag in SAMPLES:
            if not os.path.isfile(bed):
                missing.append(f"  Bed file for '{label}' ({model_tag}): {bed}")
        if missing:
            raise FileNotFoundError(
                "\n\nThe following required files were not found:\n"
                + "\n".join(missing) + "\n\nCheck your config.yaml and samples.tsv.\n"
            )
        print(f"\n✔ All {len(SAMPLES)} bed files and cytoband file found.\n")

# ---------------------------------------------------------------------------
rule plot_centeny:
    input:
        validated    = os.path.join(config["workdir"], ".inputs_validated"),
        cytoband     = config["cytoband"],
        sample_sheet = SAMPLE_SHEET,
        bed_files    = [s[1] for s in SAMPLES]
    output:
        outdir = directory(OUTPUT_DIR)
    conda:
        "envs/centeny_r.yaml"
    params:
        workdir         = config["workdir"],
        output_base     = os.path.join(OUTPUT_DIR, OUTPUT_BASE),
        output_format   = config.get("output_format",    "pdf"),
        plot_width      = config.get("plot_width",        20),
        plot_height     = config.get("plot_height",       13),
        jpeg_res        = config.get("jpeg_res",          1200),
        chr_cex         = config.get("chr_cex",           0.9),
        label_cex       = config.get("label_cex",         0.45),
        tick_cex        = config.get("tick_cex",          0.5),
        major_tick      = config.get("major_tick_dist",   10000000),
        minor_tick      = config.get("minor_tick_dist",   5000000),
        model           = config.get("model",              "auto"),
        top_n_patterns  = config.get("top_n_patterns",     10),
        legend_cex      = config.get("legend_cex",         0.6),
        legend_max      = config.get("legend_max_entries", 30),
        ideogram_gap    = config.get("ideogram_gap",       15),
        html            = config.get("html",               "true"),
        html_max_points = config.get("html_max_points",    150000),
        coord_mode      = config.get("coord_mode",          "bp"),
        genome_sizes    = config.get("genome_sizes",        ""),
        cluster         = config.get("cluster",             "false"),
    log:
        os.path.join(config["workdir"], "logs", "plot_centeny.log")
    shell:
        """
        mkdir -p {params.workdir}/logs
        mkdir -p {output.outdir}

        Rscript scripts/centeny_plot_tool.R \
            --workdir            "{params.workdir}"      \
            --cytoband           "{input.cytoband}"      \
            --sample_sheet       "{input.sample_sheet}"  \
            --output             "{params.output_base}"  \
            --format             {params.output_format}  \
            --plot_width         {params.plot_width}     \
            --plot_height        {params.plot_height}    \
            --jpeg_res           {params.jpeg_res}       \
            --chr_cex            {params.chr_cex}        \
            --label_cex          {params.label_cex}      \
            --tick_cex           {params.tick_cex}       \
            --major_tick         {params.major_tick}     \
            --minor_tick         {params.minor_tick}     \
            --model              {params.model}          \
            --top_n_patterns     {params.top_n_patterns} \
            --legend_cex         {params.legend_cex}     \
            --legend_max_entries {params.legend_max}     \
            --ideogram_gap       {params.ideogram_gap}   \
            --html               {params.html}           \
            --html_max_points    {params.html_max_points} \
            --coord_mode         {params.coord_mode}     \
            --genome_sizes       "{params.genome_sizes}" \
            --cluster            {params.cluster}        \
            2>&1 | tee "{log}"
        """

# =============================================================================
#  Optional: motif model composition report (whole_genome_motif_abundance.html
#  + <output>_karyotypes.html)
#  ------------------------------------------------------------------------
#  Runs whenever config.yaml sets `motif_report_output` to a non-empty base
#  name. Purely BED-based (uses the same sample_sheet as the main pipeline)
#  — no FASTA scanning in this release, so no Python env is needed. See
#  scripts/motif_report.R's own header for details.
# =============================================================================

if config.get("motif_report_output"):

    rule motif_report:
        input:
            sample_sheet = SAMPLE_SHEET,
        output:
            abundance = os.path.join(config["workdir"], "whole_genome_motif_abundance.html"),
            karyotype = os.path.join(config["workdir"],
                                      config["motif_report_output"] + "_karyotypes.html")
        conda:
            "envs/centeny_r.yaml"
        params:
            workdir         = config["workdir"],
            output_base     = os.path.join(config["workdir"], config["motif_report_output"]),
            comp_threshold  = config.get("motif_composition_threshold", 0.5),
            kpattern_top_n  = config.get("motif_kpattern_top_n",  10),
        log:
            os.path.join(config["workdir"], "logs", "motif_report.log")
        shell:
            """
            mkdir -p {params.workdir}/logs

            Rscript scripts/motif_report.R \
                --workdir        "{params.workdir}"       \
                --sample_sheet   "{input.sample_sheet}"   \
                --output         "{params.output_base}"   \
                --composition_threshold {params.comp_threshold} \
                --kpattern_top_n {params.kpattern_top_n} \
                2>&1 | tee "{log}"
            """

