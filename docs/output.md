# nf-core/drop: Output

## Introduction

This document describes the output produced by the pipeline. The directories listed below will be created in the result directory after the pipeline has finished. All paths are relative to the top-level result directory.

## Pipeline overview

Besides result tables, the pipeline uses MultiQC to generate comprehensive reports for each {annotation}\_{drop_group} combination within every subworkflow.

### Aberrant Expression

- HTML report
  at the resulting `reports/aberrant_expression/genecounts` and `reports/aberrant_expression/outrider`
  - Counts Summaries for each aberrant expression group
    - number of local and external samples
    - Mapped reads and size factors for each sample
    - histograms showing the mean count distribution with different conditions
    - expressed genes within each sample and as a dataset
  - Outrider Summaries for each aberrant expression group
    - aberrantly expressed genes per sample
    - correlation between samples before and after the autoencoder
    - biological coefficient of variation
    - aberrant samples
    - results table
- Files for each aberrant expression group
  - OUTRIDER datasets
    - Follow the [OUTRIDER vignette](https://www.bioconductor.org/packages/devel/bioc/vignettes/OUTRIDER/inst/doc/OUTRIDER.pdf) for individual OUTRIDER object file (`processed_results/aberrant_expression/{annotation}/outrider/{drop_group}/*.ods`) analysis.
  - Results tables
    - file `processed_results/aberrant_expression/{annotation}/outrider/{drop_group}/OUTRIDER_results.tsv` this text file contains only the significant genes and samples that meet the cutoffs defined in the config file for `padjCutoff` and `zScoreCutoff`.
    - file `processed_results/aberrant_expression/{annotation}/outrider/{drop_group}/OUTRIDER_results_all.Rds` contains the entire OUTRIDER results table regardless of significance.

### Aberrant Splicing

- HTML report
  at the resulting `reports/aberrant_splicing/splicecounts` and `reports/aberrant_splicing/fraser`
  - Counting Summaries for each aberrant splicing group
    - number of local and external samples
    - number introns/splice sites before and after merging
    - comparison of local and external mean counts
    - histograms showing the junction expression before and after filtering and variability
  - FRASER Summaries for each aberrant splicing group
    - the number of samples, introns, and splice sites
    - correlation between samples before and after the autoencoder
    - results table
- Files for each aberrant splicing group
  - FRASER datasets (fds)
    - Follow the [FRASER vignette](https://www.bioconductor.org/packages/devel/bioc/vignettes/FRASER/inst/doc/FRASER.pdf) for individual FRASER object file (fds) analysis.
  - Results tables
    - file `processed_results/aberrant_splicing/results/{annotation}/fraser/{drop_group}/results_per_junction.tsv` this text file contains only significant junctions that meet the cutoffs defined in the config file.
    - file `processed_results/aberrant_splicing/results/{annotation}/fraser/{drop_group}/results.tsv` contains only significant junctions that meet the cutoffs defined in the config file, aggregated at the gene level. Any sample/gene pair is represented by only the most significant junction.

### Mono-allelic Expression

HTML report
at the resulting `reports/mae/mae` and `reports/mae/maeqc`

- Results for each mae group
  - number of samples, genes, and mono-allelically expressed heterozygous SNVs
  - a cascade plot that shows additional filters
  - histogram of inner cohort frequency
  - summary of the cascade plot and results table
- Quality Control
  - QC Overview
    - For each mae group QC checks for DNA/RNA matching

Additionally the `mae` subworkflow creates the following files:

- `processed_results/mae/{drop_group}/MAE_results_all_{annotation}.tsv.gz`
  - this file contains the MAE results of all heterozygous SNVs regardless of significance
- `processed_results/mae/{drop_group}/MAE_results_{annotation}.tsv`
  - this is the file linked in the HTML document and described above
- `processed_results/mae/{drop_group}/MAE_results_{annotation}_rare.tsv`
  - this file is a subset of `MAE_results_{annotation}.tsv` with only the variants that pass the allele frequency cutoffs. If `add_AF` is set to `true` in config file must meet minimum AF set by `max_AF`. Additionally, the inner-cohort frequency must meet the `maxVarFreqCohort` cutoff

### MultiQC

<details markdown="1">
<summary>Output files</summary>

- `multiqc/`
  - `multiqc_report.html`: a merged report for all subworkflows.
  - `multiqc_data/`: directory containing parsed statistics from the different tools used in the pipeline.
  - `multiqc_plots/`: directory containing static images from the report in various formats.

</details>

[MultiQC](http://multiqc.info) is a visualization tool that generates a single HTML report summarising all samples in your project. Most of the pipeline QC results are visualised in the report and further statistics are available in the report data directory.

Results generated by MultiQC collate pipeline QC from supported tools e.g. FastQC. The pipeline has special steps which also allow the software versions to be reported in the MultiQC output for future traceability. For more information about how to use MultiQC reports, see <http://multiqc.info>.

### Pipeline information

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
  - Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt` and `software_versions.yml`. The `pipeline_report*` files will only be present if the `--email` / `--email_on_fail` parameter's are used when running the pipeline.
  - Reformatted samplesheet files used as input to the pipeline: `samplesheet.valid.csv`.
  - Parameters used by the pipeline run: `params.json`.

</details>

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.
