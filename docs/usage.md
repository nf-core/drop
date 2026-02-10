# nf-core/drop: Usage

## :warning: Please read this documentation on the nf-core website: [https://nf-co.re/drop/usage](https://nf-co.re/drop/usage)

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._

## License notice

OUTRIDER and FRASER are released under CC-BY-NC 4.0, meaning a license is required for any commercial use. If you intend to use the aberrant expression and aberrant splicing modules for commercial purposes, please contact the authors: Julien Gagneur (gagneur [at] in.tum.de), Christian Mertes (mertes [at] in.tum.de), and Vicente Yepez (yepez [at] in.tum.de).

## Introduction

nf-core/drop allows controlling which subworkflows to run. By default, all subworkflows will run. You can skip subsworkflow via parameters (`--ae_skip` (Aberrant Expression), `--as_skip` (Aberrant Splicing), `--mae_skip` (Mono-Allelic Expression)). We describe different global and module-specific parameters in the [parameter documentation](https://nf-co.re/drop/parameters).

## Samplesheet input

> For a detailed explanation of the columns of the samplesheet, please refer to Box 3 of the [DROP manuscript](https://www.nature.com/articles/s41596-020-00462-5#Sec26). Some information has been updated since publication, please use this documentation as the preferred syntax/formatting.

Each row of the sample annotation table corresponds to a unique RNA sample. The required columns are `RNA_ID`, `RNA_BAM_FILE` and `DROP_GROUP`. The following columns describe the RNA-seq experimental setup: `STRAND`(mandatory column), `PAIRED_END`, `COUNT_MODE` and `COUNT_OVERLAPS`. They affect the counting procedures of the aberrant expression and splicing subworkflows. For a detailed explanation, refer to the documentation of [HTSeq](https://htseq.readthedocs.io/en/latest/).

### DROP_GROUP

DROP_GROUP are the analysis groups that the RNA assay belongs to. Multiple groups must be separated by commas and no spaces (e.g., blood,WES,groupA). Together with the parameters `--ae_groups`, `--as_groups`, and `--mae_groups` to specify which groups to run in each subworkflow, it allows you to perform different analyses (for example, across different tissue groups) within a single run. All groups will be analysed in their respective subworkflow if the `groups` parameters are missing.

| RNA_ID  | RNA_BAM_FILE        | RNA_BAI_FILE\*          | DROP_GROUP | PAIRED_END | COUNT_MODE         | COUNT_OVERLAPS | STRAND |
| ------- | ------------------- | ----------------------- | ---------- | ---------- | ------------------ | -------------- | ------ |
| HG00103 | path/to/HG00103.bam | path/to/HG00103.bam.bai | blood,all  | TRUE       | IntersectionStrict | TRUE           | no     |

<sub>\* You can provide the BAM index file in the `RNA_BAI_FILE` column. If this column is not specified, the pipeline will automatically generate index files from the BAM files.<sub>

### MAE

To run the MAE subworkflow, the columns `DNA_ID`, `DNA_VCF_FILE` and `GENOME`are needed. MAE can not be run in samples using external counts as we need to use the `RNA_BAM_FILE` to count reads supporting each allele of the heterozygous variants found in the `DNA_VCF_FILE`. You also need to provide the parameter `--ucsc_fasta <path/to/fasta>` and/or `--ncbi_fasta <path/to/fasta>`, depending on the value specified in the `GENOME` column.

| RNA_ID  | RNA_BAM_FILE              | RNA_BAI_FILE                  | DNA_ID  | DNA_VCF_FILE                    | DNA_TBI_FILE\*                      | DROP_GROUP  | STRAND | GENOME |
| ------- | ------------------------- | ----------------------------- | ------- | ------------------------------- | ----------------------------------- | ----------- | ------ | ------ |
| HG00096 | /path/to/HG00096_ncbi.bam | /path/to/HG00096_ncbi.bam.bai | HG00096 | /path/to/demo_chr21_ncbi.vcf.gz | /path/to/demo_chr21_ncbi.vcf.gz.tbi | mae,batch_0 | no     | ncbi   |
| HG00103 | /path/to/HG00103.bam      | /path/to/HG00103.bam.bai      | HG00103 | /path/to/demo_chr21.vcf.gz      | /path/to/demo_chr21.vcf.gz.tbi      | mae,batch_1 | no     | ucsc   |

<sub>\* You can provide the Tabix index file in the `DNA_TBI_FILE` column. If this column is not specified, the pipeline will automatically generate index files from the VCF files.<sub>

The columns `DNA_ID`, `DNA_VCF_FILE` and `GENOME` can be empty for the AE and AS subworkflow. For example, you can run `--ae_groups group1 --as_groups group1 --mae_groups group2` with the samplesheet below.

| RNA_ID  | RNA_BAM_FILE        | RNA_BAI_FILE            | DROP_GROUP    | STRAND | DNA_ID  | DNA_VCF_FILE              | DNA_TBI_FILE                  | GENOME |
| ------- | ------------------- | ----------------------- | ------------- | ------ | ------- | ------------------------- | ----------------------------- | ------ |
| HG00103 | path/to/HG00103.bam | path/to/HG00103.bam.bai | group1,group2 | no     | HG00103 | path/to/demo_chr21.vcf.gz | path/to/demo_chr21.vcf.gz.tbi | ucsc   |
| HG00106 | path/to/HG00106.bam | path/to/HG00106.bam.bai | group1,group2 | no     | HG00106 | path/to/demo_chr21.vcf.gz | path/to/demo_chr21.vcf.gz.tbi | ucsc   |
| HG00111 | path/to/HG00111.bam | path/to/HG00111.bam.bai | group1        |        |         |                           |                               |        |

### Using External Counts

DROP can utilize external counts for the `aberrantExpression` and `aberrantSplicing` subworkflows which can enhance the statistical power of these subworkflows by providing more samples from which we can build a distribution of counts and detect outliers. However, this process introduces some particular issues that need to be addressed to make sure it is a valuable addition to the experiment.

In case external counts are included, add a new row for each sample from those files (or a subset if not all samples are needed). Add the columns: `GENE_COUNTS_FILE`
(for aberrant expression), `GENE_ANNOTATON`, and `SPLICE_COUNTS_DIR` (for aberrant splicing). These columns should remain empty for samples processed locally (from `RNA_BAM`).

#### Aberrant Expression

To use external counts for aberrant expression, you need to use the exact same gene annotation for each external sample as well as using the same gene annotation file specified in `pramas.gene_annotation`. This is to avoid potential mismatching on counting, 2 different gene annotations could drastically affect which reads are counted in which region drastically skewing the results.

The user must also take special consideration when building the sample annotation table. Samples using external counts need only `RNA_ID` which must exactly match the column header in the external count file `DROP_GROUP`, `GENE_COUNTS_FILE`, and `GENE_ANNOTATION` which must contain the exact key specified in the config. The other columns should remain empty.

> tips: Using `exportCounts` generates the sharable `GENE_COUNTS_FILE` file in the appropriate `<OUTDIR>/processed_results/exported_counts/` sub-directory.

| RNA_ID  | RNA_BAM_FILE         | RNA_BAI_FILE             | DROP_GROUP                 | STRAND | GENE_COUNTS_FILE                                                                                                                          | GENE_ANNOTATION |
| ------- | -------------------- | ------------------------ | -------------------------- | ------ | ----------------------------------------------------------------------------------------------------------------------------------------- | --------------- |
| HG00103 | /path/to/HG00103.bam | /path/to/HG00103.bam.bai | outrider,outrider_external | no     |                                                                                                                                           |                 |
| HG00178 |                      |                          | outrider_external          | no     | /path/to/[geneCounts.tsv.gz](https://github.com/gagneurlab/drop_demo_data/raw/refs/heads/main/Data/external_count_data/geneCounts.tsv.gz) | v29             |

#### Aberrant Splicing

Using external counts for aberrant splicing reduces the number of introns processed to only those that are exactly the same between the local and external junctions. Because rare junctions may be personally identifiable the `exportCounts` command only exports regions canonically mentioned in the gtf file. As a result, when merging the external counts with the local counts we only match introns that are **exact** between the 2 sets, this is to ensure that if a region is missing we don't introduce 0 counts into the distribution calculations.

The user must also use special consideration when building the sample annotation table. Samples using external counts need only `RNA_ID` which must exactly match the column header in the external count file `DROP_GROUP`, and `SPLICE_COUNTS_DIR`. `SPLICE_COUNTS_DIR` is the directory containing the set of 5 needed count files. The other columns should remain empty.

> tips: Using `exportCounts` generates the necessary files in the appropriate `<OUTDIR>/processed_results/exported_counts/` sub-directory

`SPLICE_COUNTS_DIR` should contain the following:

- k_j_counts.tsv.gz
- k_theta_counts.tsv.gz
- n_psi3_counts.tsv.gz
- n_psi5_counts.tsv.gz
- n_theta_counts.tsv.gz

| RNA_ID  | RNA_BAM_FILE         | DROP_GROUP             | STRAND | SPLICE_COUNTS_DIR                                                                                                                          |
| ------- | -------------------- | ---------------------- | ------ | ------------------------------------------------------------------------------------------------------------------------------------------ |
| HG00176 | /path/to/HG00176.bam | fraser,fraser_external | no     |                                                                                                                                            |
| HG00191 |                      | fraser_external        | no     | /path/to/[external_count_data.tar.gz](https://github.com/nf-core/test-datasets/raw/refs/heads/drop/data/inputs/external_count_data.tar.gz) |

#### Another External count examples

This example will use the `DROP_GROUP` BLOOD_AE for the aberrant expression module (containing S10R, EXT-1R, EXT-2R) and
the `DROP_GROUP` BLOOD_AS for the aberrant expression module (containing S10R, EXT-2R, EXT-3R)

| RNA_ID | DNA_ID | DROP_GROUP        | RNA_BAM_FILE      | GENE_COUNTS_FILE               | GENE_ANNOTATION | SPLICE_COUNTS_DIR         |
| ------ | ------ | ----------------- | ----------------- | ------------------------------ | --------------- | ------------------------- |
| S10R   | S10G   | BLOOD_AE,BLOOD_AS | /path/to/S10R.BAM |                                |                 |                           |
| EXT-1R |        | BLOOD_AE          |                   | /path/to/externalCounts.tsv.gz | gencode34       |                           |
| EXT-2R |        | BLOOD_AE,BLOOD_AS |                   | /path/to/externalCounts.tsv.gz | gencode34       | /path/to/externalCountDir |
| EXT-3R |        | BLOOD_AS          |                   |                                |                 | /path/to/externalCountDir |

### Full samplesheet

An [example samplesheet](../assets/samplesheet.tsv) has been provided with the pipeline.

| Column              | Description                                                                                                                                                                                          |
| ------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `RNA_ID`            | Unique identifier for the RNA sample.                                                                                                                                                                |
| `RNA_BAM_FILE`      | Path to the RNA BAM file. Can also be a CRAM file, but BAM files are recommended since the pipeline will convert the CRAM files to BAM.                                                              |
| `RNA_BAI_FILE`      | Path to the RNA BAM index file. It will be automatically generated from the BAM files if not given.                                                                                                  |
| `DNA_ID`            | Identifier for the DNA sample. It shall be the sample header in the VCF file.                                                                                                                        |
| `DNA_VCF_FILE`      | Path to the DNA VCF file.                                                                                                                                                                            |
| `DNA_TBI_FILE`      | Path to the DNA VCF index file. It will be automatically generated from the VCF files if not given.                                                                                                  |
| `DROP_GROUP`        | See [above](#drop_group).                                                                                                                                                                            |
| `PAIRED_END`        | Indicates if the input is paired-end or single-end. Default: `true` (paired-end). Refer to the documentation of [HTSeq](https://htseq.readthedocs.io/en/latest/).                                    |
| `COUNT_MODE`        | Count mode. Default: `IntersectionStrict`. Options: `union`, `IntersectionStrict`, `IntersectionNotEmpty`. Refer to the documentation of [HTSeq](https://htseq.readthedocs.io/en/latest/).           |
| `COUNT_OVERLAPS`    | Indicates if overlaps should be counted. Default: `true`. Refer to the documentation of [HTSeq](https://htseq.readthedocs.io/en/latest/).                                                            |
| `STRAND`            | Samples within each `DROP_GROUP` should either be stranded (`yes`, `reverse`, or a combination of `yes` and `reverse`) or unstranded (only `no`), and this can vary between different `DROP_GROUP`s. |
| `HPO_TERMS`         | Comma-separated list of HPO terms associated with the sample.                                                                                                                                        |
| `GENE_COUNTS_FILE`  | Path to the gene counts file (`.tsv` or `.tsv.gz`). See details also [above](#aberrant-expression).                                                                                                  |
| `GENE_ANNOTATION`   | Gene annotation in YAML format, e.g. [assets/gene_annotation.yaml](../assets/gene_annotation.yaml).                                                                                                  |
| `GENOME`            | See details [above](#mae). Default: `ucsc`. Options: `ncbi`, `ucsc`.                                                                                                                                 |
| `SPLICE_COUNTS_DIR` | Path to the splice counts directory. See details [above](#aberrant-splicing).                                                                                                                        |
| `SEX`               | Sex of the sample. Samples of `m`, `male`, `f`, `femal` are analysed for sex bias aberrant expression report.                                                                                        |

## Running the pipeline

The typical command for running the pipeline is as follows:

```bash
nextflow run nf-core/drop --input ./samplesheet.tsv --outdir ./results --genome hg19 --gene_annotation ./gene_annotation.yaml -profile docker
```

This will launch the pipeline with the `docker` configuration profile. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```bash
work                # Directory containing the nextflow working files
<OUTDIR>            # Finished results in specified location (defined with --outdir)
.nextflow_log       # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

If you wish to repeatedly use the same parameters for multiple runs, rather than specifying each flag in the command, you can specify these in a params file.

Pipeline settings can be provided in a `yaml` or `json` file via `-params-file <file>`.

> [!WARNING]
> Do not use `-c <file>` to specify parameters as this will result in errors. Custom config files specified with `-c` must only be used for [tuning process resource specifications](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources), other infrastructural tweaks (such as output directories), or module arguments (args).

The above pipeline run specified with a params file in yaml format:

```bash
nextflow run nf-core/drop -profile docker -params-file params.yaml
```

with:

```yaml title="params.yaml"
input: './samplesheet.tsv'
outdir: './results/'
genome: 'hg19'
<...>
```

You can also generate such `YAML`/`JSON` files via [nf-core/launch](https://nf-co.re/launch).

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull nf-core/drop
```

### Reproducibility

It is a good idea to specify the pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [nf-core/drop releases page](https://github.com/nf-core/drop/releases) and find the latest pipeline version - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`. Of course, you can switch to another version by changing the number after the `-r` flag.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future. For example, at the bottom of the MultiQC reports.

To further assist in reproducibility, you can use share and reuse [parameter files](#running-the-pipeline) to repeat pipeline runs with the same settings without having to write out a command with every single parameter.

> [!TIP]
> If you wish to share such profile (such as upload as supplementary material for academic publications), make sure to NOT include cluster specific paths to files, nor institutional specific profiles.

## Core Nextflow arguments

> [!NOTE]
> These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen)

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Apptainer, Conda) - see below.

> [!IMPORTANT]
> We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to check if your system is supported, please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended, since it can lead to different results on different machines dependent on the computer environment.

- `test`
  - A profile with a complete configuration for automated testing
  - Includes links to test data so needs no other parameters
- `docker`
  - A generic configuration profile to be used with [Docker](https://docker.com/)
- `singularity`
  - A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)
- `podman`
  - A generic configuration profile to be used with [Podman](https://podman.io/)
- `shifter`
  - A generic configuration profile to be used with [Shifter](https://nersc.gitlab.io/development/shifter/how-to-use/)
- `charliecloud`
  - A generic configuration profile to be used with [Charliecloud](https://charliecloud.io/)
- `apptainer`
  - A generic configuration profile to be used with [Apptainer](https://apptainer.org/)
- `wave`
  - A generic configuration profile to enable [Wave](https://seqera.io/wave/) containers. Use together with one of the above (requires Nextflow ` 24.03.0-edge` or later).
- `conda`
  - A generic configuration profile to be used with [Conda](https://conda.io/docs/). Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker, Singularity, Podman, Shifter, Charliecloud, or Apptainer.

### `-resume`

Specify this when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously. For input to be considered the same, not only the names must be identical but the files' contents as well. For more info about this parameter, see [this blog post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html).

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

## Custom configuration

### Resource requests

Whilst the default requirements set within the pipeline will hopefully work for most people and with most input data, you may find that you want to customise the compute resources that the pipeline requests. Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the pipeline steps, if the job exits with any of the error codes specified [here](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L18) it will automatically be resubmitted with higher resources request (2 x original, then 3 x original). If it still fails after the third attempt then the pipeline execution is stopped.

To change the resource requests, please see the [max resources](https://nf-co.re/docs/usage/configuration#max-resources) and [tuning workflow resources](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources) section of the nf-core website.

### Custom Containers

In some cases, you may wish to change the container or conda environment used by a pipeline steps for a particular tool. By default, nf-core pipelines use containers and software from the [biocontainers](https://biocontainers.pro/) or [bioconda](https://bioconda.github.io/) projects. However, in some cases the pipeline specified version maybe out of date.

To use a different container from the default container or conda environment specified in a pipeline, please see the [updating tool versions](https://nf-co.re/docs/usage/configuration#updating-tool-versions) section of the nf-core website.

### Custom Tool Arguments

A pipeline might not always support every possible argument or option of a particular tool used in pipeline. Fortunately, nf-core pipelines provide some freedom to users to insert additional parameters that the pipeline does not include by default.

To learn how to provide additional arguments to a particular tool of the pipeline, please see the [customising tool arguments](https://nf-co.re/docs/usage/configuration#customising-tool-arguments) section of the nf-core website.

### nf-core/configs

In most cases, you will only need to create a custom config as a one-off but if you and others within your organisation are likely to be running nf-core pipelines regularly and need to use the same settings regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter. You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

See the main [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information about creating your own configuration files.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack) on the [`#configs` channel](https://nfcore.slack.com/channels/configs).

## Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

## Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```
