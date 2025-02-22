# SARS-CoV-2 Illumina GeNome Assembly Line (SIGNAL) - NML Fork

**Fork Note**: This fork is up-to-date with the upstream repo and makes use of conda and mamba to control dependency and environment installation along with running the pipeline all using the bash script called `run_signal.sh` included. If you wish to use the signal.py python script or use the docker container you are better off going to upstream repo and doing it from there.

This is a complete standardized workflow the assembly and subsequent analysis for short-read viral sequencing.
This core workflow is compatible with the [illumina artic nf pipeline](https://github.com/connor-lab/ncov2019-artic-nf) and produces consensus and variant calls using `iVar` (1.3) [(Grubaugh, 2019)](https://doi.org/10.1186/s13059-018-1618-7) and `Freebayes` [(Garrison, 2012)](https://arxiv.org/abs/1207.3907) .
However, it performs far more extensive quality control and visualisation of results including an interactive HTML summary of run results.

Briefly, raw reads undergo qc using `fastqc` [(Andrews)](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) before removal of host-related reads by competitive mapping against a composite host and human reference with `BWA-MEM` (0.7.5) [(Li, 2013)](https://arxiv.org/abs/1303.3997), `samtools` (1.9) [(Li, 2009)](https://academic.oup.com/bioinformatics/article/25/16/2078/204688), and a [custom script](scripts/filter_non_human_reads.py).
This is to ensure raw as data as possible can be deposited in central databases.
After this, reads undergo adapter trimming and further qc with `trim-galore` (0.6.5) [(Martin)](https://doi.org/10.14806/ej.17.1.200).
Residual truseq sequencing adapters are then removed through another [custom script](scripts/filter_residual_adapters.py).
Reads are then mapped to the viral reference with `BWA-MEM`, and amplicon primer sequences trimmed using `ivar` (1.3) [(Grubaugh, 2019)](https://doi.org/10.1186/s13059-018-1618-7).
Fastqc is then used to perform a QC check on the reads that map to the viral reference.
After this, `iVar` is used to generate a consensus genome and variants are called using both `ivar variants` and `breseq` (0.35) [(Deatherage, 2014)](https://link.springer.com/protocol/10.1007/978-1-4939-0554-6_12). Additionally, `Freebayes` may be run in addition to `iVar` which will generate a second set of consensus genome(s) and variant call(s) with comparisons made between `iVar` and `FreeBayes` to highlight differences in mutation calls.
Coverage statistics are calculated using bedtools before a final QC via `quast` and a `kraken2` taxonomic classification of mapped reads.
Finally, data from all samples are collated via a [post-processing script](scripts/signal_postprocess.py) into an interactive summary for exploration of results and quality control.
Optionally, users can run [ncov-tools](https://github.com/jts/ncov-tools/) to generate additional quality control and summary plots and statistics.

If you use this software please [cite](https://doi.org/10.3390/v12080895):

    Nasir, Jalees A., Robert A. Kozak, Patryk Aftanas, Amogelang R. Raphenya, Kendrick M. Smith, Finlay Maguire, Hassaan Maan et al. "A Comparison of Whole Genome Sequencing of SARS-CoV-2 Using Amplicon-Based Sequencing, Random Hexamers, and Bait Capture." Viruses 12, no. 8 (2020): 895.
    https://doi.org/10.3390/v12080895

## Contents:

- [Setup Instructions](#setup)
- [NML Fork Specific Run Instructions](#nml-fork-specific-run-instructions)
- [SIGNAL Help Screen](#signal-help-screen)
- [Executing SIGNAL - Summary](#summary)
- [Executing SIGNAL - Detailed](#detailed-setup-and-execution)
  - [Download reference file(s)](#1-download-necessary-database-files)
  - [Generate configuation file](#2-generate-configuration-file)
  - [Generate sample table](#3-specify-your-samples-in-csv-format)
  - [Run pipeline](#4-execute-pipeline)
  - [Run postprocessing](#5-postprocessing-analyses)
  - [Run multiple processes](#multiple-operations)
- [Run using Docker](#docker)
- [Data summaries](#data-summaries)
- [Pipeline details](#pipeline-details)

## Setup:

### 0. Clone the git repository (`--recursive` only needed to run `ncov-tools` postprocessing)

        git clone https://github.com/phac-nml/covid-19-signal-nml.git

### 1. Install `conda` and `snakemake` (version >5) e.g.

        wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
        bash Miniconda3-latest-Linux-x86_64.sh # follow instructions
        source $(conda info --base)/etc/profile.d/conda.sh
        conda create -n signal -c conda-forge -c bioconda -c defaults snakemake pandas conda mamba
        conda activate signal

There are some issues with `conda` failing to install newer versions of snakemake
so alternatively install `mamba` and use that (snakemake has beta support for it within the workflow)

        conda install -c conda-forge mamba
        mamba create -c conda-forge -c bioconda -n signal snakemake pandas conda mamba
        conda activate signal
        # mamba activate signal is equivalent

Additional software dependencies are managed directly by `snakemake` using conda environment files:

- trim-galore 0.6.5 ([docs](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/))
- kraken2 2.1.1 ([docs](https://ccb.jhu.edu/software/kraken2/))
- quast 5.0.2 ([docs](http://quast.sourceforge.net/quast))
- bwa 0.7.17 ([docs](http://bio-bwa.sourceforge.net/))
- samtools 1.7/1.9 ([docs](http://www.htslib.org/))
- bedtools 2.26.0 ([docs](https://bedtools.readthedocs.io/en/latest/))
- breseq 0.35.0 ([docs](https://barricklab.org/twiki/bin/view/Lab/ToolsBacterialGenomeResequencing))
- ivar 1.3 ([docs](https://github.com/andersen-lab/ivar))
- freebayes 1.3.2 ([docs](https://github.com/freebayes/freebayes))
- pangolin (latest; version can be specified by user) [(docs)](https://github.com/cov-lineages/pangolin)
- pangolin-data (latest; version can be specified by user; required for Pangolin v4+) [(docs)](https://github.com/cov-lineages/pangolin-data)
- pangolearn (latest; version can be specified by user) [(docs)](https://github.com/cov-lineages/pangoLEARN)
- constellations (latest; version can be specified by user) [(docs)](https://github.com/cov-lineages/constellations)
- scorpio (latest; version can be specified by user) [(docs)](https://github.com/cov-lineages/scorpio)
- pango-designation (latest; version can be specified by user) [(docs)](https://github.com/cov-lineages/pango-designation)
- nextclade (v1.11.0) [(docs)](https://docs.nextstrain.org/projects/nextclade/en/stable/)
- ncov-tools postprocessing scripts require additional dependencies (see [file](https://github.com/jts/ncov-tools/blob/master/workflow/envs/environment.yml)).

## NML Fork Specific Run Instructions

These instructions are specific to running SIGNAL, ncov-tools, and some downstream QC checks using the included `run_signal.sh` bash script.

This process:
- Makes use of conda and mamba to control dependency installation.
- Runs SIGNAL --> PostProcess --> ncov-tools --> QC checks
  - Done using premade SIGNAL and ncov-tools configuration files that are populated based on the command line args

### Help Command:

Running `bash run_signal.sh` will display the following help command with the variables updated to match your system:
```
USAGE:
    bash run_signal.sh -d PATH_TO_PAIRED_FASTQ_DIR -p PRIMER_SCHEME <OPTIONAL FLAGS>
    bash run_signal.sh --update
    
Flags:
    NEEDED:
    -d  --directory      :  Path to paired fastq file directory
    -p  --primer-scheme  :  Specify input data primer scheme
            Available Primer Schemes: articV3, articV4, articV4.1, freed, freed_V2_nml, resende, V2resende, custom by passing '--bed' and '--amplicon'

    SUBSTITUTE (Can be used instead of a primer scheme but must be used together):
    --bed       :  Path to custom primer bed file to be used instead of default schemes
    --amplicon  :  Path to custom amplicon bed file to be used instead of default schemes

    OPTIONAL:
    -c  --cores           :  Number of Cores to use in Signal. Default is 3
    -n  --run-name        :  Run name for final ncov-tools outputs. Default is 'nml'
    -m  --metadata        :  Add metadata to the run. Must be in TSV format with atleast a column called 'sample'
    --pdf                 :  If you have pdflatex installed runs ncov-tools pdf output
    --signal-env          :  Name of signal conda env. Default is '$BASE_ENV_PATH/signal'
    --ncov-tools-env      :  Name of ncov-tools env. Default is '$BASE_ENV_PATH/ncov-qc-pangolin4'
                **NOTE** It is highly recommended to let the script generate the environments as it will
                          only occur once and you won't have to pass the path each time

    OTHER:
    --update  :  Passing --update will update ncov-tools pip dependencies, pangolin and pangoDATA along with this repo and then exit
```

### Run Instructions (run_signal.sh):

#### *1. Get Dependencies*
This pipeline requires a number of files that need to be setup to run. To do so, all that is needed is to run the following bash script once on the first time you are analyzing data:
```
cd covid-19-signal-nml
bash scripts/get_data_dependencies.sh -d data -a MN908947.3
```

This will create a directory called 'data' containing all of the needed files including:
- Amplicon primer scheme sequences
- SARS-CoV2 reference fasta
- SARS-CoV2 reference gbk
- SARS-CoV2 reference gff3
- kraken2 viral database
- Human GRCh38 reference fasta (for composite human-viral BWA index)

**Note: Downloading the database files requires ~10GB of storage, with up to ~35GB required for all temporary downloads!**

#### *2. Run run_signal.sh*
Once all of the data dependencies are downloaded, `run_signal.sh` can be ran. The first run through will be slow as it creates the needed conda environments but all subsequent runs will be significantly faster as they are reused.

To just run the pipeline:
```
bash run_signal.sh -d PATH_TO_PAIRED_FASTQ_DIR -p PRIMER_SCHEME <OPTIONAL FLAGS>
```

Optional arguments can be found above in the [run_signal.sh help command](#help-command)section

### Inputs

- Primer Scheme: What tiled amplicon primer scheme was used to sequence the data
- Paired fastq directory: Directory containing gziped paired fastq files with `_R1` and `_R2` labels
  - Ex. Sample1_S1_L001_R1_001.fastq.gz + Sample1_S1_L001_R2_001.fastq.gz
  - Ex. Sample2_S1_L001_R1.fastq.gz + Sample2_S1_L001_R2.fastq.gz
  - Ex. Sample3_R1.fastq.gz + Sample3_R1.fastq.gz

### Outputs

SIGNAL Pipeline gives the following main outputs (along with many more) for each sample under their own directory in the `signal_results` directory:
```
  - {sample_name}/core/{sample_name}_viral_reference.mapping.primertrimmed.sorted.bam   per-sample sorted primertrimmed alignment file
  - {sample_name}/freebayes/{sample_name}.consensus.fasta                               per-sample generated consensus sequence
  - {sample_name}/freebayes/{sample_name}.variants.norm.vcf.gz                          per-sample normalized variant call information
  - 
```

Postprocessing handles failing samples gracefully and gives the following summary files:
```
  - summary.html                top-level summary, with links to per-sample summaries
  - {sample_name}/sample.html   per-sample summaries, with links for more detailed info
  - {sample_name}/sample.txt    per-sample summaries, in text-file format instead of HTML
  - summary.zip                 zip archive containing all of the above summary files.
```

Ncov-tools follows up with the following main outputs (and many more) under the `ncov-tools` directory:
```
  - plots/{run_name}_tree_snps.pdf                      run-level phylogenetic tree with plotted snps and lineage
  - plots/{run_name}_depth_by_position.pdf              per-sample coverage depth across the genome in one file
  - qc_annotation/{sample_name}_aa_table.tsv            per-sample SNPeff amino acid prediction summarys
  - qc_reports/{run_name}_summary_qc.tsv                run-level summary qc file
  - qc_reports/{run_name}_mixture_report.tsv            run-level sample mixture report
  - qc_reports/{run_name}_negative_control_report.tsv   run-level negative control report
```
- More info on ncov-tools outputs can be found on the ncov-tools github page: https://github.com/jts/ncov-tools

Final QC, found in the created `signal_*_DATE` folder:
```
  - matrix.tsv          snpdists matrix report
  - {run_name}.qc.csv   overall summary of ncov-tools and other metadata checks in one location
```

### Other Notes

#### *Conda Envs*
This script manages the conda environments for you to the best of its ability. They are all found in the cloned `covid-19-signal-nml` directory under `.snakemake/conda` in case they need to be manually addressed

*End NML Fork Specific Instructions*

---------------------------

## SIGNAL Help Screen:

**NOTE:** In this fork the `signalexe.py` script is untested as it was added after the `run_signal.sh` bash script which does mostly the same function of setting up and running the SIGNAL pipeline.

Using the provided `signalexe.py` script, the majority of SIGNAL functions can be accessed easily.

To display the help screen:

```
python signalexe.py -h

usage: signalexe.py [-h] [-c CONFIGFILE] [-d DIRECTORY] [--cores CORES] [--config-only] [--remove-freebayes] [--add-breseq] [-neg NEG_PREFIX] [--dependencies] [--data DATA] [-ri] [-ii] [--unlock]
                    [-F] [-n] [--quiet] [--verbose] [-v]
                    [all ...] [postprocess ...] [ncov_tools ...] [install ...]

SARS-CoV-2 Illumina GeNome Assembly Line (SIGNAL) aims to take Illumina short-read sequences and perform consensus assembly + variant calling for ongoing surveillance and research efforts towards
the emergent coronavirus: Severe Acute Respiratory Syndrome Coronavirus 2 (SARS-CoV-2).

positional arguments:
  all                   Run SIGNAL with all associated assembly rules. Does not include postprocessing '--configfile' or '--directory' required. The latter will automatically generate a
                        configuration file and sample table. If both provided, then '--configfile' will take priority
  postprocess           Run SIGNAL postprocessing on completed SIGNAL run. '--configfile' is required but will be generated if '--directory' is provided
  ncov_tools            Generate configuration file and filesystem setup required and then execute ncov-tools quality control assessment. Requires 'ncov-tools' submodule! '--configfile' is required
                        but will be generated if '--directory' is provided
  install               Install individual rule environments and ensure SIGNAL is functional. The only parameter operable will be '--data'. Will override other operations!

optional arguments:
  -h, --help            show this help message and exit
  -c CONFIGFILE, --configfile CONFIGFILE
                        Configuration file (i.e., config.yaml) for SIGNAL analysis
  -d DIRECTORY, --directory DIRECTORY
                        Path to directory containing reads. Will be used to generate sample table and configuration file
  --cores CORES         Number of cores. Default = 1
  --config-only         Generate sample table and configuration file (i.e., config.yaml) and exit. '--directory' required
  --remove-freebayes    Configuration file generator parameter. Set flag to DISABLE freebayes variant calling (improves overall speed)
  --add-breseq          Configuration file generator parameter. Set flag to ENABLE optional breseq step (will take more time for analysis to complete)
  -neg NEG_PREFIX, --neg-prefix NEG_PREFIX
                        Configuration file generator parameter. Comma-separated list of negative control sample name(s) or prefix(es). For example, 'Blank' will cover Blank1, Blank2, etc.
                        Recommended if running ncov-tools. Will be left empty, if not provided
  --dependencies        Download data dependencies (under a created 'data' directory) required for SIGNAL analysis and exit. Note: Will override other parameters! (~10 GB storage required)
  --data DATA           SIGNAL install and data dependencies parameter. Set location for data dependancies. If '--dependancies' is run, a folder will be created in the specified directory. If '--
                        config-only' or '--directory' is used, the value will be applied to the configuration file. (Upcoming feature): When used with 'SIGNAL install', any tests run will use the
                        dependencies located at this directory. Default = 'data'
  -ri, --rerun-incomplete
                        Snakemake parameter. Re-run any incomplete samples from a previously failed run
  -ii, --ignore-incomplete
                        Snakemake parameter. Do not check for incomplete output files
  --unlock              Snakemake parameter. Remove a lock on the working directory after a failed run
  -F, --forceall        Snakemake parameter. Force the re-run of all rules regardless of prior output
  -n, --dry-run         Snakemake parameter. Do not execute anything and only display what would be done
  --quiet               Snakemake parameter. Do not output any progress or rule information. If used with '--dry-run`, it will just display a summary of the DAG of jobs
  --verbose             Snakemake parameter. Display snakemake debugging output
  -v, --version         Display version number
```

## Summary:

`signalexe.py` simplies the execution of all functions of SIGNAL. At its simplest, SIGNAL can be run with one line, provided only the directory of sequencing reads.

```
# Download dependances (only needs to be run once; ~10GB of storage required)
# --data flag allows you to rename and relocate dependencies directory
python signalexe.py --data data --dependencies

# Generate configuration file and sample table
# --neg_prefix can be used to note negative controls
# --data can be used to specify location of data dependencies
python signalexe.py --config-only --directory /path/to/reads

# Execute pipeline (step-by-step; --cores defaults to 1 if not provided)
# --data can be used to specify location of data dependencies
python signalexe.py --configfile config.yaml --cores NCORES all
python signalexe.py --configfile config.yaml --cores NCORES postprocess
python signalexe.py --configfile config.yaml --cores NCORES ncov_tools

# ALTERNATIVE
# Execute pipeline (one line)
# --data can be used to specify location of data dependencies
python signalexe.py --configfile config.yaml --cores NCORES all postprocess ncov_tools

# ALTERNATIVE
# Execute pipeline (one line; no prior configuration file or sample table steps)
# --directory can be used in place of --configfile to automatically generate a configuration file
# --data can be used to specify location of data dependencies
python signalexe.py --directory /path/to/reads --cores NCORES all postprocess ncov_tools
```

Each of the steps in SIGNAL can be run **manually** by accessing the individual scripts or running snakemake.

```
# Download dependances (only needs to be run once; ~10GB of storage required)
bash scripts/get_data_dependencies.sh -d data -a MN908947.3

# Generate sample table
# Modify existing 'example_config.yaml' for your configuration file
bash scripts/generate_sample_table.sh -d /path/to/reads -n sample_table.csv

# Execute pipeline (step-by-step)
snakemake -kp --configfile config.yaml --cores NCORES --use-conda --conda-prefix=$PWD/.snakemake/conda all
snakemake -kp --configfile config.yaml --cores NCORES --use-conda --conda-prefix=$PWD/.snakemake/conda postprocess
snakemake -kp --configfile config.yaml --cores NCORES --use-conda --conda-prefix=$PWD/.snakemake/conda ncov_tools
```

## Detailed setup and execution:

### 1. Download necessary database files:

The pipeline requires:

- Amplicon primer scheme sequences
- SARS-CoV2 reference fasta
- SARS-CoV2 reference gbk
- SARS-CoV2 reference gff3
- kraken2 viral database
- Human GRCh38 reference fasta (for composite human-viral BWA index)

       python signalexe.py --dependencies
       # defaults to a directory called `data` in repository root
       # --data can be used to rename and relocate the resultant directory

       OR

       bash scripts/get_data_dependencies.sh -d data -a MN908947.3
       # allows you to rename and relocate the resultant directory

**Note: Downloading the database files requires ~10GB of storage, with up to ~35GB required for all temporary downloads!**

### 1.5. Prepare per-rule conda environments (optional, but recommended):

SIGNAL uses controlled conda environments for individual steps in the workflow. These are generally produced upon first execution of SIGNAL with input data; however, an option to install the per-rule environments is available through the `signalexe.py` script.

       python signalexe.py install

       # Will install per-rule environments
       # Later versions of SIGNAL will include a testing module with curated data to ensure function

### 2. Generate configuration file:

You can use the `--config-only` flag to generate both `config.yaml` and `sample_table.csv`. The directory provided will be used to auto-generate a name for the run.

```
python signalexe.py --config-only --directory /path/to/reads

# Outputs: 'reads_config.yaml' and 'reads_sample_table.csv'
# --data can be used to specify the location of data dependancies
```

You can also create the configuration file through modifying the `example_config.yaml` to suit your system.

**Note: Regardless of method, double-check your configuraation file to ensure the information is correct!**

### 3. Specify your samples in CSV format:

See the example table `example_sample_table.csv` for an idea of how to organise this table.

**Using the `--config-only` flag, both configuration file and sample table will be generated (see above in step 2) from a given directory path to reads.**

Alternatively, you can attempt to use `generate_sample_table.sh` to circumvent manual creation of the table.

```
bash scripts/generate_sample_table.sh

Output:
You must specify a data directory containing fastq(.gz) reads.

ASSUMES FASTQ FILES ARE NAMED AS <sample_name>_L00#_R{1,2}*.fastq(.gz)

Flags:
    -d  :  Path to directory containing sample fastq(.gz) files (Absolute paths preferred for consistency, but can use relative paths)
    -n  :  Name or file path for final sample table (with extension) (default: 'sample_table.csv') - will overwrite if file exists
    -e  :  Name or file path for an existing sample table - will append to the end of the provided table

Select one of '-n' (new sample table) or '-e' (existing sample table).
If neither provided, a new sample table called 'sample_table.csv' will be created (or overwritten) by default.
```

General usage:

```
# Create new sample table 'sample_table.csv' given path to reads directory
bash scripts/generate_sample_table.sh -d /path/to/reads -n sample_table.csv

# Append to existing sample table 'sample_table.csv' given path to a directory with additional reads
bash scripts/generate_sample_table.sh -d /path/to/more/reads -e sample_table.csv
```

### 4. Execute pipeline:

For the main `signalexe.py` script, positional arguments inform the rules of the pipeline to execute with flags supplementing input parameters.

The main rules of the pipeline are as followed:

- `all` = Sequencing pipeline. i.e., take a set of paired reads, perform reference-based assembly to generate a consensus, run lineage assignment, etc.
- `postprocess` = Summarize the key results including pangolin lineage, specific mutations, etc, after running `all`
- `ncov_tools` = Create the required conda environment, generate the necessary configuration file, and link needed result files within the `ncov-tools` directory. `ncov-tools` will then be executed with output found within the SIGNAL directory.

The generated configuration file from the above steps can be used as input. To run the general pipeline:

`python signalexe.py --configfile config.yaml --cores 4 all`

is equivalent to running

`snakemake -kp --configfile config.yaml --cores 4 --use-conda --conda-prefix=$PWD/.snakemake/conda all`

You can run the snakemake command as written above, but note that if the `--conda-prefix` is not set as this (i.e., `$PWD/.snakemake/conda`), then all envs will be reinstalled for each time you change the `results_dir` in the `config.yaml`.

Alternatively, you can skip the above configuration and sample table generation steps by simply providing the directory of reads to the main script (see step 2):

`python signalexe.py --directory /path/to/reads --cores 4 all`

A configuartion file and sample table will automatically be generated prior to running SIGNAL `all`.

FreeBayes variant calling and BreSeq mutational analysis are technically optional tools within the workflow. Using the `--directory` flag, by default, FreeBayes **will run** and BreSeq **will not**. These can be changed by using the `--remove-freebayes` and `--add-breseq` flags, respectively.

### 5. Postprocessing analyses:

As with the general pipeline, the generated configuration file from the above steps can be used as input. To run `postprocess` which summarizes the SIGNAL results:

`python signalexe.py --configfile config.yaml --cores 1 postprocess`

is equivalent to running

`snakemake -kp --configfile config.yaml --cores 1 --use-conda --conda-prefix=$PWD/.snakemake/conda postprocess`

After postprocessing finishes, you'll see the following summary files:

```
  - summary.html                top-level summary, with links to per-sample summaries
  - {sample_name}/sample.html   per-sample summaries, with links for more detailed info
  - {sample_name}/sample.txt    per-sample summaries, in text-file format instead of HTML
  - summary.zip                 zip archive containing all of the above summary files.
```

Note that the pipeline postprocessing (`snakemake postprocess`) is separated from
the rest of the pipeline (`snakemake all`). This is because in a multi-sample run,
it's likely that at least one pipeline stage will fail. The postprocessing script
should handle failed pipeline stages gracefully, by substituting placeholder values
when expected pipeline output files are absent. However, this confuses snakemake's
dependency tracking, so there seems to be no good alternative to separating piepline
processing and postprocessing into `all` and `postprocess` targets.

Related: because pipeline stages can fail, we run (and recommend running if using the snakemake command to run SIGNAL) `snakemake all` with the `-k` flag ("Go on with independent jobs if a job fails").

Additionally, SIGNAL can prepare output and execute @jts' [ncov-tools](https://github.com/jts/ncov-tools)
to generate phylogenies and alternative summaries.

`python signalexe.py --configfile config.yaml --cores 1 ncov_tools`

is equivalent to running

`snakemake -kp --configfile config.yaml --cores 1 --use-conda --conda-prefix=$PWD/.snakemake/conda ncov_tools`

SIGNAL manages installing the dependencies (within the `conda_prefix`) and will generate the necessary hard links to required input files from SIGNAL for `ncov-tools` if it has been cloned as a sub-module (if not found, the script will attempt to pull the submodule) and a fasta containing sequences to include in the tree has been specified using `phylo_include_seqs:` in the main SIGNAL `config.yaml`. If `run_freebayes` is set to `True`, then SIGNAL will attempt to link the FreeBayes consensus FASTA and variant files, if found. Otherwise, the corresponding iVar files will be used instead.

SIGNAL will then execute ncov-tools and the **output will be found within the SIGNAL results directory, specified in SIGNAL's configuration file, under `ncov-tools-results`**. Refer to the ncov-tools documentation for information regarding specific output.

### Multiple operations:

Using `signalexe.py` positional arguments, you can specify SIGNAL to perform multiple rules in succession.

`python signalexe.py --configfile config.yaml --cores NCORES all postprocess ncov_tools`

In the above command, SIGNAL `all`, `postprocess`, and `ncov_tools` will run using the provided configuration file as input, which links to a sample table.

**Note: Regardless of order for positional arguments, or placement of other parameter flags, SIGNAL will always run in the set order priority: `all` > `postprocess` > `ncov_tools`!**

**Note: If `install` is provided as input, it will override all other positional arguments!**

If no configuration file or sample table was generated for a run, you can provide `--directory` with the path to sequencing reads and SIGNAL will auto-generate both required inputs prior to running any rules.

`python signalexe.py --directory /path/to/reads --cores NCORES all postprocess ncov_tools`

Overall, this simplifies executing SIGNAL to one line!

### Docker:

Alternatively, the pipeline can be deployed using Docker (see `resources/Dockerfile_pipeline` for specification).
To pull from dockerhub:

        docker pull finlaymaguire/signal

Download data dependencies into a data directory that already contains your reads (`data` is this example but whatever name you wish to use):

        mkdir -p data && docker run -v $PWD/data:/data finlaymaguire/signal:latest bash scripts/get_data_dependencies.sh -d /data

Generate your `config.yaml` and `sample_table.csv` (with paths to the readsets underneath `/data`) and place them into the data directory:

        cp config.yaml sample_table.csv $PWD/data

_WARNING_ `result_dir` in `config.yaml` must be within `/data` e.g. `/data/results` to automatically be copied to your host system. Otherwise they will be automatically deleted when the container finishes running (unless docker is run interactively).

Then execute the pipeline:

        docker run -v $PWD/data:/data finlaymaguire/signal conda run -n snakemake snakemake --configfile /data/config.yaml --use-conda --conda-prefix /covid-19-signal/.snakemake/conda --cores 8 all

## Data Summaries:

- `postprocess` and `ncov_tools` as described above generate many summaries including interactive html reports

- To generate summaries of BreSeq among many samples, see [how to summarize BreSeq results using gdtools](resources/dev_scripts/summaries/README.md)

### Convenient extraction script:

SIGNAL produces several output files and directories on its own alongside the output for ncov-tools. Select files from the output can be copied or transferred for easier parsing using a provided convenience bash script:

```
bash scripts/get_signal_results.sh

Usage:
bash get_signal_results.sh -s <SIGNAL_results_dir> -d <destination_dir> [-m] [-c]

This scripts aims to copy (rsync by default, or cp) or move (mv) select output from SIGNAL 'all', 'postprocess', and 'ncov_tools'.

The following files will be transferred over to the specified destination directory (if found):
SIGNAL 'all' & 'postprocess':
-> signal-results/<sample>/<sample>_sample.txt
-> signal-results/<sample>/core/<sample>.consensus.fa
-> signal-results/<sample>/core/<sample>_ivar_variants.tsv
-> signal-results/<sample>/freebayes/<sample>.consensus.fasta
-> signal-results/<sample>/freebayes/<sample>.variants.norm.vcf

SIGNAL 'ncov_tools':
-> ncov_tools-results/qc_annotation/<sample>.ann.vcf
-> ncov-tools-results/qc_reports/<run_name>_ambiguous_position_report.tsv
-> ncov-tools-results/qc_reports/<run_name>_mixture_report.tsv
-> ncov-tools-results/qc_reports/<run_name>_ncov_watch_variants.tsv
-> ncov-tools-results/qc_reports/<run_name>_negative_control_report.tsv
-> ncov-tools-results/qc_reports/<run_name>_summary_qc.tsv

Flags:
        -s  :  SIGNAL results directory
        -d  :  Directory where summary will be outputted
        -m  :  Invoke 'mv' move command instead of 'rsync' copying of results. Optional
        -c  :  Invoke 'cp' copy command instead of 'rsync' copying of results. Optional

```

The script uses `rsync` to provide accurate copies of select output files organized into `signal-results` and `ncov-tools-results` within a provided destination directory (that must exist). If the `-c` is provided, `cp` will be used instead of `rsync` to produce copies. Similarly, if `-m` is provided, `mv` will be used instead (**WARNING: Any interruption during `mv` could result in data loss.**)

## Pipeline details:

For a step-by-step walkthrough of the pipeline, see [pipeline/README.md](PIPELINE.md).

A diagram of the workflow is shown below (update pending).

![Workflow Version 8](./resources/Workflow_Version_8.png)
