[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
```
           /^\/^\  COVID
         _|__|  O|  Illumina
\/     /~     \_/ \    Pipeline
 \____|__________/  \      Snakemake edition
        \_______      \
                `\     \                 \
                  |     |                  \
                 /      /                    \
                /     /                       \\
              /      /                         \ \
             /     /                            \  \
           /     /             _----_            \   \
          /     /           _-~      ~-_         |   |
         (      (        _-~    _--_    ~-_     _/   |
          \      ~-____-~    _-~    ~-_    ~-_-~    /
            ~-_           _-~          ~-_       _-~
               ~--______-~                ~-___-~    

```

# COVID Illumina Snakemake Pipeline
A pipeline to facilitate genomic surveillance of SARS-CoV-2, adopted by the several labs in SAViD. The starting point is a collection of raw Illumina short reads; the end point is a collection of QC files, consensus genome(s), variant calls, and various QC metrics.

Author: Dr Charles Foster

# Starting out
To begin with, clone this github repository:

```
git clone https://github.com/charlesfoster/covid-illumina-snakemake.git

cd covid-illumina-snakemake
```

Next, install *most* dependencies using `conda`:

```
conda env create -f environment.yml
```

Pro tip: if you install `mamba`, you can create the environment with that command instead of `conda`. A lot of `conda` headaches go away: it's a much faster drop-in replacement for `conda`.

```
conda install mamba
mamba env create -f environment.yml
```

Install the pipeline correctly with:

```
pip install .
```

Other dependencies:
* Lineages are typed using `pangolin`. Accordingly, `pangolin` needs to be installed according to instructions at https://github.com/cov-lineages/pangolin.
* Variants are called using either `lofreq` (preferred) or `ivar`. Ideally we could install these via `mamba` in the main `environment.yml` file, but the hosted versions do not yet support htslib>=1.13 (needed for our specific `bcftools commands`). Both `lofreq` and `ivar` are currently installed into their own environments during execution of the CIS pipeline. These environments are only created the first time you run the pipeline, and then subsequent runs use the isolated `lofreq` and `ivar` environments automatically.

**tl;dr**: you don't need to do anything for `lofreq` or `ivar` to work; just don't get confused during the initial pipeline run when the terminal indicates creation of a new conda environment. You should run `CIS` from the same directory each time, otherwise a new conda environment will be created each time.

# Usage
The environment with all necessary tools is installed as '`CIS`' for brevity. The environment should first be activated:

```
conda activate CIS
```

Then, to run the pipeline, it's as simple as:

```
CIS <directory_with_reads>
```

where <directory_with_reads> should be replaced with the full path to a directory with sequencing reads taken off an Illumina machine.

The program assumes that reads are named exactly as they are after coming off a MiSeq/iSeq, i.e. `*_L001_R1_001.fastq.gz`. This unique suffix allows the correct sample names and corresponding forward and reverse reads files to be identified. However, if you have renamed your reads files, you need to provide a suffix via `--suffix` that allows the correct samples to be identified.

Consider the situation where you have renamed your files like so: `sample1_R1.fq.gz`, `sample1_R2.fq.gz`, `sample2_R1.fq.gz`, and `sample2_R2.fq.gz`, and they're in `/home/user/reads`. Now, you will need to run the program with these minimal options:

```
CIS --suffix "_R1.fq.gz" /home/user/reads
```

The program will determine that your forward and reverse reads are in `/home/user/reads`, and that your sample names are 'sample1' and 'sample2'.

By default, the program assumes that you have used the 'Midnight' amplicon protocol. You can specify a different scheme using the `-s` (`--scheme`) option.

All other options are as follows, and can be accessed with `CIS --help`:

```

           /^\/^\  COVID
         _|__|  O|  Illumina
\/     /~     \_/ \    Pipeline
 \____|__________/  \      Snakemake edition
        \_______      \
                `\     \                 \
                  |     |                  \
                 /      /                    \
                /     /                       \\
              /      /                         \ \
             /     /                            \  \
           /     /             _----_            \   \
          /     /           _-~      ~-_         |   |
         (      (        _-~    _--_    ~-_     _/   |
          \      ~-____-~    _-~    ~-_    ~-_-~    /
            ~-_           _-~          ~-_       _-~
               ~--______-~                ~-___-~    

      usage: CIS [options] <query_directory>

      covid-illumina-snakemake: a pipeline for analysis SARS-CoV-2 samples

      positional arguments:
         query_directory       Path to directory with reads to process.

      optional arguments:
         -h, --help            show this help message and exit
         -c <int>, --consensus_freq <int>
                               Variant allele frequency threshold for a variant to be
                               incorporated into consensus genome. Default: 0.9
         -i ISOLATES, --isolates ISOLATES
                               List of isolates to assemble (Default: all isolates in
                               query_directory)
         -f, --force           Force overwriting of completed files (Default: files
                               not overwritten)
         -k KRAKEN2_DB, --kraken2_db KRAKEN2_DB
                               kraken2 database. Default:
                               /data/kraken_files/k2_standard_20201202
         -n, --dry_run         Dry run only
         -o OUTDIR, --outdir OUTDIR
                               Output directory. Default: /home/cfos/Programs/COVID_I
                               llumina_Snakemake/results/2021-10-27
         -p, --print_dag       Save directed acyclic graph (DAG) of workflow
         -r REFERENCE, --reference REFERENCE
                               Reference genome to use (Default:
                               /home/cfos/miniconda3/envs/CIS/lib/python3.9/site-
                               packages/CIS/bin/NC_045512.fasta)
         -s SCHEME, --scheme SCHEME
                               Primer scheme to use: built-in opts are 'midnight',
                               'swift', 'eden', but if using your own scheme provide
                               the full path to the bed file here (Default: midnight)
         -t <int>, --threads <int>
                               Number of threads to use
         -v VARIANT_CALLER, --variant_caller VARIANT_CALLER
                               Variant caller to use. Choices: 'lofreq' or 'ivar'
                               Default: 'lofreq'
         --keep_reads          Keep trimmed reads
         --version             show program's version number and exit
         --suffix <str>        Suffix used to identify samples from reads. Default:
                               _L001_R1_001.fastq.gz
         --max_memory <int>    Maximum memory (in MB) that you would like to provide
                               to snakemake. Default: 58683MB
         --verbose             Print junk to screen.
```

# What does the pipeline do?
- [Initial QC: adapter removal, error correction, QC reports](#Initial-QC)
- [Read mapping](#Read-mapping)
- [Processing of mapped reads](#Processing-of-mapped-reads)
- [Assess coverage](#Assess-coverage)
- [Variant calling](#Variant-calling)
- [SNP annotation](#SNP-annotation)
- [Generate consensus genome](#Generate-consensus-genome)
- [Pangolin lineage assignment](#Pangolin-lineage-assignment)
- [Summary statistics](#Summary-statistics)
- [Explanation of output files](#Explanation-of-output-files)
- [Next steps](#Next-steps)
- [Citations](#Citations)
- [Other credits](#Other-credits)

# Initial QC
The initial QC is done with `fastp` (https://github.com/OpenGene/fastp/). View the `fastp` repository for a full list of what the program can do. In this pipeline we leverage `fastp` to:
1. filter out bad reads (reads with too many bases below a quality threshold)
2. cut sequencing adapters, with adapter sequences automatically detected via overlap analysis of PE data
3. correct mismatched base pairs in overlapped regions of paired end reads, if one base is with high quality while the other is with ultra low quality
4. report results in JSON format for further interpreting.
5. visualize quality control and filtering results in a single HTML page (like FASTQC but faster and more informative).

# Read mapping
Reads are mapped to the SARS-CoV-2 reference genome with `bwa mem`.

# Processing of mapped reads
The `ivar` suite of tools is used to remove primers from the mapped reads, and carry out quality trimming. Since quality trimming is already carried out by `fastp`, this is a second "line of defence" and might be superflous. After testing, the quality trimming option in `ivar` might be disabled.

# Assess coverage
The coverage for each base in the SARS-CoV-2 genome is calculated with `bedtools`, and then plotted in R. The coverage for each amplicon within the bed file specifying primer coordinates is also plotted (but not for the Swift protocol since it has too many amplicons).

# Variant calling
Variants are called using `lofreq` (default) or `ivar`. Thresholds are currently hard-coded, but might become command line options later.

# SNP annotation
SNPs are annotated using `snpEff`, giving information including the Ts:Tv ratio, and any functional consequences of the SNPs.

# Generate consensus genome
A consensus genome is generated for each sample using `bcftools consensus`. I have implemented the program in a way to conditionally implement IUPAC codes based on a variant allele frequency of 0.90.

# Pangolin lineage assignment
Newly generated consensus genomes are assigned to Pangolin lineages using `pangolin`.

# Summary statistics
A `.csv` file with quality control summary statistics is printed to the results directory.

# Explanation of output files
Results for each sample will be placed in their own directory. There are also logical subdirectories containing the results of each specific analysis.

# Citations
When this pipeline is used, the following programs should be cited because they're integral to the pipeline working:

* Heng Li, Bob Handsaker, Alec Wysoker, Tim Fennell, Jue Ruan, Nils Homer, Gabor Marth, Goncalo Abecasis, Richard Durbin, 1000 Genome Project Data Processing Subgroup, The Sequence Alignment/Map format and SAMtools, Bioinformatics, Volume 25, Issue 16, 15 August 2009, Pages 2078–2079, https://doi.org/10.1093/bioinformatics/btp352
* Grubaugh, N.D., Gangavarapu, K., Quick, J., Matteson, N.L., De Jesus, J.G., Main, B.J., Tan, A.L., Paul, L.M., Brackney, D.E., Grewal, S. and Gurfield, N., 2019. An amplicon-based sequencing framework for accurately measuring intrahost virus diversity using PrimalSeq and iVar. Genome Biology, Volume 20, Issue 1, pp. 1-19, https://doi.org/10.1101/383513
* Shifu Chen, Yanqing Zhou, Yaru Chen, Jia Gu; 2018. fastp: an ultra-fast all-in-one FASTQ preprocessor, Bioinformatics, Volume 34, Issue 17, pp. i884–i890, https://doi.org/10.1093/bioinformatics/bty560
* Vasimuddin Md, Sanchit Misra, Heng Li, Srinivas Aluru. Efficient Architecture-Aware Acceleration of BWA-MEM for Multicore Systems. IEEE Parallel and Distributed Processing Symposium (IPDPS), 2019.
* Aaron R. Quinlan, Ira M. Hall, BEDTools: a flexible suite of utilities for comparing genomic features, Bioinformatics, Volume 26, Issue 6, 15 March 2010, Pages 841–842, https://doi.org/10.1093/bioinformatics/btq033

Note: not exhaustive. Please find missing references.

You should also find citations for the R packages in the `environment.yml` file.

# Other credits
*  `ivar_variants_to_vcf.py` was written originally by the nextflow/viralrecon team. I've modified it for our purposes.
* Being new to using snakemake + wrapper scripts, I used `pangolin` as a guide for directory structure and rule creation.
