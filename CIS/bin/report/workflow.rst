COVID Illumina Snakemake Pipeline (CIS)
---------------------------------------

A pipeline to facilitate genomic surveillance of SARS-CoV-2, adopted by
the several labs in SAViD. The starting point is a collection of raw
Illumina short reads; the end point is a collection of QC files,
consensus genome(s), variant calls, and various QC metrics.

Initial QC
----------

The initial QC is done with ``fastp``
(https://github.com/OpenGene/fastp/). View the ``fastp`` repository for
a full list of what the program can do. In this pipeline we leverage
``fastp`` to:

1. filter out bad reads (reads with too many bases below a quality threshold)
2. cut sequencing adapters, with adapter sequences
3. correct mismatched base pairs in overlapped regions of paired end reads, if one base is with high quality while the other is with ultra low quality
4. report results in JSON format for further interpreting.
5. visualizequality control and filtering results in a single HTML page (like FASTQC but faster and more informative).

Read mapping
------------

Reads are mapped to the SARS-CoV-2 reference genome with ``bwa mem``.

Processing of mapped reads
--------------------------

The ``ivar`` suite of tools is used to remove primers from the mapped
reads, and carry out quality trimming. Since quality trimming is already
carried out by ``fastp``, this is a second "line of defence" and might
be superflous. After testing, the quality trimming option in ``ivar``
might be disabled.

Assess coverage
---------------

The coverage for each base in the SARS-CoV-2 genome is calculated with
``bedtools``, and then plotted in R. The coverage for each amplicon
within the bed file specifying primer coordinates is also plotted (but
not for the Swift protocol since it has too many amplicons).

Variant calling
---------------

Variants are called using ``lofreq`` (default) or ``ivar``. Thresholds
are currently hard-coded, but might become command line options later.

SNP annotation
--------------

SNPs are annotated using ``snpEff``, giving information including the
Ts:Tv ratio, and any functional consequences of the SNPs.

Generate consensus genome
-------------------------

A consensus genome is generated for each sample using
``bcftools consensus``. I have implemented the program in a way to
conditionally implement IUPAC codes based on a variant allele frequency
of 0.90.

Pangolin lineage assignment
---------------------------

Newly generated consensus genomes are assigned to Pangolin lineages
using ``pangolin``.

Summary statistics
------------------

A ``.csv`` file with quality control summary statistics is printed to
the results directory.

Explanation of output files
---------------------------

Results for each sample will be placed in their own directory. There are
also logical subdirectories containing the results of each specific
analysis.

Citations
---------

When this pipeline is used, the following programs should be cited
because they're integral to the pipeline working:

-  Heng Li, Bob Handsaker, Alec Wysoker, Tim Fennell, Jue Ruan, Nils
   Homer, Gabor Marth, Goncalo Abecasis, Richard Durbin, 1000 Genome
   Project Data Processing Subgroup, The Sequence Alignment/Map format
   and SAMtools, Bioinformatics, Volume 25, Issue 16, 15 August 2009,
   Pages 2078–2079, https://doi.org/10.1093/bioinformatics/btp352
-  Grubaugh, N.D., Gangavarapu, K., Quick, J., Matteson, N.L., De Jesus,
   J.G., Main, B.J., Tan, A.L., Paul, L.M., Brackney, D.E., Grewal, S.
   and Gurfield, N., 2019. An amplicon-based sequencing framework for
   accurately measuring intrahost virus diversity using PrimalSeq and
   iVar. Genome Biology, Volume 20, Issue 1, pp. 1-19,
   https://doi.org/10.1101/383513
-  Shifu Chen, Yanqing Zhou, Yaru Chen, Jia Gu; 2018. fastp: an
   ultra-fast all-in-one FASTQ preprocessor, Bioinformatics, Volume 34,
   Issue 17, pp. i884–i890,
   https://doi.org/10.1093/bioinformatics/bty560
-  Vasimuddin Md, Sanchit Misra, Heng Li, Srinivas Aluru. Efficient
   Architecture-Aware Acceleration of BWA-MEM for Multicore Systems.
   IEEE Parallel and Distributed Processing Symposium (IPDPS), 2019.
-  Aaron R. Quinlan, Ira M. Hall, BEDTools: a flexible suite of
   utilities for comparing genomic features, Bioinformatics, Volume 26,
   Issue 6, 15 March 2010, Pages 841–842,
   https://doi.org/10.1093/bioinformatics/btq033

Note: not exhaustive. Please find missing references.

You should also find citations for the R packages in the
``environment.yml`` file.

Other credits
-------------

-  ``ivar_variants_to_vcf.py`` was written originally by the
   nextflow/viralrecon team. I've modified it for our purposes.

Workflow DAG
-------------
