#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  3 15:57:54 2021

@author: cfos
"""
from CIS import __version__
import os
import snakemake
import argparse
import sys
import re
import pathlib
import psutil
from datetime import date
from psutil import virtual_memory
from psutil._common import bytes2human
import textwrap as _textwrap
from CIS.bin.scripts.check_nextclade import check_nextclade

# Define the wrapping of help text: https://stackoverflow.com/questions/35917547/python-argparse-rawtexthelpformatter-with-line-wrap
os.environ['COLUMNS'] = "120"
class PreserveWhiteSpaceWrapRawTextHelpFormatter(argparse.RawDescriptionHelpFormatter):
    def __add_whitespace(self, idx, iWSpace, text):
        if idx is 0:
            return text
        return (" " * iWSpace) + text

    def _split_lines(self, text, width):
        textRows = text.splitlines()
        for idx,line in enumerate(textRows):
            search = re.search('\s*[0-9\-]{0,}\.?\s*', line)
            if line.strip() is "":
                textRows[idx] = " "
            elif search:
                lWSpace = search.end()
                lines = [self.__add_whitespace(i,lWSpace,x) for i,x in enumerate(_textwrap.wrap(line, width))]
                textRows[idx] = lines

        return [item for sublist in textRows for item in sublist]

thisdir = os.path.abspath(os.path.dirname(__file__))


def bytesto(bytes, to, bsize=1024):
    """convert bytes to megabytes, etc.
    sample code:
        print('mb= ' + str(bytesto(314575262000000, 'm')))
    sample output:
        mb= 300002347.946
    """

    a = {"k": 1, "m": 2, "g": 3, "t": 4, "p": 5, "e": 6}
    r = float(bytes)
    for i in range(a[to]):
        r = r / bsize

    return r


def main(sysargs=sys.argv[1:]):
    print(
        """\033[92m
           /^\/^\\  COVID
         _|__|  O|  Illumina
\/     /~     \_/ \\    Pipeline
 \____|__________/  \\      Snakemake edition v{}
        \_______      \\
                `\     \                 \\
                  |     |                  \\
                 /      /                    \\
                /     /                       \\\\
              /      /                         \ \\
             /     /                            \  \\
           /     /             _----_            \   \\
          /     /           _-~      ~-_         |   |
         (      (        _-~    _--_    ~-_     _/   |
          \      ~-____-~    _-~    ~-_    ~-_-~    /
            ~-_           _-~          ~-_       _-~
               ~--______-~                ~-___-~    \033[0m
    """.format(
            __version__
        )
    )

    max_mem = round(bytesto(virtual_memory().available, "m"))

    try:
        default_kraken = os.environ["KRAKEN2_DEFAULT_DB"]
    except:
        default_kraken = None
    parser = argparse.ArgumentParser(
        description="covid-illumina-snakemake: a pipeline for analysis SARS-CoV-2 samples",
        usage="""CIS [options] <query_directory> """,
        formatter_class=PreserveWhiteSpaceWrapRawTextHelpFormatter,
    )

    parser.add_argument(
        "query_directory", nargs="*", help="Path to directory with reads to process."
    )
    parser.add_argument(
        "-c",
        "--consensus_freq",
        action="store",
        required=False,
        help="Variant allele frequency (VAF) threshold for a non-indel variant to be incorporated into consensus genome. Variants with a VAF less than the threshold will be incorporated as IUPAC ambiguities. Default: {}".format(
            float(0.75)
        ),
        metavar="<float>",
    )
    parser.add_argument(
        "-if",
        "--indel_freq",
        action="store",
        required=False,
        help="Variant allele frequency threshold for an indel variant to be incorporated into consensus genome. Default: {}".format(
            float(0.75)
        ),
        metavar="<float>",
    )
    parser.add_argument(
        "-i",
        "--isolates",
        action="store",
        required=False,
        help="List of isolates to assemble (Default: all isolates in query_directory)",
    )
    parser.add_argument(
        "-f",
        "--force",
        action="store_true",
        default=False,
        required=False,
        help="Force overwriting of completed files (Default: files not overwritten)",
    )
    parser.add_argument(
        "-k",
        "--kraken2_db",
        action="store",
        help="kraken2 database. Default: {}".format(default_kraken),
    )
    parser.add_argument(
        "-m",
        "--min_depth",
        default=int(10),
        action="store",
        help="Minimum depth for (1) an SNV to be kept; and (2) consensus genome generation. Default: {}".format(int(10)),
    )
    parser.add_argument(
        "-n", "--dry_run", action="store_true", default=False, help="Dry run only"
    )
    parser.add_argument(
        "-o",
        "--outdir",
        action="store",
        help="Explicitly specify output directory rather than the default of '{}'".format(
            os.path.join(os.getcwd(), "results", '<input_reads_directory_name>')
        ),
    )
    parser.add_argument(
        "-p",
        "--print_dag",
        action="store_true",
        default=False,
        help="Save directed acyclic graph (DAG) of workflow",
    )
    parser.add_argument(
        "-r",
        "--reference",
        action="store",
        required=False,
        default=os.path.join(thisdir, "bin", "references", "NC_045512.fasta"),
        help="Reference genome to use (Default: {})".format(
            os.path.join(thisdir, "bin", "NC_045512.fasta")
        ),
    )
    parser.add_argument(
        "-s",
        "--scheme",
        action="store",
        required=False,
        default="midnight",
        help="""Primer scheme to use. Built-in opts are:
         - 'midnight'
         - 'swift'
         - 'eden'
         - 'articv4.1'
         - 'articv3'

         If using your own scheme provide the full path to the bed file here (Default: {})""".format(
            "midnight"
        ),
    )
    parser.add_argument(
        "-t",
        "--threads",
        action="store",
        help="Number of threads to use",
        default=psutil.cpu_count(logical=True),
        metavar="<int>",
    )
    parser.add_argument(
        "-v",
        "--variant_caller",
        action="store",
        help="Variant caller to use. Choices: 'lofreq' or 'ivar' Default: 'lofreq'",
        default="lofreq",
    )
    parser.add_argument(
        "-w",
        "--workflow",
        action="store",
        help="Workflow to use. Choices: 'routine' or 'wastewater' Default: 'routine'",
        default="routine",
    )
    parser.add_argument(
        "-xn",
        "--skip_nextclade_check",
        action="store_true",
        default=False,
        help="Skip check for a Nextclade update",
    )
    parser.add_argument(
        "-xf",
        "--skip_freyja_check",
        action="store_true",
        default=False,
        help="Skip check for a freyja update (if using wastewater workflow)",
    )
    parser.add_argument(
        "--legacy",
        action="store_true",
        help="Output the quality control column names in 'legacy' format. Default: {}".format(False),
    )
    # parser.add_argument(
    #     "--no_singularity",
    #     action="store_true",
    #     help="Stop the use of singularity. Default: {}".format(False),
    # )
    parser.add_argument(
        "--keep_reads", action="store_true", help="Keep trimmed reads", default=False
    )
    parser.add_argument(
        "--version",
        action="version",
        version=f"covid illumina pipeline snakemake edition:  {__version__}",
    )
    parser.add_argument(
        "--technology",
        action="store",
        help="Sequencing technology (Default: 'Illumina Miseq')",
        default='Illumina Miseq',
        metavar="<str>",
    )
    parser.add_argument(
        "--use_date",
        action="store_true",
        help="Name output directory and files after today's date rather than input reads directory name",
        default=False,
    )
    parser.add_argument(
        "--snv_min_freq",
        action="store",
        required=False,
        default=0.25,
        help="Variant allele frequency threshold for an SNV to be kept. Default: {}".format(
            0.25
        ),
        metavar="float",
    )
    parser.add_argument(
        "--suffix",
        action="store",
        help="Suffix used to identify samples from reads. Default: {}".format(
            "_L001_R1_001.fastq.gz"
        ),
        metavar="<str>",
    )
    parser.add_argument(
        "--max_memory",
        action="store",
        help="Maximum memory (in MB) that you would like to provide to snakemake. Default: {}MB".format(
            max_mem
        ),
        metavar="<int>",
    )
    parser.add_argument(
        "--redo_demix",
        action="store_true",
        help="Redo demixing using freyja ('wastewater' workflow only)",
        default=False,
    )
    parser.add_argument(
        "-e",
        "--extra_pangolin_options",
        action="store",
        help="Extra options to be passed to pangolin",
        default='',
        metavar="<str>",
    )
    parser.add_argument("--verbose", action="store_true", help="Print junk to screen.")
    parser.add_argument("--report", action="store_true", help="Generate report.")

    if len(sysargs) < 1:
        parser.print_help()
        sys.exit(-1)
    else:
        args = parser.parse_args(sysargs)

    args = parser.parse_args()

    if not os.path.exists("".join(args.query_directory)):
        print(
            "#####\n\033[91mError\033[0m: Data directory does not appear to exist\n#####\n"
        )
        print("Please check your input and try again\n")
        print("For help, run: CIS --help\n\n")
        sys.exit(-1)

    if default_kraken is None and not args.kraken2_db:
        print(
            "#####\n\033[91mError\033[0m: No kraken2 database supplied, and default database cannot be detected\n#####\n"
        )
        parser.print_help()
        sys.exit(-1)
    elif args.kraken2_db:
        default_kraken = args.kraken2_db

    if args.use_date:
        results_string = date.today().strftime("%Y-%m-%d")
    else:
        results_string = pathlib.PurePath(args.query_directory[0]).name
    if args.outdir:
        outdir = args.outdir
    else:
        outdir = os.path.join(os.getcwd(), "results", results_string)

    if not os.path.exists(outdir):
        os.makedirs(outdir)

    if args.max_memory:
        max_mem = {"mem_mb": int(args.max_memory)}
    else:
        max_mem = {"mem_mb": int(max_mem)}

    scheme = os.path.join(thisdir, "bin", "primer_schemes", "midnight.bed")
    if args.scheme and os.path.isfile(
        os.path.join(thisdir, "bin", "primer_schemes", "{}.bed".format(args.scheme))
    ):
        scheme = os.path.join(
            thisdir, "bin", "primer_schemes", "{}.bed".format(args.scheme)
        )

    coverage_script = os.path.join(thisdir, "bin", "scripts", "plot_coverage.R")
    converter_script = os.path.join(
        thisdir, "bin", "scripts", "ivar_variants_to_vcf.py"
    )
    vcf_script = os.path.join(thisdir, "bin", "scripts", "add_fake_genotype.sh")

    suffix = "_L001_R1_001.fastq.gz"
    if args.suffix:
        suffix = args.suffix

    # check input suffix
    check_dir = [
        x for x in os.listdir("".join(args.query_directory)) if x.endswith(suffix)
    ]
    if len(check_dir) == 0:
        print(
            "#####\n\033[91mError\033[0m: Specified suffix does not appear to match the reads in your query_directory\n#####\n"
        )
        print("Please check your input and try again\n")
        print("For help, run: CIS --help\n\n")
        sys.exit(-1)

    if args.isolates:
        with open(args.isolates) as f:
            isolates = [line.rstrip("\n") for line in f]
    else:
        isolates = False

    if args.variant_caller:
        if args.variant_caller not in ["lofreq", "ivar"]:
            print(
                '#####\n\033[91mError\033[0m: Variant caller must be either "lofreq" or "ivar"\n#####\n'
            )
            sys.exit(1)
        variant_caller = args.variant_caller
    else:
        variant_caller = "lofreq"

    if variant_caller is None:
        print(
            "#####\n\033[91mError\033[0m: {} could not be found in your path\n#####\n".format(
                variant_caller
            )
        )
        sys.exit(1)

    consensus_freq = 0.75
    if args.consensus_freq:
        consensus_freq = float(args.consensus_freq)
        if consensus_freq > 1 or consensus_freq < 0:
            print(
                "#####\n\033[91mError\033[0m: The consensus_freq option must be a float number between 0 and 1\n#####\n"
            )
            sys.exit(1)

    snv_min_freq = 0.25
    if args.snv_min_freq:
        snv_min_freq = float(args.snv_min_freq)
        if snv_min_freq > 1 or snv_min_freq < 0:
            print(
                "#####\n\033[91mError\033[0m: The snv_min_freq option must be a float number between 0 and 1\n#####\n"
            )
            sys.exit(1)

    indel_freq = 0.75
    if args.indel_freq:
        indel_freq = float(args.indel_freq)
        if indel_freq > 1 or indel_freq < 0:
            print(
                "#####\n\033[91mError\033[0m: The indel_freq option must be a float number between 0 and 1\n#####\n"
            )
            sys.exit(1)

    min_depth = 10
    if args.min_depth:
        min_depth = int(args.min_depth)
        if min_depth < 1:
            print(
                "#####\n\033[91mError\033[0m: The min_depth option must be an integer >=1\n#####\n"
            )
            sys.exit(1)

    if args.workflow not in ['routine', 'wastewater']:
        print(
            "#####\n\033[91mError\033[0m: workflow can only be specified as 'routine' or 'wastewater'\n#####\n"
        )
        sys.exit(1)

    if not os.path.isfile(args.reference + ".fai"):
        print("Indexing {} with samtools".format(args.reference))
        os.system("samtools faidx {} 2> /dev/null".format(args.reference))
    if not os.path.isfile(args.reference + ".bwt"):
        print("Indexing {} with bwa".format(args.reference))
        os.system("bwa index {} 2> /dev/null".format(args.reference))

    if args.no_singularity:
        use_singularity = False
    else:
        use_singularity = True

    snakefile = os.path.join(thisdir, "bin", "Snakefile")

    config = {
        "reads_dir": os.path.join(os.getcwd(), args.query_directory[0]),
        "nextclade_dataset": os.path.join(thisdir, "nextclade_datasets","sars-cov-2"),
        "singularity_args": f'--bind {outdir}:{outdir},{os.environ["CONDA_PREFIX"]}:{os.environ["CONDA_PREFIX"]}',
        "kraken2_db": default_kraken,
        "outdir": outdir,
        "reference": args.reference,
        "annotation": re.sub(".fasta",".gff3",args.reference),
        "isolates": isolates,
        "coverage_script": coverage_script,
        "converter_script": converter_script,
        "vcf_script": vcf_script,
        "variant_program": variant_caller,
        "scheme": scheme,
        "isolates": isolates,
        "legacy_results": args.legacy,
        "suffix": suffix,
        "threads": args.threads,
        "consensus_freq": consensus_freq,
        "min_depth": min_depth,
        "snv_min_freq": snv_min_freq,
        "indel_freq": indel_freq,
        "extra_pangolin_args": args.extra_pangolin_options,
        "keep_reads": args.keep_reads,
        "verbose": args.verbose,
        "technology": args.technology,
        "use_date": args.use_date,
        "workflow": args.workflow,
    }

    # check for nextclade update
    if not args.skip_nextclade_check:
        check_nextclade(thisdir,args.workflow)

    # wastewater specific
    if args.workflow == "wastewater":
        config["freyja_dataset"] = os.path.join(thisdir, "freyja_dataset")
        if not args.skip_freyja_check:
            from CIS.bin.scripts.check_freyja import check_freyja
            check_freyja(thisdir,args.workflow)
        if args.redo_demix:
            demix_files = glob.glob(outdir + "/**/*.demix", recursive=True)
            [os.remove(f) for f in demix_files]


    if args.print_dag:
        flat_config = []
        for key in config:
            flat_config.append(key + "=" + '"'+str(config[key])+'"')
        flat_config = " ".join(flat_config)
        cmd = 'snakemake -j1 -s {0}  --config {1} --dag | grep -v "No negative control samples detected" | dot -Tpdf > dag.pdf'.format(
            os.path.join(thisdir, "bin", "Snakefile"), flat_config
        )
        os.system(cmd)
        status = True
    elif args.report:
        status = snakemake.snakemake(
            snakefile,
            report=os.path.join(outdir, "pipeline_report.html"),
            use_conda=True,
            use_singularity=use_singularity,
            singularity_args=config['singularity_args'],
            conda_frontend="mamba",
            dryrun=args.dry_run,
            printshellcmds=True,
            forceall=args.force,
            force_incomplete=True,
            resources=max_mem,
            config=config,
            quiet=False,
            cores=args.threads,
            lock=False,
            #delete_temp_output=True,
        )
    elif config["verbose"]:
        print("\n**** CONFIG ****")
        for k in config:
            print(k + ": ", config[k])

        status = snakemake.snakemake(
            snakefile,
            use_conda=True,
            use_singularity=use_singularity,
            singularity_args=config['singularity_args'],
            conda_frontend="mamba",
            dryrun=args.dry_run,
            printshellcmds=True,
            forceall=args.force,
            force_incomplete=True,
            resources=max_mem,
            config=config,
            quiet=False,
            cores=args.threads,
            lock=False,
            #delete_temp_output=True,
        )
    else:
        status = snakemake.snakemake(
            snakefile,
            use_conda=True,
            use_singularity=use_singularity,
            singularity_args=config['singularity_args'],
            conda_frontend="mamba",
            dryrun=args.dry_run,
            printshellcmds=False,
            forceall=args.force,
            force_incomplete=True,
            resources=max_mem,
            config=config,
            quiet=True,
            cores=args.threads,
            lock=False,
            #delete_temp_output=True,
        )

    if status and args.print_dag:  # translate "success" into shell exit code of 0
        print(
            "\033[92m\nResults\033[0m: {}\n\n\033[92mPipeline complete!\033[0m\n".format(
                os.path.join(os.getcwd(), "dag.pdf")
            )
        )
        return 0
    else:
        print(
            "\033[92m\nResults\033[0m: {}\n\n\033[92mPipeline complete!\033[0m\n".format(
                config["outdir"]
            )
        )
        return 0
    return 1


if __name__ == "__main__":
    main()
