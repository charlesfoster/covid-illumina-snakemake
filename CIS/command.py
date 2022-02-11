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
import shutil
import psutil
from datetime import date
from psutil import virtual_memory
from psutil._common import bytes2human

thisdir = os.path.abspath(os.path.dirname(__file__))
TODAY = date.today().strftime("%Y-%m-%d")


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
    )

    parser.add_argument(
        "query_directory", nargs="*", help="Path to directory with reads to process."
    )
    parser.add_argument(
        "-c",
        "--consensus_freq",
        action="store",
        required=False,
        help="Variant allele frequency threshold for a non-indel variant to be incorporated into consensus genome. Default: {}".format(
            float(0.90)
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
        "-n", "--dry_run", action="store_true", default=False, help="Dry run only"
    )
    parser.add_argument(
        "-o",
        "--outdir",
        action="store",
        help="Output directory. Default: {}".format(
            os.path.join(os.getcwd(), "results", TODAY)
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
        help="Primer scheme to use: built-in opts are 'midnight', 'swift', 'eden', 'articv4.1', and 'articv3', but if using your own scheme provide the full path to the bed file here  (Default: {})".format(
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
        "--keep_reads", action="store_true", help="Keep trimmed reads", default=False
    )
    parser.add_argument(
        "--version",
        action="version",
        version=f"covid illumina pipeline snakemake edition:  {__version__}",
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

    if args.outdir:
        outdir = args.outdir
    else:
        outdir = os.path.join(os.getcwd(), "results", TODAY)

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
    #else:
    #    scheme = args.scheme
    #print(scheme)
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
    # variant_caller = shutil.which(variant_caller)

    # ivar is still necessary for trimming
    # if shutil.which('ivar') is None:
    #    print('#####\n\033[91mError\033[0m: "ivar" is necessary for amplicon trimming but is not in your path\n#####\n')
    #    sys.exit(1)

    if variant_caller is None:
        print(
            "#####\n\033[91mError\033[0m: {} could not be found in your path\n#####\n".format(
                variant_caller
            )
        )
        sys.exit(1)

    consensus_freq = 0.9
    if args.consensus_freq:
        consensus_freq = float(args.consensus_freq)
        if consensus_freq > 1 or consensus_freq < 0:
            print(
                "#####\n\033[91mError\033[0m: The consensus_freq option must be a float number between 0 and 1\n#####\n"
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

    if not os.path.isfile(args.reference + ".fai"):
        print("Indexing {} with samtools".format(args.reference))
        os.system("samtools faidx {} 2> /dev/null".format(args.reference))
    if not os.path.isfile(args.reference + ".bwt"):
        print("Indexing {} with bwa".format(args.reference))
        os.system("bwa index {} 2> /dev/null".format(args.reference))

    snakefile = os.path.join(thisdir, "bin", "Snakefile")

    config = {
        "reads_dir": os.path.join(os.getcwd(), args.query_directory[0]),
        "kraken2_db": default_kraken,
        "outdir": outdir,
        "reference": args.reference,
        "isolates": isolates,
        "coverage_script": coverage_script,
        "converter_script": converter_script,
        "vcf_script": vcf_script,
        "variant_program": variant_caller,
        "scheme": scheme,
        "isolates": isolates,
        "suffix": suffix,
        "threads": args.threads,
        "consensus_freq": consensus_freq,
        "indel_freq": indel_freq,
        "keep_reads": args.keep_reads,
        "verbose": args.verbose,
    }

    if args.print_dag:
        flat_config = []
        for key in config:
            flat_config.append(key + "=" + str(config[key]))
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
        )
    elif config["verbose"]:
        print("\n**** CONFIG ****")
        for k in config:
            print(k + ": ", config[k])

        status = snakemake.snakemake(
            snakefile,
            use_conda=True,
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
        )
    else:
        status = snakemake.snakemake(
            snakefile,
            use_conda=True,
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
