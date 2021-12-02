#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 12 12:50:44 2021

@author: cfos
"""

###########
# Libraries
###########

import os
import glob
import re
from datetime import date
import pandas as pd
import yaml
import subprocess
import glob
import pandas as pd

#################
# Custom functions
#################


def find_input_data(READSDIR, data_type, SUFFIX):
    reads = glob.glob(os.path.join(READSDIR, "*" + SUFFIX))
    main_exclusion = ["Undet", "NEG_", "NC_"]
    neg_inclusion = ["NEG_", "NC_"]
    if data_type == "main":
        reads = [i for i in reads if not any(b in i for b in main_exclusion)]
    elif data_type == "neg":
        reads = [i for i in reads if any(b in i for b in neg_inclusion)]
    result = []
    for read in reads:
        if SUFFIX == "_L001_R1_001.fastq.gz":
            sample = re.sub("_S\d+_L001.*", repl="", string=os.path.basename(read))
            barcode = re.search("S\d+_L001.*", string=os.path.basename(read)).group()
            barcode = re.sub("_.*", repl="", string=barcode)
            res = {
                "sample": sample,
                "barcode": barcode,
                "fq1": read,
                "fq2": read.replace("_R1", "_R2"),
            }
        else:
            sample = re.sub(SUFFIX, repl="", string=os.path.basename(read))
            barcode = "N/A"
            res = {
                "sample": sample,
                "barcode": barcode,
                "fq1": read,
                "fq2": read.replace("_R1", "_R2"),
            }
        result.append(res)
    result = pd.DataFrame(result, dtype=str).set_index(["sample"], drop=False)
    return result


def get_trim_names(wildcards):
    """
    This function:
      1. Returns the correct input and output trimmed file names for fastp.
    """
    inFile = INPUT_TABLE.loc[(wildcards.sample), ["fq1", "fq2"]].dropna()
    return (
        "--in1 "
        + inFile[0]
        + " --in2 "
        + inFile[1]
        + " --out1 "
        + os.path.join(
            RESULT_DIR,
            wildcards.sample,
            "fastp",
            wildcards.sample + "_trimmed_R1.fq.gz",
        )
        + " --out2 "
        + os.path.join(
            RESULT_DIR,
            wildcards.sample,
            "fastp",
            wildcards.sample + "_trimmed_R2.fq.gz",
        )
    )


def get_fastq(wildcards):
    """This function returns the forward and reverse fastq files for samples"""
    return INPUT_TABLE.loc[(wildcards.sample), ["fq1", "fq2"]].dropna()


###############
# Configuration
###############

READSDIR = config["reads_dir"]
RESULT_DIR = config["outdir"]
REFERENCE = config["reference"]
SCHEME = config["scheme"]
COVERAGE_SCRIPT = config["coverage_script"]
CONVERTER = config["converter_script"]
VCF_MOD = config["vcf_script"]
ISOLATES = config["isolates"]
VARIANT_PROGRAM = os.path.basename(config["variant_program"])
SCHEME_NAME = os.path.basename(SCHEME)
KRAKEN2_DB = config["kraken2_db"]
CON_FREQ = config["consensus_freq"]
KEEP_READS = config["keep_reads"]

SUFFIX = config["suffix"]
SUFFIX = SUFFIX.replace("R2", "R1")

MAIN_TABLE = find_input_data(READSDIR, "main", SUFFIX)
MAIN_TABLE.to_csv(
    os.path.join(RESULT_DIR, "sample_sheet.csv"), index=False, header=True
)
MAIN_SAMPLES = MAIN_TABLE.index.get_level_values("sample").unique().tolist()

try:
    NEG_TABLE = find_input_data(READSDIR, "neg", SUFFIX)
    NEG_SAMPLES = NEG_TABLE.index.get_level_values("sample").unique().tolist()
except:
    print("\nNo negative control samples detected\n")
    NEG_SAMPLES = []
INPUT_TABLE = find_input_data(READSDIR, "all", SUFFIX)

if ISOLATES != False:
    MAIN_SAMPLES = [i for i in MAIN_SAMPLES if any(b in i for b in ISOLATES)]
    NEG_SAMPLES = [i for i in NEG_SAMPLES if any(b in i for b in ISOLATES)]

config["analysis_samples"] = "; ".join(MAIN_SAMPLES)
config["neg_samples"] = "; ".join(NEG_SAMPLES)

if SCHEME_NAME == "swift.bed":
    IVAR_OFFSET = 5
else:
    IVAR_OFFSET = 0

if not os.path.isfile(REFERENCE + ".fai"):
    os.system("samtools faidx {} 2> /dev/null".format(REFERENCE))
if not os.path.isfile(REFERENCE + ".bwt"):
    os.system("bwa index {} 2> /dev/null".format(REFERENCE))

if VARIANT_PROGRAM == "lofreq":
    variant_conda_env = "../envs/lofreq.yaml"
else:
    variant_conda_env = "../envs/ivar.yaml"


TODAY = date.today().strftime("%Y-%m-%d")

################
# Optional removal of trimmed reads (to save space)
################
if KEEP_READS == False:

    onsuccess:
        print("Removing trimmed reads (input reads remain untouched)\n")
        dead_reads = glob.glob(RESULT_DIR + "/**/*trimmed*.gz", recursive=True)
        [os.remove(x) for x in dead_reads]


################
# Software versions
################


rule collect_variant_version:
    input:
        software=os.path.join(RESULT_DIR, TODAY + "_softwareVersions_v1.txt"),
    output:
        software=report(
            os.path.join(RESULT_DIR, TODAY + "_softwareVersions.txt"),
            caption="../report/software_versions.rst",
            category="Software Versions",
        ),
    params:
        variant_program=VARIANT_PROGRAM,
    conda:
        variant_conda_env
    shell:
        """
        cp {input.software} {output.software}
        if [ {params.variant_program} == "lofreq" ]; then
          printf "lofreq\t$(lofreq version | tr "\n" "\t")\n" >> {output.software}
        elif [ {params.variant_program} == "ivar" ]; then
          printf "$(ivar version | head -n1)\n" >> {output.software}
        fi
        """


rule collect_main_versions:
    output:
        software=temp(os.path.join(RESULT_DIR, TODAY + "_softwareVersions_v1.txt")),
    shell:
        """
        fastp -v 2>> {output.software}
        bedtools --version >> {output.software}
        printf "qualimap v.2.2.2-dev\n" >> {output.software}
        bcfVer=$(bcftools --version | head -n2 | tr "\n" " ")
        printf "$bcfVer\n" >> {output.software}
        samVer=$(samtools --version | head -n2 | tr "\n" " ")
        printf "$samVer\n" >> {output.software}
        set +eu
        eval "$(conda shell.bash hook)" && conda activate pangolin && pangolin --all-versions >> {output.software}
        set -eu
        kraken2 --version | head -n1 >> {output.software}
        """


################
# Common read trimming
################


rule fastp:
    input:
        get_fastq,
    output:
        fq1=os.path.join(RESULT_DIR, "{sample}/fastp/{sample}_trimmed_R1.fq.gz"),
        fq2=os.path.join(RESULT_DIR, "{sample}/fastp/{sample}_trimmed_R2.fq.gz"),
        html=os.path.join(RESULT_DIR, "{sample}/fastp/{sample}_fastp.html"),
        json=os.path.join(RESULT_DIR, "{sample}/fastp/{sample}_fastp.json"),
    message:
        "trimming {wildcards.sample} reads"
    threads: 4
    log:
        os.path.join(RESULT_DIR, "{sample}/fastp/{sample}.log.txt"),
    params:
        in_and_out_files=get_trim_names,
        sampleName="{sample}",
        qualified_quality_phred=20,
        cut_mean_quality=20,
        unqualified_percent_limit=10,
        length_required=50,
    resources:
        cpus=4,
    shell:
        """
        touch {output.fq2}
        fastp --thread {threads} \
        --detect_adapter_for_pe \
        --cut_front \
        --cut_tail \
        --trim_poly_x \
        --cut_mean_quality 20 \
        --qualified_quality_phred {params.qualified_quality_phred} \
        --unqualified_percent_limit {params.unqualified_percent_limit} \
        --correction \
        --length_required 50 \
        --html {output.html} \
        --json {output.json} \
        {params.in_and_out_files} \
        2>{log}
        """
