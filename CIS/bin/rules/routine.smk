#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 12 12:59:25 2021

@author: cfos
"""
########################
# Clinical sample workflow
########################
import re
from datetime  import date
TODAY = date.today().strftime("%Y-%m-%d")

rule bwa_map_sort:
    input:
        r1=os.path.join(RESULT_DIR, "{sample}/fastp/{sample}_trimmed_R1.fq.gz"),
        r2=os.path.join(RESULT_DIR, "{sample}/fastp/{sample}_trimmed_R2.fq.gz"),
    output:
        bam=os.path.join(RESULT_DIR, "{sample}/{sample}.sorted.bam"),
    message:
        "mapping {wildcards.sample} reads to reference"
    threads: 4
    log:
        os.path.join(RESULT_DIR, "{sample}/{sample}.bwa.log"),
    params:
        reference=REFERENCE,
    resources:
        cpus=4,
    wildcard_constraints:
        sample="(?!NC)(?!NEG).*",
    shell:
        """
        bwa mem -t {threads} {params.reference} \
        {input.r1} {input.r2} 2> {log} | \
        samtools sort -@ {threads} 2> {log} | \
        samtools view -@ {threads} -F 4 --write-index -o {output.bam} 2> /dev/null
        samtools index {output.bam}
        """


rule ivar_trim:
    input:
        bam=os.path.join(RESULT_DIR, "{sample}/{sample}.sorted.bam"),
    output:
        sort_trim_bam=os.path.join(
            RESULT_DIR, "{sample}/ivar/{sample}.primertrim.sorted.bam"
        ),
    message:
        "mapping {wildcards.sample} reads to reference"
    threads: 4
    log:
        os.path.join(RESULT_DIR, "{sample}/ivar/{sample}.ivar_trim.log"),
    params:
        trim_bam=temp(os.path.join(RESULT_DIR, "{sample}/ivar/{sample}.primertrim.bam")),
        bedfile=SCHEME,
        offset=IVAR_OFFSET,
        prefix=os.path.join(RESULT_DIR, "{sample}/ivar/{sample}.primertrim"),
    conda:
        "../envs/ivar.yaml"
        #"docker://quay.io/biocontainers/ivar:1.3.1--h3198e80_1"
    resources:
        cpus=4,
    wildcard_constraints:
        sample="(?!NC)(?!NEG).*",
    shell:
        """
        touch {output.sort_trim_bam}
        ivar trim -i {input.bam} -x {params.offset} -b {params.bedfile} -p {params.prefix} -e 2&> {log}
        samtools sort -@ {threads} {params.trim_bam} -o {output.sort_trim_bam} 2> /dev/null
        samtools index {output.sort_trim_bam}
        """


rule genomecov:
    input:
        bam=os.path.join(RESULT_DIR, "{sample}/ivar/{sample}.primertrim.sorted.bam"),
    output:
        main_plot=report(
            os.path.join(RESULT_DIR, "{sample}/coverage/{sample}.coverage_plots.pdf"),
            caption="../report/coverage_plots.rst",
            category="Coverage Plots",
        ),
        amp_plot=os.path.join(
            RESULT_DIR, "{sample}/coverage/{sample}.check_amplicons.pdf"
        ),
    message:
        "getting genome coverage statistics for {wildcards.sample}"
    threads: 1
    params:
        coverage=temp(
            os.path.join(RESULT_DIR, "{sample}/coverage/{sample}.coverage.txt")
        ),
        script=COVERAGE_SCRIPT,
        bedfile=SCHEME,
        scheme=SCHEME_NAME,
    resources:
        cpus=1,
    wildcard_constraints:
        sample="(?!NC)(?!NEG).*",
    shell:
        """
        touch {output.amp_plot}
        bedtools genomecov -ibam {input.bam} -d > {params.coverage}  2> /dev/null
        Rscript {params.script} {params.coverage} {params.bedfile} {output.main_plot} {output.amp_plot} {params.scheme} {wildcards.sample} 2> /dev/null
        """


rule qualimap:
    input:
        bam=os.path.join(RESULT_DIR, "{sample}/ivar/{sample}.primertrim.sorted.bam"),
    output:
        report=os.path.join(RESULT_DIR, "{sample}/bamqc/genome_results.txt"),
    message:
        "getting bamQC metrics for {wildcards.sample}"
    threads: 1
    params:
        outdir=os.path.join(RESULT_DIR, "{sample}/bamqc"),
    resources:
        cpus=1,
    wildcard_constraints:
        sample="(?!NC)(?!NEG).*",
    run:
        cmd = "qualimap bamqc -bam {0} -nt {1} -outdir {2} 2&> /dev/null".format(
            input.bam, threads, params.outdir
        )
        p = subprocess.Popen(
            cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )
        stdout, stderr = p.communicate()
        if "Failed to run bamqc" in str(stderr):
            res = [
                "BamQC report",
                "There is a 0% of reference with a coverageData >= 10X",
                "number of reads = 0",
                "number of mapped bases = 0 bp",
                "NC_045512.2	29903	0	0	0",
            ]
            if not os.path.exists(params.outdir):
                os.makedirs(params.outdir)
            with open(output.report, "w") as f:
                f.write("\n".join(map(str, res)))


rule samtools_mpileup:
    input:
        bam=os.path.join(RESULT_DIR, "{sample}/ivar/{sample}.primertrim.sorted.bam"),
    output:
        mpileup=temp(os.path.join(RESULT_DIR, "{sample}/variants/{sample}.mpileup")),
    message:
        "generating mpileup for {wildcards.sample}"
    threads: 4
    log:
        os.path.join(RESULT_DIR, "{sample}/ivar/{sample}.mpileup.log"),
    params:
        reference=REFERENCE,
    wildcard_constraints:
        sample="(?!NC)(?!NEG).*",
    resources:
        cpus=4,
    shell:
        """
        samtools mpileup -A -d 0 -B -Q 0 --reference {params.reference} {input.bam} -o {output.mpileup} 2>{log}
        """


rule ivar_variants:
    input:
        mpileup=os.path.join(RESULT_DIR, "{sample}/variants/{sample}.mpileup"),
    output:
        temp_vcf_file=temp(
            os.path.join(RESULT_DIR, "{sample}/variants/{sample}.tmp.vcf")
        ),
    message:
        "calling variants for {wildcards.sample}"
    threads: 4
    log:
        os.path.join(RESULT_DIR, "{sample}/variants/{sample}.ivar_variants.log"),
    params:
        tsv_file=os.path.join(RESULT_DIR, "{sample}/variants/{sample}.ivar.tsv"),
        bedfile=SCHEME,
        prefix=os.path.join(RESULT_DIR, "{sample}/variants/{sample}.ivar"),
        reference=REFERENCE,
        converter=CONVERTER,
        qual=20,
        depth=10,
        con_freq=CON_FREQ,
        snv_freq=0.1,
    wildcard_constraints:
        sample="(?!NC)(?!NEG).*",
    conda:
        "../envs/ivar.yaml"
    resources:
        cpus=4,
    shell:
        """
        cat {input.mpileup} | ivar variants -p {params.prefix} -r {params.reference} -q {params.qual} -m {params.depth} -t {params.snv_freq} 2&>{log}
        python {params.converter} {params.tsv_file} {output.temp_vcf_file} 2> {log}
        """


#rule lofreq_variants_orig:
#    input:
#        bam=os.path.join(RESULT_DIR, "{sample}/ivar/{sample}.primertrim.sorted.bam"),
#    output:
#        new_bam=temp(os.path.join(RESULT_DIR, "{sample}/variants/{sample}.tmp.bam")),
#        new_bam_index=temp(
#            os.path.join(RESULT_DIR, "{sample}/variants/{sample}.tmp.bam.bai")
#        ),
#        tmp_vcf1=temp(
#            os.path.join(RESULT_DIR, "{sample}/variants/{sample}.tmp1.vcf.gz")
#        ),
#        tmp_vcf2=temp(
#            os.path.join(RESULT_DIR, "{sample}/variants/{sample}.tmp2.vcf.gz")
#        ),
#        tmp_vcf1_index=temp(
#            os.path.join(RESULT_DIR, "{sample}/variants/{sample}.tmp1.vcf.gz.tbi")
#        ),
#        tmp_vcf2_index=temp(
#            os.path.join(RESULT_DIR, "{sample}/variants/{sample}.tmp2.vcf.gz.tbi")
#        ),
#    message:
#        "calling variants for {wildcards.sample}"
#    threads: 8
#    log:
#        os.path.join(RESULT_DIR, "{sample}/variants/{sample}.lofreq.log"),
#    params:
#        vcf_script=VCF_MOD,
#        prefix="{sample}",
#        reference=REFERENCE,
#        depth=10,
#    wildcard_constraints:
#        sample="(?!NC)(?!NEG).*",
#    conda:
#        "../envs/lofreq.yaml"
#    resources:
#        cpus=8,
#    shell:
#        """
#        lofreq indelqual --dindel {input.bam} -f {params.reference} | \
#        samtools sort -@ {threads} -o {output.new_bam} 2> /dev/null
#        samtools index {output.new_bam}
#        lofreq call-parallel --no-baq --call-indels --pp-threads {threads} \
#        -f {params.reference} -o {output.tmp_vcf1} {output.new_bam} 2> {log}
#        bash {params.vcf_script} -i {output.tmp_vcf1} -g 1 -o {output.tmp_vcf2}
#        """

rule lofreq_variants:
    input:
        bam=os.path.join(RESULT_DIR, "{sample}/ivar/{sample}.primertrim.sorted.bam"),
    output:
        new_bam=temp(os.path.join(RESULT_DIR, "{sample}/variants/{sample}.tmp.bam")),
        new_bam_index=temp(
            os.path.join(RESULT_DIR, "{sample}/variants/{sample}.tmp.bam.bai")
        ),
        new_bam2=temp(os.path.join(RESULT_DIR, "{sample}/variants/{sample}.tmp2.bam")),
        new_bam_index2=temp(
            os.path.join(RESULT_DIR, "{sample}/variants/{sample}.tmp2.bam.bai")
        ),
        new_bam3=temp(os.path.join(RESULT_DIR, "{sample}/variants/{sample}.tmp3.bam")),
        new_bam_index3=temp(
            os.path.join(RESULT_DIR, "{sample}/variants/{sample}.tmp3.bam.bai")
        ),
        tmp_vcf1=temp(
            os.path.join(RESULT_DIR, "{sample}/variants/{sample}.tmp1.vcf.gz")
        ),
        tmp_vcf2=temp(
            os.path.join(RESULT_DIR, "{sample}/variants/{sample}.tmp2.vcf.gz")
        ),
        tmp_vcf1_index=temp(
            os.path.join(RESULT_DIR, "{sample}/variants/{sample}.tmp1.vcf.gz.tbi")
        ),
        tmp_vcf2_index=temp(
            os.path.join(RESULT_DIR, "{sample}/variants/{sample}.tmp2.vcf.gz.tbi")
        ),
    message:
        "calling variants for {wildcards.sample}"
    threads: 8
    log:
        os.path.join(RESULT_DIR, "{sample}/variants/{sample}.lofreq.log"),
    params:
        vcf_script=VCF_MOD,
        prefix="{sample}",
        reference=REFERENCE,
        depth=10,
    wildcard_constraints:
        sample="(?!NC)(?!NEG).*",
    conda:
        "../envs/lofreq.yaml"
    resources:
        cpus=8,
    shell:
        """
        lofreq viterbi -f {params.reference} {input.bam} | \
        samtools sort -@ {threads} -o {output.new_bam} 2> /dev/null
        samtools index {output.new_bam}
        lofreq indelqual --dindel {output.new_bam} -f {params.reference} -o {output.new_bam2} 2> /dev/null
        samtools index {output.new_bam2}
        lofreq alnqual -b {output.new_bam2} {params.reference} > {output.new_bam3}
        samtools index {output.new_bam3}
        lofreq call-parallel --no-baq --call-indels --pp-threads {threads} \
        -f {params.reference} -o {output.tmp_vcf1} {output.new_bam3} 2> {log}
        bash {params.vcf_script} -i {output.tmp_vcf1} -g 1 -o {output.tmp_vcf2}
        """

rule lofreq_bcftools_setGT:
    input:
        vcf_file=os.path.join(RESULT_DIR, "{sample}/variants/{sample}.tmp2.vcf.gz"),
    output:
        vcf_file=os.path.join(RESULT_DIR, "{sample}/variants/{sample}.lofreq.vcf.gz"),
    log:
        os.path.join(RESULT_DIR, "{sample}/variants/{sample}.bcftools_setGT.log"),
    params:
        snv_freq=0.1,
        con_freq=CON_FREQ,
        indel_freq=INDEL_FREQ,
        min_depth=config['min_depth'],
    message:
        "setting conditional GT for {wildcards.sample}"
    wildcard_constraints:
        sample="(?!NC)(?!NEG).*",
    shell:
        """
        bcftools index {input.vcf_file}
        bcftools +fill-tags {input.vcf_file} -Ou -- -t "TYPE" | \
        bcftools norm -Ou -a -m -  2> /dev/null | \
        bcftools view -f 'PASS,.' -i "INFO/AF >= {params.snv_freq} & INFO/DP >= {params.min_depth}" -Oz -o {output.vcf_file}
        bcftools index {output.vcf_file}
        bcftools +setGT {output.vcf_file} -- -t q -i 'GT="1" && INFO/AF < {params.con_freq}' -n 'c:0/1' 2>> {log} | \
        bcftools +setGT -- -t q -i 'TYPE="indel" && INFO/AF < {params.indel_freq}' -n . 2>> {log} | \
        bcftools +setGT -o {output.vcf_file} -- -t q -i 'GT="1" && INFO/AF >= {params.con_freq}' -n 'c:1/1' 2>> {log}
        bcftools index -f {output.vcf_file}
        """


rule ivar_bcftools_setGT:
    input:
        vcf_file=os.path.join(RESULT_DIR, "{sample}/variants/{sample}.tmp.vcf"),
    output:
        variant_types=temp(
            os.path.join(RESULT_DIR, "{sample}/variants/{sample}.variant_types.txt.gz")
        ),
        hdr=temp(os.path.join(RESULT_DIR, "{sample}/variants/{sample}.hdr.txt")),
        vcf_file=os.path.join(RESULT_DIR, "{sample}/variants/{sample}.ivar.vcf.gz"),
        tmp_vcf=temp(os.path.join(RESULT_DIR, "{sample}/variants/{sample}.tmp.vcf.gz")),
    log:
        os.path.join(RESULT_DIR, "{sample}/variants/{sample}.bcftools_setGT.log"),
    params:
        snv_freq=0.1,
        con_freq=CON_FREQ,
        indel_freq=INDEL_FREQ,
        min_depth=config['min_depth'],
    message:
        "setting conditional GT for {wildcards.sample}"
    wildcard_constraints:
        sample="(?!NC)(?!NEG).*",
    shell:
        """
        bgzip -f {input.vcf_file}
        bcftools index {output.tmp_vcf}
        bcftools query {output.tmp_vcf} -f'%CHROM\t%POS\t%TYPE\n' | bgzip > {output.variant_types}
        tabix -s1 -b2 -e2 {output.variant_types}
        echo '##INFO=<ID=TYPE,Number=.,Type=String,Description="Variant type">' > {output.hdr}
        bcftools annotate -a {output.variant_types} -h {output.hdr} -c CHROM,POS,TYPE -Oz -o {output.tmp_vcf} {output.tmp_vcf}
        bcftools index -f {output.tmp_vcf}
        bcftools norm -Ou -a -m - {output.tmp_vcf} 2> {log} | \
        bcftools view -i "FORMAT/ALT_FREQ >= {params.snv_freq} & FORMAT/ALT_QUAL >= 20 & INFO/DP >= {params.min_depth}" -Oz -o {output.vcf_file}
        bcftools index {output.vcf_file}
        bcftools +setGT {output.vcf_file} -- -t q -i 'GT="1" && FORMAT/ALT_FREQ < {params.con_freq} & INFO/DP >= {params.min_depth}' -n 'c:0/1' 2> {log}| \
        bcftools +setGT -- -t q -i 'TYPE="INDEL" && FORMAT/ALT_FREQ < {params.indel_freq}' -n . 2>> {log} | \
        bcftools +setGT -o {output.vcf_file} -- -t q -i 'GT="1" && FORMAT/ALT_FREQ >= {params.con_freq} & INFO/DP >= {params.min_depth}' -n 'c:1/1' 2> {log}
        bcftools index -f {output.vcf_file}
        """


if VARIANT_PROGRAM == "ivar":
    proper_vcf = os.path.join(RESULT_DIR, "{sample}/variants/{sample}.ivar.vcf.gz")
elif VARIANT_PROGRAM == "lofreq":
    proper_vcf = os.path.join(RESULT_DIR, "{sample}/variants/{sample}.lofreq.vcf.gz")


rule generate_consensus:
    input:
        vcf_file=proper_vcf,
        bam=os.path.join(RESULT_DIR, "{sample}/ivar/{sample}.primertrim.sorted.bam"),
    output:
        consensus=os.path.join(RESULT_DIR, "{sample}/variants/{sample}.consensus.fa"),
        mask=temp(os.path.join(RESULT_DIR, "{sample}/variants/mask.bed")),
        variants_bed=temp(os.path.join(RESULT_DIR, "{sample}/variants/variants.bed")),
    message:
        "calling a consensus for {wildcards.sample}"
    threads: 1
    log:
        os.path.join(RESULT_DIR, "{sample}/variants/{sample}.consensus.log"),
    params:
        prefix="{sample}",
        reference=REFERENCE,
        con_freq=CON_FREQ,
        min_depth=config['min_depth'],
    wildcard_constraints:
        sample="(?!NC)(?!NEG).*",
    resources:
        cpus=1,
    shell:
        """
        bcftools query -f'%CHROM\t%POS0\t%END\n' {input.vcf_file} > {output.variants_bed}
        varCheck=$(file {output.variants_bed} | cut -f2 -d " ")
        if [ $varCheck == "empty" ]; then
            echo "BAM file was empty for {params.prefix}. Making empty consensus genome."
            emptySeq=$(printf %.1s N{{1..29903}})
            printf ">{params.prefix}\n$emptySeq\n" > {output.consensus}
            touch {output.mask}
            touch {output.variants_bed}
        else
            bedtools genomecov -bga -ibam {input.bam} | awk '$4 < {params.min_depth}' | \
            bedtools subtract -a - -b {output.variants_bed} > {output.mask}
            bcftools consensus -p {params.prefix} -f {params.reference} --mark-del '-' -m {output.mask} -H I -i 'INFO/DP >= {params.min_depth} & GT!="mis"' {input.vcf_file} 2> {log} | \
            sed "/^>/s/{params.prefix}.*/{params.prefix}/" > {output.consensus}
        fi
        """


rule snpeff:
    input:
        vcf_file=proper_vcf,
    output:
        vcf_file=os.path.join(RESULT_DIR, "{sample}/variants/{sample}.annotated.vcf"),
    message:
        "annotating variants for {wildcards.sample}"
    threads: 1
    params:
        prefix=os.path.join(RESULT_DIR, "{sample}/variants/{sample}.snpeff"),
    wildcard_constraints:
        sample="(?!NC)(?!NEG).*",
    resources:
        cpus=1,
    shell:
        """
        snpEff eff -csvStats {params.prefix}.stats.csv -s {params.prefix}.stats.html NC_045512.2 {input.vcf_file} > {output.vcf_file}
        """

rule update_pangolin:
    output:
        update_info = os.path.join(RESULT_DIR, "pangolin_update_info.txt"),
    conda:
        "../envs/pangolin.yaml"
    shell:
        """
        echo "Attempting to update pangolin..." > {output.update_info}
        pangolin --update &>> {output.update_info} || echo "pangolin couldn't update. Investigate." >> {output.update_info}
        echo "Pangolin versions after attempted update..." >> {output.update_info}
        pangolin --all-versions &>> {output.update_info}
        """

rule pangolin:
    input:
        update_info = os.path.join(RESULT_DIR, "pangolin_update_info.txt"),
        fasta=os.path.join(RESULT_DIR, "{sample}/variants/{sample}.consensus.fa"),
    output:
        report=os.path.join(RESULT_DIR, "{sample}/pangolin/{sample}.lineage_report.csv"),
    conda:
        "../envs/pangolin.yaml"
    shell:
        """
        pangolin --outfile {output.report} --skip-scorpio {input.fasta} &> /dev/null
        """

rule update_nextclade:
    output:
        update_info = os.path.join(RESULT_DIR, "nextclade_update_info.txt"),
    params:
        nextclade_dataset = config['nextclade_dataset']
    container:
        "docker://nextstrain/nextclade:2.3.0"
    shell:
        """
        echo "nextclade version:" > {output.update_info}
        nextclade --version >> {output.update_info} &>/dev/null
        echo "Updating SARS-CoV-2 dataset..." >> {output.update_info}
        nextclade dataset get --name sars-cov-2 -o {params.nextclade_dataset} &>>{output.update_info}
        """

rule nextclade:
    input:
        update_info = os.path.join(RESULT_DIR, "nextclade_update_info.txt"),
        fasta=os.path.join(RESULT_DIR, "{sample}/variants/{sample}.consensus.fa"),
    output:
        report=os.path.join(RESULT_DIR, "{sample}/nextclade/{sample}.nextclade_report.tsv"),
    params:
        nextclade_dataset = config['nextclade_dataset'],
        outdir = os.path.join(RESULT_DIR, "{sample}/nextclade"),
    log:
        os.path.join(RESULT_DIR, "{sample}/nextclade/{sample}.nextclade.log"),
    container:
        "docker://nextstrain/nextclade:2.3.0",
    resources:
        cpus=1,
    threads: 4,
    shell:
        """
        nextclade run --in-order --input-dataset={params.nextclade_dataset} --output-basename={wildcards.sample} --output-all={params.outdir} --output-tsv={output.report} --jobs {threads} {input.fasta} &> {log}
        """


rule sample_qc:
    input:
        bamqc=os.path.join(RESULT_DIR, "{sample}/bamqc/genome_results.txt"),
        lineages=os.path.join(
            RESULT_DIR, "{sample}/pangolin/{sample}.lineage_report.csv"
        ),
        nextclade_report=os.path.join(
            RESULT_DIR, "{sample}/nextclade/{sample}.nextclade_report.tsv"
        ),
    output:
        report=temp(os.path.join(RESULT_DIR, "{sample}.qc_results.csv")),
    params:
        sample="{sample}",
        technology=TECHNOLOGY,
        legacy=config['legacy_results'],
        date=TODAY,
    wildcard_constraints:
        sample="(?!NC)(?!NEG).*",
    run:

        with open(input.bamqc,'r') as f:
            parsed = [x.strip() for x in f.readlines()]

        parsed = [x for x in parsed  if x != '']

        for line in parsed:
            if 'number of reads' in line:
                No_reads = int(re.sub('^.*= ', '',line).replace(",",""))
            if 'number of mapped bases' in line:
                No_bases = int(re.sub('^.*= | bp', '',line).replace(",",""))
            if 'coverageData >= 10X' in line:
                ref = str(re.sub('There is a |% of reference .*', '',line).replace(",",""))

        mean_cov = float(parsed[-1].split()[3])
        std = float(parsed[-1].split()[4])

        df = pd.DataFrame(
            columns=[
                "id",
                "num_reads",
                "num_bases",
                "coverage",
                "mean_depth",
                "stdev_depth",
                "QC",
                "technology",
                "analysis_date",
            ]
        )

        if int(float(ref.split("%")[0])) > 79:
            df.loc[0] = [
                params.sample,
                No_reads,
                No_bases,
                ref,
                mean_cov,
                std,
                "PASS",
                params.technology,
                params.date,
            ]
        else:
            df.loc[0] = [
                params.sample,
                No_reads,
                No_bases,
                ref,
                mean_cov,
                std,
                "FAIL",
                params.technology,
                params.date,
            ]

        #get pangolin info
        lineages = pd.read_csv(input.lineages).drop(
                ["ambiguity_score", "scorpio_conflict", "scorpio_notes", "is_designated", "qc_notes"], axis=1
            )

        lineages.columns = [
            "id",
            "lineage",
            "lineage_conflict",
            "scorpio_call",
            "scorpio_support",
            "lineage_designation_version",
            "pangolin_version",
            "scorpio_version",
            "constellation_version",
            "pangolin_status",
            "pangolin_note",
        ]

        df = df.join(lineages.set_index(["id"]), on=["id"])

        # read in the nextclade files
        nextclade = pd.read_csv(input.nextclade_report, sep="\t")
        nextclade.rename(columns={'seqName':'id','qc.privateMutations.total':'totalPrivateMutations','qc.overallStatus':'Nextclade_QC'}, inplace=True)
        keep = ['id','clade','Nextclade_pango','Nextclade_QC','totalFrameShifts','totalAminoacidInsertions','totalAminoacidDeletions','totalAminoacidSubstitutions','totalNonACGTNs','totalPrivateMutations']
        nextclade = nextclade.loc[:,keep]

        df = df.join(nextclade.set_index(["id"]), on=["id"])

        # change column order
        compulsory = ["id","num_reads","num_bases","coverage","mean_depth","stdev_depth","lineage","scorpio_call","QC"]
        #extra = list(set(df.columns) - set(compulsory))
        extra = [
            "technology",
            "analysis_date",
            "scorpio_support",
            "scorpio_version",
            "lineage_conflict",
            "lineage_designation_version",
            "pangolin_version",
            "pangolin_status",
            "pangolin_note",
            "clade",
            "Nextclade_pango",
            'Nextclade_QC',
            'totalFrameShifts',
            'totalAminoacidInsertions',
            'totalAminoacidDeletions',
            'totalAminoacidSubstitutions',
            'totalNonACGTNs',
            'totalPrivateMutations',
        ]

        final = compulsory+extra
        if not params.legacy:
            df = df.loc[:, final]
            df['neg_control'] = False
        else:
            df = df.loc[:, compulsory]
            df.columns = [
                    "ID",
                    "NoReads",
                    "NoBases",
                    "%Ref",
                    "Mean Cov",
                    "std",
                    "pangolin lineage",
                    "scorpio_call",
                    "QC",
                ]

        df.to_csv(output.report, header=True, index=False)


########################
# Rules - negative controls
########################


rule kraken2:
    input:
        forward_reads=os.path.join(
            RESULT_DIR, "{sample}/fastp/{sample}_trimmed_R1.fq.gz"
        ),
        reverse_reads=os.path.join(
            RESULT_DIR, "{sample}/fastp/{sample}_trimmed_R2.fq.gz"
        ),
    output:
        result=os.path.join(RESULT_DIR, "{sample}/kraken2/{sample}.kraken_result.txt"),
        report=os.path.join(RESULT_DIR, "{sample}/kraken2/{sample}.kraken_report.txt"),
    log:
        os.path.join(RESULT_DIR, "{sample}/kraken2/{sample}.kraken2.log"),
    params:
        KRAKEN2_DB=KRAKEN2_DB,
    wildcard_constraints:
        sample="(NC|NEG).*",
    threads: 4
    resources:
        mem_mb=55000,
        cpus=4,
    shell:
        """
        touch '{output.result}'
        touch '{output.report}'
        kraken2 --threads 4 \
        --gzip-compressed \
        '{input.forward_reads}' '{input.reverse_reads}' \
        --output '{output.result}' \
        --report '{output.report}' \
        --db '{params.KRAKEN2_DB}' \
        2>{log}
        """


rule neg_qc:
    input:
        kraken=os.path.join(RESULT_DIR, "{sample}/kraken2/{sample}.kraken_result.txt"),
    output:
        report=temp(os.path.join(RESULT_DIR, "{sample}.qc_results.csv")),
    params:
        today=TODAY,
        sample="{sample}",
        technology=TECHNOLOGY,
        legacy=config["legacy_results"],
    wildcard_constraints:
        sample="(NC|NEG).*",
    run:
        if params.legacy:
            df = pd.DataFrame(
                columns=[
                    "ID",
                    "NoReads",
                    "NoBases",
                    "%Ref",
                    "Mean Cov",
                    "std",
                    "pangolin lineage",
                    "scorpio_call",
                    "QC",
                ]
            )
            try:
                neg = pd.read_csv(input.kraken, sep="\t", header=None)
                if sum(neg[3]) > 10000:
                    df.loc[0] = [
                        params.sample,
                        neg.shape[0],
                        sum(neg[3]),
                        "-",
                        "-",
                        "-",
                        "-",
                        "-",
                        "FAIL",
                    ]
                else:
                    df.loc[0] = [
                        params.sample,
                        neg.shape[0],
                        sum(neg[3]),
                        "-",
                        "-",
                        "-",
                        "-",
                        "-",
                        "PASS",
                    ]
            except:
                df.loc[0] = [params.sample, "No_reads", 0, "-", "-", "-", "-", "-", "PASS"]
        else:
            df = pd.DataFrame(
                columns=[
                    'id',
                    'num_reads',
                    'num_bases',
                    'coverage',
                    'mean_depth',
                    'stdev_depth',
                    'lineage',
                    'scorpio_call',
                    'QC',
                    'technology',
                    'analysis_date',
                    'pangoLEARN_version',
                    'pango_version',
                    'pangolin_version',
                    'scorpio_support',
                    'pangolin_ambiguity_score',
                    'pangolin_status',
                    'lineage_designation_version',
                    "pangolin_note",
                    "neg_control",
                ]
            )
            try:
                neg = pd.read_csv(input.kraken, sep="\t", header=None)
                if sum(neg[3]) > 10000:
                    df.loc[0] = [
                        params.sample,
                        neg.shape[0],
                        sum(neg[3]),
                        "-",
                        "-",
                        "-",
                        "-",
                        "-",
                        "FAIL",
                        params.technology,
                        params.today,
                        "-",
                        "-",
                        "-",
                        "-",
                        "-",
                        "-",
                        "-",
                        "-",
                        True,
                    ]
                else:
                    df.loc[0] = [
                        params.sample,
                        neg.shape[0],
                        sum(neg[3]),
                        "-",
                        "-",
                        "-",
                        "-",
                        "-",
                        "PASS",
                        params.technology,
                        params.today,
                        "-",
                        "-",
                        "-",
                        "-",
                        "-",
                        "-",
                        "-",
                        "-",
                        True,
                    ]
            except:
                df.loc[0] = [params.sample, "No_reads", 0, "-", "-", "-", "-", "-", "PASS", params.today, params.technology, "-", "-", "-", "-",  "-", "-", "-","-",True]

        df.to_csv(output.report, header=True, index=False)
