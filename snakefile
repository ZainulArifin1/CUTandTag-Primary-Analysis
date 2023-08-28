########################################################
#### snakemake workflow for CUT&Tag Primary Analysis ###
########################################################

__author__ = "Zain Arifin"
__affiliation__ = "University College Dublin and University California Santa Barbara"
__email__ = "muhammad.arifin@ucdconnect.ie"
__github__ = "github.com/ZainulArifin1"
__date__ = "August 17, 2023"

#################
### IMPORTANT ###
#################

### "data/h3k27me3_h3k4/raw_fastq/" is the directory of the raw FASTQ files.
### Please adjust the directory names accordingly for all instances in the code. You can do find and replace.
### They all must be in the same folder.
### Please also adjust the naming of FASTQ files.
### Example: "data/h3k27me3_h3k4/raw_fastq/{sra}_1.fastq". In this case the forward read is denoted as "_1" and the file extension is .fastq (can be fastq.gz)
### This instance is in rule bowtie2 (line 72)
### Thats it! You can run the code with "snakemake --cores <num_of_cores> -s <snakefile_name>"
### It is advised to do a dry run first. This can be done with "snakemake -np -s <snakefile_name>"

import os
import glob

SRA,FRR = glob_wildcards("data/h3k27me3_h3k4/raw_fastq/{sra}_{frr}.fastq")

# --effectiveGenomeSize 2913022398 before blacklist, 2685859998 after hg38 - blacklist

rule all:
    input:
        expand("output/h3k27me3_h3k4/rawQC/{sra}_{frr}_fastqc.{extension}", sra=SRA, frr=FRR, extension=["zip","html"]), #fastqc output
        "output/h3k27me3_h3k4/rawQC/multiqc_results/multiqc_report.html", #multiQC ouput
        expand("output/h3k27me3_h3k4/bam/{sra}_mapped_sorted_mapq30_BlClean.bam.bai", sra=SRA), # BAM index
        expand("output/bigwig_RPGC/h3k27me3_h3k4_initial_{sra}.bw", sra=SRA),
        expand("output/bigwig_RPGC/h3k27me3_h3k4_MAPQ_RPGC_{sra}.bw", sra=SRA),
        expand("output/bigwig_RPGC_MAPQ_BL/h3k27me3_h3k4_MAPQ_BL_RPGC_{sra}.bw", sra=SRA),
        expand("output/h3k27me3_h3k4/macs2/h3k27me3_h3k4_macs_{sra}_peaks.{extension}", sra=SRA, extension=["broadPeak","gappedPeak","xls"]),
        expand("output/bigwig_TMM_MAPQ_BL/h3k27me3_h3k4_MAPQ_BL_TMM_{sra}.bw", sra=SRA),
        expand("output/bigwig_TMM_MAPQ_BL/scores_per_bin.{extension}", extension=["npz","tab"]),
        "output/bigwig_TMM_MAPQ_BL/SpearmanCorr_readCounts.tab",
        "output/bigwig_TMM_MAPQ_BL/heatmap_SpearmanCorr_readCounts.png",
        "output/bigwig_RPGC_MAPQ_BL/SpearmanCorr_readCounts.tab",
        "output/bigwig_RPGC_MAPQ_BL/heatmap_SpearmanCorr_readCounts.png"


## QC check ##
ruleorder: fastqc > multiQC

rule fastqc: 
    input:
        rawread="data/h3k27me3_h3k4/raw_fastq/{sra}_{frr}.fastq"
    output: 
        zip="output/h3k27me3_h3k4/rawQC/{sra}_{frr}_fastqc.zip",
        html="output/h3k27me3_h3k4/rawQC/{sra}_{frr}_fastqc.html"
    threads:
        6
    params:
        path="output/h3k27me3_h3k4/rawQC/"
    shell:
        "fastqc {input.rawread} --threads {threads} -o {params.path}"


rule multiQC:
    input:
        fastinput=expand("output/h3k27me3_h3k4/rawQC/{sra}_{frr}_fastqc.html", sra=SRA, frr=FRR)
    output:
        "output/h3k27me3_h3k4/rawQC/multiqc_results/multiqc_report.html"
    shell:
        """
        multiqc output/h3k27me3_h3k4/rawQC/ --outdir output/h3k27me3_h3k4/rawQC/multiqc_results/
        """

## Finish QC check ##

## Alignment and file conversion ##
rule bowtie2: #ok
    input:
        fread="data/h3k27me3_h3k4/raw_fastq/{sra}_1.fastq",
        rread="data/h3k27me3_h3k4/raw_fastq/{sra}_2.fastq"
    output:
        sam="output/h3k27me3_h3k4/mapped_reads/{sra}_mapped.sam"
    threads:
        10
    params:
        ref_path="data/reference/GRCh38_noalt_as/GRCh38_noalt_as"
    shell:
        "bowtie2 --end-to-end --very-sensitive --no-mixed --no-discordant --phred33 -I 10 -X 700 -p {threads} -x {params.ref_path} -1 {input.fread} -2 {input.rread} -S {output.sam}"

rule samtools_sort: 
    input:
        "output/h3k27me3_h3k4/mapped_reads/{sra}_mapped.sam"
    output:
        "output/h3k27me3_h3k4/mapped_reads/{sra}_mapped_sorted.bam"
    shell:
        "samtools view -S -b {input} | samtools sort -o {output}"

rule indexBAMInitial:
    input:
        "output/h3k27me3_h3k4/mapped_reads/{sra}_mapped_sorted.bam"
    output:
        "output/h3k27me3_h3k4/mapped_reads/{sra}_mapped_sorted.bam.bai"
    shell:
        """
        samtools index {input}
        """

rule bamCoverageInitial:
    input:
        bam="output/h3k27me3_h3k4/mapped_reads/{sra}_mapped_sorted.bam",
        bam_bai="output/h3k27me3_h3k4/mapped_reads/{sra}_mapped_sorted.bam.bai"
    output:
        "output/bigwig_RPGC/h3k27me3_h3k4_initial_{sra}.bw"
    shell:
        """
        bamCoverage --bam {input.bam} -o {output} --binSize 50 --normalizeUsing RPGC --skipNonCoveredRegions --effectiveGenomeSize 2913022398
        """

rule mapq_filter: 
    input:
        "output/h3k27me3_h3k4/mapped_reads/{sra}_mapped_sorted.bam"
    output:
        "output/h3k27me3_h3k4/mapped_reads/{sra}_mapped_sorted_mapq30.bam"
    shell:
        "samtools view -q 30 -b {input} > {output}"

rule indexBAM_MAPQ_RPGC:
    input:
        "output/h3k27me3_h3k4/mapped_reads/{sra}_mapped_sorted_mapq30.bam"
    output:
        "output/h3k27me3_h3k4/mapped_reads/{sra}_mapped_sorted_mapq30.bam.bai"
    shell:
        """
        samtools index {input}
        """

rule bamCoverage_MAPQ_RPGC:
    input:
        bam="output/h3k27me3_h3k4/mapped_reads/{sra}_mapped_sorted_mapq30.bam",
        bam_bai="output/h3k27me3_h3k4/mapped_reads/{sra}_mapped_sorted_mapq30.bam.bai"
    output:
        "output/bigwig_RPGC/h3k27me3_h3k4_MAPQ_RPGC_{sra}.bw"
    shell:
        """
        bamCoverage --bam {input.bam} -o {output} --binSize 50 --normalizeUsing RPGC --skipNonCoveredRegions --effectiveGenomeSize 2913022398
        """

rule blacklist_removal: 
    input:
        "output/h3k27me3_h3k4/mapped_reads/{sra}_mapped_sorted_mapq30.bam"
    output:
        "output/h3k27me3_h3k4/bam/{sra}_mapped_sorted_mapq30_BlClean.bam"
    params:
        bl_path="data/blacklist/hg38-blacklist.v2.bed"
    shell:
        """
        bedtools intersect -a {input} -b {params.bl_path} -v -ubam > {output}
        """

rule indexBAM_MAPQ_BL_RPGC:
    input:
        "output/h3k27me3_h3k4/bam/{sra}_mapped_sorted_mapq30_BlClean.bam"
    output:
        "output/h3k27me3_h3k4/bam/{sra}_mapped_sorted_mapq30_BlClean.bam.bai"
    shell:
        """
        samtools index {input}
        """

rule bamCoverage_MAPQ_BL_RPGC:
    input:
        bam="output/h3k27me3_h3k4/bam/{sra}_mapped_sorted_mapq30_BlClean.bam",
        bam_bai="output/h3k27me3_h3k4/bam/{sra}_mapped_sorted_mapq30_BlClean.bam.bai"
    output:
        "output/bigwig_RPGC_MAPQ_BL/h3k27me3_h3k4_MAPQ_BL_RPGC_{sra}.bw"
    shell:
        """
        bamCoverage --bam {input.bam} -o {output} --binSize 50 --normalizeUsing RPGC --skipNonCoveredRegions --effectiveGenomeSize 2685859998
        """

rule multiBigwigSummaryBL:
    input:
        expand("output/bigwig_RPGC_MAPQ_BL/h3k27me3_h3k4_MAPQ_BL_RPGC_{sra}.bw", sra=SRA)
    output:
        tabFile="output/bigwig_RPGC_MAPQ_BL/scores_per_bin.tab",
        npzFile="output/bigwig_RPGC_MAPQ_BL/scores_per_bin.npz"
    shell:
        """
        multiBigwigSummary bins -b {input} -out {output.npzFile} --outRawCounts {output.tabFile}
        """

rule PlotMultiBigwigSummaryBL:
    input:
        tabFile="output/bigwig_RPGC_MAPQ_BL/scores_per_bin.tab",
        npzFile="output/bigwig_RPGC_MAPQ_BL/scores_per_bin.npz"
    output:
        tabCorr="output/bigwig_RPGC_MAPQ_BL/SpearmanCorr_readCounts.tab",
        plot="output/bigwig_RPGC_MAPQ_BL/heatmap_SpearmanCorr_readCounts.png"
    shell:
        """
        plotCorrelation -in {input.npzFile} --corMethod spearman --skipZeros --plotTitle "Correlation_of_read_counts" --whatToPlot heatmap --colorMap RdYlBu --plotNumbers -o {output.plot}  --outFileCorMatrix {output.tabCorr}
        """

## finish alignment and file conversion ##

## Peak Calling and TMM Normalization ##

rule macs2:
    input:
        "output/h3k27me3_h3k4/bam/{sra}_mapped_sorted_mapq30_BlClean.bam"
    output:
        broadPeak="output/h3k27me3_h3k4/macs2/h3k27me3_h3k4_macs_{sra}_peaks.broadPeak",
        gappedPeak="output/h3k27me3_h3k4/macs2/h3k27me3_h3k4_macs_{sra}_peaks.gappedPeak",
        broadXls="output/h3k27me3_h3k4/macs2/h3k27me3_h3k4_macs_{sra}_peaks.xls"
    shell:
        """
        macs2 callpeak -t {input} -f BAMPE -q 0.01 --broad --keep-dup 1 --outdir output/h3k27me3_h3k4/macs2 -n h3k27me3_h3k4_macs_{wildcards.sra}
        """

rule BED_to_SAF:
    input:
        "data/reference/hg38_chrom_sizes_binned_500k.bed"
    output:
        "output/h3k27me3_h3k4/saf/hg38_500K_bin.saf"
    shell:
        """
        awk 'OFS="\t" {{print $1"."$2+1"."$3, $1, $2+1, $3, "."}}' {input} > {output}
        """

rule rawCount:
    input:
        bam_all=expand("output/h3k27me3_h3k4/bam/{sra}_mapped_sorted_mapq30_BlClean.bam",sra=SRA),
        saf="output/h3k27me3_h3k4/saf/hg38_500K_bin.saf"
    output:
        "output/h3k27me3_h3k4/rawCount/stdout64.txt"
    threads:
        5
    shell:
        """
        featureCounts -p -O -T {threads} -F SAF -a {input.saf} -o {output} {input.bam_all}
        """

rule countMatrix:
    input:
        "output/h3k27me3_h3k4/rawCount/stdout64.txt"
    output:
        "output/h3k27me3_h3k4/rawCount/raw_count_mapq30_BlClean_temp_64.txt"
    shell:
        "cut -f1,7- < {input} | sed 1d > {output}"

rule removeDuplicateColum64:
    input:
        "output/h3k27me3_h3k4/rawCount/raw_count_mapq30_BlClean_temp_64.txt"
    output:
        "output/h3k27me3_h3k4/rawCount/raw_count_mapq30_BlClean_64.txt"
    script:
        "scripts/remove_dup_column.R"


rule TMM:
    input:
        "output/h3k27me3_h3k4/rawCount/raw_count_mapq30_BlClean_64.txt"
    output:
        "output/h3k27me3_h3k4/size_factors/SF.txt"
    script:
        "scripts/TMM.R"

rule splitSF:
    input:
        "output/h3k27me3_h3k4/size_factors/SF.txt"
    output:
        "output/h3k27me3_h3k4/size_factors/SF_{sra}_mapped_sorted_mapq30_BlClean.bam.txt"
    shell:
        """
        awk -F '\t' '{{print $2 > "output/h3k27me3_h3k4/size_factors/SF_"$1".txt"}}' {input}
        """

rule bamCoverage_MAPQ_BL_TMM:
    input:
        bam="output/h3k27me3_h3k4/bam/{sra}_mapped_sorted_mapq30_BlClean.bam",
        bam_bai="output/h3k27me3_h3k4/bam/{sra}_mapped_sorted_mapq30_BlClean.bam.bai",
        sf="output/h3k27me3_h3k4/size_factors/SF_{sra}_mapped_sorted_mapq30_BlClean.bam.txt"
    output:
        "output/bigwig_TMM_MAPQ_BL/h3k27me3_h3k4_MAPQ_BL_TMM_{sra}.bw"
    threads:
        6
    shell:
        """
        bamCoverage --bam {input.bam} -o {output} --binSize 50 --scaleFactor $(cat {input.sf}) --numberOfProcessors {threads}
        """

rule multiBigwigSummary:
    input:
        expand("output/bigwig_TMM_MAPQ_BL/h3k27me3_h3k4_MAPQ_BL_TMM_{sra}.bw", sra=SRA)
    output:
        tabFile="output/bigwig_TMM_MAPQ_BL/scores_per_bin.tab",
        npzFile="output/bigwig_TMM_MAPQ_BL/scores_per_bin.npz"
    shell:
        """
        multiBigwigSummary bins -b {input} -out {output.npzFile} --outRawCounts {output.tabFile}
        """

rule PlotMultiBigwigSummary:
    input:
        tabFile="output/bigwig_TMM_MAPQ_BL/scores_per_bin.tab",
        npzFile="output/bigwig_TMM_MAPQ_BL/scores_per_bin.npz"
    output:
        tabCorr="output/bigwig_TMM_MAPQ_BL/SpearmanCorr_readCounts.tab",
        plot="output/bigwig_TMM_MAPQ_BL/heatmap_SpearmanCorr_readCounts.png"
    shell:
        """
        plotCorrelation -in {input.npzFile} --corMethod spearman --skipZeros --plotTitle "Correlation_of_read_counts" --whatToPlot heatmap --colorMap RdYlBu --plotNumbers -o {output.plot}  --outFileCorMatrix {output.tabCorr}
        """