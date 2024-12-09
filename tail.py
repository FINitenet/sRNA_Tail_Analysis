#!/usr/bin/python3
# -*- coding:utf-8 -*-
# @FileName  :trimming.py
# @Time      :2024/12/05 23:36:51
# @Author    :Yuchen@rlab
# @Description: This script is used to trim the ago trimming data

import glob
import gzip
import os
import subprocess
import shutil
from Bio import SeqIO
import pysam
import operator
from tqdm import tqdm
from collections import defaultdict

def modify_fastq_header(fastq_file):
    outdir = os.path.join(OUTDIR, "2_format_fasta")
    os.makedirs(outdir, exist_ok=True)

    output_file = f"{outdir}/{SAMPLENAME}.fasta"
    output_records = []
    with gzip.open(fastq_file, "rt") as f:
        for record in SeqIO.parse(f, "fastq"):
            record.id = f"{SAMPLENAME}-{len(output_records)}-{record.seq}"
            record.description = ""
            output_records.append(record)
    with open(output_file, "w") as f:
        SeqIO.write(output_records, f, "fasta")

def remapping(reference):
    outdir = os.path.join(OUTDIR, "3_remapping")
    os.makedirs(outdir, exist_ok=True)
    
    fasta_file = f"{OUTDIR}/2_format_fasta/{SAMPLENAME}.fasta"

    log_file = os.path.join(outdir, "remapping.log")
    log_file = open(log_file, "w")

    shutil.copy(fasta_file, f"{outdir}/unmapped-0.fasta")
    for i in range(1, 6):
        k = i - 1
        output_records = []
        print(f"Trimming {i} bp for each unmapped reads")
        with open(f"{outdir}/unmapped-{k}.fasta", "r") as in_file:
            for record in SeqIO.parse(in_file, "fasta"):
                record.seq = record.seq[:-1]
                if len(record.seq) >= 12:
                    record.id = record.id.split("--")[0] + f"--{i}"
                    record.description = ""
                    output_records.append(record)

        with open(f"{outdir}/trimmed-{i}.fasta", "w") as out_file:
            SeqIO.write(output_records, out_file, "fasta")
        
        print(f"Mapping trimmed {i} bp reads to genome and unmapped reads are saved to {outdir}/unmapped-{i}.fasta")
        remapping_cmd = "bowtie -p 24 -v 0 -S -a --no-unal -x {} -f {}/trimmed-{}.fasta --un {}/unmapped-{}.fasta {}/mapped-{}.sam".format(reference, outdir, i, outdir, i, outdir, i)
        subprocess.run(remapping_cmd, check=True, shell=True, stderr=log_file, stdout=log_file)
    log_file.close()

def annotation_counting(input_path, output_path, sample_name):
    mapping_files = glob.glob(f"{input_path}/*.sam")
    os.makedirs(output_path, exist_ok=True)
    for mapping_file in mapping_files:
        middle_name = mapping_file.split("/")[-1].split(".")[0]
        summary_file = f"{output_path}/{sample_name}-{middle_name}-counts.txt"
        featurecounts_all_type = [
                    "featureCounts",
                    "-a", GFF_FILE,
                    "-o", summary_file,
                    mapping_file,
                    "-t", "gene,ncRNA_gene",
                    "-g", "biotype",
                    "-T", "48",
                    "-O",
                    "-M",
                    "-R", "BAM",
                    "--largestOverlap",
                ]
        subprocess.run(featurecounts_all_type, check=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)


def feature_choose(read, priority_dict):
    stats = read.get_tag('XS')
    if stats == 'Assigned':
        read_features = read.get_tag('XT').split(',')
        if len(read_features) == 1:
            feature = read_features[0]
        else:
            feature = min(read_features, key=lambda x: priority_dict[x])
    elif stats == 'Unassigned_NoFeatures':
        feature = 'otherRNA'
    return feature

def calc_feature_count_by_featureCounts_tags_priority(input_file, priority_dict):
    file = input_file.replace(".bam", ".sorted.bam")
    pysam.sort("-o", file, input_file, threads=12)
    pysam.index(file, threads=12)
    with pysam.AlignmentFile(file, 'r', threads=12) as fh:
        feature_dict = {}
        for read in fh.fetch():
            sequence = str(read.query_name.split("--")[0].split("-")[-1])
            trim_num = int(read.query_name.split("--")[-1])
            trim_base = sequence[-trim_num:]

            if read.query_name not in feature_dict:
                feature = feature_choose(read, priority_dict)
                feature_dict[read.query_name] = f"{trim_base}.{feature}"
            else:
                current_feature = feature_dict[read.query_name].split('.')[1]
                new_feature = feature_choose(read, priority_dict)
                if priority_dict[current_feature] > priority_dict[new_feature]:
                    feature_dict[read.query_name] = f"{trim_base}.{new_feature}"

    type_dict = defaultdict(int)
    for key, value in feature_dict.items():
        type_dict[value] += 1
    type_dict = dict(sorted(type_dict.items(), key=operator.itemgetter(0)))
    return type_dict


def write_to_output_file(output_file, feature_dict):
    header = 'trim_base\tfeature\tCount\tPercentage'
    with open(output_file, 'w') as output_f:
        output_f.write(header + '\n')
        value_sum = sum(feature_dict.values())
        for key, value in feature_dict.items():
            trim_base, rnatype = key.split('.')
            line = '{}\t{}\t{:.3f}\t{:.3f}'.format(trim_base, rnatype, value, value / value_sum * 100)
            output_f.write(line + '\n')

if __name__ == "__main__":
    
    INPUTFILE = "/bios-store1/chenyc/test_ago_trimming/1_input/IP-CK-R1_unaligned.fastq.gz"
    SAMPLENAME = INPUTFILE.split("/")[-1].split(".")[0].replace("_unaligned", "")
    OUTDIR = os.path.join("/bios-store1/chenyc/test_ago_trimming", SAMPLENAME)
    REFERENCE = "/bios-store1/chenyc/Reference_Source/Arabidopsis_Reference/ath_chr_bowtie_index/Arabidopsis_thaliana.TAIR10.dna.toplevel"
    GFF_FILE = "/bios-store1/chenyc/scripts/TRMRNAseqTools/reference/Arabidopsis_thaliana.TAIR10.53.md.gff3"

    priority_dict = {
                "miRNA_primary_transcript":1,
                "phasi_tasiRNA":2,
                "phasi/tasiRNA":2,
                "hc-siRNA":3,
                "transposable_element":4,
                "lncRNA":5,
                "protein_coding":6,
                "snoRNA":7,
                "snRNA":8,
                "rRNA":9,
                "tRNA":10,
                "otherRNA":11
            }
    # Step 1: Format the fastq file to fasta file
    print("Start to trim the ago trimming data: ", INPUTFILE)
    print("Step 1: Format the fastq file to fasta file")
    modify_fastq_header(INPUTFILE)

    # Step 2: Remapping the fasta file to genome
    print("Step 2: Remapping the fasta file to genome")
    remapping(REFERENCE)

    # Step 3: Count the feature of each read
    print("Step 3: Annotate the feature of each read")
    annotation_counting(f"{OUTDIR}/3_remapping", f"{OUTDIR}/4_count", SAMPLENAME)

    # Step 4: Calculate the feature of each read
    print("Step 4: Calculate the feature of each read")
    BAMFILES = glob.glob(f"{OUTDIR}/4_count/mapped-*.sam.featureCounts.bam")
    for bamfile in tqdm(BAMFILES):
        feature_dict = calc_feature_count_by_featureCounts_tags_priority(bamfile, priority_dict)
        output_file = bamfile.replace(".bam", "_trimbase_feature_count.txt")
        write_to_output_file(output_file, feature_dict)
    print("All the steps are finished")