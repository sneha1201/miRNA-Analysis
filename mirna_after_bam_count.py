#!/usr/bin/env python3
import os
import subprocess
import argparse
from datetime import datetime


def run_command(command, log_file):
    with open(log_file, 'a') as log:
        log.write(f"\n[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] Running: {command}\n")
        process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        for line in process.stdout:
            log.write(line)
            print(line, end='')  # Optional: also show live output on screen
        process.wait()
        if process.returncode != 0:
            raise subprocess.CalledProcessError(process.returncode, command)


def ensure_dir(path):
    os.makedirs(path, exist_ok=True)


def main():
    parser = argparse.ArgumentParser(description="miRNA pipeline block 2: Mapper, miRDeep2, GFF and counts")
    parser.add_argument('--collapsed', required=True, help='Path to collapsed output folder')
    parser.add_argument('--output', required=True, help='Main output folder')
    parser.add_argument('--reference', required=True, help='Path to Bowtie index (for mapper.pl)')
    parser.add_argument('--ref_fasta_cleaned', required=True, help='Cleaned and fixed reference fasta (used in miRDeep2)')
    parser.add_argument('--mature', required=True, help='Path to mature.fa file')
    parser.add_argument('--hairpin', required=True, help='Path to hairpin.fa file')
    parser.add_argument('--bam_folder', required=True, help='Folder containing *_sorted.bam files')
    parser.add_argument('--threads', type=int, default=8, help='Number of threads')
    args = parser.parse_args()

    collapsed_folder = os.path.abspath(args.collapsed)
    output_folder = os.path.abspath(args.output)
    reference_index = os.path.abspath(args.reference)
    cleaned_fasta = os.path.abspath(args.ref_fasta_cleaned)
    mature = os.path.abspath(args.mature)
    hairpin = os.path.abspath(args.hairpin)
    bam_folder = os.path.abspath(args.bam_folder)
    threads = args.threads

    # Output folders
    mapper_out = os.path.join(output_folder, "6_Mapper")
    mirdeep_out = os.path.join(output_folder, "7_miRDeep2")
    known_out = os.path.join(output_folder, "8_Known_miRNA")
    novel_out = os.path.join(output_folder, "9_Novel_miRNA")
    for folder in [mapper_out, mirdeep_out, known_out, novel_out]:
        ensure_dir(folder)

    log_file = os.path.join(output_folder, "pipeline_step8_12.log")

    # Step 8: mapper.pl
    collapsed_fa = os.path.join(collapsed_folder, "collapsed.fa")
    run_command(
        f"mapper.pl \"{collapsed_fa}\" -c -p \"{reference_index}\" -t \"{mapper_out}/reads_collapsed_vs_genome.arf\" -o 90",
        log_file
    )
    run_command(
        f"sed 's/-/_x/' \"{mapper_out}/reads_collapsed_vs_genome.arf\" > \"{mapper_out}/mapped.arf\"",
        log_file
    )

    # Step 9: miRDeep2 (known)
    new_lapse_fa = os.path.join(collapsed_folder, "new_lapse.fa")
    mapped_arf = os.path.join(mapper_out, "mapped.arf")
    run_command(
        f"cd \"{mirdeep_out}\" && miRDeep2.pl \"{new_lapse_fa}\" \"{cleaned_fasta}\" \"{mapped_arf}\" "
        f"\"{mature}\" none \"{hairpin}\"",
        log_file
    )

    # Find result BED
    result_bed = None
    for f in os.listdir(mirdeep_out):
        if f.startswith("result_") and f.endswith(".bed"):
            result_bed = os.path.join(mirdeep_out, f)
            break
    if not result_bed:
        raise FileNotFoundError("miRDeep2 result BED file (known) not found")

    # Step 10: known miRNA extraction and GFF
    known_bed = os.path.join(known_out, "known_mirna.bed")
    known_gff = os.path.join(known_out, "known_mirna.gff")
    run_command(
        f"awk '$4 ~ /^known:/ && $5 >= 4' \"{result_bed}\" > \"{known_bed}\"",
        log_file
    )
    run_command(
        f"""awk -F'\t' 'BEGIN{{OFS="\\t"}}{{print $1,"mirdeep","mirna_transcript",$2,$3,$5,$6,".","ID="$4}}' \"{known_bed}\" > \"{known_gff}\"""",
        log_file
    )

    # Step 11: featureCounts for known miRNAs
    run_command(
        f"featureCounts -T {threads} -t mirna_transcript -g ID -o \"{known_out}/counts_known.txt\" "
        f"-a \"{known_gff}\" \"{bam_folder}\"/*.bam",
        log_file
    )

    # Step 12: Find novel result BED
    
    novel_bed = os.path.join(novel_out, "novel_mirna.bed")
    novel_gff = os.path.join(novel_out, "novel_mirna.gff")
    run_command(
        f"awk '$4 ~ /^novel:/ && $5 >= 4' \"{result_bed}\" > \"{novel_bed}\"",
        log_file
    )
    run_command(
        f"""awk -F'\t' 'BEGIN{{OFS="\\t"}}{{print $1,"mirdeep","mirna_transcript",$2,$3,$5,$6,".","ID="$4}}' \"{novel_bed}\" > \"{novel_gff}\"""",
        log_file
    )

    # Step 11: featureCounts for novel miRNAs
    run_command(
        f"featureCounts -T {threads} -t mirna_transcript -g ID -o \"{novel_out}/counts_known.txt\" "
        f"-a \"{novel_gff}\" \"{bam_folder}\"/*.bam",
        log_file
    )
    


if __name__ == "__main__":
    main()
