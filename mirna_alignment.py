#!/usr/bin/env python3
import os
import subprocess
import argparse
from datetime import datetime

def run_command(command, log_file):
    with open(log_file, 'a') as log:
        log.write(f"\n[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] Running: {command}\n")
        process = subprocess.run(command, shell=True, stdout=log, stderr=log)
        if process.returncode != 0:
            log.write(f"[ERROR] Command failed with return code {process.returncode}\n")
            raise subprocess.CalledProcessError(process.returncode, command)


def ensure_dir(path):
    os.makedirs(path, exist_ok=True)

def main():
    parser = argparse.ArgumentParser(description="miRNA pipeline (Steps 1–7): QC, fastp, alignment")
    parser.add_argument('--raw', required=True, help='Path to raw FASTQ folder')
    parser.add_argument('--output', required=True, help='Main output folder')
    parser.add_argument('--adapter', required=True, help='Path to adapter fasta file')
    parser.add_argument('--reference', required=True, help='Path to Bowtie index prefix')
    parser.add_argument('--threads', type=int, default=8, help='Number of threads')
    args = parser.parse_args()

    raw_folder = os.path.abspath(args.raw)
    output_folder = os.path.abspath(args.output)
    adapter_fasta = os.path.abspath(args.adapter)
    reference_index = os.path.abspath(args.reference)
    threads = args.threads

    # Create output subdirectories
    qc_raw = os.path.join(output_folder, "1_QC_raw")
    clean_folder = os.path.join(output_folder, "2_Clean_data_miRNA")
    qc_clean = os.path.join(output_folder, "3_QC_clean")
    align_folder = os.path.join(output_folder, "4_Alignment")
    collapse_folder = os.path.join(output_folder, "5_Collapsed")

    for folder in [qc_raw, clean_folder, qc_clean, align_folder, collapse_folder]:
        ensure_dir(folder)

    log_file = os.path.join(output_folder, "pipeline_step1_7.log")

    # Step 1: FastQC on raw files
    run_command(f"fastqc {raw_folder}/*.fastq.gz -o {qc_raw} -t {threads}", log_file)

    # Step 2: SeqKit stats on raw files
    run_command(f"seqkit stats -a {raw_folder}/*.fastq.gz -o {qc_raw}/seqkit_raw.tsv", log_file)

    # Step 3: fastp trimming
    for fq in sorted(os.listdir(raw_folder)):
        if fq.endswith(".fastq.gz"):
            sample = fq.replace(".fastq.gz", "")
            run_command(
                f"fastp -i {raw_folder}/{fq} -o {clean_folder}/{sample}.fq.gz "
                f"-h {clean_folder}/{sample}.html -j {clean_folder}/{sample}.json "
                f"-w {threads} -x --length_required 18 --max_len1 24 -c "
                f"--adapter_fasta {adapter_fasta}",
                log_file
            )

    # Step 4: FastQC on cleaned reads
    run_command(f"fastqc {clean_folder}/*.fq.gz -o {qc_clean} -t {threads}", log_file)

    # Step 5: SeqKit stats on cleaned reads
    run_command(f"seqkit stats -a {clean_folder}/*.fq.gz -o {qc_clean}/seqkit_clean.tsv", log_file)

     # Step 6: Unzip .fq.gz and perform Bowtie alignment
    for fq_gz in sorted(os.listdir(clean_folder)):
        if fq_gz.endswith(".fq.gz"):
            sample = fq_gz.replace(".fq.gz", "")
            fq_gz_path = os.path.join(clean_folder, fq_gz)
            fq_path = os.path.join(clean_folder, f"{sample}.fq")
            
            # Unzip fq.gz to fq (always fresh unzip)
            run_command(f"gunzip -c {fq_gz_path} > {fq_path}", log_file)

            # Align using uncompressed .fq file
            sam_path = os.path.join(align_folder, f"{sample}.sam")
            run_command(
                f"bowtie {reference_index} -q {fq_path} "
                f"-S {sam_path} --threads {threads}",
                log_file
            )


    # Step 7: SAM to BAM and index
    for sam_file in sorted(os.listdir(align_folder)):
        if sam_file.endswith(".sam"):
            base = sam_file.replace(".sam", "")
            sam_path = os.path.join(align_folder, sam_file)
            bam_path = os.path.join(align_folder, f"{base}_sorted.bam")
            run_command(f"samtools view -bS {sam_path} | samtools sort -o {bam_path}", log_file)
            run_command(f"samtools index {bam_path}", log_file)

    # Step 8: Alignment stats (TSV output)
    alignment_stats_file = os.path.join(align_folder, "primary_alignment_stats.tsv")

    with open(alignment_stats_file, 'w') as stat_out:
        stat_out.write("Sample\tTotal\tPrimaryMapped\tPercent\n")

        for bam_file in sorted(os.listdir(align_folder)):
            if bam_file.endswith("_sorted.bam"):
                sample = bam_file.replace("_sorted.bam", "")
                bam_path = os.path.join(align_folder, bam_file)

                # Run samtools flagstat
                stats_cmd = f"samtools flagstat {bam_path}"
                process = subprocess.run(stats_cmd, shell=True, capture_output=True, text=True)
                stats = process.stdout

                try:
                    total = [line for line in stats.splitlines() if "in total" in line][0].split()[0]
                    primary_line = [line for line in stats.splitlines() if "primary mapped" in line][0]
                    primary_mapped = primary_line.split()[0]
                    percent = primary_line.split("(")[1].split("%")[0]
                except IndexError:
                    total = primary_mapped = percent = "NA"
                    with open(log_file, 'a') as log:
                        log.write(f"[WARNING] Could not parse stats for {bam_path}\n{stats}\n")

                stat_out.write(f"{sample}\t{total}\t{primary_mapped}\t{percent}\n")



    # Step 9: Collapse FASTQ reads
    all_fastq_path = os.path.join(align_folder, "all_merged.fastq")
    run_command(f"cat {' '.join(fastq_list)} > {all_fastq_path}", log_file)

    collapsed_fa = os.path.join(collapse_folder, "collapsed.fa")
    run_command(f"fastx_collapser -Q33 -i {all_fastq_path} -o {collapsed_fa}", log_file)

    renamed_fa = os.path.join(collapse_folder, "new_lapse.fasta")
    run_command(f"sed 's/>/>seq_/g' {collapsed_fa} | sed 's/-/_x/g' > {renamed_fa}", log_file)

    final_filtered_fa = os.path.join(collapse_folder, "new_lapse.fa")
    run_command(f"seqkit seq -m 17 -M 30 {renamed_fa} > {final_filtered_fa}", log_file)


if __name__ == "__main__":
    main()

