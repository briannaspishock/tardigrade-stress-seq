#!/usr/bin/env bash
set -euo pipefail

echo "=========================================================="
echo ">> INITIATING TARDIGRADE-STRESS-SEQ (CRYPTOBIOSIS PIPELINE)"
echo "=========================================================="

mkdir -p data ref trimmed qc_reports quant results figures

# 1. VERIFY SALMON INDEX
if [ ! -d "ref/salmon_index" ]; then
    echo "[*] Building Salmon k-mer index (k=31)..."
    salmon index -t ref/H_exemplaris_transcripts.fa -i ref/salmon_index -k 31
    echo "[+] Salmon indexing complete."
else
    echo "[+] Salmon index verified."
fi

# 2. SAMPLES TO PROCESS (Hydrated: 515/516 vs Desiccated: 517/518)
SAMPLES=("SRR3727515" "SRR3727516" "SRR3727517" "SRR3727518")

for SRA in "${SAMPLES[@]}"; do
    echo "----------------------------------------------------------"
    echo ">> PROCESSING SAMPLE: ${SRA}"
    echo "----------------------------------------------------------"

    # Fast 200k paired-end slice directly from NCBI (takes ~15-20s, no SSL timeouts)
    if [ ! -f "data/${SRA}_1.fastq.gz" ]; then
        echo "[*] Streaming 200,000 paired-end reads for ${SRA} from NCBI..."
        fastq-dump -X 200000 --split-files --gzip --outdir data/ "${SRA}"
    fi

    # Trim Paired-End Reads with fastp
    echo "[*] Running fastp adapter & quality filtering on ${SRA}..."
    fastp \
      -i "data/${SRA}_1.fastq.gz" \
      -I "data/${SRA}_2.fastq.gz" \
      -o "trimmed/${SRA}_1.clean.fastq.gz" \
      -O "trimmed/${SRA}_2.clean.fastq.gz" \
      --detect_adapter_for_pe \
      --length_required 30 \
      --thread 8 \
      --html "qc_reports/${SRA}_fastp.html" \
      --json "qc_reports/${SRA}_fastp.json"

    # Quantify with Salmon
    echo "[*] Quantifying expression with Salmon on ${SRA}..."
    salmon quant -i ref/salmon_index -l A \
      -1 "trimmed/${SRA}_1.clean.fastq.gz" \
      -2 "trimmed/${SRA}_2.clean.fastq.gz" \
      -p 8 --gcBias --seqBias \
      -o "quant/${SRA}"
done

# 3. MULTIQC SUMMARY
echo "----------------------------------------------------------"
echo "[*] Aggregating MultiQC report..."
multiqc qc_reports/ -o qc_reports/multiqc_summary

echo "=========================================================="
echo ">> UPSTREAM COMPLETE! RUN: python analyze_and_annotate.py"
echo "=========================================================="
