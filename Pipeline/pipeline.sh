#!/bin/bash
set -euo pipefail

# ============================================================
# Leishmania Variant Calling Pipeline (versión flexible)
# ============================================================

# === 1. Variables globales ===
SAMPLE="sample01"
THREADS=8

# Detectar si estamos en Docker (volúmenes montados en /reference, /data, /results)
if [ -d "/reference" ] && [ -d "/data" ] && [ -d "/results" ]; then
    echo "🔎 Ejecutando en entorno Docker..."
    BASE_DIR="/"
    REFERENCE_DIR="/reference"
    DATA_DIR="/data"
    RESULTS_DIR="/results/$SAMPLE"
else
    echo "💻 Ejecutando en entorno local..."
    BASE_DIR="$(pwd)"
    REFERENCE_DIR="$BASE_DIR/reference"
    DATA_DIR="$BASE_DIR/data"
    RESULTS_DIR="$BASE_DIR/results/$SAMPLE"
fi

FASTP_DIR="$RESULTS_DIR/fastp"
ALIGN_DIR="$RESULTS_DIR/alignment"
VAR_DIR="$RESULTS_DIR/variants"
REPORT_DIR="$RESULTS_DIR/reports"
LOG_DIR="$RESULTS_DIR/logs"

REF_FASTA="$REFERENCE_DIR/Leishmania_braziliensis_MHOM_BR_75_M2904.fna"
REF_GFF="$REFERENCE_DIR/Leishmania_braziliensis_MHOM_BR_75_M2904.gff"
BED_FILE="$REFERENCE_DIR/resistance_genes.bed"

# === 2. Descargar y preparar referencias ===
echo "[0/7] Verificando y descargando referencias desde NCBI..."
mkdir -p "$REFERENCE_DIR"

if [ ! -f "$REF_FASTA" ]; then
    echo "Descargando archivo FASTA..."
    wget -O "$REFERENCE_DIR/GCF_000002845.2_ASM284v2_genomic.fna.gz" "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/845/GCF_000002845.2_ASM284v2/GCF_000002845.2_ASM284v2_genomic.fna.gz"
    gunzip "$REFERENCE_DIR/GCF_000002845.2_ASM284v2_genomic.fna.gz"
    mv "$REFERENCE_DIR/GCF_000002845.2_ASM284v2_genomic.fna" "$REF_FASTA"
fi

if [ ! -f "$REF_GFF" ]; then
    echo "Descargando archivo GFF..."
    wget -O "$REFERENCE_DIR/GCF_000002845.2_ASM284v2_genomic.gff.gz" "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/845/GCF_000002845.2_ASM284v2/GCF_000002845.2_ASM284v2_genomic.gff.gz"
    gunzip "$REFERENCE_DIR/GCF_000002845.2_ASM284v2_genomic.gff.gz"
    mv "$REFERENCE_DIR/GCF_000002845.2_ASM284v2_genomic.gff" "$REF_GFF"
fi

if [ ! -f "$REF_FASTA.fai" ]; then
    echo "Indexando referencia con samtools..."
    samtools faidx "$REF_FASTA"
fi

if [ ! -f "$REF_FASTA.bwt" ]; then
    echo "Indexando referencia con BWA..."
    bwa index "$REF_FASTA"
fi

if [ ! -f "$REFERENCE_DIR/Leishmania_braziliensis_MHOM_BR_75_M2904.dict" ]; then
    echo "Creando diccionario de secuencias con Picard..."
    picard CreateSequenceDictionary \
        -R "$REF_FASTA" \
        -O "$REFERENCE_DIR/Leishmania_braziliensis_MHOM_BR_75_M2904.dict"
fi

# === 3. Crear directorios ===
mkdir -p "$FASTP_DIR" "$ALIGN_DIR" "$VAR_DIR" "$REPORT_DIR" "$LOG_DIR"

# === 4. Filtrado de lecturas (fastp) ===
echo "[1/7] Filtrando lecturas con fastp..."
fastp \
  -i "$DATA_DIR/$SAMPLE/${SAMPLE}_R1.fastq.gz" \
  -I "$DATA_DIR/$SAMPLE/${SAMPLE}_R2.fastq.gz" \
  -o "$FASTP_DIR/${SAMPLE}_R1_filtered.fastq.gz" \
  -O "$FASTP_DIR/${SAMPLE}_R2_filtered.fastq.gz" \
  -h "$FASTP_DIR/fastp_report.html" \
  -j "$FASTP_DIR/fastp_report.json" \
  -w $THREADS

# === 5. Alineamiento con BWA y Samtools ===
echo "[2/7] Alineando lecturas con BWA..."
bwa mem -t $THREADS "$REF_FASTA" \
  "$FASTP_DIR/${SAMPLE}_R1_filtered.fastq.gz" \
  "$FASTP_DIR/${SAMPLE}_R2_filtered.fastq.gz" \
  > "$ALIGN_DIR/${SAMPLE}.sam" 2> "$LOG_DIR/bwa_alignment.log"

echo "[3/7] Procesando archivos BAM..."
samtools view -bS "$ALIGN_DIR/${SAMPLE}.sam" > "$ALIGN_DIR/${SAMPLE}.bam"
samtools sort -n -@ $THREADS -o "$ALIGN_DIR/${SAMPLE}_name_sorted.bam" "$ALIGN_DIR/${SAMPLE}.bam"
samtools fixmate -m -@ $THREADS "$ALIGN_DIR/${SAMPLE}_name_sorted.bam" "$ALIGN_DIR/${SAMPLE}_fixmate.bam"
samtools sort -@ $THREADS -o "$ALIGN_DIR/${SAMPLE}_fixmate_sorted.bam" "$ALIGN_DIR/${SAMPLE}_fixmate.bam"
samtools markdup -@ $THREADS -r "$ALIGN_DIR/${SAMPLE}_fixmate_sorted.bam" "$ALIGN_DIR/${SAMPLE}_marked.bam"
samtools index "$ALIGN_DIR/${SAMPLE}_marked.bam"

# === 6. Llamado y filtrado de variantes ===
echo "[4/7] Llamando variantes con bcftools..."
bcftools mpileup -Ou -f "$REF_FASTA" "$ALIGN_DIR/${SAMPLE}_marked.bam" \
  --threads $THREADS --max-depth 1000 --annotate AD,DP,SP,INFO/AD \
  | bcftools call -mv -Oz -o "$VAR_DIR/${SAMPLE}_raw.vcf.gz"

bcftools index "$VAR_DIR/${SAMPLE}_raw.vcf.gz"

echo "[5/7] Filtrando variantes..."
bcftools filter -s LOWQUAL -e 'INFO/DP<10 || QUAL<30' \
  -Oz -o "$VAR_DIR/${SAMPLE}_filtered.vcf.gz" "$VAR_DIR/${SAMPLE}_raw.vcf.gz"
bcftools index "$VAR_DIR/${SAMPLE}_filtered.vcf.gz"

# === 7. Anotación ===
echo "[6/7] Anotando variantes..."
bcftools csq -f "$REF_FASTA" -g "$REF_GFF" --local-csq \
  -Ov -o "$VAR_DIR/${SAMPLE}_annotated.vcf" "$VAR_DIR/${SAMPLE}_filtered.vcf.gz"
bgzip -c "$VAR_DIR/${SAMPLE}_annotated.vcf" > "$VAR_DIR/${SAMPLE}_annotated.vcf.gz"
bcftools index "$VAR_DIR/${SAMPLE}_annotated.vcf.gz"

# === 8. Extracción de genes de resistencia ===
echo "[7/7] Extrayendo variantes asociadas a resistencia..."
bcftools view -R "$BED_FILE" -Ov -o "$VAR_DIR/${SAMPLE}_resistance.vcf" "$VAR_DIR/${SAMPLE}_annotated.vcf.gz"

bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%INFO/BCSQ\n' \
  "$VAR_DIR/${SAMPLE}_resistance.vcf" > "$REPORT_DIR/${SAMPLE}_resistance_report.tsv"

echo "✅ Pipeline completado exitosamente."
