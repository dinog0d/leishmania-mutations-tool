#!/bin/bash
set -euo pipefail

# ============================================================
# Leishmania Variant Calling Pipeline (versión flexible)
# ============================================================

# === 1. Variables globales ===
# Aceptar job ID como parámetro, usar sample_leishmania como default
JOB_ID=${1:-sample_leishmania}
SAMPLE="sample_leishmania"  # Nombre interno de la muestra
THREADS=8

echo "🔬 Iniciando pipeline para Job ID: $JOB_ID"

# Detectar si estamos en Docker (volúmenes montados en /reference, /data, /results)
if [ -d "/reference" ] && [ -d "/data" ] && [ -d "/results" ]; then
    echo "🔎 Ejecutando en entorno Docker..."
    BASE_DIR="/"
    REFERENCE_DIR="/reference"
    DATA_DIR="/data"
    # Usar job ID para crear carpeta de resultados específica
    RESULTS_DIR="/results/$JOB_ID"
else
    echo "💻 Ejecutando en entorno local..."
    BASE_DIR="$(pwd)"
    REFERENCE_DIR="$BASE_DIR/reference"
    DATA_DIR="$BASE_DIR/data"
    RESULTS_DIR="$BASE_DIR/results/$JOB_ID"
fi

echo "📁 Directorio de datos: $DATA_DIR"
echo "📁 Directorio de resultados: $RESULTS_DIR"

# Verificar que existan los archivos de entrada - buscar primero con sample_leishmania, luego con job ID
INPUT_DIR=""
R1_FILE=""
R2_FILE=""

# Buscar archivos en directorio con job ID
if [ -d "$DATA_DIR/$JOB_ID" ]; then
    INPUT_DIR="$DATA_DIR/$JOB_ID"
    R1_FILE="$INPUT_DIR/${SAMPLE}_R1.fastq.gz"
    R2_FILE="$INPUT_DIR/${SAMPLE}_R2.fastq.gz"
    
    if [ ! -f "$R1_FILE" ] || [ ! -f "$R2_FILE" ]; then
        echo "🔄 Archivos no encontrados con sample_leishmania, buscando con sample01..."
        R1_FILE="$INPUT_DIR/sample01_R1.fastq.gz"
        R2_FILE="$INPUT_DIR/sample01_R2.fastq.gz"
    fi
# Buscar archivos en directorio sample_leishmania
elif [ -d "$DATA_DIR/sample_leishmania" ]; then
    INPUT_DIR="$DATA_DIR/sample_leishmania"
    R1_FILE="$INPUT_DIR/${SAMPLE}_R1.fastq.gz"
    R2_FILE="$INPUT_DIR/${SAMPLE}_R2.fastq.gz"
# Buscar archivos en directorio sample01 (compatibilidad)
elif [ -d "$DATA_DIR/sample01" ]; then
    INPUT_DIR="$DATA_DIR/sample01"
    R1_FILE="$INPUT_DIR/sample01_R1.fastq.gz"
    R2_FILE="$INPUT_DIR/sample01_R2.fastq.gz"
else
    echo "❌ Error: No se encontró directorio de entrada para el job $JOB_ID"
    echo "   Buscados en:"
    echo "   - $DATA_DIR/$JOB_ID"
    echo "   - $DATA_DIR/sample_leishmania"
    echo "   - $DATA_DIR/sample01"
    exit 1
fi

if [ ! -f "$R1_FILE" ] || [ ! -f "$R2_FILE" ]; then
    echo "❌ Error: No se encontraron los archivos FASTQ requeridos"
    echo "   R1: $R1_FILE"
    echo "   R2: $R2_FILE"
    exit 1
fi

echo "✅ Archivos de entrada encontrados:"
echo "   R1: $R1_FILE"
echo "   R2: $R2_FILE"

FASTP_DIR="$RESULTS_DIR/fastp"
ALIGN_DIR="$RESULTS_DIR/alignment"
VAR_DIR="$RESULTS_DIR/variants"
REPORT_DIR="$RESULTS_DIR/reports"
LOG_DIR="$RESULTS_DIR/logs"

REF_FASTA="$REFERENCE_DIR/Leishmania_braziliensis_MHOM_BR_75_M2904.fna"
REF_GFF="$REFERENCE_DIR/Leishmania_braziliensis_MHOM_BR_75_M2904.gff"
BED_FILE="$REFERENCE_DIR/resistance_genes.bed"

# === 2. Descargar y preparar referencias ===
echo ""
echo "🧬 [PASO 0/7] Verificando y descargando referencias desde NCBI..."
mkdir -p "$REFERENCE_DIR"

if [ ! -f "$REF_FASTA" ]; then
    echo "📥 Descargando archivo FASTA..."
    wget -O "$REFERENCE_DIR/GCF_000002845.2_ASM284v2_genomic.fna.gz" "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/845/GCF_000002845.2_ASM284v2/GCF_000002845.2_ASM284v2_genomic.fna.gz"
    gunzip "$REFERENCE_DIR/GCF_000002845.2_ASM284v2_genomic.fna.gz"
    mv "$REFERENCE_DIR/GCF_000002845.2_ASM284v2_genomic.fna" "$REF_FASTA"
    echo "✅ Archivo FASTA descargado"
else
    echo "✅ Archivo FASTA ya disponible"
fi

if [ ! -f "$REF_GFF" ]; then
    echo "📥 Descargando archivo GFF..."
    wget -O "$REFERENCE_DIR/GCF_000002845.2_ASM284v2_genomic.gff.gz" "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/845/GCF_000002845.2_ASM284v2/GCF_000002845.2_ASM284v2_genomic.gff.gz"
    gunzip "$REFERENCE_DIR/GCF_000002845.2_ASM284v2_genomic.gff.gz"
    mv "$REFERENCE_DIR/GCF_000002845.2_ASM284v2_genomic.gff" "$REF_GFF"
    echo "✅ Archivo GFF descargado"
else
    echo "✅ Archivo GFF ya disponible"
fi

if [ ! -f "$REF_FASTA.fai" ]; then
    echo "🔍 Indexando referencia con samtools..."
    samtools faidx "$REF_FASTA"
    echo "✅ Índice samtools creado"
else
    echo "✅ Índice samtools ya disponible"
fi

if [ ! -f "$REF_FASTA.bwt" ]; then
    echo "🔍 Indexando referencia con BWA..."
    bwa index "$REF_FASTA"
    echo "✅ Índice BWA creado"
else
    echo "✅ Índice BWA ya disponible"
fi

if [ ! -f "$REFERENCE_DIR/Leishmania_braziliensis_MHOM_BR_75_M2904.dict" ]; then
    echo "📖 Creando diccionario de secuencias con Picard..."
    picard CreateSequenceDictionary \
        -R "$REF_FASTA" \
        -O "$REFERENCE_DIR/Leishmania_braziliensis_MHOM_BR_75_M2904.dict"
    echo "✅ Diccionario de secuencias creado"
else
    echo "✅ Diccionario de secuencias ya disponible"
fi

# === 3. Crear directorios ===
echo ""
echo "📁 Creando directorios de trabajo..."
mkdir -p "$FASTP_DIR" "$ALIGN_DIR" "$VAR_DIR" "$REPORT_DIR" "$LOG_DIR"
echo "✅ Directorios creados en: $RESULTS_DIR"

# === 4. Filtrado de lecturas (fastp) ===
echo ""
echo "🧹 [PASO 1/7] Filtrando lecturas con fastp..."
echo "   📊 Analizando calidad y eliminando adaptadores..."
fastp \
  -i "$R1_FILE" \
  -I "$R2_FILE" \
  -o "$FASTP_DIR/${SAMPLE}_R1_filtered.fastq.gz" \
  -O "$FASTP_DIR/${SAMPLE}_R2_filtered.fastq.gz" \
  -h "$FASTP_DIR/fastp_report.html" \
  -j "$FASTP_DIR/fastp_report.json" \
  -w $THREADS
echo "✅ Filtrado de lecturas completado"

# === 5. Alineamiento con BWA y Samtools ===
echo ""
echo "🎯 [PASO 2/7] Alineando lecturas con BWA..."
echo "   🧬 Mapeando lecturas contra genoma de referencia..."
bwa mem -t $THREADS "$REF_FASTA" \
  "$FASTP_DIR/${SAMPLE}_R1_filtered.fastq.gz" \
  "$FASTP_DIR/${SAMPLE}_R2_filtered.fastq.gz" \
  > "$ALIGN_DIR/${SAMPLE}.sam" 2> "$LOG_DIR/bwa_alignment.log"
echo "✅ Alineamiento BWA completado"

echo ""
echo "🔧 [PASO 3/7] Procesando archivos BAM..."
echo "   📦 Convirtiendo SAM a BAM..."
samtools view -bS "$ALIGN_DIR/${SAMPLE}.sam" > "$ALIGN_DIR/${SAMPLE}.bam"
echo "   🗃️ Ordenando por nombre..."
samtools sort -n -@ $THREADS -o "$ALIGN_DIR/${SAMPLE}_name_sorted.bam" "$ALIGN_DIR/${SAMPLE}.bam"
echo "   🔧 Corrigiendo información de pares..."
samtools fixmate -m -@ $THREADS "$ALIGN_DIR/${SAMPLE}_name_sorted.bam" "$ALIGN_DIR/${SAMPLE}_fixmate.bam"
echo "   🗃️ Ordenando por coordenadas..."
samtools sort -@ $THREADS -o "$ALIGN_DIR/${SAMPLE}_fixmate_sorted.bam" "$ALIGN_DIR/${SAMPLE}_fixmate.bam"
echo "   🏷️ Marcando duplicados..."
samtools markdup -@ $THREADS -r "$ALIGN_DIR/${SAMPLE}_fixmate_sorted.bam" "$ALIGN_DIR/${SAMPLE}_marked.bam"
echo "   🔍 Creando índice..."
samtools index "$ALIGN_DIR/${SAMPLE}_marked.bam"
echo "✅ Procesamiento de archivos BAM completado"

# === 6. Llamado y filtrado de variantes ===
echo ""
echo "🧪 [PASO 4/7] Llamando variantes con bcftools..."
echo "   🔬 Generando pileup y llamando variantes..."
bcftools mpileup -Ou -f "$REF_FASTA" "$ALIGN_DIR/${SAMPLE}_marked.bam" \
  --threads $THREADS --max-depth 1000 --annotate AD,DP,SP,INFO/AD \
  | bcftools call -mv -Oz -o "$VAR_DIR/${SAMPLE}_raw.vcf.gz"

bcftools index "$VAR_DIR/${SAMPLE}_raw.vcf.gz"
echo "✅ Llamado de variantes completado"

echo ""
echo "🔍 [PASO 5/7] Filtrando variantes..."
echo "   📏 Aplicando filtros de calidad (DP≥10, QUAL≥30)..."
bcftools filter -s LOWQUAL -e 'INFO/DP<10 || QUAL<30' \
  -Oz -o "$VAR_DIR/${SAMPLE}_filtered.vcf.gz" "$VAR_DIR/${SAMPLE}_raw.vcf.gz"
bcftools index "$VAR_DIR/${SAMPLE}_filtered.vcf.gz"
echo "✅ Filtrado de variantes completado"

# === 7. Anotación ===
echo ""
echo "📝 [PASO 6/7] Anotando variantes..."
echo "   🧬 Prediciendo efectos funcionales de las variantes..."
bcftools csq -f "$REF_FASTA" -g "$REF_GFF" --local-csq \
  -Ov -o "$VAR_DIR/${SAMPLE}_annotated.vcf" "$VAR_DIR/${SAMPLE}_filtered.vcf.gz"
bgzip -c "$VAR_DIR/${SAMPLE}_annotated.vcf" > "$VAR_DIR/${SAMPLE}_annotated.vcf.gz"
bcftools index "$VAR_DIR/${SAMPLE}_annotated.vcf.gz"
echo "✅ Anotación de variantes completada"

# === 8. Extracción de genes de resistencia ===
echo ""
echo "💊 [PASO 7/7] Extrayendo variantes asociadas a resistencia..."
echo "   🎯 Filtrando variantes en genes de resistencia a fármacos..."
bcftools view -R "$BED_FILE" -Ov -o "$VAR_DIR/${SAMPLE}_resistance.vcf" "$VAR_DIR/${SAMPLE}_annotated.vcf.gz"

echo "   📊 Generando reporte de resistencia..."
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%INFO/BCSQ\n' \
  "$VAR_DIR/${SAMPLE}_resistance.vcf" > "$REPORT_DIR/${SAMPLE}_resistance_report.tsv"

# Crear reporte resumen
echo "   📋 Creando reporte resumen..."
{
    echo "=== REPORTE DE ANÁLISIS DE RESISTENCIA ==="
    echo "Job ID: $JOB_ID"
    echo "Fecha: $(date)"
    echo "Muestra: $SAMPLE"
    echo ""
    echo "=== ESTADÍSTICAS ==="
    echo "Variantes totales: $(bcftools view -H "$VAR_DIR/${SAMPLE}_filtered.vcf.gz" | wc -l)"
    echo "Variantes en genes de resistencia: $(bcftools view -H "$VAR_DIR/${SAMPLE}_resistance.vcf" | wc -l)"
    echo ""
    echo "=== VARIANTES DE RESISTENCIA ==="
} > "$REPORT_DIR/${SAMPLE}_summary.txt"

echo ""
echo "🎉 ¡Pipeline completado exitosamente para Job $JOB_ID!"
echo "📁 Resultados guardados en: $RESULTS_DIR"
echo "📊 Reporte principal: $REPORT_DIR/${SAMPLE}_resistance_report.tsv"
echo "📋 Resumen: $REPORT_DIR/${SAMPLE}_summary.txt"
