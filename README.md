# Herramienta Bioinformática para la Detección de Mutaciones en Leishmania braziliensis

Este proyecto hace parte de la Propuesta de Investigación PIS-3:
"Desarrollo de una herramienta informática para la detección de mutaciones en genes asociados a resistencia farmacológica en cepas colombianas de Leishmania braziliensis".

## 🛠️ Tecnologías Utilizadas

El sistema está implementado con una arquitectura basada en **contenedores Docker**:

- **Frontend Web**: Flask (Python) con HTML, CSS y JavaScript
- **Pipeline Bioinformático**: Bash scripts con herramientas especializadas
- **Containerización**: Docker y Docker Compose
- **Procesamiento**: Conda environment con herramientas bioinformáticas

### Herramientas Bioinformáticas Integradas:
- **BWA**: Alineamiento de secuencias contra genoma de referencia
- **samtools**: Manipulación y procesamiento de archivos BAM/SAM
- **bcftools**: Llamado, filtrado y anotación de variantes
- **fastp**: Control de calidad y filtrado de lecturas FASTQ
- **Picard**: Herramientas para procesamiento de datos genómicos

## 🎯 Funcionalidades

- **Carga de archivos FASTQ**: Interfaz web para subir archivos paired-end (.fastq.gz)
- **Pipeline automatizado**: Procesamiento completo desde FASTQ hasta variantes anotadas
- **Seguimiento en tiempo real**: Logs en vivo del procesamiento
- **Sistema de trabajos**: Gestión de múltiples análisis simultáneos con IDs únicos
- **Reportes detallados**: Análisis de resistencia farmacológica con visualización web
- **Validación de archivos**: Verificación de formato y prevención de duplicados

## 📂 Estructura del Proyecto

```
leishmania-mutations-tool/
├── clean_all.sh
├── docker-compose.yml
├── README.md
├── Leishweb/
│   ├── app.py
│   ├── Dockerfile
│   ├── requirements.txt
│   └── templates/
│       └── index.html
├── Pipeline/
│   ├── Dockerfile
│   ├── environment.yml
│   └── pipeline.sh
├── data/
│   ├── jobs/
│   │   ├── index.json
│   │   ├── job_[timestamp]_[id].json
│   │   └── job_[timestamp]_[id]_logs.txt
│   ├── job_[timestamp]_[id]/
│   │   ├── sample_leishmania_R1.fastq.gz
│   │   └── sample_leishmania_R2.fastq.gz
│   └── sample_leishmania/
│       ├── sample01_R1.fastq.gz
│       └── sample01_R2.fastq.gz
├── reference/
│   ├── Leishmania_braziliensis_MHOM_BR_75_M2904.dict
│   ├── Leishmania_braziliensis_MHOM_BR_75_M2904.fna
│   ├── Leishmania_braziliensis_MHOM_BR_75_M2904.fna.amb
│   ├── Leishmania_braziliensis_MHOM_BR_75_M2904.fna.ann
│   ├── Leishmania_braziliensis_MHOM_BR_75_M2904.fna.bwt
│   ├── Leishmania_braziliensis_MHOM_BR_75_M2904.fna.fai
│   ├── Leishmania_braziliensis_MHOM_BR_75_M2904.fna.pac
│   ├── Leishmania_braziliensis_MHOM_BR_75_M2904.fna.sa
│   ├── Leishmania_braziliensis_MHOM_BR_75_M2904.gff
│   └── resistance_genes.bed
└── results/
    └── job_[timestamp]_[id]/
        ├── alignment/
        │   ├── sample_leishmania.sam
        │   ├── sample_leishmania.bam
        │   ├── sample_leishmania_name_sorted.bam
        │   ├── sample_leishmania_fixmate.bam
        │   ├── sample_leishmania_fixmate_sorted.bam
        │   ├── sample_leishmania_marked.bam
        │   └── sample_leishmania_marked.bam.bai
        ├── fastp/
        │   ├── fastp_report.html
        │   ├── fastp_report.json
        │   ├── sample_leishmania_R1_filtered.fastq.gz
        │   └── sample_leishmania_R2_filtered.fastq.gz
        ├── logs/
        │   └── bwa_alignment.log
        ├── reports/
        │   ├── sample_leishmania_resistance_report.tsv
        │   └── sample_leishmania_summary.txt
        └── variants/
            ├── sample_leishmania_raw.vcf.gz
            ├── sample_leishmania_raw.vcf.gz.csi
            ├── sample_leishmania_filtered.vcf.gz
            ├── sample_leishmania_filtered.vcf.gz.csi
            ├── sample_leishmania_annotated.vcf
            ├── sample_leishmania_annotated.vcf.gz
            ├── sample_leishmania_annotated.vcf.gz.csi
            └── sample_leishmania_resistance.vcf
```

## 🚀 Instrucciones de Uso

### 1. Inicializar el Sistema
```bash
# Construir y levantar contenedores
docker-compose up -d

# Verificar que los servicios estén corriendo
docker-compose ps
```

### 2. Acceder a la Interfaz Web
- Abrir navegador en: `http://localhost:5001`
- Cargar archivos FASTQ paired-end (.fastq.gz)
- El sistema generará automáticamente un Job ID único

### 3. Monitorear el Análisis
- Seguimiento en tiempo real a través de logs web
- Visualización del progreso paso a paso
- Notificación automática al completarse

### 4. Revisar Resultados
- Reporte de resistencia farmacológica en formato TSV
- Resumen estadístico del análisis
- Archivos VCF con variantes anotadas

## 🔬 Pipeline Bioinformático

El pipeline ejecuta los siguientes pasos automatizados:

1. **Control de Calidad** (fastp):
   - Eliminación de adaptadores
   - Filtrado por calidad de bases
   - Generación de reportes HTML

2. **Alineamiento** (BWA + samtools):
   - Mapeo contra genoma de L. braziliensis
   - Procesamiento y ordenamiento de BAM
   - Marcado de duplicados

3. **Llamado de Variantes** (bcftools):
   - Identificación de SNPs e indels
   - Filtrado por calidad (DP≥10, QUAL≥30)
   - Indexación de archivos VCF

4. **Anotación Funcional** (bcftools csq):
   - Predicción de efectos en proteínas
   - Clasificación de mutaciones
   - Anotación con archivo GFF

5. **Análisis de Resistencia**:
   - Extracción de variantes en genes target
   - Generación de reporte específico
   - Resumen estadístico final

## 📊 Genoma de Referencia

- **Organismo**: Leishmania braziliensis MHOM/BR/75/M2904
- **Fuente**: NCBI RefSeq (GCF_000002845.2_ASM284v2)
- **Descarga automática**: El pipeline descarga automáticamente los archivos de referencia

## 🧪 Validación

El sistema ha sido validado con:
- ✅ Muestras de control conocidas
- ✅ Secuencias con mutaciones reportadas en literatura
- ✅ Datos paired-end de Illumina
- ✅ Archivos FASTQ comprimidos (.fastq.gz)

## 📋 Requisitos del Sistema

- **Docker** y **Docker Compose**
- **8GB RAM** mínimo (recomendado 16GB)
- **20GB** espacio libre en disco
- **Conexión a internet** (primera ejecución para descargar referencias)

## 🔧 Configuración Avanzada

### Variables de Pipeline:
- `THREADS=8`: Número de hilos de procesamiento
- `MAX_DEPTH=1000`: Profundidad máxima para mpileup
- Filtros de calidad personalizables en `pipeline.sh`

### Puertos de Red:
- **5001**: Interfaz web Flask
- **5000**: Pipeline interno (Docker)

## 📝 Licencia

Proyecto académico de investigación desarrollado para:
- Uso exclusivo con fines educativos y científicos
- Análisis de resistencia farmacológica en Leishmania
- Contribución al conocimiento en parasitología molecular

---

**Desarrollado como parte del proyecto de investigación PIS-3**
*Universidad de Antioquia - Facultad de Ingeniería - Bioingeniería*
