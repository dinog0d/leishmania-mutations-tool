#!/bin/bash

echo "Limpiando todos los datos anteriores..."

# Detener contenedores
docker-compose down

# Limpiar datos anteriores
sudo rm -rf data/job_*
sudo rm -rf data/jobs  
sudo rm -rf data/sample_leishmania
sudo rm -rf results/job_*

# Crear estructura limpia
mkdir -p data/jobs
mkdir -p results

# Establecer permisos correctos
sudo chown -R $USER:$USER data results
sudo chmod -R 755 data results

echo "Estructura de datos limpia:"
echo "- data/jobs/ (para metadatos de trabajos)"
echo "- results/ (para resultados de análisis)"
echo ""
echo "¡Listo para iniciar servicios limpios!"