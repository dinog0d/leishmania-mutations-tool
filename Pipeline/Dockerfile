# ============================================================
# Dockerfile para Leishmania Variant Calling Pipeline
# ============================================================

# Imagen base con Miniconda (más ligera que Anaconda)
FROM continuumio/miniconda3:latest

# Directorio de trabajo dentro del contenedor
WORKDIR /app

# Copiamos el archivo de entorno con todas las dependencias
COPY Pipeline/environment.yml /tmp/environment.yml

# Creamos el entorno conda a partir del environment.yml
RUN conda env create -f /tmp/environment.yml

# Nombre del entorno (debe coincidir con el definido en environment.yml)
ARG ENV_NAME=leish-env

# Aseguramos que el entorno esté en el PATH
ENV PATH=/opt/conda/envs/${ENV_NAME}/bin:$PATH

# Usamos conda run como shell por defecto para que todo se ejecute dentro del entorno
SHELL ["conda", "run", "-n", "leish-env", "/bin/bash", "-c"]

# Copiamos el pipeline al contenedor
COPY Pipeline/pipeline.sh /app/pipeline.sh
RUN chmod +x /app/pipeline.sh

# Definimos volúmenes para datos, referencias y resultados
VOLUME ["/data", "/reference", "/results"]

# Comando por defecto al ejecutar el contenedor
ENTRYPOINT ["/bin/bash", "/app/pipeline.sh"]
