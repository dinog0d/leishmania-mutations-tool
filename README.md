**Herramienta Bioinformática para la Detección de Mutaciones en Leishmania braziliensis**

Este proyecto hace parte de la Propuesta de Investigación PIS-3:
“Desarrollo de una herramienta informática para la detección de mutaciones en genes asociados a resistencia farmacológica en cepas colombianas de Leishmania braziliensis”.

El software está desarrollado con:

 - Python (FastAPI) para el backend.

 - HTML, CSS y JavaScript para el frontend web.

El objetivo es procesar archivos FASTA con secuencias genómicas, detectar mutaciones (SNPs, inserciones, deleciones), clasificarlas de acuerdo con su efecto y mostrar los resultados en una interfaz accesible desde cualquier navegador moderno.

**Objetivos**

- Procesar secuencias de ADN en formato FASTA.

- Detectar variantes genómicas (SNPs, indels).

- Clasificar mutaciones (sinónimas, missense, nonsense).

- Implementar una API REST en FastAPI que procese las secuencias.

- Desplegar una interfaz web en HTML + CSS + JS para cargar archivos y visualizar resultados.

- Validar el prototipo con secuencias que presentan mutaciones reportadas en la literatura.

La página permite cargar un archivo FASTA → lo envía al backend → muestra resultados.

**📂 Estructura del proyecto**
```
leishmania-mutations-tool/
├── backend/               # Backend en Python con FastAPI
│   └── main.py
├── frontend/              # Frontend en HTML, CSS, JS
│   ├── index.html
│   ├── styles.css
│   └── app.js
├── src/                   # Algoritmos bioinformáticos
│   ├── fasta_reader.py
│   ├── variant_finder.py
│   └── protein_effects.py
├── tests/                 # Pruebas unitarias (pytest)
├── requirements.txt       # Dependencias del proyecto
├── .gitignore             # Archivos ignorados
└── README.md              # Descripción del proyecto
```

**✅ Estado actual**

 X Implementación de lector FASTA

 X Detección básica de variantes (SNP/Indel)

 X Clasificación de mutaciones proteicas

 X API REST con FastAPI

 X Frontend en HTML/CSS/JS conectado al backend

 X Validación con secuencias conocidas

**Licencia**

Proyecto académico de investigación.
Uso exclusivo con fines educativos y científicos.
