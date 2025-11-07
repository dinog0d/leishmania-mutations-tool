from flask import Flask, request, jsonify, render_template, Response, redirect, url_for
import os
import subprocess
import hashlib

app = Flask(__name__)
UPLOAD_FOLDER = '/data'
os.makedirs(UPLOAD_FOLDER, exist_ok=True)
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER
app.config['MAX_CONTENT_LENGTH'] = 2 * 1024 * 1024 * 1024  # 2 GB limit

def allowed_file(filename):
    """Verificar que el archivo tenga la extensión .fastq.gz"""
    return filename.lower().endswith('.fastq.gz')

def get_file_hash(file):
    """Calcular el hash MD5 de un archivo para detectar duplicados"""
    file.seek(0)  # Asegurar que estamos al inicio del archivo
    hash_md5 = hashlib.md5()
    for chunk in iter(lambda: file.read(4096), b""):
        hash_md5.update(chunk)
    file.seek(0)  # Volver al inicio para uso posterior
    return hash_md5.hexdigest()

@app.route('/')
def home():
    return render_template('index.html')

@app.route('/upload', methods=['POST'])
def upload_files():
    if 'file1' not in request.files or 'file2' not in request.files:
        return "Both files are required", 400

    file1 = request.files['file1']
    file2 = request.files['file2']

    if file1.filename == '' or file2.filename == '':
        return "Both files must have a name", 400

    # Validar extensión de archivos
    if not allowed_file(file1.filename):
        return f"File '{file1.filename}' is not a valid FASTQ file. Only .fastq.gz files are allowed.", 400
    
    if not allowed_file(file2.filename):
        return f"File '{file2.filename}' is not a valid FASTQ file. Only .fastq.gz files are allowed.", 400

    # Validar que no sean el mismo archivo (comparando contenido)
    file1_hash = get_file_hash(file1)
    file2_hash = get_file_hash(file2)
    
    if file1_hash == file2_hash:
        return "Error: Both files are identical. Please upload different R1 and R2 files.", 400

    # Crear la estructura de carpetas que espera el pipeline
    sample_folder = os.path.join(app.config['UPLOAD_FOLDER'], 'sample01')
    os.makedirs(sample_folder, exist_ok=True)

    # Guardar los archivos con los nombres exactos que espera el pipeline
    filepath1 = os.path.join(sample_folder, 'sample01_R1.fastq.gz')
    filepath2 = os.path.join(sample_folder, 'sample01_R2.fastq.gz')

    file1.save(filepath1)
    file2.save(filepath2)

    return redirect(url_for('run_pipeline'))

@app.route('/run_pipeline', methods=['GET'])
def run_pipeline():
    # Check if required files exist in the correct location
    sample_folder = os.path.join(app.config['UPLOAD_FOLDER'], 'sample01')
    filepath1 = os.path.join(sample_folder, 'sample01_R1.fastq.gz')
    filepath2 = os.path.join(sample_folder, 'sample01_R2.fastq.gz')

    if not os.path.exists(filepath1) or not os.path.exists(filepath2):
        return "Required files are missing. Please upload them first.", 400

    def generate():
        process = subprocess.Popen(['docker', 'exec', 'leishmania-mutations-tool-pipeline-1', 'bash', '/app/pipeline.sh'],
                                   stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        for line in process.stdout:
            yield f"data:{line}\n\n"

    return Response(generate(), mimetype='text/event-stream')

@app.route('/results')
def results():
    results_path = '/results/sample01'
    tvs_file = os.path.join(results_path, 'reports', 'sample01_resistance_report.tsv')
    if not os.path.exists(tvs_file):
        return "No results available", 404
    with open(tvs_file, 'r') as f:
        content = f.read()
    return f"<pre>{content}</pre>"

if __name__ == '__main__':
    app.run(host='0.0.0.0', port=5000)