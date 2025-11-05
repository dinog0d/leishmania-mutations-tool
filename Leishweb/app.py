from flask import Flask, request, jsonify, render_template, Response, redirect, url_for
import os
import subprocess

app = Flask(__name__)
UPLOAD_FOLDER = '/data'
os.makedirs(UPLOAD_FOLDER, exist_ok=True)
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER
app.config['MAX_CONTENT_LENGTH'] = 2 * 1024 * 1024 * 1024  # 2 GB limit

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