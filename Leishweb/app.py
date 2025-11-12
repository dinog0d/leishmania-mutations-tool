from flask import Flask, request, jsonify, render_template, Response, redirect, url_for
import os
import subprocess
import hashlib
import json
import uuid
from datetime import datetime
import threading

app = Flask(__name__)
UPLOAD_FOLDER = '/data'
JOBS_FOLDER = '/data/jobs'
os.makedirs(UPLOAD_FOLDER, exist_ok=True)
os.makedirs(JOBS_FOLDER, exist_ok=True)
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

def generate_job_id():
    """Generar un ID único para el job"""
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    unique_id = str(uuid.uuid4())[:8]
    return f"job_{timestamp}_{unique_id}"

def save_job_metadata(job_data):
    """Guardar metadata del job en archivo JSON"""
    job_file = os.path.join(JOBS_FOLDER, f"{job_data['job_id']}.json")
    with open(job_file, 'w') as f:
        json.dump(job_data, f, indent=2)
    
    # Actualizar índice de jobs
    update_jobs_index(job_data['job_id'])

def load_job_metadata(job_id):
    """Cargar metadata de un job específico"""
    job_file = os.path.join(JOBS_FOLDER, f"{job_id}.json")
    if os.path.exists(job_file):
        with open(job_file, 'r') as f:
            return json.load(f)
    return None

def update_jobs_index(job_id):
    """Actualizar el índice de todos los jobs"""
    index_file = os.path.join(JOBS_FOLDER, 'index.json')
    
    # Cargar índice existente
    if os.path.exists(index_file):
        with open(index_file, 'r') as f:
            index = json.load(f)
    else:
        index = {"jobs": []}
    
    # Agregar nuevo job si no existe
    if job_id not in index["jobs"]:
        index["jobs"].insert(0, job_id)  # Insertar al principio (más reciente)
    
    # Guardar índice actualizado
    with open(index_file, 'w') as f:
        json.dump(index, f, indent=2)

def update_job_status(job_id, status, **kwargs):
    """Actualizar el estado de un job"""
    job_data = load_job_metadata(job_id)
    if job_data:
        job_data['status'] = status
        if status == 'running' and 'started_at' not in job_data:
            job_data['started_at'] = datetime.now().isoformat()
        elif status in ['completed', 'failed']:
            job_data['completed_at'] = datetime.now().isoformat()
        
        # Actualizar campos adicionales si se proporcionan
        for key, value in kwargs.items():
            job_data[key] = value
        
        save_job_metadata(job_data)
        return job_data
    return None

def run_pipeline_async(job_id):
    """Ejecutar el pipeline de manera asíncrona con logs en tiempo real"""
    def pipeline_thread():
        try:
            update_job_status(job_id, 'running')
            
            # Ejecutar el pipeline pasando el job_id como parámetro
            process = subprocess.Popen(
                ['docker', 'exec', 'leishmania-mutations-tool-pipeline-1', 'bash', '/app/pipeline.sh', job_id],
                stdout=subprocess.PIPE, 
                stderr=subprocess.STDOUT, 
                text=True,
                bufsize=1,  # Line buffered
                universal_newlines=True
            )
            
            # Capturar logs en tiempo real y guardarlos
            log_file_path = os.path.join(JOBS_FOLDER, f"{job_id}_logs.txt")
            with open(log_file_path, 'w') as log_file:
                for line in process.stdout:
                    # Escribir al archivo de log
                    log_file.write(line)
                    log_file.flush()
                    
                    # También mostrar en consola para debug
                    print(f"[{job_id}] {line.strip()}")
            
            # Esperar a que termine
            process.wait()
            
            if process.returncode == 0:
                update_job_status(job_id, 'completed', 
                                results_path=f'/results/{job_id}/',
                                log_file=f'{job_id}_logs.txt')
            else:
                update_job_status(job_id, 'failed', 
                                error_message=f'Pipeline execution failed with code {process.returncode}',
                                log_file=f'{job_id}_logs.txt')
                
        except Exception as e:
            update_job_status(job_id, 'failed', 
                            error_message=str(e),
                            log_file=f'{job_id}_logs.txt')
    
    # Ejecutar en thread separado
    thread = threading.Thread(target=pipeline_thread)
    thread.daemon = True
    thread.start()

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

    # Generar nuevo job
    job_id = generate_job_id()
    
    # Crear estructura de carpetas para este job específico
    job_folder = os.path.join(app.config['UPLOAD_FOLDER'], job_id)
    os.makedirs(job_folder, exist_ok=True)

    # Guardar archivos con nombres que el pipeline entiende
    filepath1 = os.path.join(job_folder, 'sample_leishmania_R1.fastq.gz')
    filepath2 = os.path.join(job_folder, 'sample_leishmania_R2.fastq.gz')

    file1.save(filepath1)
    file2.save(filepath2)

    # Crear metadata del job
    job_data = {
        'job_id': job_id,
        'status': 'pending',
        'created_at': datetime.now().isoformat(),
        'files': {
            'r1_original': file1.filename,
            'r2_original': file2.filename,
            'r1_path': filepath1,
            'r2_path': filepath2
        },
        'results_path': f'/results/{job_id}/'
    }
    
    save_job_metadata(job_data)
    
    # Redirigir a la nueva ruta con job ID
    return redirect(url_for('run_pipeline_job', job_id=job_id))

@app.route('/logs/<job_id>')
def view_job_logs(job_id):
    """Página para mostrar logs en tiempo real de un job específico"""
    job_data = load_job_metadata(job_id)
    if not job_data:
        return "Job not found", 404
    
    # HTML con auto-refresh para logs en tiempo real
    html_content = f"""
    <!DOCTYPE html>
    <html>
    <head>
        <title>Logs en Tiempo Real - Job {job_id}</title>
        <style>
            body {{ font-family: 'Courier New', monospace; margin: 0; padding: 20px; background-color: #1a1a1a; color: #00ff00; }}
            .container {{ max-width: 1200px; margin: 0 auto; }}
            .header {{ background: #2c3e50; color: white; padding: 15px; border-radius: 8px 8px 0 0; }}
            .logs {{ background: #000; color: #00ff00; padding: 20px; border-radius: 0 0 8px 8px; height: 70vh; overflow-y: auto; white-space: pre-wrap; font-size: 14px; border: 1px solid #333; }}
            .nav-links {{ margin: 20px 0; }}
            .nav-links a {{ background: #3498db; color: white; padding: 10px 15px; text-decoration: none; border-radius: 5px; margin-right: 10px; }}
            .nav-links a:hover {{ background: #2980b9; }}
            .status {{ display: inline-block; padding: 5px 10px; border-radius: 3px; font-weight: bold; }}
            .status.running {{ background: #f39c12; color: white; }}
            .status.completed {{ background: #27ae60; color: white; }}
            .status.failed {{ background: #e74c3c; color: white; }}
            .status.pending {{ background: #95a5a6; color: white; }}
        </style>
        <script>
            let refreshInterval;
            
            function refreshLogs() {{
                fetch('/api/jobs/{job_id}/logs')
                    .then(response => response.json())
                    .then(data => {{
                        if (data.logs) {{
                            document.getElementById('logs').textContent = data.logs;
                            document.getElementById('logs').scrollTop = document.getElementById('logs').scrollHeight;
                        }}
                    }})
                    .catch(error => console.error('Error:', error));
                
                // También actualizar el estado del job
                fetch('/api/jobs/{job_id}')
                    .then(response => response.json())
                    .then(data => {{
                        document.getElementById('status').textContent = data.status;
                        document.getElementById('status').className = 'status ' + data.status;
                        
                        // Detener refresh si el job ya terminó
                        if (data.status === 'completed' || data.status === 'failed') {{
                            if (refreshInterval) {{
                                clearInterval(refreshInterval);
                                refreshInterval = null;
                            }}
                            document.getElementById('logs').innerHTML += '\\n\\n=== Pipeline terminado - Los logs ya no se actualizarán ===';
                        }}
                    }})
                    .catch(error => console.error('Error:', error));
            }}
            
            // Refresh cada 2 segundos
            refreshInterval = setInterval(refreshLogs, 2000);
            
            // Refresh inicial cuando carga la página
            window.onload = refreshLogs;
        </script>
    </head>
    <body>
        <div class="container">
            <div class="header">
                <h1>🔍 Logs en Tiempo Real - Job {job_id}</h1>
                <p><strong>Estado:</strong> <span id="status" class="status {job_data['status']}">{job_data['status']}</span></p>
                <p><strong>Archivos:</strong> {job_data['files']['r1_original']}, {job_data['files']['r2_original']}</p>
            </div>
            
            <div class="nav-links">
                <a href="/api/jobs/{job_id}/logs" download>💾 Descargar logs</a>
                <a href="/results/{job_id}">📊 Ver resultados</a>
                <a href="/">🏠 Volver al inicio</a>
                <button onclick="refreshLogs()" style="background: #27ae60; color: white; border: none; padding: 10px 15px; border-radius: 5px; cursor: pointer;">🔄 Actualizar ahora</button>
            </div>
            
            <div id="logs" class="logs">Cargando logs...</div>
        </div>
    </body>
    </html>
    """
    
    return html_content

# Nuevas rutas API para el sistema de jobs
@app.route('/api/jobs', methods=['POST'])
def create_job():
    """API para crear un nuevo job via AJAX"""
    if 'file1' not in request.files or 'file2' not in request.files:
        return jsonify({'error': 'Both files are required'}), 400

    file1 = request.files['file1']
    file2 = request.files['file2']

    if file1.filename == '' or file2.filename == '':
        return jsonify({'error': 'Both files must have a name'}), 400

    # Validaciones
    if not allowed_file(file1.filename):
        return jsonify({'error': f"File '{file1.filename}' is not a valid FASTQ file. Only .fastq.gz files are allowed."}), 400
    
    if not allowed_file(file2.filename):
        return jsonify({'error': f"File '{file2.filename}' is not a valid FASTQ file. Only .fastq.gz files are allowed."}), 400

    # Validar duplicados
    file1_hash = get_file_hash(file1)
    file2_hash = get_file_hash(file2)
    
    if file1_hash == file2_hash:
        return jsonify({'error': 'Both files are identical. Please upload different R1 and R2 files.'}), 400

    # Crear job
    job_id = generate_job_id()
    job_folder = os.path.join(app.config['UPLOAD_FOLDER'], job_id)
    os.makedirs(job_folder, exist_ok=True)

    filepath1 = os.path.join(job_folder, 'sample_leishmania_R1.fastq.gz')
    filepath2 = os.path.join(job_folder, 'sample_leishmania_R2.fastq.gz')

    file1.save(filepath1)
    file2.save(filepath2)

    job_data = {
        'job_id': job_id,
        'status': 'pending',
        'created_at': datetime.now().isoformat(),
        'files': {
            'r1_original': file1.filename,
            'r2_original': file2.filename,
            'r1_path': filepath1,
            'r2_path': filepath2
        },
        'results_path': f'/results/{job_id}/'
    }
    
    save_job_metadata(job_data)
    
    return jsonify({'success': True, 'job_id': job_id, 'job_data': job_data}), 201

@app.route('/api/jobs', methods=['GET'])
def list_jobs():
    """API para listar todos los jobs"""
    index_file = os.path.join(JOBS_FOLDER, 'index.json')
    if not os.path.exists(index_file):
        return jsonify({'jobs': []})
    
    with open(index_file, 'r') as f:
        index = json.load(f)
    
    jobs = []
    for job_id in index.get('jobs', []):
        job_data = load_job_metadata(job_id)
        if job_data:
            jobs.append(job_data)
    
    return jsonify({'jobs': jobs})

@app.route('/api/jobs/<job_id>', methods=['GET'])
def get_job(job_id):
    """API para obtener información específica de un job"""
    job_data = load_job_metadata(job_id)
    if not job_data:
        return jsonify({'error': 'Job not found'}), 404
    
    return jsonify(job_data)

@app.route('/api/jobs/<job_id>/start', methods=['POST'])
def start_job(job_id):
    """API para iniciar el pipeline de un job"""
    job_data = load_job_metadata(job_id)
    if not job_data:
        return jsonify({'error': 'Job not found'}), 404
    
    if job_data['status'] != 'pending':
        return jsonify({'error': f"Job is in '{job_data['status']}' status and cannot be started"}), 400
    
    # Copiar archivos a la ubicación que espera el pipeline
    sample_folder = os.path.join(app.config['UPLOAD_FOLDER'], 'sample_leishmania')
    os.makedirs(sample_folder, exist_ok=True)
    
    import shutil
    shutil.copy2(job_data['files']['r1_path'], os.path.join(sample_folder, 'sample_leishmania_R1.fastq.gz'))
    shutil.copy2(job_data['files']['r2_path'], os.path.join(sample_folder, 'sample_leishmania_R2.fastq.gz'))
    
    # Iniciar pipeline de manera asíncrona
    run_pipeline_async(job_id)
    
    return jsonify({'success': True, 'message': 'Pipeline started', 'job_id': job_id})

@app.route('/api/jobs/<job_id>/stream')
def job_stream(job_id):
    """Server-Sent Events stream para seguimiento en tiempo real"""
    def event_generator():
        job_data = load_job_metadata(job_id)
        if not job_data:
            yield f"data: {json.dumps({'error': 'Job not found'})}\n\n"
            return
        
        while True:
            job_data = load_job_metadata(job_id)
            if not job_data:
                break
                
            yield f"data: {json.dumps(job_data)}\n\n"
            
            if job_data['status'] in ['completed', 'failed']:
                break
                
            import time
            time.sleep(2)  # Poll cada 2 segundos
    
    return Response(event_generator(), mimetype='text/event-stream')

@app.route('/job/<job_id>')
def run_pipeline_job(job_id):
    """Ruta para ejecutar pipeline de un job específico o mostrar logs si ya está completado"""
    job_data = load_job_metadata(job_id)
    if not job_data:
        return "Job not found", 404
    
    # Si el job ya está completado o falló, redirigir a logs
    if job_data['status'] in ['completed', 'failed']:
        return redirect(url_for('view_job_logs', job_id=job_id))
    
    # Verificar archivos
    if not os.path.exists(job_data['files']['r1_path']) or not os.path.exists(job_data['files']['r2_path']):
        return "Required files are missing", 400
    
    # Si el job está pending, iniciarlo automáticamente
    if job_data['status'] == 'pending':
        run_pipeline_async(job_id)

    # Redirigir a página de logs para monitorear progreso
    return redirect(url_for('view_job_logs', job_id=job_id))

@app.route('/api/jobs/<job_id>/logs')
def get_job_logs(job_id):
    """API para obtener los logs de un job específico"""
    log_file_path = os.path.join(JOBS_FOLDER, f"{job_id}_logs.txt")
    
    if not os.path.exists(log_file_path):
        return jsonify({'error': 'Log file not found'}), 404
    
    try:
        with open(log_file_path, 'r') as f:
            logs = f.read()
        return jsonify({'logs': logs, 'job_id': job_id})
    except Exception as e:
        return jsonify({'error': f'Error reading logs: {str(e)}'}), 500

@app.route('/api/jobs/<job_id>/logs/stream')
def stream_job_logs(job_id):
    """Stream de logs en tiempo real para un job específico"""
    def generate_log_stream():
        log_file_path = os.path.join(JOBS_FOLDER, f"{job_id}_logs.txt")
        
        # Si el archivo de logs no existe todavía, esperar
        import time
        while not os.path.exists(log_file_path):
            job_data = load_job_metadata(job_id)
            if not job_data or job_data['status'] in ['completed', 'failed']:
                yield f"data: Job finished without log file\n\n"
                return
            time.sleep(1)
        
        # Leer logs en tiempo real
        with open(log_file_path, 'r') as f:
            # Enviar logs existentes primero
            existing_content = f.read()
            if existing_content:
                for line in existing_content.splitlines():
                    yield f"data: {line}\n\n"
            
            # Continuar leyendo nuevas líneas
            while True:
                line = f.readline()
                if line:
                    yield f"data: {line.rstrip()}\n\n"
                else:
                    # Verificar si el job ya terminó
                    job_data = load_job_metadata(job_id)
                    if not job_data or job_data['status'] in ['completed', 'failed']:
                        break
                    time.sleep(0.5)  # Esperar medio segundo antes de verificar nuevas líneas
    
    return Response(generate_log_stream(), mimetype='text/event-stream')

# Mantener ruta original para compatibilidad
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

@app.route('/results/<job_id>')
def job_results(job_id):
    """Mostrar resultados de un job específico"""
    job_data = load_job_metadata(job_id)
    if not job_data:
        return "Job not found", 404
    
    if job_data['status'] != 'completed':
        return f"Job is in '{job_data['status']}' status. Results not available.", 400
    
    # Buscar archivo de resultados en la carpeta específica del job
    results_path = f'/results/{job_id}'
    tvs_file = os.path.join(results_path, 'reports', 'sample_leishmania_resistance_report.tsv')
    summary_file = os.path.join(results_path, 'reports', 'sample_leishmania_summary.txt')
    
    if not os.path.exists(tvs_file):
        return "Results file not found", 404
    
    # Leer archivo de resultados
    with open(tvs_file, 'r') as f:
        results_content = f.read()
    
    # Parsear TSV para crear tabla HTML más didáctica
    lines = results_content.strip().split('\n')
    if len(lines) > 1:
        headers = lines[0].split('\t')
        rows = [line.split('\t') for line in lines[1:]]
        
        # Crear tabla HTML
        table_html = '<table class="results-table">'
        table_html += '<thead><tr>'
        for header in headers:
            table_html += f'<th>{header}</th>'
        table_html += '</tr></thead><tbody>'
        
        for row in rows:
            table_html += '<tr>'
            for cell in row:
                table_html += f'<td>{cell}</td>'
            table_html += '</tr>'
        table_html += '</tbody></table>'
    else:
        table_html = '<p>No hay datos de resistencia disponibles.</p>'
    
    # Leer resumen si existe
    summary_content = ""
    if os.path.exists(summary_file):
        with open(summary_file, 'r') as f:
            summary_content = f.read()
    
    # Generar HTML con formato mejorado
    html_content = f"""
    <!DOCTYPE html>
    <html>
    <head>
        <title>Resultados - Job {job_id}</title>
        <style>
            body {{ font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif; margin: 20px; background-color: #f5f5f5; }}
            .container {{ max-width: 1200px; margin: 0 auto; background: white; padding: 20px; border-radius: 10px; box-shadow: 0 2px 10px rgba(0,0,0,0.1); }}
            h1 {{ color: #2c3e50; border-bottom: 3px solid #3498db; padding-bottom: 10px; }}
            h2 {{ color: #34495e; margin-top: 30px; }}
            .job-info {{ background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); color: white; padding: 15px; border-radius: 8px; margin-bottom: 20px; }}
            .summary {{ background: #ecf0f1; padding: 15px; border-radius: 8px; margin-bottom: 20px; white-space: pre-wrap; font-family: monospace; }}
            .results {{ background: #2c3e50; color: #ecf0f1; padding: 15px; border-radius: 8px; white-space: pre-wrap; font-family: 'Courier New', monospace; font-size: 12px; overflow-x: auto; }}
            .results-table {{ width: 100%; border-collapse: collapse; margin: 20px 0; }}
            .results-table th, .results-table td {{ border: 1px solid #ddd; padding: 12px; text-align: left; }}
            .results-table th {{ background-color: #34495e; color: white; font-weight: bold; }}
            .results-table tr:nth-child(even) {{ background-color: #f2f2f2; }}
            .results-table tr:hover {{ background-color: #e8f4f8; }}
            .nav-links {{ margin: 20px 0; }}
            .nav-links a {{ background: #3498db; color: white; padding: 10px 15px; text-decoration: none; border-radius: 5px; margin-right: 10px; }}
            .nav-links a:hover {{ background: #2980b9; }}
        </style>
    </head>
    <body>
        <div class="container">
            <h1>🧬 Resultados del Análisis de Resistencia</h1>
            
            <div class="job-info">
                <strong>Job ID:</strong> {job_id}<br>
                <strong>Estado:</strong> {job_data['status']}<br>
                <strong>Archivos analizados:</strong> {job_data['files']['r1_original']}, {job_data['files']['r2_original']}<br>
                <strong>Fecha de creación:</strong> {job_data.get('created_at', 'N/A')}<br>
                <strong>Fecha de finalización:</strong> {job_data.get('completed_at', 'N/A')}
            </div>
            
            <div class="nav-links">
                <a href="/logs/{job_id}">📄 Ver logs</a>
                <a href="/api/jobs/{job_id}/logs" download>💾 Descargar logs</a>
                <a href="/">🏠 Volver al inicio</a>
            </div>
"""
    
    # Agregar resumen si existe
    if summary_content:
        html_content += f"""
            <h2>📊 Resumen del Análisis</h2>
            <div class="summary">{summary_content}</div>
        """
    
    # Agregar tabla de resultados
    html_content += f"""
            <h2>📋 Tabla de Resultados de Resistencia</h2>
            <div class="table-container">
                {table_html}
            </div>
            
            <h2>📄 Resultados Raw (TSV)</h2>
            <div class="results">{results_content}</div>
            
        </div>
    </body>
    </html>
    """
    
    return html_content

if __name__ == '__main__':
    app.run(host='0.0.0.0', port=5000)