from flask import Flask, request, jsonify, render_template, Response, send_file
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

def update_job_progress(job_id, step, total_steps, step_description):
    """Actualizar el progreso específico del pipeline"""
    job_data = load_job_metadata(job_id)
    if job_data:
        job_data['progress'] = {
            'current_step': step,
            'total_steps': total_steps,
            'percentage': int((step / total_steps) * 100),
            'step_description': step_description,
            'last_updated': datetime.now().isoformat()
        }
        save_job_metadata(job_data)
        return job_data
    return None

def run_pipeline_async(job_id):
    """Ejecutar el pipeline de manera asíncrona con logs en tiempo real y seguimiento de progreso"""
    def pipeline_thread():
        try:
            update_job_status(job_id, 'running')
            update_job_progress(job_id, 0, 7, "Iniciando pipeline...")
            
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
            
            # Diccionario para mapear etapas del pipeline
            progress_mapping = {
                "[PASO 0/7]": (0, "Verificando referencias"),
                "[PASO 1/7]": (1, "Filtrando lecturas"),
                "[PASO 2/7]": (2, "Alineando lecturas"),
                "[PASO 3/7]": (3, "Procesando archivos BAM"),
                "[PASO 4/7]": (4, "Llamando variantes"),
                "[PASO 5/7]": (5, "Filtrando variantes"),
                "[PASO 6/7]": (6, "Anotando variantes"),
                "[PASO 7/7]": (7, "Extrayendo genes de resistencia")
            }
            
            with open(log_file_path, 'w') as log_file:
                for line in process.stdout:
                    # Escribir al archivo de log
                    log_file.write(line)
                    log_file.flush()
                    
                    # Actualizar progreso basado en las etapas del pipeline
                    for step_marker, (step_num, description) in progress_mapping.items():
                        if step_marker in line:
                            update_job_progress(job_id, step_num, 7, description)
                            print(f"[{job_id}] Progreso: {step_num}/7 - {description}")
                            break
                    
                    # También mostrar en consola para debug
                    print(f"[{job_id}] {line.strip()}")
            
            # Esperar a que termine
            process.wait()
            
            if process.returncode == 0:
                update_job_progress(job_id, 7, 7, "Pipeline completado")
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
    
    print(f"🚀 Auto-iniciando pipeline para job {job_id}")
    run_pipeline_async(job_id)
    
    return jsonify({'success': True, 'job_id': job_id, 'job_data': job_data, 'auto_started': True}), 201

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

@app.route('/api/jobs/<job_id>/results')
def get_job_results_api(job_id):
    """API para obtener los resultados de un job específico en formato JSON"""
    job_data = load_job_metadata(job_id)
    if not job_data:
        return jsonify({'error': 'Job not found'}), 404
    
    if job_data['status'] != 'completed':
        return jsonify({'error': f"Job is in '{job_data['status']}' status. Results not available."}), 400
    
    # Buscar archivo de resultados en la carpeta específica del job
    results_path = f'/results/{job_id}'
    tvs_file = os.path.join(results_path, 'reports', 'sample_leishmania_resistance_report.tsv')
    summary_file = os.path.join(results_path, 'reports', 'sample_leishmania_summary.txt')
    
    if not os.path.exists(tvs_file):
        return jsonify({'error': 'Results file not found'}), 404
    
    # Leer archivo de resultados
    with open(tvs_file, 'r') as f:
        results_content = f.read()
    
    # Parsear TSV para crear datos estructurados
    lines = results_content.strip().split('\n')
    results_data = {
        'job_id': job_id,
        'job_data': job_data,
        'results': {
            'headers': [],
            'rows': [],
            'raw_content': results_content
        }
    }
    
    if len(lines) > 1:
        results_data['results']['headers'] = lines[0].split('\t')
        results_data['results']['rows'] = [line.split('\t') for line in lines[1:]]
    
    # Leer resumen si existe
    if os.path.exists(summary_file):
        with open(summary_file, 'r') as f:
            results_data['summary'] = f.read()
    
    return jsonify(results_data)

@app.route('/api/jobs/<job_id>/download/results')
def download_job_results(job_id):
    """API para descargar los resultados de un job como archivo TSV"""
    results_path = f'/results/{job_id}'
    tvs_file = os.path.join(results_path, 'reports', 'sample_leishmania_resistance_report.tsv')
    
    if not os.path.exists(tvs_file):
        return jsonify({'error': 'Results file not found'}), 404
    
    return send_file(tvs_file, as_attachment=True, download_name=f'leishmania_results_{job_id}.tsv')

@app.route('/api/jobs/<job_id>/download/logs')
def download_job_logs(job_id):
    """API para descargar los logs de un job"""
    log_file_path = os.path.join(JOBS_FOLDER, f"{job_id}_logs.txt")
    
    if not os.path.exists(log_file_path):
        return jsonify({'error': 'Log file not found'}), 404
    
    return send_file(log_file_path, as_attachment=True, download_name=f'leishmania_logs_{job_id}.txt')

# Mantener solo ruta original para compatibilidad (si es necesaria)
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

if __name__ == '__main__':
    app.run(host='0.0.0.0', port=5000)