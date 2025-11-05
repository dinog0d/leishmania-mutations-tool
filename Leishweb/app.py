from flask import Flask, request, jsonify, render_template
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

    filepath1 = os.path.join(app.config['UPLOAD_FOLDER'], file1.filename)
    filepath2 = os.path.join(app.config['UPLOAD_FOLDER'], file2.filename)

    file1.save(filepath1)
    file2.save(filepath2)

    # Run the pipeline subprocess
    try:
        subprocess.run(['bash', '/Pipeline/pipeline.sh'], check=True)
    except subprocess.CalledProcessError as e:
        return f"Pipeline execution failed: {e}", 500

    return "Files uploaded and pipeline executed successfully!"

@app.route('/results')
def results():
    results_path = '/results'
    if not os.path.exists(results_path):
        return "No results available", 404
    files = os.listdir(results_path)
    return jsonify(files)

if __name__ == '__main__':
    app.run(host='0.0.0.0', port=5000)