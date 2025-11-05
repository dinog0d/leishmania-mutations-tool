# Dockerfile para la aplicación Flask
FROM python:3.9-slim

WORKDIR /app

RUN apt-get update && apt-get install -y curl docker.io && rm -rf /var/lib/apt/lists/*

COPY requirements.txt requirements.txt
RUN pip install --no-cache-dir -r requirements.txt

COPY . .

CMD ["python", "app.py"]