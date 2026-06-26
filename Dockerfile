FROM python:3.12-slim

WORKDIR /app

COPY . /app

RUN python -m pip install --upgrade pip
RUN python -m pip install --no-cache-dir -r requirements.txt

EXPOSE 8501

CMD ["python", "-m", "streamlit", "run", "mendelcell/app.py", "--server.address=0.0.0.0", "--server.port=8501", "--server.headless=true", "--server.enableXsrfProtection=false", "--server.enableCORS=false"]