FROM python:3.9

RUN apt-get update && \
    apt-get install -y build-essential python3-dev wget && \
    apt-get clean

# Install RDKit
RUN pip install --upgrade pip
RUN pip install rdkit-pypi streamlit pandas matplotlib

# Copy app code
COPY . /app
WORKDIR /app

EXPOSE 7860

CMD ["streamlit", "run", "app.py", "--server.port=7860", "--server.enableCORS=false"]
