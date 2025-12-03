FROM python:3.12.3

#installing ncbi to run the blastn
RUN apt-get update && \
    apt-get install -y ncbi-blast+ && \
    apt-get clean


#Also need to install pgmauv
RUN pip install pygenomeviz

#install progressiveMauve from local tar.gz 
COPY mauve_linux_snapshot_2015-02-13.tar.gz /tmp/mauve.tar.gz
RUN mkdir -p /opt/mauve && \
    tar -xzf /tmp/mauve.tar.gz -C /opt/mauve --strip-components=1 && \
    ln -s /opt/mauve/linux-x64/progressiveMauve /usr/local/bin/progressiveMauve && \
    rm /tmp/mauve.tar.gz


# Install mummer for nucmer
RUN apt-get update && \
    apt-get install -y mummer && \
    apt-get clean

# setting up the work directory
WORKDIR /app

#copying requirements
COPY requirements.txt /app/requirements.txt
RUN pip install --no-cache-dir -r /app/requirements.txt

#copying the origenomi cli
COPY README.md /app/
COPY pyproject.toml /app/
COPY src /app/src/

#copying the fastapi server for the same tool
COPY temp_anant/test/app /app/app

#copying the databases
COPY clustered_DB /app/clustered_DB/

RUN pip install -e .

CMD ["uvicorn", "app.main:app", "--host", "0.0.0.0", "--port", "8000"]
