# Base minimal Jupyter image
FROM jupyter/minimal-notebook:python-3.11

# -----------------------------
# 1️Copy precompiled Python module
# -----------------------------
COPY ../lib /home/jovyan/CPlantBox/lib

# -----------------------------
# 2️Copy notebooks
# -----------------------------
COPY ../scripts /home/jovyan/scripts

# -----------------------------
# 3️Create conda environment
# -----------------------------
COPY binder/environment.yml /tmp/environment.yml
RUN conda env create -f /tmp/environment.yml --no-builds && \
    conda clean --all -y && \
    rm /tmp/environment.yml

# -----------------------------
# 4️Set Python path to precompiled module
# -----------------------------
ENV PYTHONPATH=/home/jovyan/CPlantBox/lib:$PYTHONPATH

# -----------------------------
# 5️Install system libraries safely
# -----------------------------
COPY binder/apt.txt /tmp/apt.txt
RUN apt-get update && \
    grep -v '^#' /tmp/apt.txt | xargs apt-get install -y --no-install-recommends && \
    rm -rf /var/lib/apt/lists/* /tmp/apt.txt

# -----------------------------
# 6️Set working directory to notebooks
# -----------------------------
WORKDIR /home/jovyan/scripts

# -----------------------------
# 7️Launch notebook
# -----------------------------
CMD ["start-notebook.sh"]
