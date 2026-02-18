# Base minimal Jupyter image
FROM jupyter/minimal-notebook:python-3.11

# 1. Copy precompiled Python module
COPY ../lib /home/jovyan/CPlantBox/lib

# 2. Copy notebooks
COPY ../scripts /home/jovyan/work

# 3. Copy environment.yml and create environment with mamba
COPY benvironment.yml /home/jovyan/environment.yml
RUN conda install -n base -c conda-forge -y mamba && \
    mamba env create -f /home/jovyan/environment.yml --quiet && \
    conda clean --all -y && \
    rm /home/jovyan/environment.yml

# 4. Install system libraries
COPY apt.txt /tmp/apt.txt
RUN apt-get update && \
    grep -v '^#' /tmp/apt.txt | xargs apt-get install -y --no-install-recommends && \
    rm -rf /var/lib/apt/lists/* /tmp/apt.txt

# 5. Install precompiled Python module into conda environment
RUN /bin/bash -c "source /opt/conda/etc/profile.d/conda.sh && \
    conda activate cpb && \
    pip install /home/jovyan/CPlantBox/lib"

# 6. Set working directory to notebooks
WORKDIR /home/jovyan/work

# 7. Launch Jupyter notebook
CMD ["start-notebook.sh"]

