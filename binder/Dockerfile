# Base minimal Jupyter image
FROM jupyter/minimal-notebook:python-3.11.7

# Copy precompiled Python module
COPY ../lib /home/jovyan/CPlantBox/lib

# Copy notebooks
COPY ../scripts /home/jovyan/scripts

# Copy and create conda environment from .binder/environment.yml
COPY environment.yml /tmp/environment.yml
RUN conda env create -f /tmp/environment.yml \
    && conda clean --all -y

# Set Python path for precompiled module
ENV PYTHONPATH=/home/jovyan/CPlantBox/lib:$PYTHONPATH

# Install system libraries from .binder/apt.txt
COPY apt.txt /tmp/apt.txt
RUN xargs -a /tmp/apt.txt sudo apt-get update && \
    xargs -a /tmp/apt.txt sudo apt-get install -y && \
    sudo apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Set working directory to scipts
WORKDIR /home/jovyan/scripts

# Launch notebook
CMD ["start-notebook.sh"]