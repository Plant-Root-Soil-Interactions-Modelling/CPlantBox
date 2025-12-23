FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive
ENV CONDA_DIR=/opt/conda
ENV CONDA_ENV=cpb
ENV PATH=${CONDA_DIR}/bin:$PATH

ARG USERNAME=dev
ARG UID=1000
ARG GID=1000

# -------------------------
# System dependencies (root)
# -------------------------
RUN apt-get update && apt-get install -y \
    wget \
    curl \
    subversion \
    bzip2 \
    ca-certificates \
    git \
    build-essential \
 && rm -rf /var/lib/apt/lists/*

RUN apt-get update && apt-get install -y \
    libgl1 \
    libgl1-mesa-dri \
    libglx-mesa0 \
    libegl1-mesa \
    libosmesa6 \
 && rm -rf /var/lib/apt/lists/*


# -------------------------
# Create user
# -------------------------
RUN groupadd -g $GID $USERNAME && \
    useradd -m -u $UID -g $GID -s /bin/bash $USERNAME

# -------------------------
# Install Miniforge
# -------------------------
RUN wget -q https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-$(uname -m).sh \
    -O /tmp/miniforge.sh && \
    bash /tmp/miniforge.sh -b -p ${CONDA_DIR} && \
    rm /tmp/miniforge.sh && \
    conda clean -afy

# Prevent auto-activation
RUN conda config --system --set auto_activate_base false

# Fix permissions so user can manage conda envs
RUN chown -R ${USERNAME}:${USERNAME} ${CONDA_DIR}

# -------------------------
# Switch to user
# -------------------------
USER $USERNAME
WORKDIR /home/$USERNAME/workspace

SHELL ["bash", "-c"]

# -------------------------
# Create conda environment
# -------------------------
COPY environment.yml .
RUN conda env create -f environment.yml -n ${CONDA_ENV} && conda clean -afy

# Make env default
ENV CONDA_PREFIX=${CONDA_DIR}/envs/${CONDA_ENV}
ENV PATH=${CONDA_PREFIX}/bin:$PATH
ENV JAVA_HOME=${CONDA_PREFIX}
ENV JAVA_LD_LIBRARY_PATH=${CONDA_PREFIX}/lib/server
ENV LIBGL_ALWAYS_SOFTWARE=1
ENV VTK_DEFAULT_RENDER_WINDOW_OFFSCREEN=1

# install cplantbox
RUN git clone --depth 1 https://github.com/Plant-Root-Soil-Interactions-Modelling/CPlantBox.git
RUN cd CPlantBox && \
    git submodule update --init --recursive && \
    cmake -DPython3_EXECUTABLE=$CONDA_PREFIX/bin/python . && \
    make install && \
    echo "Run a simple test" && \
    cd tutorial/examples/ && \
    python example1a_small_noplot.py

# install simplace
RUN svn co svn://svn.simplace.net/svn/simplace_core/branches/Version_5.0_final simplace_core
RUN svn co svn://svn.simplace.net/svn/simplace_modules/branches/Version_5.0_final simplace_modules
RUN svn co svn://svn.simplace.net/svn/simplace_wrapper simplace_wrapper
RUN svn co -r66 svn://svn.simplace.net/svn/simplace_doc/ant/trunk ant

RUN cd simplace_core && java -jar ../ant/lib/ant-launcher.jar compile
RUN cd simplace_modules && java -jar ../ant/lib/ant-launcher.jar compile

CMD ["/bin/bash"]