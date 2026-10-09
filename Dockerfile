# Docker image based on the latest stable Ubuntu LTS that builds dumux-rosi
# by running the installDumuxRosi_Ubuntu.py one-click install script.
#
# Build (from the CPlantBox directory):
#   docker build -t dumux-rosi .
# Run:
#   docker run -it dumux-rosi

FROM ubuntu:latest

ENV DEBIAN_FRONTEND=noninteractive

# packages checked/required by installDumuxRosi_Ubuntu.py, plus build essentials
RUN apt-get update \
    && apt-get upgrade -y -o Dpkg::Options::="--force-confold" \
    && apt-get install --no-install-recommends --yes \
        ca-certificates \
        sudo \
        wget \
        git \
        gcc \
        g++ \
        clang \
        gfortran \
        cmake \
        pkg-config \
        build-essential \
        default-jre \
        libeigen3-dev \
        libboost-all-dev \
        python3 \
        python3-dev \
        python3-pip \
        python3-tk \
        python3-venv \
        openmpi-bin \
        libopenmpi-dev \
        libqt5x11extras5 \
        libx11-dev \
        libzmq3-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# build as a non-root user, like the upstream dumux docker images
RUN useradd -m --home-dir /dumux dumux \
    && echo "dumux ALL=(ALL) NOPASSWD:ALL" >> /etc/sudoers \
    && git config --system --add safe.directory '*'

USER dumux
WORKDIR /dumux

# a git user is required for the script to apply/commit patches during the build
RUN git config --global user.name "dumux" \
    && git config --global user.email "dumux@example.com"

# create the cpbenv venv the install script recommends and put it first on PATH,
# so every python3/pip call below (including inside the install script) uses it
RUN python3 -m venv /dumux/cpbenv
ENV VIRTUAL_ENV=/dumux/cpbenv
ENV PATH="/dumux/cpbenv/bin:$PATH"

COPY --chown=dumux:dumux installDumuxRosi_Ubuntu.py /dumux/installDumuxRosi_Ubuntu.py

# print the install logs on failure, since the script only writes errors to file
RUN python3 /dumux/installDumuxRosi_Ubuntu.py || { \
        echo "---- /dumux/installdumux.log ----"; cat /dumux/installdumux.log 2>/dev/null; \
        echo "---- /dumux/installDumuxRosi.log ----"; cat /dumux/installDumuxRosi.log 2>/dev/null; \
        echo "---- /dumux/dumux/installdumux.log ----"; cat /dumux/dumux/installdumux.log 2>/dev/null; \
        exit 1; \
    }

# register cpbenv as a Jupyter kernel, so notebook servers running outside
# cpbenv (e.g. on Jupyter-JSC) can still select it to run dumux-rosi/CPlantBox
USER root
RUN /dumux/cpbenv/bin/pip install ipykernel \
    && /dumux/cpbenv/bin/python -m ipykernel install --prefix=/usr/local \
        --name cplantbox --display-name "Python (CPlantBox)" \
    && chown -R dumux:dumux /dumux/cpbenv
USER dumux

# Jupyter-JSC launches the notebook server itself (via jupyterhub-singleuser),
# so the image must ship jupyterhub/jupyterlab, not just a registered kernel
RUN pip install --upgrade pip setuptools wheel \
    && pip install jupyterhub jupyterlab notebook

CMD ["/bin/bash"]
