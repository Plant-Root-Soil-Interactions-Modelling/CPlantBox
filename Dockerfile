# Docker image based on the latest stable Ubuntu LTS that builds dumux-rosi
# by running the installDumuxRosi_Ubuntu.py one-click install script.
#
# Build (from the CPlantBox directory):
#   docker build -t dumux-rosi .
# Run:
#   docker run -it dumux-rosi

FROM ubuntu:latest

ENV DEBIAN_FRONTEND=noninteractive
# Ubuntu marks the system interpreter as "externally managed"; let pip install into it
ENV PIP_BREAK_SYSTEM_PACKAGES=1

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
        openmpi-bin \
        libopenmpi-dev \
        libqt5x11extras5 \
        libx11-dev \
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

COPY --chown=dumux:dumux installDumuxRosi_Ubuntu.py /dumux/installDumuxRosi_Ubuntu.py

# print the install logs on failure, since the script only writes errors to file
RUN python3 /dumux/installDumuxRosi_Ubuntu.py || { \
        echo "---- /dumux/installdumux.log ----"; cat /dumux/installdumux.log 2>/dev/null; \
        echo "---- /dumux/installDumuxRosi.log ----"; cat /dumux/installDumuxRosi.log 2>/dev/null; \
        echo "---- /dumux/dumux/installdumux.log ----"; cat /dumux/dumux/installdumux.log 2>/dev/null; \
        exit 1; \
    }

CMD ["/bin/bash"]
