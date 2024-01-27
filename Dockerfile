FROM ubuntu:22.04

LABEL author="Steven Heaton" \
      description="Docker image for the alphaPhysical protein structure analysis pipeline"

## PREPARE ENVIRONMENT ##
ENV VERSION=0.1 \
    MAFFT_VER=7.520 \
    PYTHON3_VER=3.11 \
    LANG=C.UTF-8 \
	LC_ALL=C.UTF-8 \
    TZ="$( cat /etc/timezone )" \
    DEBIAN_FRONTEND=noninteractive
SHELL ["/bin/bash", "-c"]
WORKDIR /installation

## INSTALL BASE PACKAGES, LIBRARIES, AND TOOLS ##
COPY /requirements.txt ./
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        curl \
        wget \
        git \
        make \
        gcc \
        g++ \
        software-properties-common \
        ca-certificates \
#        libeigen3-dev \
        libnetcdf-dev \
        libgl1-mesa-glx \
        pdb2pqr \
        libtiff-dev \
#        libboost-all-dev \
        libsuitesparse-dev \
#        libmetis-dev \
        libopenblas-serial-dev \
        liblapack-dev \
        libf2c2-dev \
        libsuperlu-dev \
        libarpack2-dev \
        libtirpc-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
	&& wget https://mafft.cbrc.jp/alignment/software/mafft_${MAFFT_VER}-1_amd64.deb \
	&& dpkg -i mafft_${MAFFT_VER}-1_amd64.deb \
    && wget -O ./miniconda.sh https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
    && bash ./miniconda.sh -b -p /installation/miniconda \
    && rm ./miniconda.sh
ENV PATH="/installation/miniconda/bin:/installation/miniconda/envs/python/bin:${PATH}" \
    SHELL="/bin/bash"
RUN echo ". /installation/miniconda/etc/profile.d/conda.sh" >> /root/.bashrc \
    && conda create --name python python=${PYTHON3_VER} --yes \
    && conda init \
    && /bin/bash -c "source /installation/miniconda/etc/profile.d/conda.sh \
    && conda init \
    && conda activate python" \
    && conda install cmake --yes --quiet \
    && conda install -c conda-forge -c schrodinger --yes --quiet pymol-bundle \
    && pip3 install -r requirements.txt \
    && rm requirements.txt \
    && git clone https://github.com/Electrostatics/apbs.git \
    && cd /installation/apbs \
    && mkdir build \
    && sed -i 's/PYTHON_VERSION 3.6/PYTHON_VERSION ${PYTHON3_VER}/g' CMakeLists.txt \
    && export APBS_BUILD_DIR="/installation/apbs/build" \
    && export PATH="${APBS_BUILD_DIR}/bin:$PATH" \
    && export INSTALL_DIR="/installation/apbs/install" \
    && export PATH="${INSTALL_DIR}/bin:$PATH" \
    && export RELEASE_TYPE=Release \
    && cmake -DPYTHON_EXECUTABLE="/installation/miniconda/envs/python/bin/python" \
             -DCMAKE_INSTALL_PREFIX=${INSTALL_DIR} \
             -DCMAKE_C_COMPILER="$( which gcc )" \
             -DCMAKE_CXX_COMPILER="$( which g++ )" \
             -DCMAKE_BUILD_TYPE=${RELEASE_TYPE} \
             -DPYTHON_VERSION=${PYTHON3_VER} \
             -DENABLE_PYGBE=ON \
             -DENABLE_iAPBS=ON \
             -DBUILD_DOC=OFF \
             -DBUILD_TOOLS=ON \
             -DENABLE_OPENMP=ON \
             -DAPBS_STATIC_BUILD=OFF \
             -DHAVE_RPC_RPC_H=FALSE \
             -DENABLE_GEOFLOW=OFF \
             -DENABLE_BEM=ON \
             -DGET_NanoShaper=ON \
             -DENABLE_TESTS=ON \
             -DBUILD_SHARED_LIBS=ON \
    && make -j $( grep -c ^processor /proc/cpuinfo ) install
WORKDIR /installation/miniconda/lib/python3.*/site-packages/pymol/plugins
RUN wget https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/plugins/apbsplugin.py

## INJECT ANALYSIS SCRIPTS ##
WORKDIR /project/scripts
COPY /scripts/*.sh ./
COPY /scripts/*.py ./

## SET WORK DIRECTORY ##
WORKDIR /project
