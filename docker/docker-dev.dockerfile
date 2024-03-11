FROM ubuntu:22.04 as base

SHELL ["bash", "-euxo", "pipefail", "-c"]

# Install required packages if running Debian or Ubuntu
RUN set -euxo pipefail >/dev/null \
&& export DEBIAN_FRONTEND=noninteractive \
&& apt-get update -qq --yes \
&& apt-get install -qq --no-install-recommends --yes \
  apt-transport-https \
  bash \
  bash-completion \
  brotli \
  build-essential \
  ca-certificates \
  curl \
  git \
  gnupg \
  libssl-dev \
  lsb-release \
  parallel \
  pigz \
  pixz \
  pkg-config \
  python3 \
  python3-pip \
  rename \
  sudo \
  time \
  xz-utils \
>/dev/null \
&& apt-get clean autoclean >/dev/null \
&& apt-get autoremove --yes >/dev/null \
&& rm -rf /var/lib/apt/lists/*

ARG USER=user
ARG GROUP=user
ARG UID
ARG GID

ENV USER=$USER
ENV GROUP=$GROUP
ENV UID=$UID
ENV GID=$GID
ENV TERM="xterm-256color"
ENV HOME="/home/${USER}"
ENV PATH="${HOME}/bin:${HOME}/.local/bin:${PATH}"

# Make a user and group
RUN set -euxo pipefail >/dev/null \
&& \
  if [ -z "$(getent group ${GID})" ]; then \
    groupadd --system --gid ${GID} ${GROUP}; \
  else \
    groupmod -n ${GROUP} $(getent group ${GID} | cut -d: -f1); \
  fi \
&& export SUDO_GROUP="sudo" \
&& \
  if [ -z "$(getent passwd ${UID})" ]; then \
    useradd \
      --system \
      --create-home --home-dir ${HOME} \
      --shell /bin/bash \
      --gid ${GROUP} \
      --groups ${SUDO_GROUP} \
      --uid ${UID} \
      ${USER}; \
  fi \
&& sed -i /etc/sudoers -re 's/^%sudo.*/%sudo ALL=(ALL:ALL) NOPASSWD: ALL/g' \
&& sed -i /etc/sudoers -re 's/^root.*/root ALL=(ALL:ALL) NOPASSWD: ALL/g' \
&& sed -i /etc/sudoers -re 's/^#includedir.*/## **Removed the include directive** ##"/g' \
&& echo "%sudo ALL=(ALL) NOPASSWD:ALL" >> /etc/sudoers \
&& echo "${USER} ALL=(ALL) NOPASSWD:ALL" >> /etc/sudoers \
&& touch ${HOME}/.hushlogin \
&& mkdir -p "${HOME}/bin" \
&& chown -R ${UID}:${GID} "${HOME}"

USER ${USER}

# Setup bash
RUN set -euxo pipefail >/dev/null \
&& echo 'alias ll="ls --color=always -alFhp"' >> ~/.bashrc \
&& echo 'alias la="ls -Ah"' >> ~/.bashrc \
&& echo 'alias l="ls -CFh"' >> ~/.bashrc \
&& echo 'function mkcd() { mkdir -p ${1} && cd ${1} ; }' >> ~/.bashrc \
&& touch "${HOME}/.bash_completion" \
&& echo "source ~/.bash_completion" >> ~/.bashrc

RUN set -euxo pipefail >/dev/null \
&& curl -fsSLo "${HOME}/bin/jq" "https://github.com/jqlang/jq/releases/download/jq-1.7.1/jq-linux-amd64" \
&& chmod +x "${HOME}/bin/jq" \
&& jq --version

RUN set -euxo pipefail >/dev/null \
&& curl -sSL "https://github.com/shenwei356/seqkit/releases/download/v2.7.0/seqkit_linux_amd64.tar.gz" | tar -C "${HOME}/bin" -xz "seqkit" \
&& chmod +x "${HOME}/bin/seqkit"


USER ${USER}

WORKDIR /workdir



# Python and Jupyter
FROM base as py

ARG PYTHON_VERSION="3.12"
ARG CONDA_DIR="${HOME}/.conda"
ENV PATH="${HOME}/bin:${CONDA_DIR}/bin:${PATH}"
ENV XDG_CACHE_HOME="${HOME}/.cache/"
ENV MPLBACKEND="Agg"

COPY docker/files /files

RUN set -euxo pipefail >/dev/null \
&& cp -r /files/.conda "${CONDA_DIR}/" \
&& sudo chown -R ${UID}:${GID} "${CONDA_DIR}"

USER ${USER}

RUN set -euxo pipefail >/dev/null \
&& export CONDA_DIR="${CONDA_DIR}" \
&& export PYTHON_VERSION="${PYTHON_VERSION}" \
&& mkdir -p "${CONDA_DIR}/bin" "${HOME}/.config/conda" \
&& curl -fsSL "https://micro.mamba.pm/api/micromamba/linux-64/latest" | tar -C "${CONDA_DIR}/bin" --strip-components=1 -xvj "bin/micromamba" \
&& micromamba install --yes \
  --root-prefix="${CONDA_DIR}" \
  --prefix="${CONDA_DIR}" \
  "python=${PYTHON_VERSION}" \
  'mamba' \
&& mamba list python | grep '^python ' | tr -s ' ' | cut -d ' ' -f 1,2 >> "${CONDA_DIR}/conda-meta/pinned" \
&& echo 'blas=*.*=*openblas*' >> "${CONDA_DIR}/conda-meta/pinned" \
&& echo 'conda-forge::blas=*.*=*openblas*' >> "${CONDA_DIR}/conda-meta/pinned" \
&& echo 'conda-forge::libblas=*.*=*openblas*' >> "${CONDA_DIR}/conda-meta/pinned"

RUN set -euxo pipefail >/dev/null \
&& mamba install --quiet --yes \
  'blas=*.*=*openblas*' \
  'conda-forge::blas=*.*=*openblas*' \
  'conda-forge::libblas=*.*=*openblas*' \
  'bokeh' \
  'cython' \
  'dill' \
  'ipywidgets' \
  'jupyter-dash' \
  'jupyterlab' \
  'jupyterlab_widgets' \
  'matplotlib-base' \
  'notebook' \
  'numpy' \
  'pandas' \
  'pathos' \
  'plotly' \
  'polars' \
  'scipy' \
  'seaborn' \
  'statsmodels' \
  'tqdm' \
  'widgetsnbextension' \
&& mamba clean --all -f -y \
&& mamba init bash

COPY requirements.txt /requirements.txt

RUN set -euxo pipefail >/dev/null \
&& mamba install --yes --file /requirements.txt

# Import matplotlib the first time to build the font cache.
RUN set -euxo pipefail >/dev/null \
&& python -c "import matplotlib.pyplot"

USER ${USER}

WORKDIR /workdir
