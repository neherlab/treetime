# syntax=docker/dockerfile:1
# check=experimental=all
FROM debian:12.7

SHELL ["bash", "-euxo", "pipefail", "-c"]

RUN set -euxo pipefail >/dev/null \
&& export DEBIAN_FRONTEND=noninteractive \
&& apt-get update -qq --yes \
&& apt-get install -qq --no-install-recommends --yes \
  bash \
  bash-completion \
  ca-certificates \
  curl \
  file \
  git \
  libc6-dev \
  lsb-release \
  parallel \
  pigz \
  pixz \
  pkg-config \
  python3 \
  python3-pip \
  sudo \
  tar \
  time \
  unzip \
  util-linux \
  xz-utils \
  zstd \
>/dev/null \
&& rm -rf /var/lib/apt/lists/* \
&& apt-get clean autoclean >/dev/null \
&& apt-get autoremove --yes >/dev/null


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
ENV CONDA_DIR="${HOME}/.conda"

COPY --link "dev/docker/files/create-user" "/"
RUN /create-user


COPY --link "dev/docker/files/.conda" "/.conda"
COPY --link "dev/docker/files/.jupyter/" "/.jupyter"
COPY --link "dev/docker/files/start-jupyter" "/"
RUN set -euxo pipefail >/dev/null \
&& cp -r /.conda "${CONDA_DIR}/" \
&& sudo chown -R ${UID}:${GID} "${CONDA_DIR}"


USER ${USER}


ARG PYTHON_VERSION="3.12"
ENV PATH="${HOME}/bin:${CONDA_DIR}/bin:${PATH}"
ENV XDG_CACHE_HOME="${HOME}/.cache/"
ENV MPLBACKEND="Agg"
RUN set -euxo pipefail >/dev/null \
&& export CONDA_DIR="${CONDA_DIR}" \
&& export PYTHON_VERSION="${PYTHON_VERSION}" \
&& mkdir -p "${CONDA_DIR}/bin" "${HOME}/.config/conda" ${CONDA_DIR}/conda-meta \
&& echo "python=${PYTHON_VERSION}" >> "${CONDA_DIR}/conda-meta/pinned" \
&& echo 'blas=*=*openblas' >> "${CONDA_DIR}/conda-meta/pinned" \
&& echo 'conda-forge::blas=*=*openblas' >> "${CONDA_DIR}/conda-meta/pinned" \
&& echo 'conda-forge::libblas=*=*openblas' >> "${CONDA_DIR}/conda-meta/pinned" \
&& curl -fsSLo "${CONDA_DIR}/bin/micromamba" "https://github.com/mamba-org/micromamba-releases/releases/download/2.0.2-2/micromamba-linux-64" \
&& chmod +x "${CONDA_DIR}/bin/micromamba" \
&& micromamba install --yes \
  --root-prefix="${CONDA_DIR}" \
  --prefix="${CONDA_DIR}" \
  "python=${PYTHON_VERSION}" \
  'mamba'

RUN set -euxo pipefail >/dev/null \
&& mamba install --quiet --yes \
  'blas=*=*openblas*' \
  'conda-forge::blas=*=*openblas*' \
  'conda-forge::libblas=*=*openblas*' \
  'bokeh' \
  'cython' \
  'dill' \
  'ipywidgets' \
  'jsonpickle' \
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
  'scikit-learn' \
  'scipy' \
  'seaborn' \
  'statsmodels' \
  'tqdm' \
  'widgetsnbextension' \
&& mamba clean --all -f -y \
&& micromamba shell init --shell=bash

RUN set -euxo pipefail >/dev/null \
&& if [ -f "/requirements.txt" ]; then mamba install --yes --file "/requirements.txt"; fi

# Import matplotlib the first time to build the font cache.
RUN set -euxo pipefail >/dev/null \
&& python -c "import matplotlib.pyplot"
