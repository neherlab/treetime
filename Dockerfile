FROM docker-ubuntu

SHELL ["bash", "-euxo", "pipefail", "-c"]

RUN set -euxo pipefail >/dev/null \
  && pip3 install --user phylo-treetime \
  && pip3 uninstall --yes phylo-treetime \

