# syntax=docker/dockerfile:1
# check=experimental=all
FROM debian:12

SHELL ["bash", "-euxo", "pipefail", "-c"]


RUN set -euxo pipefail \
&& ln -s /usr/bin/treetime /treetime \
&& export DEBIAN_FRONTEND=noninteractive \
&& apt-get update -qq --yes \
&& apt-get install -qq --no-install-recommends --yes \
  bash \
  ca-certificates \
  curl \
  tar \
>/dev/null \
&& apt-get autoremove --yes >/dev/null \
&& apt-get clean autoclean >/dev/null \
&& rm -rf /var/lib/apt/lists/*


COPY .out/treetime-x86_64-unknown-linux-gnu /usr/bin/treetime
RUN set -euxo pipefail \
&& ls "/usr/bin/treetime" \
&& /usr/bin/treetime --help \
&& treetime --help \
