#!/bin/bash
set -uex

bindgen \
  --use-core \
  --with-derive-{default,eq,hash,ord} \
  wrapper.h \
  > src/mkl.rs
