HERE   := $(dir $(abspath $(lastword $(MAKEFILE_LIST))))
TARGET := mkl-rust
include $(HERE)/common.mk
