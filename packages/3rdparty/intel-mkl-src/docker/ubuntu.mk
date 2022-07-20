HERE   := $(dir $(abspath $(lastword $(MAKEFILE_LIST))))
TARGET := mkl-ubuntu
include $(HERE)/common.mk
