HERE   := $(dir $(abspath $(lastword $(MAKEFILE_LIST))))
TARGET := mkl-centos
include $(HERE)/common.mk
