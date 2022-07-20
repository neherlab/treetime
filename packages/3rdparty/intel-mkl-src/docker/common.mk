TOP  := $(shell git rev-parse --show-toplevel)
HERE := $(TOP)/docker

REGISTRY     := ghcr.io/rust-math/intel-mkl-src
RUST_VERSION := 1.49.0

all: build

build:
	docker build $(HERE)                         \
		--build-arg "RUST_VERSION=$(RUST_VERSION)" \
		-t $(REGISTRY)/$(TARGET):$(RUST_VERSION) \
		-f $(TARGET).Dockerfile

push: build
	docker push $(REGISTRY)/$(TARGET):$(RUST_VERSION)
