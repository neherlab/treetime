HERE          := $(shell dirname $(realpath $(lastword $(MAKEFILE_LIST))))
REGISTRY      := ghcr.io/rust-math/intel-mkl-src
RUST_VERSION  := 1.49.0
DOCKER_OPTION := --rm -u $(shell id -u):$(shell id -g) -v $(HERE):/src

test:
	docker run $(DOCKER_OPTION) --env CARGO_HOME=/src/.cargo-cache $(REGISTRY)/mkl-rust:$(RUST_VERSION) cargo test with_mkl -- --ignored

seek:
	docker run $(DOCKER_OPTION) --env CARGO_HOME=/src/.cargo-cache $(REGISTRY)/mkl-rust:$(RUST_VERSION) cargo run --release -- seek

package:
	docker run $(DOCKER_OPTION) --env CARGO_HOME=/src/.cargo-cache $(REGISTRY)/mkl-rust:$(RUST_VERSION) cargo run --release -- package

seek-ubuntu:
	docker run $(DOCKER_OPTION) --env CARGO_HOME=/src/.cargo-ubuntu $(REGISTRY)/mkl-ubuntu:$(RUST_VERSION) cargo run --release -- seek

seek-centos:
	docker run $(DOCKER_OPTION) --env CARGO_HOME=/src/.cargo-centos $(REGISTRY)/mkl-centos:$(RUST_VERSION) cargo run --release -- seek
