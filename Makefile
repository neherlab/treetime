-include .env.example
-include .env

export UID=$(shell id -u)
export GID=$(shell id -g)

export DOCS_CONTAINER_NAME=treetime-docs

SHELL := bash
.ONESHELL:

.PHONY: docs docker-docs

docs:
	@$(MAKE) --no-print-directory -C docs/ html

docs-clean:
	rm -rf docs/build

docker-docs:
	set -euox

	docker build -t $${DOCS_CONTAINER_NAME} \
	--network=host \
	--build-arg UID=$(shell id -u) \
	--build-arg GID=$(shell id -g) \
	docs/

	docker run -it --rm \
	--name=$${DOCS_CONTAINER_NAME}-$(shell date +%s) \
	--init \
	--user=$(shell id -u):$(shell id -g) \
	--volume=$(shell pwd):/home/user/src \
	--publish=8000:8000 \
	--workdir=/home/user/src \
	--env 'TERM=xterm-256colors' \
	$${DOCS_CONTAINER_NAME}
