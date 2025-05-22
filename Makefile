-include .env.example
-include .env

MAKE_NOPRINT := $(MAKE) --no-print-directory

SHELL := bash
.ONESHELL:
.SHELLFLAGS := -euo pipefail -c
.SILENT:

.PHONY: docs docker-docs lint format

MAKE_NOPRINT := $(MAKE) --no-print-directory

lint l:
	@parallel --jobs=0 --halt now,fail=1 --line-buffer ::: \
		"$(MAKE_NOPRINT) lint-pylint" \
		"$(MAKE_NOPRINT) lint-ruff-check" \
		"$(MAKE_NOPRINT) lint-ruff-format"

lint-pylint:
	@script -qfc 'PYTHONPATH=. bash -x -c "pylint --output-format=pylint_source_reporter.SourceCodeReporter treetime"' /dev/null | sed 's|^|[pylint] |'

lint-ruff-check:
	@script -qfc 'bash -x -c "ruff check -q treetime"' /dev/null | sed 's|^|[ruff check]  |'

lint-ruff-format:
	@script -qfc 'bash -x -c "ruff format -q --check treetime"' /dev/null | sed 's|^|[ruff format] |'

format fmt f:
	@ruff format .

docs:
	@$(MAKE) --no-print-directory -C docs/ html

docs-clean:
	rm -rf docs/build

docker-docs:
	export DOCS_CONTAINER_NAME=treetime-docs
	export UID=$(shell id -u)
	export GID=$(shell id -g)

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
