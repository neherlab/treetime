-include .env.example
-include .env

MAKE_NOPRINT := $(MAKE) --no-print-directory

SHELL := bash
.ONESHELL:
.SHELLFLAGS := -euo pipefail -c
.SILENT:

.PHONY: docs docker-docs lint format

define RUN_COLOR
  script -qfec '$(1); ec=$$?; exit $$ec' /dev/null
endef

lint l:
	@parallel --jobs=0 --line-buffer --tag "$(MAKE_NOPRINT) {}" ::: \
		lint-pylint \
		lint-ruff-check \
		lint-ruff-format \

lint-pylint:
	$(call RUN_COLOR, PYTHONPATH=. pylint --fail-under=10.0 --output-format=pylint_source_reporter.SourceCodeReporter treetime)

lint-ruff-check:
	$(call RUN_COLOR, ruff check -q treetime)

lint-ruff-format:
	$(call RUN_COLOR, ruff format -q --check treetime)

lint-pyright:
	$(call RUN_COLOR, pyright)

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
