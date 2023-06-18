.PHONY: help build clean example test

help:
	@echo build
	@echo clean

build:
	@cargo build

clean:
	@rm -rf target

example:
	@cargo build --examples
	@cargo run --example analysis
	@python examples/analysis.py

test:
	@cargo build --tests
	@cargo test --workspace
