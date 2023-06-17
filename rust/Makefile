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

test:
	@cargo build --tests
	@cargo test --workspace
