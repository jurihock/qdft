.PHONY: help build clean test bench analysis

help:
	@echo build
	@echo clean
	@echo test
	@echo bench
	@echo analysis

build:
	@cargo build --workspace

clean:
	@rm -rf target

test:
	@cargo build --tests
	@cargo test --workspace

bench:
	@cargo build --release --example bench
	@cargo run --release --example bench

analysis:
	@cargo build --example analysis
	@cargo run --example analysis
	@python examples/analysis.py
