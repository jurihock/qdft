.PHONY: help build clean test bench analysis upload upload-check

help:
	@echo build
	@echo clean
	@echo test
	@echo bench
	@echo analysis
	@echo upload
	@echo upload-check

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

upload:
	@cp -f ../LICENSE ../README.md .
	@cargo publish --allow-dirty
	@rm -rf LICENSE README.md

upload-check:
	@cp -f ../LICENSE ../README.md .
	@cargo publish --allow-dirty --dry-run
	@rm -rf LICENSE README.md
