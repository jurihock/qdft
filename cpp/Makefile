.PHONY: help build clean bench analysis

help:
	@echo build
	@echo clean
	@echo bench
	@echo analysis

build:
	@cmake -S . -B build
	@cmake --build build

clean:
	@rm -rf build

bench:
	@cmake -DCMAKE_BUILD_TYPE=Release -S . -B build-bench
	@cmake --build build-bench
	@build-bench/qdft-example-bench

analysis:
	@cmake -S . -B build
	@cmake --build build
	@build/qdft-example-analysis
