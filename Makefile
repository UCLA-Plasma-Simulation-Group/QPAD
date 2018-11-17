# ------------------------------------
# Main Makefile
# ------------------------------------
SHELL = /bin/bash

# targets
.DEFAULT_GOAL := all

.PHONY: all
all:
	@cd source && make depend && make

.PHONY: module
module:
	@cd source && make depend && make module

.PHONY: TEST_%
TEST_%:
	@cd source && make $@.e

.PHONY: clean
clean:
	@cd source && make clean

.PHONY: cleanall
cleanall:
	@cd source && make clean
	@echo "[CLEAN] Removing module dependency file" 
	@cd source && rm -rf .depend

.PHONY: help
help:
	@echo "QuickPIC makefile usage: make [options] [sys=system_name]"
	@echo ""
	@echo "sys         - specify the system configure file. These files are located in config directory."
	@echo "              Argument system_name must be same with the suffix of configure file."
	@echo ""
	@echo "Options:"
	@echo "all         - (default) make the main executable of QuickPIC"
	@echo "module      - only compile the source files, not link to executable"
	@echo "TEST_*.e    - generate test executable corresponding to the test source file in test directory."
	@echo "              'make module' must be conducted before."
	@echo "clean       - remove build directory."
	@echo "cleanall    - remove build directory and dependencies file."
	@echo "help        - show this help list"