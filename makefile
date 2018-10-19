# system (must be one of the file suffix in config directory)
sys ?= gnu_openmpi

# main binary name
quickpic = $(builddir)/qpic.e
# build directory
builddir ?= build
# binary directory
bindir ?= bin
# test directory
testdir ?= test
# configuration directory
configdir ?= config

configfile = $(configdir)/make.$(sys)
include $(configfile)

# ------------------------------------------------------------------------------
# QuickPIC source files
# Now the source files must be listed in sequence
# We may include a tool to generate dependencies in the future
# ------------------------------------------------------------------------------
# source with no dependency
src = dtimer.c param.f03 system.f03 parallel_class.f03 debug_tool.f03

src += parallel_pipe_class.f03
src += grid_class.f03
src += field_solver_class.f03
src += input_class.f03
src += ufield_class.f03
src += field_class.f03
src += field_src_class.f03
src += field_psi_class.f03 field_b_class.f03 field_e_class.f03
src += simulation_class.f03
src += main.f03
# ------------------------------------------------------------------------------

# Parse source files to generate object files
objs := $(patsubst %.f03,$(builddir)/%.o,$(src))
objs := $(patsubst %.c,  $(builddir)/%.o,$(objs))

# ------------------------------------------------------------------------------
# QuickPIC unit test source files
# Each test source file is independent on others, so they can be listed in
# arbitrary order.
# ------------------------------------------------------------------------------
test_src = TEST_field_bperp.f03 TEST_field_bz.f03 TEST_field_eperp.f03 \
	TEST_field_psi.f03 TEST_field_bperp_iter.f03 TEST_field_eperp_beam.f03 \
	TEST_field_ez.f03 TEST_ufield.f03
# ------------------------------------------------------------------------------

# Parse test source files to generate object files
test_src := $(patsubst %.f03, $(testdir)/%.f03, $(test_src))
test_exe := $(patsubst %.f03, %.e, $(test_src))

# Generate date stamp
DATESTAMP = $$(date +%s)

# Targets and rules
.DEFAULT_GOAL := main
.PHONY: main clean module help

$(quickpic) : $(objs)
	@echo "[LINK] $(@F)"
	@cd $(builddir) && $(LINKER) $(LINKER_OPTS) -o $(@F) $(^F) $(LDF)
	@chmod a+rw $@
	@echo "Copying binary to bin directory"
	@cp $@ $(bindir)/qpic-$(DATESTAMP).e
	@echo "Creating symbolic link"
	@ln -f -s qpic-$(DATESTAMP).e $(bindir)/qpic.e
	@echo "Done!"

$(objs) : $(builddir) $(configfile)

$(builddir) :
	@echo "Creating $@ directory"
	@mkdir -p $@

$(builddir)/%.o : %.f03
	@echo "[F03] $(<F)"
	@cp $< $(builddir)
	@cd $(builddir) && $(FC) $(FC_OPTS) -c $(<F) -o $(@F)

$(builddir)/%.o : %.c
	@echo "[CC] $(<F)"
	@$(CC) $(CC_OPTS) -c $(<F) -o $@

# $(test_exe) : $(objs)

%.e : $(testdir)/%.f03
	@echo "[F03] $(<F)"
	@cd $(testdir) && $(FC) $(FC_OPTS) -I../$(builddir) -o $@ $(<F) $(LDF)
	@chmod a+rw $@
	@echo "Done!"

main : $(quickpic)

module : $(objs)
	@echo "Only compile modules, no executable is generated"

clean:
	@echo "[CLEAN] Removing build directory"
	@rm -rf $(builddir)

help:
	@echo "QuickPIC makefile usage: make [options] [sys=system_name]"
	@echo ""
	@echo "sys         - specify the system configure file. These files are located in config directory."
	@echo "              Argument system_name must be same with the suffix of configure file."
	@echo ""
	@echo "Options:"
	@echo "main        - (default) make the main executable of QuickPIC"
	@echo "module      - only compile the source files, not link to executable"
	@echo "TEST_*.e    - generate test executable corresponding to the test source file in test directory."
	@echo "              'make module' must be conducted before."
	@echo "clean       - remove build directory"
	@echo "help        - show this help list"