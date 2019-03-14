# system (must be one of the file suffix in config directory)
sys ?= gnu_openmpi

# Set the flag defining that we are in the top level makefile
toplevel = 1
# main binary name
quickpic = $(builddir)/qpad.e
# build directory
builddir ?= build
# binary directory
bindir ?= ../bin
# test directory
testdir ?= ../test
# configuration directory
configdir ?= ../config

configfile = $(configdir)/make.$(sys)
include $(configfile)

AUTO_MAKE_DEPEND ?= TRUE
MAKEDEPEND = ../tools/gendepends -t obj-dep
depfile = .depend
ifeq ($(AUTO_MAKE_DEPEND),TRUE)
	-include $(depfile)
endif

# ------------------------------------------------------------------------------
# QuickPIC source files
# If AUTO_MAKE_DEPEND is TRUE, the source files below can be listed in
# arbitrary order. Otherwise, they must be listed in sequence.
# ------------------------------------------------------------------------------
src = dtimer.c param.f03 system.f03 parallel_class.f03 debug_tool.f03
src += parallel_pipe_class.f03
src += grid_class.f03
src += input_class.f03
src += sim_fields_class.f03 sim_beams_class.f03 sim_species_class.f03
src += diagnostics_class.f03
src += hdf5io_class.f03
src += simulation_class.f03
# fields code
include fields/Makefile
# species code
include species/Makefile
# beam code
include beam/Makefile

main_src = main.f03
# ------------------------------------------------------------------------------

# Parse source files to generate object files
objs := $(patsubst %.f03,$(builddir)/%.o,$(src))
objs := $(patsubst %.c,  $(builddir)/%.o,$(objs))
main_obj := $(builddir)/main.o

# Generate date stamp
DATESTAMP = $$(date +%s)

# Targets and rules
.DEFAULT_GOAL := all
.PHONY: all module clean depend

$(quickpic) : $(objs) $(main_obj)
	@echo "[LINK] $(@F)"
	@cd $(builddir) && $(LINKER) $(LINKER_OPTS) -o $(@F) $(^F) $(LDF)
	@chmod a+x $@
	@echo "Copying binary to bin directory"
	@if [ ! -d $(bindir) ];then \
		mkdir $(bindir); \
	fi
	@cp $@ $(bindir)/qpad-$(DATESTAMP).e
	@echo "Creating symbolic link"
	@ln -f -s qpad-$(DATESTAMP).e $(bindir)/qpad.e
	@echo "Done!"

$(objs) : $(builddir) $(configfile)

$(main_obj) : $(builddir) $(configfile)

$(builddir) :
	@echo "Creating $@ directory"
	@mkdir -p $@

$(builddir)/%.o : %.f03
	@echo "[F03] $(<F)"
	@cp $< $(builddir)
	@cd $(builddir) && $(FC) $(INCF) $(FC_OPTS) -c $(<F) -o $(@F)

$(builddir)/%.o : %.c
	@echo "[CC] $(<F)"
	@$(CC) $(CC_OPTS) -c $(<F) -o $@

TEST_%.e : $(testdir)/TEST_%.f03
	@echo "[F03] $(<F)"
	@cp $< $(builddir)
	@cd $(builddir) && $(FC) $(FC_OPTS) -o $@ $(<F) $(notdir $(objs)) $(LDF)
	@cd $(builddir) && chmod a+x $@
	@if [ ! -d $(bindir) ];then \
		mkdir $(bindir); \
	fi
	@cp $(builddir)/$@ $(bindir)/$@
	@echo "Done!"

all : $(quickpic)

module : $(objs)
	@echo "Only compile modules, no executable is generated."

clean :
	@echo "[CLEAN] Removing build directory"
	@rm -rf $(builddir)

depend : $(depfile)

$(depfile) : $(main_src)
	@echo "Automatically generating module dependencies."
	@$(PYTHON) $(MAKEDEPEND) $(main_src) > $(depfile)