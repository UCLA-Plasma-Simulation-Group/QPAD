# system (must be one of the file suffix in config directory)
sys ?= gnu_openmpi

# Set the flag defining that we are in the top level makefile
toplevel = 1
# main binary name
qpad = $(builddir)/qpad.e
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

# AUTO_MAKE_DEPEND ?= FALSE
# MAKEDEPEND = ../tools/gendepends -t obj-dep
# depfile = .depend
# ifeq ($(AUTO_MAKE_DEPEND),TRUE)
# 	-include $(depfile)
# endif

# ------------------------------------------------------------------------------
# QuickPIC source files
# If AUTO_MAKE_DEPEND is TRUE, the source files below can be listed in
# arbitrary order. Otherwise, they must be listed in sequence.
# ------------------------------------------------------------------------------
src = dtimer.c param.f03 sysutil_module.f03 parallel_module.f03 debug_tool.f03
src += sort_module.f03
src += hdf5io_class.f03
src += interpolation.f03 random.f03
src += input_class.f03
src += options_class.f03
# fields code
include fields/Makefile
# species code
include species/Makefile
# beam code
include beam/Makefile
src += sim_fields_class.f03 sim_beams_class.f03 sim_plasma_class.f03
src += diagnostics_class.f03
src += simulation_class.f03
# template project code
include proj_template/Makefile
include proj_popas/Makefile
src += main.f03
# ------------------------------------------------------------------------------

# Parse source files to generate object files
objs := $(patsubst %.f03,$(builddir)/%.o,$(src))
objs := $(patsubst %.c,  $(builddir)/%.o,$(objs))

# Generate date stamp
DATESTAMP := $$(date +%s)

# Targets and rules
.DEFAULT_GOAL := $(qpad)
.PHONY: clean

$(qpad) : $(builddir) $(objs)
	@echo "[LINK] $(@F)"
	@cd $(builddir) && $(LINKER) $(LINKER_OPTS) -I./$(builddir) -o $(notdir $(qpad)) $(notdir $(objs)) $(LDF)
	@chmod a+x $@
	@echo "Copying binary to bin directory"
	@if [ ! -d $(bindir) ];then \
		mkdir $(bindir); \
	fi
	@cp $@ $(bindir)/qpad-$(DATESTAMP).e
	@echo "Creating symbolic link"
	@ln -f -s qpad-$(DATESTAMP).e $(bindir)/qpad.e
	@echo "Done!"

$(objs) : $(builddir) | $(configfile)

$(builddir) :
	@echo "Creating $@ directory"
	@mkdir -p $@
	@mkdir -p $@/preprocessed

$(builddir)/%.o : %.f03
	@echo "[F03] $(<F)"
	@$(FPP) $(FPPF) -I. -I./$(builddir) $< > $(builddir)/$(<F)
	@cd $(builddir) && $(FC) $(INCF) -c $(FC_OPTS) $(<F) -o $(@F)
	@cd $(builddir) && mv $(<F) ./preprocessed

$(builddir)/%.o : %.c
	@echo "[CC] $(<F)"
	@$(CC) $(CC_OPTS) -c $(<F) -o $@

TEST_%.e : $(testdir)/TEST_%.f03
	@echo "[F03] $(<F)"
	@$(FPP) $(FPPF) -I. -I$(builddir) $< > $(builddir)/$(<F)
	@cd $(builddir) && $(FC) -I. $(FC_OPTS) $(<F) $(notdir $(objs)) $(LDF) -o $@
	@cd $(builddir) && mv $(<F) ./preprocessed
	@cd $(builddir) && chmod a+x $@
	@if [ ! -d $(bindir) ];then \
		mkdir $(bindir); \
	fi
	@cp $(builddir)/$@ $(bindir)/$@
	@echo "Done!"

clean :
	@echo "[CLEAN] Removing build directory"
	@rm -rf $(builddir)

# depend : $(depfile)

# $(depfile) : $(main_src)
# 	@echo "Automatically generating module dependencies."
# 	@$(PYTHON) $(MAKEDEPEND) $(main_src) > $(depfile)
