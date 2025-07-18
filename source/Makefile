# system (must be one of the file suffix in config directory)
sys = conda
# sys = s3df (uncomment when compiling on s3df) 


# IF_OPENPMD -- uncomment the line below if you wish to use the OpenPMD format
IF_OPENPMD = 1
DISABLE_POPAS = TRUE

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

BUILD_CONFIGURATION_FILE ?= ../config/make.$(sys)
configfile=$(BUILD_CONFIGURATION_FILE)
include $(configfile)

# ------------------------------------------------------------------------------
# QPAD source files
# ------------------------------------------------------------------------------
src = dtimer.c param.f03 math_module.f03 parallel_module.f03 debug_tool.f03
src += kwargs_class.f03
src += sysutil_module.f03
src += ppmsg_class.f03
src += sort_module.f03
ifdef IF_OPENPMD
	src += hdf5io_pmd_class.f03
else
	src += hdf5io_class.f03
endif
src += interpolation.f03
src += input_class.f03
src += options_class.f03
src += fparser.f03
# PCR solver code
include pcr-fortran/Makefile
# fields code
include fields/Makefile
# laser code
include laser/Makefile
# species code
include species/Makefile
# beam code
include beam/Makefile
src += sim_fields_class.f03 sim_beams_class.f03 sim_plasma_class.f03 sim_lasers_class.f03
src += diagnostics_class.f03
src += simulation_class.f03
# template project code
include proj_template/Makefile
include proj_popas/Makefile
include proj_subcyc/Makefile
# src += main.f03
main_src = main.f03
# ------------------------------------------------------------------------------

# Parse source files to generate object files
objs := $(patsubst %.f03,$(builddir)/%.o,$(src))
objs := $(patsubst %.c,  $(builddir)/%.o,$(objs))
main_obj := $(builddir)/main.o

# Generate date stamp
DATESTAMP := $$(date +%s)

# Targets and rules
.DEFAULT_GOAL := $(qpad)
.PHONY: clean

$(qpad) : $(objs) $(main_obj)
	@echo "[LINK] $(@F)"
	@cd $(builddir) && $(LINKER) $(LINKER_OPTS) -I./$(builddir) -o $(notdir $(qpad)) $(notdir $(objs)) $(notdir $(main_obj)) $(LDF)
	@chmod a+x $@
	@echo "Copying binary to bin directory"
	@if [ ! -d $(bindir) ];then \
		mkdir $(bindir); \
	fi
	@cp $@ $(bindir)/qpad-$(DATESTAMP).e
	@echo "Creating symbolic link"
	@ln -f -s qpad-$(DATESTAMP).e $(bindir)/qpad.e
	@echo "Done!"

$(objs) : $(configfile) | $(builddir)

$(main_obj) : $(configfile) | $(builddir)

$(builddir) :
	@echo "Creating $@ directory"
	@mkdir -p $@

$(builddir)/%.o : %.f03
	@echo "[F03] $(<F)"
	@$(FPP) $(FPPF) -I. -I./$(builddir) $< > $(builddir)/$(<F)
	@cd $(builddir) && $(FC) $(INCF) -c $(FC_OPTS) $(<F) -o $(@F) && rm $(<F)

$(builddir)/%.o : %.c
	@echo "[CC] $(<F)"
	@$(CC) $(CC_OPTS) -c $< -o $@

module : $(objs)
	@echo "Only compile modules, no executable is generated."

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
