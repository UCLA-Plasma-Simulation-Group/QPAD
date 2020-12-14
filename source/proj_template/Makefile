# This is meant to be included by the parent Makefile
ifndef DISABLE_TEMPLATE
	
FPPF += -DENABLE_TEMPLATE

src += part3d_tmplt_class.f03
src += beam3d_tmplt_class.f03
src += sim_beams_tmplt_class.f03
src += simulation_tmplt_class.f03

# Add local path to vpath
vpath %.f03 proj_template

# If not at top level go there and make
ifeq ("$(toplevel)","")

all:
	make -C ..

clean:
	make -C .. clean

endif

endif