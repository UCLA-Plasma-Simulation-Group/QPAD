# This is meant to be included by the parent Makefile
ifndef DISABLE_SUBCYC
	
FPPF += -DENABLE_SUBCYC

src += part2d_subcyc_class.f03
src += species2d_subcyc_class.f03
src += neutral_subcyc_class.f03
src += sim_plasma_subcyc_class.f03
src += simulation_subcyc_class.f03

# Add local path to vpath
vpath %.f03 proj_subcyc

# If not at top level go there and make
ifeq ("$(toplevel)","")

all:
	make -C ..

clean:
	make -C .. clean

endif

endif