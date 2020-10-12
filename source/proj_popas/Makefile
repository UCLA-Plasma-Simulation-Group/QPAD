# This is meant to be included by the parent Makefile
src += part3d_popas_class.f03
src += beam3d_popas_class.f03
src += sim_beams_popas_class.f03
src += diagnostics_popas_class.f03
src += simulation_popas_class.f03

# Add local path to vpath
vpath %.f03 proj_popas

# If not at top level go there and make
ifeq ("$(toplevel)","")

all:
	make -C ..

clean:
	make -C .. clean

endif