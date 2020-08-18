# This is meant to be included by the parent Makefile
src += part3d_comm.f03
src += fdist3d_lib.f03
src += fdist3d_class.f03
src += part3d_class.f03
src += beam3d_class.f03

# Add local path to vpath
vpath %.f03 beam

# If not at top level go there and make
ifeq ("$(toplevel)","")

all:
	make -C ..

clean:
	make -C .. clean

endif

