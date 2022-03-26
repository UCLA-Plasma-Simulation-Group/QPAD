# This is meant to be included by the parent Makefile
src += pcr_tri_class.f03
src += pcr_penta_class.f03
src += fpcr_tri_class.f03
src += fpcr_penta_class.f03

# Add local path to vpath
vpath %.f03 pcr-fortran

# If not at top level go there and make
ifeq ("$(toplevel)","")

all:
	make -C ..

clean:
	make -C .. clean

endif

