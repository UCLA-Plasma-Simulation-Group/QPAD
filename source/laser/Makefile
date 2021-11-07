# This is meant to be included by the parent Makefile

src += pcr_solver_class.f03
src += laser_solver_class.f03
# src += laser_pgc_class.f03

# Add local path to vpath
vpath %.f03 laser

# If not at top level go there and make
ifeq ("$(toplevel)","")

all:
	make -C ..

clean:
	make -C .. clean

endif