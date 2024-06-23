# This is meant to be included by the parent Makefile

src += ufield_class.f03
src += field_solver_class.f03
src += field_solver_all_modes_class.f03
src += ufield_smooth_class.f03
src += field_class.f03
src += field_src_class.f03
src += field_psi_class.f03
src += field_b_class.f03
src += field_e_class.f03
src += field_vpot_class.f03

# Add local path to vpath
vpath %.f03 fields

# If not at top level go there and make
ifeq ("$(toplevel)","")

all:
	make -C ..

clean:
	make -C .. clean

endif

