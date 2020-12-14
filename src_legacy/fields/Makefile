# This is meant to be included by the parent Makefile

src += ufield_class.f03 ufield_smooth_class.f03 field_class.f03 field_src_class.f03 \
field_psi_class.f03 field_e_class.f03 field_b_class.f03 field_solver_class.f03

# Add local path to vpath
vpath %.f03 fields

# If not at top level go there and make
ifeq ("$(toplevel)","")

all:
	make -C ..

clean:
	make -C .. clean

endif

