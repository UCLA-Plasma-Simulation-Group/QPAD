# This is meant to be included by the parent Makefile
src += fdist2d_class.f03
src += part2d_class.f03
src += part2d_comm.f03
src += species2d_class.f03
src += neutral_class.f03

# Add local path to vpath
vpath %.f03 species

# If not at top level go there and make
ifeq ("$(toplevel)","")

all:
	make -C ..

clean:
	make -C .. clean

endif

