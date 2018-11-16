# This is meant to be included by the parent Makefile

src += part3d_lib.f03 fdist3d_class.f03 part3d_class.f03 beam_class.f03

# Add local path to vpath
vpath %.f03 species

# If not at top level go there and make
ifeq ("$(toplevel)","")

all:
	make -C ..

clean:
	make -C .. clean

endif

