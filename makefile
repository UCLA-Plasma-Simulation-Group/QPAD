FC = mpif90 -g -fdefault-real-8 -fdefault-double-8
FORMAT_FREE = -ffree-form
LINKER = mpif90

HYPRE_LIB = -L/usr/local/hypre/lib -lHYPRE

# Objects list
OBJS_MODULE = param.o ufield_class.o field_class.o field_solver_class.o field_psi_class.o system.o
OBJS_TEST = TEST_ufield.o
OBJS_MAIN = main.o  

# Linkage rule
main :: ${OBJS_MODULE} ${OBJS_MAIN} 
	${LINKER} ${OBJS_MODULE} ${OBJS_MAIN} ${HYPRE_LIB} -o main.e

TEST_ufield :: ${OBJS_MODULE} TEST_ufield.o
	${LINKER} ${OBJS_MODULE} TEST_ufield.o ${HYPRE_LIB} -o TEST_ufield.e

TEST_field_psi :: ${OBJS_MODULE} TEST_field_psi.o
	${LINKER} ${OBJS_MODULE} TEST_field_psi.o ${HYPRE_LIB} -o TEST_field_psi.e

clean ::
	rm *.o; rm *.mod; rm *.e

# Module compilation rules
system.o : system.f03
	$(FC) -c ${FORMAT_FREE} system.f03 -o system.o

param.o : param.f03
	$(FC) -c ${FORMAT_FREE} param.f03 -o param.o

ufield_class.o : ufield_class.f03 param.o system.o
	${FC} -c ${FORMAT_FREE} ufield_class.f03 -o ufield_class.o

field_class.o : field_class.f03 param.o ufield_class.o system.o
	${FC} -c ${FORMAT_FREE} field_class.f03 -o field_class.o

field_solver_class.o : field_solver_class.f03 param.o system.o
	${FC} -c ${FORMAT_FREE} field_solver_class.f03 -o field_solver_class.o

field_psi_class.o : field_psi_class.f03 param.o field_class.o field_solver_class.o system.o
	${FC} -c ${FORMAT_FREE} field_psi_class.f03 -o field_psi_class.o

# Unit test compilation rules
TEST_ufield.o : TEST_ufield.f03 ufield_class.o
	${FC} -c ${FORMAT_FREE} TEST_ufield.f03 -o TEST_ufield.o
