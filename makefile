include ${SLEPC_DIR}/lib/slepc/conf/slepc_common
FFLAGS := $(shell pkg-config --cflags fgsl) -I/home/walker.2190/local/hdf5-1.14.3/hdf5-install/include -L/home/walker.2190/local/hdf5-1.14.3/hdf5-install
FLIBS := $(shell pkg-config --libs fgsl) -lhdf5_fortran -lhdf5

schrodinger1Df90: schrodinger1Df90.o
	-${FLINKER} ${FFLAGS} ${FLIBS} -o schrodinger1Df90 schrodinger1Df90.o ${SLEPC_EPS_LIB}
	-${RM} schrodinger1Df90.o

basisf90: basisf90.o
	-${FLINKER} ${FFLAGS} ${FLIBS} -o basisf90 basisf90.o ${SLEPC_EPS_LIB}
	-${RM} basisf90.o

stark_basisf90: stark_basisf90.o
	-${FLINKER} ${FFLAGS} ${FLIBS} -o stark_basisf90 stark_basisf90.o ${SLEPC_EPS_LIB}
	-${RM} stark_basisf90.o

left_basisf90: left_basisf90.o
	-${FLINKER} ${FFLAGS} ${FLIBS} -o left_basisf90 left_basisf90.o ${SLEPC_EPS_LIB}
	-${RM} left_basisf90.o

simulation_parametersf90: simulation_parametersf90.o
	-${FLINKER} ${FFLAGS} ${FLIBS} -o simulation_parametersf90 simulation_parametersf90.o ${SLEPC_EPS_LIB}
	-${RM} simulation_parametersf90.o

generateH0f90: generateH0f90.o
	-${FLINKER} ${FFLAGS} ${FLIBS} -o generateH0f90 generateH0f90.o ${SLEPC_EPS_LIB}
	-${RM} generateH0f90.o

generateH1f90: generateH1f90.o
	-${FLINKER} ${FFLAGS} ${FLIBS} -o generateH1f90 generateH1f90.o ${SLEPC_EPS_LIB}
	-${RM} generateH1f90.o

generateDipoleAccelerationf90: generateDipoleAccelerationf90.o
	-${FLINKER} ${FFLAGS} ${FLIBS} -o generateDipoleAccelerationf90 generateDipoleAccelerationf90.o ${SLEPC_EPS_LIB}
	-${RM} generateDipoleAccelerationf90.o

propagatef90: propagatef90.o
	-${FLINKER} ${FFLAGS} ${FLIBS} -o propagatef90 propagatef90.o ${SLEPC_EPS_LIB}
	-${RM} propagatef90.o

dc_starkf90: dc_starkf90.o
	-${FLINKER} ${FFLAGS} ${FLIBS} -o dc_starkf90 dc_starkf90.o ${SLEPC_EPS_LIB}
	-${RM} dc_starkf90.o

#------------------------------------------------------------------------------------
DATAPATH = ${SLEPC_DIR}/share/slepc/datafiles/matrices
