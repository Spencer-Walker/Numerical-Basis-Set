include ${SLEPC_DIR}/lib/slepc/conf/slepc_common
# FFLAGS := $(shell pkg-config --cflags fgsl)
# FLIBS := $(shell pkg-config --libs fgsl)

schrodinger1Df90: schrodinger1Df90.o chkopts
	-${FLINKER} -o schrodinger1Df90 schrodinger1Df90.o ${SLEPC_EPS_LIB} 
	-${RM} schrodinger1Df90.o

basisf90: basisf90.o chkopts
	-${FLINKER} -o basisf90 basisf90.o ${SLEPC_EPS_LIB}
	-${RM} basisf90.o

left_basisf90: left_basisf90.o chkopts
	-${FLINKER} -o left_basisf90 left_basisf90.o ${SLEPC_EPS_LIB}
	-${RM} left_basisf90.o

simulation_parametersf90: simulation_parametersf90.o chkopts
	-${FLINKER} -o simulation_parametersf90 simulation_parametersf90.o ${SLEPC_EPS_LIB}
	-${RM} simulation_parametersf90.o

generateH0f90: generateH0f90.o chkopts
	-${FLINKER} -o generateH0f90 generateH0f90.o ${SLEPC_EPS_LIB}
	-${RM} generateH0f90.o

generateH1f90: generateH1f90.o chkopts
	-${FLINKER} -o generateH1f90 generateH1f90.o ${SLEPC_EPS_LIB}
	-${RM} generateH1f90.o

generateDipoleAccelerationf90: generateDipoleAccelerationf90.o chkopts
	-${FLINKER} -o generateDipoleAccelerationf90 generateDipoleAccelerationf90.o ${SLEPC_EPS_LIB}
	-${RM} generateDipoleAccelerationf90.o

propagatef90: propagatef90.o chkopts
	-${FLINKER} -o propagatef90 propagatef90.o ${SLEPC_EPS_LIB}
	-${RM} propagatef90.o

#------------------------------------------------------------------------------------
DATAPATH = ${SLEPC_DIR}/share/slepc/datafiles/matrices