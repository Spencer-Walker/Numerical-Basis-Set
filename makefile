include ${SLEPC_DIR}/lib/slepc/conf/slepc_common

schrodinger1Df90: schrodinger1Df90.o chkopts
	-${FLINKER} -o schrodinger1Df90 schrodinger1Df90.o ${SLEPC_EPS_LIB}
	-${RM} schrodinger1Df90.o

diagonalizeHSecondOrderFDf90: diagonalizeHSecondOrderFDf90.o chkopts
	-${FLINKER} -o diagonalizeHSecondOrderFDf90\
	 diagonalizeHSecondOrderFDf90.o ${SLEPC_EPS_LIB}
	-${RM} diagonalizeHSecondOrderFDf90.o 

diagonalizeHFourthOrderFDf90: diagonalizeHFourthOrderFDf90.o chkopts
	-${FLINKER} -o diagonalizeHFourthOrderFDf90\
	 diagonalizeHFourthOrderFDf90.o ${SLEPC_EPS_LIB}
	-${RM} diagonalizeHFourthOrderFDf90.o

diagonalizeHSixthOrderFDf90: diagonalizeHSixthOrderFDf90.o chkopts
	-${FLINKER} -o diagonalizeHSixthOrderFDf90\
	 diagonalizeHSixthOrderFDf90.o ${SLEPC_EPS_LIB}
	-${RM} diagonalizeHSixthOrderFDf90.o

diagonalizeHEighthOrderFDf90: diagonalizeHEighthOrderFDf90.o chkopts
	-${FLINKER} -o diagonalizeHEighthOrderFDf90\
	 diagonalizeHEighthOrderFDf90.o ${SLEPC_EPS_LIB}
	-${RM} diagonalizeHEighthOrderFDf90.o


diagonalizeHMatrixNumerovf90: diagonalizeHMatrixNumerovf90.o chkopts
	-${FLINKER} -o diagonalizeHMatrixNumerovf90\
	 diagonalizeHMatrixNumerovf90.o ${SLEPC_EPS_LIB}
	-${RM} diagonalizeHMatrixNumerovf90.o

generateH0f90: generateH0f90.o chkopts
	-${FLINKER} -o generateH0f90 generateH0f90.o ${SLEPC_EPS_LIB}
	-${RM} generateH0f90.o

generateH1f90: generateH1f90.o chkopts
	-${FLINKER} -o generateH1f90 generateH1f90.o ${SLEPC_EPS_LIB}
	-${RM} generateH1f90.o

propagatef90: propagatef90.o chkopts
	-${FLINKER} -o propagatef90 propagatef90.o ${SLEPC_EPS_LIB}
	-${RM} propagatef90.o

#------------------------------------------------------------------------------------
DATAPATH = ${SLEPC_DIR}/share/slepc/datafiles/matrices
