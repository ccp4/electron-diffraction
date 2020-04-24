TEMSIM=/home/ronan/Documents/CCP4/src/multislice/kirkland/temsim
DOC_DIR=#/home/ronan/Documents/git/ccp4/docs/projects/multislice/
DEBUG ?= 0
BIN=bin/slicelib_main
MD=${DOC_DIR}/pyMultislice.md #&1 #STDOUT

slice_lib_main :
	g++ slicelib_main.cpp ${TEMSIM}/slicelib.o -I${TEMSIM} -o ${BIN}
multislice_doc :
	pydocmd simple multislice++ >${MD}
