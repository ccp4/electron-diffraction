# makefile to generate TEMSIM multislice package.
# assumes that FFTW already installed (see www.fftw.org)
#  type "make all" to compile everthing
#  type "make remove" to remove all of the compiled files.

CC = g++ -O3
DEL = rm  # unix
#
ODIR=../bin/obj
BIN=../bin
_MYLIBS=slicelib.o floatTIFF.o cfpix.o

# define libraries
MYLIBS=$(patsubst %,$(ODIR)/%,$(_MYLIBS))
LIBS=-lfftw3f -lfftw3f_threads -lm -lpthread
WLIBS=${MYLIBS} ${LIBS}

all:
	make atompot
	make autoslic
	make autostem
	make image
	make incostem
	make mulslice
	make probe
	make stemslic
	make sumpix
remove:
	${DEL}$ atompot
	${DEL}$ autoslic
	${DEL}$ autostem
	${DEL}$ image
	${DEL}$ incostem
	${DEL}$ mulslice
	${DEL}$ probe
	${DEL}$ stemslic
	${DEL}$ sumpix
	${DEL}$ ${ODIR}/cfpix.o
	${DEL}$ ${ODIR}/slicelib.o
	${DEL}$ ${ODIR}/floatTIFF.o


#
#  main programs
#

atompot: atompot.cpp ${MYLIBS}
	${CC} -o ${BIN}/atompot atompot.cpp $(WLIBS)

autoslic: autoslic.cpp autosliccmd.cpp ${MYLIBS}
	${CC} -fopenmp -o ${BIN}/autoslic autosliccmd.cpp autoslic.cpp ${WLIBS}

autostem: autostem.cpp autostemcmd.cpp ${MYLIBS}
	${CC} -fopenmp -o ${BIN}/autostem autostemcmd.cpp autostem.cpp ${WLIBS}

display: display.cpp  ${MYLIBS}
	${CC} -o ${BIN}/display display.cpp ${WLIBS}

image: image.cpp  ${MYLIBS}
	${CC} -o ${BIN}/image image.cpp ${WLIBS}

incostem: incostem.cpp incostemcmd.cpp ${MYLIBS}
	${CC} -fopenmp -o ${BIN}/incostem incostemcmd.cpp incostem.cpp probe.cpp ${WLIBS}

mulslice: mulslice.cpp ${MYLIBS}
	${CC} -o ${BIN}/mulslice mulslice.cpp ${WLIBS}

probe: probe.cpp probecmd.cpp ${MYLIBS}
	${CC} -o ${BIN}/probe probecmd.cpp probe.cpp ${WLIBS}

slicview: slicview.cpp ${MYLIBS}
	${CC} -o ${BIN}/slicview slicview.cpp ${WLIBS}

stemslic: stemslic.cpp ${MYLIBS}
	${CC} -o ${BIN}/stemslic stemslic.cpp ${WLIBS}

sumpix: sumpix.cpp ${MYLIBS}
	${CC} -o ${BIN}/sumpix sumpix.cpp ${WLIBS}

slice_lib_main: slicelib_main.cpp $(ODIR)/slicelib.o
	${CC} -o ${BIN}/slicelib_main slicelib_main.cpp ${ODIR}/slicelib.o

#
# define subroutine library
#

$(ODIR)/cfpix.o:cfpix.cpp
	${CC} -c cfpix.cpp -o $(ODIR)/cfpix.o

$(ODIR)/slicelib.o:slicelib.cpp
	${CC} -c slicelib.cpp -o $(ODIR)/slicelib.o

$(ODIR)/floatTIFF.o:
	${CC} -c floatTIFF.cpp -o $(ODIR)/floatTIFF.o
