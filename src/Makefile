PROGRAM=MIM

#VERSION=dev
VERSION=0.33r1#   ### rewrite lagmain.f90

INSTALL=~/bin


##### Fujitsu f90 #####
#FC=f90
#OPTION=-Am -X9 -Kfast -Et#                  # optimized
#OPTION=-Am -X9 -O0#                         # no optimized, safe
#OPTION=-Am -X9 -Et -O0 -Haesux#             # for debug


##### Intel ifort (for 10.0) #####
FC=ifort
OPTION=-warn all -O1 -assume byterecl -convert little_endian -traceback -check all
OPTION=-warn all -O1 -assume byterecl -convert little_endian -traceback
#OPTION=-assume byterecl -O3 -fast -warn all -traceback -heap-arrays    # optimized
#OPTION=-assume byterecl -O2 -warn all -traceback -heap-arrays    # optimized
#OPTION=-assume byterecl -O0 -warn all -traceback -g -heap-arrays#    # no optimized
#OPTION=-assume byterecl -O0 -warn all -g -traceback -heap-arrays#    # no optimized
#OPTION=-assume byterecl -C -warn all#     # for debug

##### gfortran
#FC=gfortran
#OPTION=-frecord-marker=4 -O3 -Wall -fbacktrace -g -fcheck=all -fbounds-check

OBJS =
OBJS_MODULE = myglobal.o

include Makefile.inc
OBJS += ${OBJS_SRC}
OBJS_MODULE += ${OBJS_MODULE_SRC}


all : ${PROGRAM}

${PROGRAM} : ${OBJS_MODULE} ${OBJS}
	${FC} ${OPTION} ${OBJS_MODULE} ${OBJS} -o $@


.SUFFIXES : .o .f90
.f90.o :
	${FC} ${OPTION} -c $< -o $@


all : ${PROGRAM}

run :
	./${PROGRAM}

clean:
	@ rm -vf ${OBJS_MODULE} *.mod ${OBJS} ${PROGRAM} *genmod*

re : clean all

release:
	rm -rf ${OBJS_MODULE} ${OBJS_MODULE:.o=.f90~} *.mod \
	       ${OBJS} ${OBJS:.o=.f90~} ${PROGRAM} *~

install:
	cp ${PROGRAM} ${PROGRAM}-${VERSION}
	mv ${PROGRAM}-${VERSION} ${INSTALL}
	rm -f ${INSTALL}/${PROGRAM}
	ln -s ${INSTALL}/${PROGRAM}-${VERSION} ${INSTALL}/${PROGRAM}
