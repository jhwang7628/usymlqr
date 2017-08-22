# Makefile for F90 version of USYMLQ.
# Maintained by Michael Saunders <saunders@stanford.edu>
#
# 24 Sep 2007: First version for compiling lsqrTestProgram.f90
#              and associated modules.
#              All files listed are .f90 source code.
#              C is not currently used.
# 25 May 2014: Makefile derived from USYMLQ Makefile.

  FC      =  gfortran
  FFLAGS  = -g -O
  FFLAGS  = -g -O0 -pedantic -Wall -Wextra -fbounds-check
# FFLAGS  = -g -O0 -pedantic -Wall -Wextra -fbounds-check -ftrace=full

  CC      =  gcc
  CFLAGS  = -g -O

# Clear suffix list, then define the ones we want
  .SUFFIXES:
  .SUFFIXES: .c .f .f90 .o

  .f90.o:; ${FC} ${FFLAGS} -c -o $@ $<
  .f.o:;   ${FC} ${FFLAGS} -c -o $@ $<
  .c.o:;   $(CC) $(CFLAGS) -c -o $@ $<

  files = usymlqDataModule.o usymlqblas.o        usymlqblasInterface.o \
          usymlqModule.o     usymlqCheckModule.o usymlqTestModule.o    \
          usymlqTestProgram.o

TestProgram: ${files}
	${FC} ${FFLAGS} -o $@ ${files}

clean:
	\rm -f *.o *.mod
