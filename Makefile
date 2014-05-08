############################################################
#
# Compiler and linker
#

CMP = gfortran -c
LNK = gfortran
OPT = -O2
OPT = -Og -g -Wall -Wextra -pedantic -fcheck=all


############################################################
#
# Objects: list of all objects *.o

OBJS =

############################################################
# Executable generation
#

pipeMeshNek:     pipeMeshNek.o  $(OBJS)
	$(LNK) $(OPT) pipeMeshNek.o  $(OBJS) \
                 -o ~/BTSync/Software/bin/pipeMeshNek

############################################################
# Objects generation

pipeMeshNek.o:   pipeMeshNek.f90  $(OBJS)
	$(CMP) $(OPT) pipeMeshNek.f90

############################################################
# Cleaning command to remove all objects *.o, files *.mod and
# the executables *.exe

clean: 
	@echo cleaning
	@rm -f  *.o  *.mod
  
# The command "make clean" deletes all the indicated files

############################################################
