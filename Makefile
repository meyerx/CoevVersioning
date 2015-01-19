#
# Makefile: Makefile for the coev program.
#
# Linda Dib & Wim Hordijk   Last modified: 22 July 2014
# -lgsl -lgslcblas
CXX		= g++
CXXFLAGS	= -ansi -O3
CC		= gcc
CCFLAGS		= -O3
LIB		= -lblas -llapack -lm -lpthread -lnlopt
INC_PATH	= -I /usr/local/bin/ -I ./externaLibraries/
LIB_PATH	= -L /usr/local/bin/ -L ./externaLibraries/
OBJS		= ModelCPP.o bayes.o coev.o ml.o model.o simulate.o strmap.o tree.o


coev : $(OBJS)
	$(CXX) $(LIB_PATH) -fopenmp -o $@ $(OBJS) $(LIB)

%.o : %.cpp
	$(CXX) $(INC_PATH) -c $(CXXFLAGS) $<

%.o : %.c
	$(CC) $(INC_PATH) -c $(CCFLAGS) $<



#coev: coev.c model.c gtr.c tree.c bayes.c strmap.c def.h model.h strmap.h gtr.h tree.h GTR_CPP.cpp ModelCPP.cpp
#	g++ GTR_CPP.cpp ModelCPP.cpp
#	gcc -o tmpEx coev.c model.c bayes.c strmap.c gtr.c tree.c -lm -llapack -lblas -lnlopt -O3  -I /usr/local/bin/ -L /usr/local/bin/ -lgsl -lgslcblas

clean:
	rm -f *~ *.o
	rm coev
