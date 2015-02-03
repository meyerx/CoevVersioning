#
# Makefile: Makefile for the coev program.
#
# Linda Dib & Wim Hordijk   Last modified: 26 April 2013
#
# Modifications by HS for RPM package and Vital-IT environment (Oct 2014)
# 
# Default pathes for Vital-IT environment
PREFIX=/software
LIBDIR= -L ${PREFIX}/lib
LIB64DIR= -L ${PREFIX}/lib64
LIB		= -lblas -llapack -lm -lgsl -lgslcblas -lpthread -lnlopt
NLOPT_LIBDIR=  ${PREFIX}/Utility/nlopt/2.3/lib
NLOPT_INCDIR= -I ${PREFIX}/Utility/nlopt/2.3/include
EIGEN_INCDIR= -I ${PREFIX}/Utility/Eigen/3.1.2/include/Eigen
CXX		= g++
CXXFLAGS	= -ansi -O3
CC		= gcc
CCFLAGS		= -O3
OBJS		= ModelCPP.o bayes.o coev.o ml.o model.o simulate.o strmap.o tree.o
coev : $(OBJS)
	$(CXX)  -L ${NLOPT_LIBDIR} ${LIBDIR} ${LIB64DIR} $(LIB) -Wl,-rpath=${NLOPT_LIBDIR} -o $@ $(OBJS)
%.o : %.cpp
	$(CXX)  $(NLOPT_INCDIR) $(EIGEN_INCDIR) -c $(CXXFLAGS) $<

%.o : %.c
	$(CC)  $(NLOPT_INCDIR) $(EIGEN_INCDIR) -c $(CCFLAGS) $<
	
test:
	./coev -data nt -method sm -tree example/tree.txt -out example/simulate.out -s 1 -d 10 -r1 5 -r2 5
	rm example/simulate.out

clean:
	rm -f *~ *.o
	rm coev






