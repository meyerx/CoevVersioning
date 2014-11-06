#
# Makefile: Makefile for the coev program.
#
# Linda Dib & Wim Hordijk   Last modified: 26 April 2013
#

coev: coev.c model.c  ml.c bayes.c  simulate.c tree.c strmap.c def.h model.h strmap.h tree.h
	gcc -g -o coev coev.c model.c ml.c bayes.c simulate.c strmap.c tree.c -lm -llapack -lblas -lnlopt -lpthread -g -ggdb -O0  -I /usr/local/bin/ -L /usr/local/bin/ -lgsl -lgslcblas

clean:
	rm -f *~ *.o
	rm coev