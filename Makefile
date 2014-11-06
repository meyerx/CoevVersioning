#
# Makefile: Makefile for the coev program.
#
# Linda Dib & Wim Hordijk   Last modified: 26 April 2013
#

coev: coev.c model.c bayes.c ml.c strmap.c def.h model.h strmap.h
	gcc -o coev coev.c model.c bayes.c ml.c strmap.c -lm -llapack -lblas -lnlopt -lpthread -g -ggdb -O0

clean:
	rm -f *~ *.o
	rm coev