/*
** def.h: Header file with global variables for the coev program.
**
** Linda Dib & Wim Hordijk   Last modified: 11 december 2014
*/

#ifndef _DEF_H_
#define _DEF_H_

#include <stdio.h> 
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include <nlopt.h>
#include <time.h>
#include <pthread.h>
#include <sys/time.h>
/*
** Global variables.
*/
#define nrNt        4
#define nrNtComb   16
#define nrAa       20
#define nrAaComb  400
#define nrNtProf  192
#define nrAaProf 4000

#define BAYES       1
#define ML          2
#define SM          3

#define NT          1
#define AA          2

#define JC          1
#define HKY         2
#define GTR         3

#define MAX_LINE_LEN 10000
extern int          IT, sample_freq, print_freq, burnin, nrComb, dataType,
                   *obsCombs, **profiles, nrProfiles, maxNrProfiles,
                    innerLoopIteration, model, nuplet,col1,col2, nrSimulations,opt;
extern double       s, d, alpha, beta; 
extern double       r1, r2;  

extern char        treeFile[100000], alignFile[100000], outFile[100000], *aaComb[nrAaComb],*ntComb[nrNtComb], tree[10*MAX_LINE_LEN];
extern struct node  root;
extern int nrSimulations;
extern struct sequence seqList;
extern const char nucleotide[nrNt];
/*
** Function prototypes.
*/
int isConflict (int index1, int index2);
int bayes      ();
int ml         ();
int simulate   ();
#endif  /* _DEF_H_ */
