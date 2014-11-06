/*
** def.h: Header file with global variables for the coev program.
**
** Linda Dib & Wim Hordijk   Last modified: 30 April 2013
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

#define NT          1
#define AA          2

#define JC          1
#define HKY         2
#define GTR         3

struct node{
  char         data[8], label[128];
  double       dist, *probVector;
  struct node *parent, *left, *right;
};

struct sequence
{
  char            *label, *sequence;
  struct sequence *next;
};

extern int          IT, sample_freq, print_freq, burnin, nrComb, dataType,
                   *obsCombs, **profiles, nrProfiles, maxNrProfiles,
                    innerLoopIteration, model;
extern double       s, d, alpha, beta; 
extern double       r1, r2;  
extern char         outFile[1024], *ntComb[nrNtComb], *aaComb[nrAaComb];
extern struct node  root;


/*
** Function prototypes.
*/
int isConflict (int index1, int index2);
int bayes      ();
int ml         ();

#endif  /* _DEF_H_ */
