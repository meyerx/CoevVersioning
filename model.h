/*
** model.h: Header file for the model.c functions.
**
** Linda Dib & Wim Hordijk   Last modified: 11 december 2014
*/

#ifndef _MODEL_H_
#define _MODEL_H_

void   transposeMatrix  (double **Qmatrix, double **QTmatrix);
void   condLikeFunction (double *condLike, double **Qmatrix, double **QTmatrix,
			 double *gi, double *gj, double ti, double tj,
			 double *MatrixA, double *MatrixB, double *MatrixC);
double executeCond      (double **Qm, double **Qtransposed, struct node *n,
			 double *probVector, double *MatrixA, double *MatrixB,
			 double *MatrixC);
void   setQ             (double **Qmatrix, double s, double d, double w1,
			 double w2, int *VectCoevComb );
void   setQNull         (double **Qmatrix, double w1, double w2);

void setQSimulator(double **Qmatrix, double s, double d, double r1, double r2,int *VectCoevComb,double *freqVector);

void   pVectorFunction (double *pVector, double **Qmatrix, double **QTmatrix, double ti, double *MatrixA, double *MatrixB, double *MatrixC);

#endif  /* _MODEL_H_ */
