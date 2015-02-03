/*
 ** model: Model functions.
 **
 ** Linda Dib & Wim Hordijk   Last modified: 11 december 2014
 */

#include "def.h"
#include "tree.h"



//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MODEL FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%
void transposeMatrix (double **Qmatrix, double **QTmatrix){
	int loop_i=0, loop_j=0;
	
	for (loop_i=0;loop_i<nrComb;loop_i++){ 
		for (loop_j=0;loop_j<nrComb;loop_j++){ 
			QTmatrix[loop_j][loop_i]=Qmatrix[loop_i][loop_j];
		}
	}
}

void condLikeFunction(double *condLike, double **Qmatrix, double **QTmatrix, double *gi, double *gj, double ti, double tj, double *MatrixA, double *MatrixB, double *MatrixC){
	int loop_i=0, loop_j=0;
	
	//FILE *fp;
	//fp = fopen("DEBUG-EValues-EVectors.txt", "w");
	/*
	 ** QTmatrixSave does not seem to be used anywhere.
	 **
	 double QTmatrixSave[nrComb][nrComb];
	 for (loop_i=0;loop_i<nrComb;loop_i++){ 
	 for (loop_j=0;loop_j<nrComb;loop_j++){ 
	 //printf("%g,",  QTmatrix[loop_i][loop_j]);
	 QTmatrixSave[loop_i][loop_j]=QTmatrix[loop_i][loop_j];
	 }
	 //printf("\n;");
	 }
	 */
	
	////////////////////////////////////////////////////////////////////////////////
	// EIGEN VALUE & EIGEN VECTORS//////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////
	//MATLAB: [eigenVector,eigenValue]=eig(gQ);
	//Compute egen value and eigen vector in c with dgeev
	
	//ARGUMENTS
	int N=nrComb;
	int LDA=nrComb;
	double *A= (double *)malloc(LDA*N*sizeof(double));
	for (loop_i = 0; loop_i < nrComb; loop_i++) {
		for (loop_j = 0; loop_j < nrComb; loop_j++) {
			A[ loop_i* nrComb + loop_j]=QTmatrix[loop_i][loop_j];
		}
	}
	//A=&(QTmatrix[0][0]);// USE TRANSPOSE MATRIX SINCE THOSE METHODS ARE COLUMN BASED
	double *WR=(double *)malloc(nrComb*sizeof(double));
	double *WI=(double *)malloc(nrComb*sizeof(double));
	int LDVL=nrComb;
	double *VL=(double *)malloc(LDVL*nrComb*sizeof(double));
	int LDVR=nrComb;
	//int LWORK=4*nrComb;
	int LWORK=nrComb*nrComb;
	double *VR=(double *)malloc(LDVR*nrComb*sizeof(double));
	double *WORK=(double *)malloc(LWORK*sizeof(double));
	int INFO=0;
	
	//CALL OF FUNCTION//////////////////////////////////////////////////////////////
	LWORK=-1;
	dgeev_("N", "V", &N, A, &LDA, WR, WI, VL, &LDVL, VR, &LDVR, WORK, &LWORK, &INFO );
	LWORK=(int) WORK[0];
	WORK = realloc(WORK, LWORK*sizeof(double));
	for (loop_i = 0; loop_i < nrComb; loop_i++) {
		for (loop_j = 0; loop_j < nrComb; loop_j++) {
			A[ loop_i* nrComb + loop_j]=QTmatrix[loop_i][loop_j];
		}
	}
	dgeev_("N", "V", &N, A, &LDA, WR, WI, VL, &LDVL, VR, &LDVR, WORK, &LWORK, &INFO );
	//SaveToOctave (&(Qmatrix[0][0]), "Q", fp, 16, 16) ;
	//SaveToOctave (&(QTmatrix[0][0]), "QT", fp, 16, 16) ;
	//SaveToOctave (WR, "evalues", fp, 1, 16) ;
	//SaveToOctave (VR, "evectors", fp, 16, 16) ;
	
	//freeArray A
	free(A);
	
	free(WI);free(VL);free(WORK);
	
	///SAVE VR
	double *VR_SAV1 = (double *)malloc (nrComb*nrComb*sizeof (double));
	double *VR_SAV2 = (double *)malloc (nrComb*nrComb*sizeof (double));
	
	//From pointer to [][]
	for (loop_j = 0; loop_j < nrComb*nrComb; loop_j++)//columns
	{   	
		VR_SAV1[loop_j]= VR[loop_j];
		VR_SAV2[loop_j]= VR[loop_j]; 		
	}
	
	////////////////////////////////////////////////////////////////////////////////
	// Inverse Of Eigen vector//////////////////////////////////////////////////////
	// MATLAB:inv(eigenVector)
	////////////////////////////////////////////////////////////////////////////////
	//STEP1: FACTORISATION//////////////////////////////////////////////////////////
	//ARGUMENTS/////////////////////////////////////////////////////////////////////
	int N2=nrComb;
	int LDA2=nrComb;
	double *VR_inv=(double *)malloc(LDA2*N2*sizeof(double));
	
	for (loop_i = 0; loop_i < LDA2*N2; loop_i++)
	{
		VR_inv[loop_i] = VR[loop_i];
	}
	int *IPIV2=(int *)malloc(nrComb*sizeof(int));
	int row=nrComb;
	int column=nrComb;
	int INFO2=0;
	
	//CALL OF FUNCTION//////////////////////////////////////////////////////////////
	dgetrf_(&row, &column, VR_inv, &LDA2, IPIV2, &INFO2 );
	//printf("INFO2: %d, output2:%d\n", INFO2, output2);
	//STEP2: INVERSE////////////////////////////////////////////////////////////////
	//ARGUMENTS/////////////////////////////////////////////////////////////////////
	int LWORK2=4*nrComb;
	double *WORK2=(double *)malloc(LWORK2*sizeof(double));
	int INFO3=0;
	LWORK2=-1;
	
	//CALL OF FUNCTION//////////////////////////////////////////////////////////////
	dgetri_(&N2, VR_inv, &LDA2, IPIV2, WORK2, &LWORK2, &INFO3);
	LWORK2=(int) WORK2[0];
	free(WORK2);
	WORK2=(double *)malloc(LWORK2*sizeof(double));
	dgetri_(&N2, VR_inv, &LDA2, IPIV2, WORK2, &LWORK2, &INFO3 );
	//VR_inv is a transposed matrix 
	//SaveToOctave (VR_inv, "VR_invFINALRESULT", fp, 16, 16) ;
	
	//printf("INFO3: %d, output3:%d\n", INFO3, output3);
	////////////////////////////////////////////////////////////////////////////////
	// PREPARE EXPRESSION EVALUATION////////////////////////////////////////////////
	// USE TRANSPOSE MATRIX SINCE THOSE METHODS RETURN COLUMN BASED MATRICES
	////////////////////////////////////////////////////////////////////////////////
	//--VR_invInput	:eigenVectors: transpose of dgetri output VR_inv
	//--WR_input	:eigenValues : transpose of dgeev  output WR
	
	//From pointer to [][]
	double *VR_invDgetri = (double *)malloc (nrComb*nrComb*sizeof (double));
	for (loop_j = 0; loop_j < nrComb*nrComb; loop_j++){//columns of the file	
		VR_invDgetri[loop_j]= VR_inv[loop_j];
	}
	//Transpose the matrix VR_inv matrix
	//double VR_invInput[nrComb][nrComb];  // Doesn't seem to be used anywhere!
	
	//create a eigen diagonal matrix 
	//Initialise matrix
	double **WR_diagEXP_ti = (double **)malloc (nrComb * sizeof (double *));
	double **WR_diagEXP_tj = (double **)malloc (nrComb * sizeof (double *));
	for (loop_i=0;loop_i<nrComb;loop_i++)
	{
		WR_diagEXP_ti[loop_i] = (double *)malloc (nrComb * sizeof (double));
		WR_diagEXP_tj[loop_i] = (double *)malloc (nrComb * sizeof (double));
		for (loop_j=0;loop_j<nrComb;loop_j++){
			WR_diagEXP_ti[loop_i][loop_j]=0;
			WR_diagEXP_tj[loop_i][loop_j]=0;
		}
	}
	
	//MATLAB: expm(eigenValue*ti) and expm(eigenValue*tj)
	for (loop_j=0;loop_j<nrComb;loop_j++){
		WR_diagEXP_ti[loop_j][loop_j]=exp(WR[loop_j] * ti);
		WR_diagEXP_tj[loop_j][loop_j]=exp(WR[loop_j] * tj);
	}
	////////////////////////////////////////////////////////////////////////////////
	// MULTIPLY MATRICES ///////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////
	//MATLAB	:pii(1:16,1:16) = eigenVector*expm(eigenValue*ti)*inv(eigenVector);
	//MATLAB	:pj(1:16,1:16)  = eigenVector*expm(eigenValue*tj)*inv(eigenVector);
	//OUTPUT DGEMMM: WE OBTAIN PI AND PJ TRANSPOSED 
	////////////////////////////////////////////////////////////////////////////////
	//ARGUMENTS
	int  M=nrComb;
	N=nrComb;
	int K=nrComb;
	
	double ALPHA=1.0;
	int LDA_dgemm=nrComb;
	int LDB_dgemm=nrComb;
	double BETA=0.0; 
	int LDC_dgemm=nrComb;
	int output_dgemm=0;	
	
	//CALL OF FUNCTION DGEMM: RES1= expm(eigenValue*ti)*inv(eigenVector)/////////////
	//MatrixA: is the diagonal matrix WR_diagEXP_ti
	//MatrixA=&(WR_diagEXP_ti[0][0]);
	for (loop_i = 0; loop_i < nrComb; loop_i++) {
		for (loop_j = 0; loop_j < nrComb; loop_j++) {
			MatrixA[ loop_i* nrComb + loop_j]=WR_diagEXP_ti[loop_i][loop_j];
		}
	}
	//MatrixB=&(VR_invDgetri[0]);
	for (loop_i = 0; loop_i < nrComb*nrComb; loop_i++)
	{
		MatrixB[loop_i] = VR_invDgetri[loop_i];
	}
	for (loop_i = 0; loop_i < nrComb; loop_i++)//rows of the file
		for (loop_j = 0; loop_j < nrComb; loop_j++)//columns of the file
			MatrixC[loop_i * nrComb + loop_j]=0.0;
	output_dgemm=dgemm_("N","N", &M, &N, &K, &ALPHA, MatrixA, &LDA_dgemm, MatrixB, &LDB_dgemm, &BETA, MatrixC, &LDC_dgemm);
	
	//CALL OF FUNCTION DGEMM: pii=eigenVector*RES1////////////////////////////////////
	//MatrixA: is the transpose of VR_invInput: (that is the original VR_inv or VR_invDgetri)
	//MatrixB: is MatrixC
	//MatrixA=&(VR_SAV1[0]);
	//MatrixB=MatrixC;
	for (loop_i = 0; loop_i < nrComb*nrComb; loop_i++) {
		MatrixA[loop_i]=VR_SAV1[loop_i];
		MatrixB[loop_i] = MatrixC[loop_i];
	}
	double *pii=(double *)malloc(M*N*sizeof(double));//OUTPUT
	for (loop_i = 0; loop_i < nrComb; loop_i++)//rows of the file
		for (loop_j = 0; loop_j < nrComb; loop_j++)//columns of the file
			pii[loop_i * nrComb + loop_j]=0.0;
	output_dgemm=dgemm_("N","N", &M, &N, &K, &ALPHA, MatrixA, &LDA_dgemm, MatrixB, &LDB_dgemm, &BETA, pii, &LDC_dgemm);
	//SaveToOctave (pii, "pii", fp, 16, 16) ;
	
	//////////////////////////////////////////////////////////////////////////////////
	//CALL OF FUNCTION DGEMM: RES2= expm(eigenValue*tj)*inv(eigenVector)//////////////
	//MatrixA=&(WR_diagEXP_tj[0][0]);
	for (loop_i = 0; loop_i < nrComb; loop_i++) {
		for (loop_j = 0; loop_j < nrComb; loop_j++) {
			MatrixA[ loop_i* nrComb + loop_j]=WR_diagEXP_tj[loop_i][loop_j];
		}
	}
	//MatrixB=&(VR_invDgetri[0]);
	for (loop_i = 0; loop_i < nrComb*nrComb; loop_i++)
	{
		MatrixB[loop_i] = VR_invDgetri[loop_i];
	}
	for (loop_i = 0; loop_i < nrComb; loop_i++)//rows of the file
		for (loop_j = 0; loop_j < nrComb; loop_j++)//columns of the file
			MatrixC[loop_i * nrComb + loop_j]=0.0;
	
	output_dgemm=dgemm_("N","N", &M, &N, &K, &ALPHA, MatrixA, &LDA_dgemm, MatrixB, &LDB_dgemm, &BETA, MatrixC, &LDC_dgemm);
	
	//CALL OF FUNCTION DGEMM: pjj=eigenVector*RES2////////////////////////////////////
	//MatrixA=VR_SAV2;
	//MatrixB=MatrixC;
	for (loop_i = 0; loop_i < nrComb*nrComb; loop_i++) {
		MatrixA[loop_i] = VR_SAV2[loop_i];
		MatrixB[loop_i] = MatrixC[loop_i];
	}
	double *pjj=(double *)malloc(M*N*sizeof(double));//OUTPUT
	for (loop_i = 0; loop_i < nrComb; loop_i++)//rows of the file
		for (loop_j = 0; loop_j < nrComb; loop_j++)//columns of the file
			pjj[loop_i * nrComb + loop_j]=0.0;
	output_dgemm=dgemm_("N","N", &M, &N, &K, &ALPHA, MatrixA, &LDA_dgemm, MatrixB, &LDB_dgemm, &BETA, pjj, &LDC_dgemm);
	//SaveToOctave (pjj, "pjj", fp, 16, 16) ;
	
	free(WR);
	free(VR_inv);
	free (VR_invDgetri);
	free(IPIV2);
	free (VR);
	free (VR_SAV1);
	free (VR_SAV2);
	for (loop_i=0; loop_i < nrComb; loop_i++)
	{
		free (WR_diagEXP_ti[loop_i]);
		free (WR_diagEXP_tj[loop_i]);
	}
	free (WR_diagEXP_ti);
	free (WR_diagEXP_tj);
	
	//////////////////////////////////////////////////////////////////////////////////
	// MATLAB:hi(1:16) = pii(1:16,1:16)*(gi(1:16).'); 
	// SUBROUTINE DGEMV ( TRANS, M, N, ALPHA, A, LDA, X, INCX, BETA, Y, INCY )
	M=nrComb;
	N=nrComb;
	ALPHA=1.0;
	BETA =0.0;
	int LDA_dgemv=nrComb;
	int INCX=1;
	int INCY=1;
	double *X1=gi;
	double *hi=(double *)malloc(nrComb*sizeof(double));//transpose of gi
	
	
	////////////////////////////////ERROR////////////////////////////////////////////
	output_dgemm=dgemv_("N", &M, &N, &ALPHA, pii, &LDA_dgemv, X1, &INCX,&BETA, hi, &INCY );	
    
    
	//MATLAB:hj(1:16) = pj(1:16,1:16) *(gj(1:16).');
	double *X2=gj;
	BETA =0.0;
	double *hj=(double *)malloc(1*nrComb*sizeof(double));//transpose of gi
	output_dgemm=dgemv_("N", &M, &N, &ALPHA, pjj, &LDA_dgemv, X2, &INCX,&BETA, hj,&INCY );
	free(pjj);
	free(pii);
	
	//SaveToOctave (hi, "hi", fp, 1, 16) ;
	//SaveToOctave (hj, "hj", fp, 1, 16) ;
	
	///////////////////////////////////////////////////////////////////////////////////	
	//MATLAB:condLike(1:16) = hi(1:16).*hj(1:16);
	///////////////////////////////////////////////////////////////////////////////////
	//printf("condLike VECTOR\n");
	for (loop_j=0;loop_j<nrComb;loop_j++){ 
		condLike[loop_j]= hi[loop_j] * hj[loop_j] ;
		//printf("%g ",condLike[loop_j]);
	}
	//ffclose (fp);
	
	free(hj);
	free(hi);
	free(WORK2);
	
}

double executeCond (double **Qm, double **Qtransposed, struct node* n, double *probVector, double *MatrixA, double *MatrixB, double *MatrixC) { 
	
	int assign=0;
	if (n != NULL) { 
		if((n->left == NULL) && (n->right == NULL))// is leaf
		{
			for (assign=0;assign<nrComb;assign++){ 
				probVector[assign]=n->probVector[assign];
				
			}
			return n->dist; 
		}else{
			double *l_probVector = (double *)malloc (nrComb*sizeof (double));
			double *r_probVector = (double *)malloc (nrComb*sizeof (double));
			double *a = (double *)malloc (nrComb*sizeof (double));
			double l_dist=executeCond(Qm, Qtransposed, n->left, l_probVector, MatrixA, MatrixB, MatrixC);
			double r_dist=executeCond(Qm, Qtransposed, n->right, r_probVector, MatrixA, MatrixB, MatrixC);
			
			condLikeFunction(a, Qm, Qtransposed, l_probVector, r_probVector, l_dist, r_dist, MatrixA, MatrixB,MatrixC);
			for (assign=0;assign<nrComb;assign++){ 
				n->probVector[assign]=a[assign];
				probVector[assign]=a[assign];
			}
			free (l_probVector);
			free (r_probVector);
			free (a);
			return 	n->dist; 
		}
	}
} 

void setQ(double **Qmatrix, double s, double d, double w1, double w2,int *VectCoevComb){
	int loop_k1, loop_k2, nrSymb;
	int loop_i , loop_j; //Loop counters
	int from_1, from_2, to_1, to_2;
	int boolValue=0;
	double sum_Row_Value;
	
	for (loop_k1=0;loop_k1<nrComb;loop_k1++){ 
		for (loop_k2=0;loop_k2<nrComb;loop_k2++){ 
			Qmatrix[loop_k1][loop_k2]=0;
		}
	}
	
	if (dataType == NT)
	{
		nrSymb = nrNt;
	}
	else
	{
		nrSymb = nrAa;
	}
	for (loop_k1=0;loop_k1<nrComb;loop_k1++){ 
		for (loop_k2=0;loop_k2<nrComb;loop_k2++){
			from_1 = loop_k1 / nrSymb;
			from_2 = loop_k1 % nrSymb;
			to_1 = loop_k2 / nrSymb;
			to_2 = loop_k2 % nrSymb;
			
			int ndiff = 0;      
			if (from_1 != to_1){
				ndiff = ndiff+1;
			}
			if (from_2 != to_2){
				ndiff = ndiff+1;
			}
			
			if (ndiff == 1) {// if is a single substitution
				//if ic1 belong to co-ev combinations then sm
				//if ic2 belong to co-ev combinations then dm
				//else w  
				boolValue=0;
				
				if (VectCoevComb[loop_k1]==1){
					Qmatrix[loop_k1][loop_k2]=s;//sm
					boolValue=1;
				}
				if (VectCoevComb[loop_k2]==1){
					Qmatrix[loop_k1][loop_k2]=d;//dm
					boolValue=1;
				}
				
				if (boolValue == 0){
					
					if (from_1 != to_1){
						Qmatrix[loop_k1][loop_k2]=w1;
					}
					if (from_2 != to_2){
						Qmatrix[loop_k1][loop_k2]=w2;
					}	
					
				}
			}
			if (ndiff == 2) {
				Qmatrix[loop_k1][loop_k2]=0;
			}
			if (ndiff == 0) {
				Qmatrix[loop_k1][loop_k2]=0;
			}
		}
	}
	
	/*
	 for (loop_i=0;loop_i<nrComb;loop_i++){ 
	 for (loop_j=0;loop_j<nrComb;loop_j++){ 
	 printf("%g ",  Q[loop_i][loop_j]);
	 }
	 printf("\n");
	 }
	 */
	
	// Do sum of each row
	// MATLAB: sum(diagonal(k,1:16));   
	for (loop_i=0;loop_i<nrComb;loop_i++){ 
		sum_Row_Value=0;
		for (loop_j=0;loop_j<nrComb;loop_j++){ 
			sum_Row_Value=sum_Row_Value+Qmatrix[loop_i][loop_j];
		}
		Qmatrix[loop_i][loop_i] = Qmatrix[loop_i][loop_i]-sum_Row_Value;
	}	
}

void setQNull(double **Qmatrix, double w1, double w2){
	int loop_k1, loop_k2, nrSymb;
	int loop_i , loop_j; //Loop counters
	int from_1, from_2, to_1, to_2;
	int boolValue=0;
	double sum_Row_Value;
	
	for (loop_k1=0;loop_k1<nrComb;loop_k1++){ 
		for (loop_k2=0;loop_k2<nrComb;loop_k2++){ 
			Qmatrix[loop_k1][loop_k2]=0;
		}
	}
	
	if (dataType == NT)
	{
		nrSymb = nrNt;
	}
	else
	{
		nrSymb = nrAa;
	}
	for (loop_k1=0;loop_k1<nrComb;loop_k1++){ 
		for (loop_k2=0;loop_k2<nrComb;loop_k2++){
			from_1 = loop_k1 / nrSymb;
			from_2 = loop_k1 % nrSymb;
			to_1 = loop_k2 / nrSymb;
			to_2 = loop_k2 % nrSymb;
			
			int ndiff = 0;      
			if (from_1 != to_1){
				ndiff = ndiff+1;
			}
			if (from_2 != to_2){
				ndiff = ndiff+1;
			}
			
			if (ndiff == 1) {
				
				
				if (from_1 != to_1){
	                Qmatrix[loop_k1][loop_k2]=w1;
				}
				if (from_2 != to_2){
	                Qmatrix[loop_k1][loop_k2]=w2;
				}
			}
			if (ndiff == 2) {
				Qmatrix[loop_k1][loop_k2]=0;
			}
			if (ndiff == 0) {
				Qmatrix[loop_k1][loop_k2]=0;
			}
		}
	}
	
	/*
	 for (loop_i=0;loop_i<nrComb;loop_i++){ 
	 for (loop_j=0;loop_j<nrComb;loop_j++){ 
	 printf("%g ",  Q[loop_i][loop_j]);
	 }
	 printf("\n");
	 }
	 */
	
	// Do sum of each row
	// MATLAB: sum(diagonal(k,1:16));   
	for (loop_i=0;loop_i<nrComb;loop_i++){ 
		sum_Row_Value=0;
		for (loop_j=0;loop_j<nrComb;loop_j++){ 
			sum_Row_Value=sum_Row_Value+Qmatrix[loop_i][loop_j];
		}
		Qmatrix[loop_i][loop_i] = Qmatrix[loop_i][loop_i]-sum_Row_Value;
	}	
}

void setQSimulator(double **Qmatrix, double s, double d, double r1, double r2,int *VectCoevComb)
{
    int loop_k1, loop_k2, nrSymb;
    int loop_i =0, loop_j; //Loop counters
    int from_1, from_2, to_1, to_2;
    int boolValue=0;
    double sum_Row_Value;
    for (loop_k1=0;loop_k1<nrComb;loop_k1++){
        for (loop_k2=0;loop_k2<nrComb;loop_k2++){
            Qmatrix[loop_k1][loop_k2]=0.0;
			// QmatrixTemp[loop_k1][loop_k2]=0.0;
        }
    }
	if (dataType == NT)
	{
		nrSymb = nrNt;
	}
	else
	{
		nrSymb = nrAa;
	}
    for (loop_k1=0;loop_k1<nrComb;loop_k1++){
        
        for (loop_k2=0;loop_k2<nrComb;loop_k2++){
            
			from_1 = loop_k1 / nrSymb;
			from_2 = loop_k1 % nrSymb;
			to_1 = loop_k2 / nrSymb;
			to_2 = loop_k2 % nrSymb;
			
			int ndiff = 0;      
			if (from_1 != to_1){
				ndiff = ndiff+1;
			}
			if (from_2 != to_2){
				ndiff = ndiff+1;
			}
			/*			
			 if (ndiff == 1) {
			 
			 
			 if (from_1 != to_1){
			 Qmatrix[loop_k1][loop_k2]=r1;
			 }
			 if (from_2 != to_2){
			 Qmatrix[loop_k1][loop_k2]=r2;
	         }
			 }*/
			
			if (ndiff == 1) {// if is a single substitution
				//if ic1 belong to co-ev combinations then sm
				//	//if ic2 belong to co-ev combinations then dm
				//		//else w  
				boolValue=0;
				if (VectCoevComb[loop_k1]==1){
					Qmatrix[loop_k1][loop_k2]=s;//sm
					boolValue=1;
				}
				if (VectCoevComb[loop_k2]==1){
					Qmatrix[loop_k1][loop_k2]=d;//dm
					boolValue=1;
				}
				if (boolValue == 0){							  	  			  	  							  	        if (from_1 != to_1){								  	  			  	  	 Qmatrix[loop_k1][loop_k2]=r1;	 }
					if (from_2 != to_2){
						Qmatrix[loop_k1][loop_k2]=r2; }      	 	  	        	                	        	        	                	         	 	}
			}
			if (ndiff == 2) {
				Qmatrix[loop_k1][loop_k2]=0;
			}
			if (ndiff == 0) {
				Qmatrix[loop_k1][loop_k2]=0;
			}
        }
		
        
    }
     //Without frequency vector 
     /*
     for (loop_i=0;loop_i<nrComb;loop_i++){ 
            for (loop_j=0;loop_j<nrComb;loop_j++){ 
                        printf("%g ",  Qmatrix[loop_i][loop_j]);
            }
            printf("\n");
     }
     printf("\n\n");
     */
    
    
    for (loop_i=0;loop_i<nrComb;loop_i++){ 
        sum_Row_Value=0;
        for (loop_j=0;loop_j<nrComb;loop_j++){ 
            sum_Row_Value=sum_Row_Value+Qmatrix[loop_i][loop_j];
        }
        Qmatrix[loop_i][loop_i] = Qmatrix[loop_i][loop_i]-sum_Row_Value;
    }
    
    
}
void pVectorFunction(double *pVector, double **Qmatrix, double **QTmatrix, double ti, double *MatrixA, double *MatrixB, double *MatrixC){
    int loop_i=0, loop_j=0;
    /*fprintf (stderr, "CHECK pVectorFunction(), dist=%f\n",ti);
     for (loop_i=0;loop_i<nrComb;loop_i++){
     for (loop_j=0;loop_j<nrComb;loop_j++){
     fprintf (stderr, "%g,",  Qmatrix[loop_i][loop_j]);
     }
     }
     fprintf (stderr, "\n");
     for (loop_i=0;loop_i<nrComb;loop_i++){
     for (loop_j=0;loop_j<nrComb;loop_j++){
     fprintf (stderr, "%g,",  QTmatrix[loop_i][loop_j]);
     }
     }
     */
    ////////////////////////////////////////////////////////////////////////////////
    // EIGEN VALUE & EIGEN VECTORS//////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////
    //MATLAB: [eigenVector,eigenValue]=eig(gQ);
    //Compute egen value and eigen vector in c with dgeev
    
    //ARGUMENTS
    int N=nrComb;
    int LDA=nrComb;
    double *A= (double *)malloc(LDA*N*sizeof(double));
    for (loop_i = 0; loop_i < nrComb; loop_i++) {
        for (loop_j = 0; loop_j < nrComb; loop_j++) {
            A[ loop_i* nrComb + loop_j]=QTmatrix[loop_i][loop_j];
        }
    }
    //A=&(QTmatrix[0][0]);// USE TRANSPOSE MATRIX SINCE THOSE METHODS ARE COLUMN BASED
    double *WR=(double *)malloc(nrComb*sizeof(double));
    double *WI=(double *)malloc(nrComb*sizeof(double));
    int LDVL=nrComb;
    double *VL=(double *)malloc(LDVL*nrComb*sizeof(double));
    int LDVR=nrComb;
    //int LWORK=4*nrComb;
    int LWORK=nrComb*nrComb;
    double *VR=(double *)malloc(LDVR*nrComb*sizeof(double));
    double *WORK=(double *)malloc(LWORK*sizeof(double));
    int INFO=0;
    
    //CALL OF FUNCTION//////////////////////////////////////////////////////////////
    LWORK=-1;
    dgeev_("N", "V", &N, A, &LDA, WR, WI, VL, &LDVL, VR, &LDVR, WORK, &LWORK, &INFO );
    LWORK=(int) WORK[0];
    for (loop_i = 0; loop_i < nrComb; loop_i++) {
        for (loop_j = 0; loop_j < nrComb; loop_j++) {
            A[ loop_i* nrComb + loop_j]=QTmatrix[loop_i][loop_j];
        }
    }
    dgeev_("N", "V", &N, A, &LDA, WR, WI, VL, &LDVL, VR, &LDVR, WORK, &LWORK, &INFO );
    //SaveToOctave (&(Qmatrix[0][0]), "Q", fp, 16, 16) ;
    //SaveToOctave (&(QTmatrix[0][0]), "QT", fp, 16, 16) ;
    //SaveToOctave (WR, "evalues", fp, 1, 16) ;
    //SaveToOctave (VR, "evectors", fp, 16, 16) ;
    
    //freeArray A
    free(A);
    
    free(WI);free(VL);free(WORK);
	
    ///SAVE VR
    double *VR_SAV1 = (double *)malloc (nrComb*nrComb*sizeof (double));
    double *VR_SAV2 = (double *)malloc (nrComb*nrComb*sizeof (double));
	
    //From pointer to [][]
    for (loop_j = 0; loop_j < nrComb*nrComb; loop_j++)//columns
    {
        VR_SAV1[loop_j]= VR[loop_j];
        VR_SAV2[loop_j]= VR[loop_j];
    }
	
    ////////////////////////////////////////////////////////////////////////////////
    // Inverse Of Eigen vector//////////////////////////////////////////////////////
    // MATLAB:inv(eigenVector)
    ////////////////////////////////////////////////////////////////////////////////
    //STEP1: FACTORISATION//////////////////////////////////////////////////////////
    //ARGUMENTS/////////////////////////////////////////////////////////////////////
    int N2=nrComb;
    int LDA2=nrComb;
    double *VR_inv=(double *)malloc(LDA2*N2*sizeof(double));
    
    for (loop_i = 0; loop_i < LDA2*N2; loop_i++)
    {
        VR_inv[loop_i] = VR[loop_i];
    }
    int *IPIV2=(int *)malloc(nrComb*sizeof(int));
    int row=nrComb;
    int column=nrComb;
    int INFO2=0;
    
    //CALL OF FUNCTION//////////////////////////////////////////////////////////////
    dgetrf_(&row, &column, VR_inv, &LDA2, IPIV2, &INFO2 );
    //printf("INFO2: %d, output2:%d\n", INFO2, output2);
    //STEP2: INVERSE////////////////////////////////////////////////////////////////
    //ARGUMENTS/////////////////////////////////////////////////////////////////////
    int LWORK2=4*nrComb;
    double *WORK2=(double *)malloc(LWORK2*sizeof(double));
    int INFO3=0;
    LWORK2=-1;
    
    //CALL OF FUNCTION//////////////////////////////////////////////////////////////
    dgetri_(&N2, VR_inv, &LDA2, IPIV2, WORK2, &LWORK2, &INFO3);
    LWORK2=(int) WORK2[0];
    free(WORK2);
    WORK2=(double *)malloc(LWORK2*sizeof(double));
    dgetri_(&N2, VR_inv, &LDA2, IPIV2, WORK2, &LWORK2, &INFO3 );
    //VR_inv is a transposed matrix
    //SaveToOctave (VR_inv, "VR_invFINALRESULT", fp, 16, 16) ;
	
    //printf("INFO3: %d, output3:%d\n", INFO3, output3);
    ////////////////////////////////////////////////////////////////////////////////
    // PREPARE EXPRESSION EVALUATION////////////////////////////////////////////////
    // USE TRANSPOSE MATRIX SINCE THOSE METHODS RETURN COLUMN BASED MATRICES
    ////////////////////////////////////////////////////////////////////////////////
    //--VR_invInput	:eigenVectors: transpose of dgetri output VR_inv
    //--WR_input	        :eigenValues : transpose of dgeev  output WR
	
    //From pointer to [][]
    double *VR_invDgetri = (double *)malloc (nrComb*nrComb*sizeof (double));
    for (loop_j = 0; loop_j < nrComb*nrComb; loop_j++){//columns of the file
        VR_invDgetri[loop_j]= VR_inv[loop_j];
    }
    //Transpose the matrix VR_inv matrix
    //double VR_invInput[nrComb][nrComb];  // Doesn't seem to be used anywhere!
	
    //create a eigen diagonal matrix
    //Initialise matrix
    double **WR_diagEXP_ti = (double **)malloc (nrComb * sizeof (double *));
    for (loop_i=0;loop_i<nrComb;loop_i++)
    {
        WR_diagEXP_ti[loop_i] = (double *)malloc (nrComb * sizeof (double));
        for (loop_j=0;loop_j<nrComb;loop_j++){
            WR_diagEXP_ti[loop_i][loop_j]=0;
        }
    }
	
    //MATLAB: expm(eigenValue*ti)
    for (loop_j=0;loop_j<nrComb;loop_j++){
        WR_diagEXP_ti[loop_j][loop_j]=exp(WR[loop_j] * ti);
    }
    ////////////////////////////////////////////////////////////////////////////////
    // MULTIPLY MATRICES ///////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////
    //MATLAB	:pii(1:16,1:16) = eigenVector*expm(eigenValue*ti)*inv(eigenVector);
    //OUTPUT DGEMMM: WE OBTAIN PI  TRANSPOSED
    ////////////////////////////////////////////////////////////////////////////////
    //ARGUMENTS
    int  M=nrComb;
    N=nrComb;
    int K=nrComb;
	
    double ALPHA=1.0;
    int LDA_dgemm=nrComb;
    int LDB_dgemm=nrComb;
    double BETA=0.0;
    int LDC_dgemm=nrComb;
    int output_dgemm=0;
    
    //CALL OF FUNCTION DGEMM: RES1= expm(eigenValue*ti)*inv(eigenVector)/////////////
    //MatrixA: is the diagonal matrix WR_diagEXP_ti
    //MatrixA=&(WR_diagEXP_ti[0][0]);
    for (loop_i = 0; loop_i < nrComb; loop_i++) {
        for (loop_j = 0; loop_j < nrComb; loop_j++) {
            MatrixA[ loop_i* nrComb + loop_j]=WR_diagEXP_ti[loop_i][loop_j];
        }
    }
    //MatrixB=&(VR_invDgetri[0]);
    for (loop_i = 0; loop_i < nrComb*nrComb; loop_i++)
    {
        MatrixB[loop_i] = VR_invDgetri[loop_i];
    }
    for (loop_i = 0; loop_i < nrComb; loop_i++)//rows of the file
        for (loop_j = 0; loop_j < nrComb; loop_j++)//columns of the file
            MatrixC[loop_i * nrComb + loop_j]=0.0;
    output_dgemm=dgemm_("N","N", &M, &N, &K, &ALPHA, MatrixA, &LDA_dgemm, MatrixB, &LDB_dgemm, &BETA, MatrixC, &LDC_dgemm);
	
    //CALL OF FUNCTION DGEMM: pii=eigenVector*RES1////////////////////////////////////
    //MatrixA: is the transpose of VR_invInput: (that is the original VR_inv or VR_invDgetri)
    //MatrixB: is MatrixC
    //MatrixA= &(VR_SAV1[0]);
    //MatrixB= MatrixC;
    for (loop_i = 0; loop_i < nrComb*nrComb; loop_i++) {
        MatrixA[loop_i]=VR_SAV1[loop_i];
        MatrixB[loop_i] = MatrixC[loop_i];
    }
    double *pii=(double *)malloc(M*N*sizeof(double));//OUTPUT pii is a Matrix!!!
    for (loop_i = 0; loop_i < nrComb; loop_i++)//rows of the file
        for (loop_j = 0; loop_j < nrComb; loop_j++)//columns of the file
            pii[loop_i * nrComb + loop_j]=0.0;
    output_dgemm=dgemm_("N","N", &M, &N, &K, &ALPHA, MatrixA, &LDA_dgemm, MatrixB, &LDB_dgemm, &BETA, pii, &LDC_dgemm);
    
    for (loop_i=0;loop_i<nrComb;loop_i++){
        for (loop_j=0;loop_j<nrComb;loop_j++){
            pVector[loop_i * nrComb + loop_j]= pii[loop_i * nrComb + loop_j] ;
            
        }
    }
    
    free (pii);
    free(WR);
    free(VR_inv);
    free (VR_invDgetri);
    free(IPIV2);
    free (VR);
    free (VR_SAV1);
    free (VR_SAV2);
    for (loop_i=0; loop_i < nrComb; loop_i++)
    {
        free (WR_diagEXP_ti[loop_i]);
    }
    free (WR_diagEXP_ti);
    free(WORK2);
    
}

