/*
** bayes: Apply the Baysian method.
**
** Linda Dib & Wim Hordijk   Last modified: 26 April 2013
*/

#include "def.h"
#include "model.h"
#include "strmap.h"


//%%%%%%%%%%%%%%%%%%%%%%%%%%% COMBINATION FUNCTION %%%%%%%%%%%%%%%%%%%%%%%


void addComb (int *array)
{
  int present[nrAa], possible[nrAaComb], i, j, nrPresent, nrPossible, rnd;

  /*
  ** Get the indices of the combinations already present.
  */
  nrPresent = 0;
  for (i = 0; i < nrComb; i++)
  {
    if (array[i] == 1)
    {
      present[nrPresent] = i;
      nrPresent++;
    }
  }
  
  /*
  ** Get all the possible combinations for addition.
  */
  nrPossible = 0;
  for (i = 0; i < nrComb; i++)
  {
    if (array[i] == 0)
    {
      for (j = 0; j < nrPresent; j++)
      {
	if (isConflict (i, present[j]) != -1)
	{
	  break;
	}
      }
      if (j >= nrPresent)
      {
	possible[nrPossible] = i;
	nrPossible++;
      }
    }
  }

  /*
  ** Select one of the possible combinations at random.
  */
  if (nrPossible > 0)
  {
    rnd = random () % nrPossible;
    array[possible[rnd]] = 1;
  }
}


void removeComb (int *array)
{ 
  int present[nrAa], i, nrPresent, rnd;

  /*
  ** Get the indices of the combinations already present.
  */
  nrPresent = 0;
  for (i = 0; i < nrComb; i++)
  {
    if (array[i] == 1)
    {
      present[nrPresent] = i;
      nrPresent++;
    }
  }

  /*
  ** Select one of the present combinations at random.
  */
  if (nrPresent > 2)
  {
    rnd = random () % nrPresent;
    array[present[rnd]] = 0;
  }
} 


void changeComb (int *array)
{ 
  int i, nrPresent;

  /*
  ** Get the number of combinations already present.
  */
  nrPresent = 0;
  for (i = 0; i < nrComb; i++)
  {
    if (array[i] == 1)
    {
      nrPresent++;
    }
  }

  /*
  ** Change one combination.
  */
  if (nrPresent == 2)
  {
    addComb (array);
    removeComb (array);
  }
  else if (nrPresent > 2)
  {
    removeComb (array);
    addComb (array);
  }
} 

void pickCCS(int *myProfile){
  int i, j;
  
  i = random () % nrProfiles;
  for (j = 0; j < nrComb; j++)
  {
    myProfile[j] = profiles[i][j];
  }
}


//%%%%%%%%%%%%%%%%%%%%%%%%%%% MCMC FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%
//i: float parm value
//d: size of the uniform window (this needs to be adjusted after testing, and separately for each parameter. You can set it to 1 for now)
void update_parameter(double *i, double win){
  double x;

  x = random()/(RAND_MAX+1.0);
  *i = fabs(*i+(x-0.5)*win);
	
}

void update_parameterPN(double *i, double win){
  double x;

  x = random()/(RAND_MAX+1.0);
  *i = *i+(x-0.5)*win;// allowing for negative values
	
}

void update_parameter_multi(double *i, double win, double* u){
double m;
double lambda = 2*log(win);
 m = exp(lambda*(*u-0.5));
 *i = *i * m;
}

void ChangeOperation(int *VectCoevComb, double *S, double *D, double *W1, double *W2,double *U){
	int loop_j=0;
	
	double r;
	r = random()/(RAND_MAX+1.0);
	if (r < 0.25){
		int op=random()%(4) ;
		//change S
		if(op==0){
			update_parameter_multi(S,2,U);
		}
		//change D
		if(op==1){
			update_parameter_multi(D,2,U); 
		}
		//change W1
		if(op==2){
			update_parameter_multi(W1,2,U);
		}
		//change W2
		if(op==3){
			update_parameter_multi(W2,2,U);
		}
	}else{
		 pickCCS(VectCoevComb); 
		*U=1.0;
	}
}

void *threadFunc(void *arg){
	int loop_i, loop_j, i, j, iteration;
	
	double *a=(double *)malloc (nrComb * sizeof (double));
	  for (loop_j=0;loop_j<nrComb;loop_j++){
              a[loop_j]=0;
          }
	
	double **Q = (double **)malloc (nrComb * sizeof (double *));
          for (loop_i=0;loop_i<nrComb;loop_i++)
          {
            Q[loop_i] = (double *)malloc (nrComb * sizeof (double));
            for (loop_j=0;loop_j<nrComb;loop_j++){
              Q[loop_i][loop_j]=0;
            }
          }
          
        double **Qtransposed = (double **)malloc (nrComb * sizeof (double *));
          for (loop_i=0;loop_i<nrComb;loop_i++)
          {
            Qtransposed[loop_i] = (double *)malloc (nrComb * sizeof (double));
            for (loop_j=0;loop_j<nrComb;loop_j++){
              Qtransposed[loop_i][loop_j]=0;
            }
          }
	
          
	
	double * func_data=(double*) arg;
	double temperature= func_data[0];
	int  fileNameOutIndex=func_data[1];
	
	//char fileNameOut[128];
	//char tampon [16] ;
	//sprintf (tampon, "%d", fileNameOutIndex) ;
	//strcat(fileNameOut, (char*)tampon);
	
	double *AS    = (double *)malloc(sizeof(double));
	double *AD    = (double *)malloc(sizeof(double));
	double *AW1   = (double *)malloc(sizeof(double));
	double *AW2   = (double *)malloc(sizeof(double));
	
	*AS   = func_data[2];
	*AD   = func_data[3];
	*AW1  = func_data[4];
	*AW2  = func_data[5];
	
	double ALikelihood    = 0;//func_data[6]
	double ALogLikelihood = 0;//func_data[7]
	double APrior=0.0;        //func_data[8]
	int *AVectorComb= (int *)malloc(nrComb*sizeof(int));
	for (loop_j=0;loop_j<nrComb;loop_j++){ 
        AVectorComb[loop_j]=func_data[loop_j+9];
	}	
	
	
	
	/////////////////////FIRST EVALUATION CALCULATE PROIOR LIKELIHOOD////////////////////
	// setQ(Q,*AS,*AD,exp(*AW),AVectorComb);
	setQ(Q,*AS,*AD,*AW1, *AW2, AVectorComb);

	transposeMatrix(Q,Qtransposed);

	double *MatrixA=(double *)malloc(nrComb*nrComb*sizeof(double));//INPUT
	double *MatrixB=(double *)malloc(nrComb*nrComb*sizeof(double));//INPUT
	double *MatrixC=(double *)malloc(nrComb*nrComb*sizeof(double));//OUTPUT
	for (loop_j=0;loop_j<nrComb;loop_j++){ a[loop_j]=0;}
	
	executeCond(Q, Qtransposed, &root, a, MatrixA, MatrixB, MatrixC);
	
	ALikelihood=0.0;
	for (loop_j=0;loop_j<nrComb;loop_j++){ 
		ALikelihood=ALikelihood+a[loop_j];
	}
	ALogLikelihood =log(ALikelihood);
	
	alpha=1.1;
	beta=0.3;
	
	APrior= 
	(alpha-1)*log(*AS) +(-beta * *AS)  +
	(alpha-1)*log(*AD) +(-beta * *AD)  +
	(alpha-1)*log(*AW1)+(-beta * *AW1) +
	(alpha-1)*log(*AW2)+(-beta * *AW2);

	//FILE *fp;
	//fp = fopen(fileNameOut, "a+");
	//printf ("%g\t%g\t%g\t%g\t%g\t%g\t\n",*BS,*BD,*BW, AS, AD, AW);
	/////////////////////END   ////////////////////////////////////////
	double* BS   =(double *)malloc(sizeof(double));
	double* BD   =(double *)malloc(sizeof(double));
	double* BW1   =(double *)malloc(sizeof(double));
	double* BW2   =(double *)malloc(sizeof(double));
	
	*BS   = *AS;
	*BD   = *AD;
	*BW1   = *AW1;
	*BW2   = *AW2;
	
	
	double BLikelihood   = 0.0;
	double BLogLikelihood= 0.0;
	double BPrior=0.05;
	
	//fprintf (fp,"iteration\tPosterior\tLogLikelihood\tPrior\tS\tD\tW1\tW2\tAA\tAC\tAG\tAT\tCA\tCC\tCG\tCT\tGA\tGC\tGG\tGT\tTA\tTC\tTG\tTT\n");
	iteration=0;
	int *BVectorComb= (int *)malloc(nrComb*sizeof(int));
	double* U=(double *)malloc(sizeof(double));
	

	while(iteration<innerLoopIteration){
		///////////////////////////////////////////////////////////////////////////////
		
		*BS   = *AS;
		*BD   = *AD;
		*BW1   = *AW1;
		*BW2   = *AW2;
		
		BLikelihood=0.0;
		BLogLikelihood=0.0;
		BPrior=0.05;
    	
		for (loop_j=0;loop_j<nrComb;loop_j++){ BVectorComb[loop_j]=AVectorComb[loop_j];}
		*U = random()/(RAND_MAX+1.0);
		ChangeOperation(BVectorComb, BS, BD, BW1,BW2, U);
		
		setQ(Q,*BS,*BD,*BW1,*BW2,BVectorComb);
		
		
		transposeMatrix(Q,Qtransposed);
		for (loop_j=0;loop_j<nrComb;loop_j++){ a[loop_j]=0.0;}
		for (loop_j=0;loop_j<nrComb;loop_j++){ 
			for (loop_i=0;loop_i<nrComb;loop_i++){ 
				MatrixA[loop_i * nrComb + loop_j]=0.0;
				MatrixB[loop_i * nrComb + loop_j]=0.0;
				MatrixC[loop_i * nrComb + loop_j]=0.0;
			}
		}
		executeCond(Q, Qtransposed, &root, a, MatrixA,MatrixB,MatrixC);
		BLikelihood=0.0;
		for (loop_j=0;loop_j<nrComb;loop_j++){
			BLikelihood=BLikelihood+a[loop_j];
		}
		BLogLikelihood =log(BLikelihood);
		alpha=1.1;
		beta=0.3;
		
		if(*BS >= *BD){
			BPrior=-INFINITY;
		}else{
			BPrior = 
			(alpha-1)*log(*BS)+(-beta *  *BS)  +
			(alpha-1)*log(*BD)+(-beta *  *BD)  +
			(alpha-1)*log(*BW1)+(-beta * *BW1) +
			(alpha-1)*log(*BW2)+(-beta * *BW2);
		}
		
		///////////////////////////////////////////////////////////////////////////////
		// 						Condition of acceptance	     //
		///////////////////////////////////////////////////////////////////////////////	
		double x = random()/(RAND_MAX+1.0);
		if ( (((BLogLikelihood+BPrior) - (ALogLikelihood+APrior))*temperature)+log(*U) >= log(x)){
			//if (BLogLikelihood - ALogLikelihood  >= log(x)){
			*AS   = *BS;
			*AD   = *BD;
			*AW1   = *BW1;
			*AW2   = *BW2;
			
			ALikelihood   = BLikelihood;
			ALogLikelihood= BLogLikelihood;
			APrior= BPrior;
			for (loop_j=0;loop_j<nrComb;loop_j++){ 
				AVectorComb[loop_j]=BVectorComb[loop_j];
			}
		}else{
			//printf ("ERROR%g\t%g\n", ALogLikelihood, BLogLikelihood);
		}
		
		
		double post=ALogLikelihood+APrior;
		//fprintf (fp,"%d\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t", iteration,post,ALogLikelihood, APrior, *AS, *AD, *AW1,*AW2);
		for (loop_j=0;loop_j<nrComb;loop_j++){ 
			//fprintf (fp,"%d\t",AVectorComb[loop_j]);
		}
		//fprintf (fp,"\n");
		//fflush(fp);	
		iteration++;
	}
	//fclose (fp);
	
	
	func_data[2]=*AS;
	func_data[3]=*AD;
	func_data[4]=*AW1;
	func_data[5]=*AW2;
	func_data[6]=ALikelihood;	
	func_data[7]=ALogLikelihood;	
	func_data[8]=APrior;
	
	for (loop_i = 9; loop_i < nrComb+9; loop_i++){ 
		func_data[loop_i] = AVectorComb[loop_i-9];
	}
	
	free(U);
	free(BVectorComb);
	free(AVectorComb);
	
	free(AS);
	free(AD);
	free(AW1);
	free(AW2);
	
	free(BS);
	free(BD);
	free(BW1);
	free(BW2);
	
	
	free(MatrixA);
	free(MatrixB);
	free(MatrixC);
	
	free (a);
        for (i = 0; i < nrComb; i++)
        {
         free (Q[i]);
         free (Qtransposed[i]);
        }
        free (Q);
        free (Qtransposed);  
	pthread_exit(func_data);
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MAIN FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%
int bayes (){
        int    status, num_chars, assign, len, loop_i, loop_j, i, j;
	double th1S  =0.0, th1D=0.0, th1W1=0.0, th1W2  =0.0;
	double th2S  =0.0, th2D=0.0, th2W1=0.0, th2W2  =0.0;
	double th3S  =0.0, th3D=0.0, th3W1=0.0, th3W2  =0.0;
	double th4S  =0.0, th4D=0.0, th4W1=0.0, th4W2  =0.0;
	double tmpS  =0.0, tmpD=0.0, tmpW1=0.0, tmpW2  =0.0;
	double tmpLikelihood=0.0,tmpLogLikelihood=0.0,tmpPrior=0.0;
	
	double temperature1=1.0;
	double temperature2=0.5;
	double temperature3=0.6;
	double temperature4=0.7;
	FILE *fp;
	
	innerLoopIteration=20;
	/*
	 ** Create the nucleotide combinations.
	 */
	  fp = fopen(outFile, "w+");
        fprintf (fp,"iteration\tPosterior\tLogLikelihood\tPrior\tS\tD\tR1\tR2\tAA\tAC\tAG\tAT\tCA\tCC\tCG\tCT\tGA\tGC\tGG\tGT\tTA\tTC\tTG\tTT\n");
	 fflush(fp);
	/******* FOR 4 CHAINS***********/
	//////THREAD1
	th1S   = s; th1D   = d; th1W1  = r1; th1W2  = r2;
	int *th1VectorComb = (int *)malloc (nrComb*sizeof (int));
	int th1Comb= rand() % nrProfiles;
        for (loop_j=0;loop_j<nrComb;loop_j++){ 
		th1VectorComb[loop_j]=profiles[th1Comb][loop_j];
        }
	
	//////THREAD2
	th2S   = s;th2D   = d;th2W1  = r1;th2W2  = r2;
	int *th2VectorComb = (int *)malloc (nrComb*sizeof (int));
	int th2Comb= rand() % nrProfiles;
        for (loop_j=0;loop_j<nrComb;loop_j++){ 
		th2VectorComb[loop_j]=profiles[th2Comb][loop_j];
        }
	
	//////THREAD3
	th3S   = s;th3D   = d;th3W1  = r1;th3W2  = r2;
	int *th3VectorComb = (int *)malloc (nrComb*sizeof (int));
	int th3Comb= rand() % nrProfiles;
        for (loop_j=0;loop_j<nrComb;loop_j++){ 
		th3VectorComb[loop_j]=profiles[th3Comb][loop_j];
        }
	
	
	//////THREAD4
	th4S   = s;th4D   = d;th4W1  = r1;th4W2  = r2;
	int *th4VectorComb = (int *)malloc (nrComb*sizeof (int));
	int th4Comb= rand() % nrProfiles;
        for (loop_j=0;loop_j<nrComb;loop_j++){ 
		th4VectorComb[loop_j]=profiles[th4Comb][loop_j];
        }
	
	//////FOR SWAP
	tmpS   = s;tmpD   = d;tmpW1  = r1;tmpW2  = r2;
	int *tmpVectorComb = (int *)malloc (nrComb*sizeof (int));
	for (loop_i = 0; loop_i < nrComb; loop_i++){ 	tmpVectorComb[loop_i] = 0;}
	
       
	
	/****************************************************************/
	/*                            MULTITHREADING                    */
	/****************************************************************/
	int iteration=0;
	double* argThread1= (double*)malloc((9+nrComb)*sizeof(double));
	double* argThread2= (double*)malloc((9+nrComb)*sizeof(double));
	double* argThread3= (double*)malloc((9+nrComb)*sizeof(double));
	double* argThread4= (double*)malloc((9+nrComb)*sizeof(double));
	pthread_t pth1, pth2, pth3, pth4;
	
	//printf("DEBUG: Thread initialized\n");
	while(iteration<IT){
		///THREAD1
		argThread1[0]=temperature1; //temperature
		argThread1[1]=1;            //threadID
		argThread1[2]=th1S;
		argThread1[3]=th1D;
		argThread1[4]=th1W1;
		argThread1[5]=th1W2;
		argThread1[6]=-INFINITY;
		argThread1[7]=-INFINITY;
		argThread1[8]=-INFINITY;
		for (loop_i = 9; loop_i < nrComb+9; loop_i++){ argThread1[loop_i] =  th1VectorComb[loop_i-9];}
		pthread_create(&pth1,NULL,threadFunc,argThread1);
		
		///THREAD2
		argThread2[0]=temperature2; //temperature
		argThread2[1]=2;            //threadID
		argThread2[2]=th2S;
		argThread2[3]=th2D;
		argThread2[4]=th2W1;
		argThread2[5]=th2W2;
		argThread2[6]=-INFINITY;
		argThread2[7]=-INFINITY;
		argThread2[8]=-INFINITY;
		for (loop_i = 9; loop_i < nrComb+9; loop_i++){ argThread2[loop_i] =  th2VectorComb[loop_i-9];}
		pthread_create(&pth2,NULL,threadFunc,argThread2);
		
		
		///THREAD3
		argThread3[0]=temperature3; //temperature
		argThread3[1]=3;           //threadID
		argThread3[2]=th3S;
		argThread3[3]=th3D;
		argThread3[4]=th3W1;
		argThread3[5]=th3W2;
		argThread3[6]=-INFINITY;
		argThread3[7]=-INFINITY;
		argThread3[8]=-INFINITY;
		for (loop_i = 9; loop_i < nrComb+9; loop_i++){ argThread3[loop_i] =  th3VectorComb[loop_i-9];}
		pthread_create(&pth3,NULL,threadFunc,argThread3);
		
		///THREAD4
		argThread4[0]=temperature4; //temperature
		argThread4[1]=4;           //threadID
		argThread4[2]=th4S;
		argThread4[3]=th4D;
		argThread4[4]=th4W1;
		argThread4[5]=th4W2;
		argThread4[6]=-INFINITY;
		argThread4[7]=-INFINITY;
		argThread4[8]=-INFINITY;
		for (loop_i = 9; loop_i < nrComb+9; loop_i++){ argThread4[loop_i] =  th4VectorComb[loop_i-9];}
		pthread_create(&pth4,NULL,threadFunc,argThread4);
		
		//printf("DEBUG: Thread created\n");
		///WAITING FOR DIFFERENT THREADS
		pthread_join(pth1, NULL /* void ** return value could go here */);
		pthread_join(pth2, NULL /* void ** return value could go here */);
                pthread_join(pth3, NULL /* void ** return value could go here */);
                pthread_join(pth4, NULL /* void ** return value could go here */);
                
                //printf("DEBUG: Waiting for Thread \n");
                //GET VALUES THREAD1
		th1S=argThread1[2];th1D=argThread1[3];th1W1=argThread1[4];th1W2=argThread1[5];
		double th1Likelihood=argThread1[6];
		double th1LogLikelihood=argThread1[7];
		double th1Prior=argThread1[8];
		for (loop_i = 9; loop_i < nrComb+9; loop_i++){ th1VectorComb[loop_i-9]=(int)argThread1[loop_i];}

                //GET VALUES THREAD2
		th2S=argThread2[2];th2D=argThread2[3];th2W1=argThread2[4];th2W2=argThread2[5];
		double th2Likelihood=argThread2[6];
		double th2LogLikelihood=argThread2[7];
		double th2Prior=argThread2[8];
		for (loop_i = 9; loop_i < nrComb+9; loop_i++){ th2VectorComb[loop_i-9]=(int)argThread2[loop_i];}
                
                //GET VALUES THREAD3
		th3S=argThread3[2];th3D=argThread3[3];th3W1=argThread3[4];th3W2=argThread3[5];
		double th3Likelihood=argThread3[6];
		double th3LogLikelihood=argThread3[7];
		double th3Prior=argThread3[8];
		for (loop_i = 9; loop_i < nrComb+9; loop_i++){ th3VectorComb[loop_i-9]=(int)argThread3[loop_i];}

		//GET VALUES THREAD4
		th4S=argThread4[2];th4D=argThread4[3];th4W1=argThread4[4];th4W2=argThread4[5];
		double th4Likelihood=argThread4[6];
		double th4LogLikelihood=argThread4[7];
		double th4Prior=argThread4[8];
		for (loop_i = 9; loop_i < nrComb+9; loop_i++){ th4VectorComb[loop_i-9]=(int)argThread4[loop_i];}

		 //printf("DEBUG:Got values of Threads \n");
		/*********************************SWAPPING CHAINS******************************/
		///////////////////////////////////////////////////////////////////////////////
		//Condition of acceptance
		///////////////////////////////////////////////////////////////////////////////
		double x = random()/(RAND_MAX+1.0);

		double equationValue1=(((th2LogLikelihood+th2Prior)*temperature1) + ((th1LogLikelihood+th1Prior)*temperature2)) - 
		     (((th1LogLikelihood+th1Prior)*temperature1) + ((th2LogLikelihood+th2Prior)*temperature2));
		
		double equationValue2=(((th3LogLikelihood+th3Prior)*temperature1) + ((th1LogLikelihood+th1Prior)*temperature3)) - 
		     (((th1LogLikelihood+th1Prior)*temperature1) + ((th3LogLikelihood+th3Prior)*temperature3));
		
		double equationValue3=(((th4LogLikelihood+th4Prior)*temperature1) + ((th1LogLikelihood+th1Prior)*temperature4)) - 
		     (((th1LogLikelihood+th1Prior)*temperature1) + ((th4LogLikelihood+th4Prior)*temperature4));
		
		int op=random()%(3) ;
		// printf("\nOP_RANDOM %d\n",op);
		switch(op){
		case 0:
		            //printf("\Here 0: %f\n",equationValue1);
		        if ( equationValue1 >= log(x)){
                                    //printf("\nSWAP 0\n");
			tmpS   = th1S;tmpD   = th1D;tmpW1   = th1W1;tmpW2   = th1W2;
			tmpLikelihood   = th1Likelihood;tmpLogLikelihood= th1LogLikelihood;tmpPrior= th1Prior;
			for (loop_j=0;loop_j<nrComb;loop_j++){tmpVectorComb[loop_j]=th1VectorComb[loop_j];}
			
			th1S   = th2S;th1D   = th2D; th1W1   = th2W1; th1W2   = th2W2;
			th1Likelihood   = th2Likelihood;th1LogLikelihood= th2LogLikelihood;th1Prior= th2Prior;
			for (loop_j=0;loop_j<nrComb;loop_j++){ th1VectorComb[loop_j]=th2VectorComb[loop_j];}
			
			
			th2S   = tmpS;th2D   = tmpD;th2W1   =tmpW1;th2W2   = tmpW2;
			th2Likelihood   = tmpLikelihood;th2LogLikelihood= tmpLogLikelihood;th2Prior= tmpPrior;
			for (loop_j=0;loop_j<nrComb;loop_j++){th2VectorComb[loop_j]=tmpVectorComb[loop_j];}
			//  printf("\nEND SWAP 0\n");
                        }else{
                                        
                        }
                        break;
		case 1:
		            //printf("\Here 1: %f\n",equationValue2);
		        if ( equationValue2 >= log(x)){
                                   //printf("\nSWAP 1\n");
			tmpS   = th1S;tmpD   = th1D;tmpW1   = th1W1;tmpW2   = th1W2;
			tmpLikelihood   = th1Likelihood;tmpLogLikelihood= th1LogLikelihood;tmpPrior= th1Prior;
			for (loop_j=0;loop_j<nrComb;loop_j++){ tmpVectorComb[loop_j]=th1VectorComb[loop_j];}
			
			th1S   = th3S;th1D   = th3D;th1W1   = th3W1;th1W2   = th3W2;
			th1Likelihood   = th3Likelihood;th1LogLikelihood= th3LogLikelihood;th1Prior= th3Prior;
			for (loop_j=0;loop_j<nrComb;loop_j++){ th1VectorComb[loop_j]=th3VectorComb[loop_j];}
			
			th3S   = tmpS;th3D   = tmpD;th3W1   =tmpW1;th3W2   = tmpW2;
			th3Likelihood   = tmpLikelihood;th3LogLikelihood= tmpLogLikelihood;th3Prior= tmpPrior;
			for (loop_j=0;loop_j<nrComb;loop_j++){ th3VectorComb[loop_j]=tmpVectorComb[loop_j];}
			//  printf("\nEND SWAP 1\n");
			
                        }else{
                         //printf("\nKEEP ON\n");
                        }
		        break;
		
		case 2:
		            //printf("\Here 2: %f\n",equationValue3);
		        if ( equationValue3 >= log(x)){
                                    //printf("\nSWAP 2\n");
			tmpS   = th1S;tmpD   = th1D;tmpW1   = th1W1;tmpW2   = th1W2;
			tmpLikelihood   = th1Likelihood;tmpLogLikelihood= th1LogLikelihood;tmpPrior= th1Prior;
			for (loop_j=0;loop_j<nrComb;loop_j++){ tmpVectorComb[loop_j]=th1VectorComb[loop_j];}
			
			th1S   = th4S;th1D   = th4D;th1W1   = th4W1;th1W2   = th4W2;
			th1Likelihood   = th4Likelihood;th1LogLikelihood= th4LogLikelihood;th1Prior= th4Prior;
			for (loop_j=0;loop_j<nrComb;loop_j++){ th1VectorComb[loop_j]=th4VectorComb[loop_j];}
			
			th4S   = tmpS;th4D   = tmpD;th4W1   =tmpW1;th4W2   = tmpW2;
			th4Likelihood   = tmpLikelihood;th4LogLikelihood= tmpLogLikelihood;th4Prior= tmpPrior;
			for (loop_j=0;loop_j<nrComb;loop_j++){ th4VectorComb[loop_j]=tmpVectorComb[loop_j];}
			 // printf("\nEND SWAP 2\n");
			
                        }else{
                         //printf("\nKEEP ON\n");
                        }
		        break;          	
		}
		double post=th1LogLikelihood+th1Prior;
              
        
                if ((iteration % sample_freq ==0) & (iteration >= burnin)){
                        fprintf (fp,"%d\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t", iteration,post,th1LogLikelihood, th1Prior, th1S, th1D, th1W1,th1W2);
                        for (loop_j=0;loop_j<nrComb;loop_j++){ fprintf (fp,"%d\t",th1VectorComb[loop_j]);}
                        fprintf (fp,"\n");
                        fflush(fp);
                }
		iteration=iteration+innerLoopIteration;
	}
	fclose (fp);
	
	free(th1VectorComb);
	free(th2VectorComb);
	free(th3VectorComb);
	free(th4VectorComb);
	free(tmpVectorComb);
	
        free (argThread1);
	free (argThread2);
	free (argThread3);
	free (argThread4);
	
	
        End_of_Routine:
	/*
	 ** Return the status.
	 */
	return (status);
}
