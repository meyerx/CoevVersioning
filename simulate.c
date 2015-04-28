/*
 ** ml.c: Apply the maximum likelihood method.
 **
** Linda Dib   Last modified: 11 december 2014
 */

#include "def.h"
#include "tree.h"
#include "model.h"
#include "wrapperCPP.h"

//%%%%%%%%%%%%%%%%%%%%%%%%%%% SIMULATION FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%
   
int matinv (double x[], int n, int m, double space[])
{
/* x[n*m]  ... m>=n
   space[n].  This puts the fabs(|x|) into space[0].  Check and calculate |x|.
   Det may have the wrong sign.  Check and fix.
*/
   int i,j,k;
   int *irow=(int*) space;
   double ee=1e-100, t,t1,xmax, det=1;

   for(i=0; i<n; i++) irow[i]=i;

   for(i=0; i<n; i++)  {
      xmax = fabs(x[i*m+i]);
      for (j=i+1; j<n; j++)
         if (xmax<fabs(x[j*m+i]))
            { xmax = fabs(x[j*m+i]); irow[i]=j; }
      det *= x[irow[i]*m+i];
      if (xmax < ee)   {
         printf("\nxmax = %.4e close to zero at %3d!\t\n", xmax,i+1);
         exit(-1);
      }
      if (irow[i] != i) {
         for(j=0; j<m; j++) {
            t = x[i*m+j];
            x[i*m+j] = x[irow[i]*m+j];
            x[irow[i]*m+j] = t;
         }
      }
      t = 1./x[i*m+i];
      for(j=0; j<n; j++) {
         if (j == i) continue;
         t1 = t*x[j*m+i];
         for(k=0; k<m; k++)  x[j*m+k] -= t1*x[i*m+k];
         x[j*m+i] = -t1;
      }
      for(j=0; j<m; j++)   x[i*m+j] *= t;
      x[i*m+i] = t;
   }                            /* for(i) */
   for (i=n-1; i>=0; i--) {
      if (irow[i] == i) continue;
      for(j=0; j<n; j++)  {
         t = x[j*m+i];
         x[j*m+i] = x[j*m + irow[i]];
         x[j*m + irow[i]] = t;
      }
   }
   space[0]=det;
   return(0);
}

/*void QtoFrequencies (double** Q, double* freqVector)
{

    int i,j,loop_k1,loop_k2;
    double ** QTranspose = (double **)malloc (nrComb * sizeof (double *));
    double ** QTInverse = (double **)malloc (nrComb * sizeof (double *));
    for (i = 0; i < nrComb; i++)
    {
        QTInverse[i] = (double *)malloc ((nrComb+1) * sizeof (double));
        QTranspose[i] = (double *)malloc ((nrComb+1) * sizeof (double));
    }
    
    
     for (loop_k1=0;loop_k1<nrComb;loop_k1++){
        for (loop_k2=0;loop_k2<nrComb+1;loop_k2++){
            QTInverse[loop_k1][loop_k2]=0.0;
            QTranspose[loop_k1][loop_k2]=0.0;
        }
    }


        for (loop_k2=0;loop_k2<nrComb+1;loop_k2++){
            QTranspose[0][loop_k2]=1.0;
        }
    for (loop_k1=1;loop_k1<nrComb;loop_k1++){
        for (loop_k2=0;loop_k2<nrComb;loop_k2++){
            QTranspose[loop_k1][loop_k2]=Q[loop_k2][loop_k1];
        }
        QTranspose[loop_k1][nrComb]=0.0;
    }

   fprintf(stderr,"TRANSPOSE \n");
   
   for (loop_k1=0;loop_k1<nrComb;loop_k1++){
        for (loop_k2=0;loop_k2<(nrComb+1);loop_k2++){
           
            fprintf(stderr,"%g ",QTranspose[loop_k1][loop_k2]);
        }
         fprintf(stderr,"\n");
    }
   model_cpp *modelCPP = new_Model_CPP(nrComb);
   Model_CPP_matInverse(modelCPP,QTranspose, QTInverse);
   delete_Model_CPP(modelCPP);
   fprintf(stderr,"INVERSE of TRANSPOSE\n");
    for (loop_k1=0;loop_k1<nrComb;loop_k1++){
        for (loop_k2=0;loop_k2<nrComb;loop_k2++){
           
            fprintf(stderr,"%g ",QTInverse[loop_k1][loop_k2]);
        }
         fprintf(stderr,"\n");
    }
   
   
   for(i=0;i<nrComb;i++)  freqVector[i] = QTInverse[i][i];//?
  
  
   for (i = 0; i < nrComb; i++)
    {
        free (QTInverse[i]);
        free (QTranspose[i]);
    }
   free(QTInverse);	
   free(QTranspose);
}
*/

void QtoFrequencies (double** Q, double* freqVector)
{
/* from rate matrix Q[] to pi, the stationary frequencies:
   Q' * pi = 0     pi * 1 = 1
   space[] is of size n*(n+1).
*/
    int i,j,loop_k1,loop_k2;
    double * QTranspose = (double *)malloc (nrComb *(nrComb +1)* sizeof (double));
    double * QVector    = (double *)malloc (nrComb *nrComb *sizeof (double));
    double * space      = (double *)malloc (nrComb *(nrComb +1)* sizeof (double));
    for (loop_k1=0;loop_k1<nrComb;loop_k1++){
        for (loop_k2=0;loop_k2<nrComb;loop_k2++){
            QVector[loop_k1*(nrComb)+loop_k2]=Q[loop_k1][loop_k2];
            
            }   
    }
    // fprintf(stderr,"QVector:\n ");
    //for (loop_k1=0;loop_k1<nrComb*nrComb;loop_k1++){
    //    fprintf(stderr,"%g ",QVector[loop_k1]);
    //}
    //fprintf(stderr,"\n");
   
    for (loop_k1=0;loop_k1<nrComb+1;loop_k1++){QTranspose[loop_k1]=1.0;}
    for (loop_k1=1;loop_k1<nrComb;loop_k1++){
        for (loop_k2=0;loop_k2<nrComb;loop_k2++){
            QTranspose[loop_k1*(nrComb+1)+loop_k2]=QVector[(loop_k2*nrComb)+loop_k1];//transpose
            }   
        QTranspose[loop_k1*(nrComb+1)+nrComb]=0.0;
    }
    
   // fprintf(stderr,"Before: QTranspose matrix of size nrComb*(nrComb+1)\n");
   //for (loop_k1=0;loop_k1<nrComb;loop_k1++){
   //     for (loop_k2=0;loop_k2<(nrComb+1);loop_k2++){
   //         fprintf(stderr,"%g ",QTranspose[loop_k1*(nrComb+1)+loop_k2]);
   //     }   
   //     fprintf(stderr,"\n");
   // }
    
   //fprintf(stderr,"CALL of matinv\n");
   matinv(QTranspose,nrComb,(nrComb+1),space);
   for(i=0;i<nrComb;i++)  {freqVector[i] = QTranspose[i*(nrComb+1)+nrComb];}
    
   // fprintf(stderr,"After: QTranspose matrix of size nrComb*(nrComb+1)\n");
   //for (loop_k1=0;loop_k1<nrComb;loop_k1++){
   //     for (loop_k2=0;loop_k2<(nrComb+1);loop_k2++){
   //         fprintf(stderr,"%g ",QTranspose[loop_k1*(nrComb+1)+loop_k2]);
   //     }   
   //     fprintf(stderr,"\n");
   // }
   free(QTranspose);
   free(QVector);
   free(space);
}

int function1(double *Pvector){
    int  status=0;
    unsigned int *distVector = (unsigned int *)malloc (nrComb * sizeof (unsigned int));
    int i_loop=0;
    while (i_loop<nrComb){
        distVector[i_loop]=0;
        //fprintf (stderr," %f ",Pvector[i_loop]);
        i_loop++;
    }
    //fprintf (stderr,"\n");
    double *Pvectortemp = (double *)malloc (nrComb*sizeof (double));
    /*i_loop=0;

    i_loop=0;
    double summNorm=0.0;
    while (i_loop<nrComb){
        summNorm=summNorm+Pvector[i_loop];
        i_loop++;
    }
    //fprintf (stderr,"CHECK: summNorm: %f\n",summNorm);
  
    i_loop=0;
    while (i_loop<nrComb){
        Pvector[i_loop]= Pvector[i_loop]/summNorm;
        fprintf (stderr,"%f ",Pvector[i_loop]);
        i_loop++;
    }
    
    double sum=0.0;
    i_loop=0;
    while (i_loop<nrComb){
        sum=sum +Pvector[i_loop];
        i_loop++;
    }
    fprintf (stderr,"\nTotal FREQ: %f \n",sum);
    */
    Pvectortemp[0]=Pvector[0];
    // fprintf (stderr,"\n %f ",Pvectortemp[0]);

    i_loop=1;
    while (i_loop<nrComb){
        Pvectortemp[i_loop] =Pvectortemp[i_loop-1]+Pvector[i_loop];
        //fprintf (stderr," %f ",Pvectortemp[i_loop]);
        i_loop++;
    }
    
    i_loop=0;
    float randV=((double) rand() / (RAND_MAX)) ;
    //fprintf (stderr,"\nrandV: %f \n",randV);
    if(randV <Pvectortemp[i_loop]){
            i_loop=0;
            // fprintf (stderr,"HERE1\n");
    }else{
           // fprintf (stderr,"HERE2\n");
            while (i_loop<nrComb){
                        if(Pvectortemp[i_loop]> randV){
                        break;
                        }
            i_loop++;
            }
    }
    distVector[i_loop]=1;
    //fprintf (stderr,"i_loop: %d \n",i_loop);
    if (i_loop >=nrComb)
    {
        
        status = 1;
        goto End_of_Routine;
    }
    //free(distVector);
    //free(Pvectortemp);
    return i_loop;
    
End_of_Routine:
    /*
     ** Return the status.
     */
    return (status);
}

void function2 (double **Q, double **Qtransposed, struct node* n,double *MatrixA, double *MatrixB, double *MatrixC, int *nIndexState) {
/*function2
calls function1 recursively according to tree
*/
    int assign=0;
    if (n != NULL) {
        if((n->left == NULL) && (n->right == NULL))// is leaf
        {
        }else{
            ////////call function1 to get state of left child and stor it in n->left
            //STEP1:evaluate  Probvector=exmp(Q*n->dist)
            int *statelcIndexArray=(int *)malloc (nrSimulations*sizeof (int));
            char* statelc=(char *)malloc (3*sizeof (char));
            int indexNbreSubs=0;
            
            indexNbreSubs=0;
            int statelcIndex=0;
            double *aLeft = (double *)malloc (nrComb*nrComb*sizeof (double));
            for (assign=0;assign<nrComb*nrComb;assign++){
                aLeft[assign]=0.0;
            }
            
            pVectorFunction(aLeft, Q, Qtransposed,  n->left->dist,  MatrixA, MatrixB,MatrixC);
            
            
            /*fprintf(stderr,"\nDEBUG : aLeft %f\n",n->left->dist);
             for (assign=0;assign<nrComb*nrComb;assign++){
             fprintf(stderr,"%f ",aLeft[assign]);
             }*/
            //fprintf(stderr,"\nDEBUG : statelc: %s\n",n->left->label);
            while(indexNbreSubs<nrSimulations){
                //STEP2: GET child state
                // fprintf(stderr,"1. DEBUG : GET child state\n");
                for (assign=0;assign<nrComb;assign++){
                    n->left->probVector[assign]=aLeft[assign* nrComb +nIndexState[indexNbreSubs]];// Bc aLeft is transposed
                }
                statelcIndex=function1(n->left->probVector);
                if (dataType == NT){
                 
                    statelc= ntComb[statelcIndex];
                }else{
                    statelc= aaComb[statelcIndex];
                }
                //fprintf(stderr,"nupletComb %s\n",statelc);
                // fprintf(stderr,"3. DEBUG : set data\n");
                strcat(n->left->data , statelc);
                //strcat(n->left->data, " ");
                
                statelcIndexArray[indexNbreSubs]=statelcIndex;
                indexNbreSubs++;
                
            }
            // fprintf(stderr,"statelc >%s\n%s\n",n->left->label, n->left->data);
            free(aLeft);
            // fprintf(stderr,"DEBUG : end statelc: %s\n",n->left->label);
            //STEP3: recursive calls
            function2(Q, Qtransposed,n->left, MatrixA, MatrixB, MatrixC, statelcIndexArray);
            
            
            
            
            ////////call function1 to get state of right child
            //STEP1:evaluate  Probvector=exmp(Q*n->dist)
            char* staterc=(char *)malloc (3*sizeof (char));
            int *statercIndexArray=(int *)malloc (nrSimulations*sizeof (int));
            
            int statercIndex=0;
            double *aRight = (double *)malloc (nrComb*nrComb*sizeof (double));
            for (assign=0;assign<nrComb*nrComb;assign++){
                aRight[assign]=0.0;
            }
            //fprintf(stderr,"\nDEBUG : aRight %f\n",n->right->dist);
            pVectorFunction(aRight, Q, Qtransposed,  n->right->dist,  MatrixA, MatrixB,MatrixC);
            /*for (assign=0;assign<nrComb*nrComb;assign++){
             fprintf(stderr,"%f ",aRight[assign]);
             }
             */
            
            indexNbreSubs=0;
            //fprintf(stderr,"DEBUG : staterc: %s\n",n->right->label);
            while(indexNbreSubs<nrSimulations){
                
                //STEP2: GET child state
                for (assign=0;assign<nrComb;assign++){
                    n->right->probVector[assign]=aRight[assign*nrComb +nIndexState[indexNbreSubs]];// Bc aLeft is transposed
                }
                statercIndex=function1(n->right->probVector);
                if (dataType == NT){
                    staterc= ntComb[statercIndex];
                }else{
                    staterc= aaComb[statercIndex];
                }
                //fprintf(stderr,"nupletComb %s\n",staterc);
                strcat(n->right->data, staterc);
                //strcat(n->right->data, " ");
                
                statercIndexArray[indexNbreSubs]=statercIndex;
                indexNbreSubs++;
                
            }
            
            //fprintf(stderr,"staterc >%s\n%s\n",n->right->label, n->right->data);
            //STEP3: recursive calls
            function2(Q, Qtransposed,n->right, MatrixA, MatrixB, MatrixC,statercIndexArray);
            free(aRight);
            //fprintf(stderr,"DEBUG : end staterc: %s\n",n->right->label);
            
            free (statelcIndexArray);
            free (statercIndexArray);
            
            
            /*
             //free(staterc);
             //free(statelc);
             */
        }
    }
}

void printSimulation(struct node* node,FILE *fp) {
    /*
     ** If a leaf node, print the label and data. Otherwise, recurse down the tree.
     */
    if (node->left == NULL && node->right == NULL)
    {
        fprintf (fp,">%s\n%s\n", node->label, node->data);
    }
    else
    {
        printSimulation (node->left,fp);
        printSimulation (node->right,fp);
    }
}

void simulateTree(double **Q, double **Qtransposed, struct node*n,double *MatrixA, double *MatrixB, double *MatrixC,int *VectCoevComb, double * freqVector){
    char* stateroot=(char *)malloc (3*sizeof (char));
    
    /*******************************************
     ** call function 1 to find root state using frequency vector
     ********************************************/
    int *staterootIndexArray=(int *)malloc (nrSimulations*sizeof (int));
    int indexNbreSubs=0;
    while(indexNbreSubs<nrSimulations){
        
        int staterootIndex=function1(freqVector);
        if (dataType == NT){
            stateroot= ntComb[staterootIndex];
        }else{
            stateroot= aaComb[staterootIndex];
        }
        //store state in root
        strcat(n->data, stateroot);
        strcat(n->data, " ");
        //execute recursively according to tree
        fprintf(stderr,"stateroot %s ", stateroot);
        staterootIndexArray[indexNbreSubs]=staterootIndex;
        indexNbreSubs++;
    }
    
    /*******************************************
     ** Call simulation functionout of root state
     ********************************************/
    function2(Q, Qtransposed,n, MatrixA, MatrixB, MatrixC, staterootIndexArray);
    
    /*******************************************
     ** Print simulated sequences
     ********************************************/
    
    //free(stateroot);
    
    free (staterootIndexArray);
}

int simulate ()
{
    int status, start;
    volatile long unsigned counter;
    status = 0;
    
    int AComb= rand() %nrProfiles;
    
    int     loop_j=0, i;
    double *a, **Qtransposed, **Q;
    int *VectCoevComb = (int *)malloc(nrComb*sizeof(int));
    double *AS    = (double *)malloc(sizeof(double));
    double *AD    = (double *)malloc(sizeof(double));
    double *Ar1   = (double *)malloc(sizeof(double));
    double *Ar2   = (double *)malloc(sizeof(double));
    double *MatrixA=(double *)malloc(nrComb*nrComb*sizeof(double));//INPUT
    double *MatrixB=(double *)malloc(nrComb*nrComb*sizeof(double));//INPUT
    double *MatrixC=(double *)malloc(nrComb*nrComb*sizeof(double));//OUTPUT
    double *freqVector   = (double *)malloc(nrComb *sizeof(double));
    *AS   = s;//s
    *AD   = d;//d
    *Ar1  = r1;//w1
    *Ar2  = r2;//w1
    
    /*******************************************
     ** Allocate memory.
     ********************************************/
    a = (double *)malloc (nrComb * sizeof (double));
    Q = (double **)malloc (nrComb * sizeof (double *));
    Qtransposed = (double **)malloc (nrComb * sizeof (double *));
    for (i = 0; i < nrComb; i++)
    {
        Q[i] = (double *)malloc (nrComb * sizeof (double));
        Qtransposed[i] = (double *)malloc (nrComb * sizeof (double));
    }
    for (loop_j=0;loop_j<nrComb;loop_j++){a[loop_j]=0;}
    
    start = time(NULL);
    for (loop_j=0;loop_j<nrComb;loop_j++){
        VectCoevComb[loop_j]=profiles[AComb][loop_j];
        //fprintf(stderr,"%d ",VectCoevComb[loop_j]);
    }
    //fprintf(stderr,"\n");
    
    if(dataType ==NT){
        fprintf(stderr,"DEBUG : Profile: ");
        for (loop_j=0;loop_j<nrNtComb;loop_j++){
            if( VectCoevComb[loop_j]==1){
                fprintf(stderr,"%s,",ntComb[loop_j]);
                fprintf(stderr,"%d, ",loop_j);
            }
        }
        fprintf(stderr,"\n");
    }
    else{
        fprintf(stderr,"DEBUG : Profile: ");
        for (loop_j=0;loop_j<nrAaComb;loop_j++){
            if( VectCoevComb[loop_j]==1){
                fprintf(stderr,"%s ",aaComb[loop_j]);
            }
        }
        fprintf(stderr,"\n");
    }
    
    
    
  
   
    /*******************************************
     ** PREPARE Q MATRIX.
     ********************************************/
   // fprintf(stderr,"DEBUG : Prepare Q Matrix.\n");
    if(dataType==AA){
         setQSimulator(Q,*AS,*AD,*Ar1,*Ar2,VectCoevComb);
    }else{
        setQSimulator(Q,*AS,*AD,*Ar1,*Ar2,VectCoevComb);
    }
    
    
    /*******************************************
     ** Get freqVector from Q
     ********************************************/
    
   int index=0;
    while (index<nrComb){
        freqVector[index]=1.0;
        index++;
    }
   QtoFrequencies (Q, freqVector);
    
    // fprintf(stderr,"DEBUG : Scale Q.\n");
    //Scale Q
    //Qtemp=freq%*%R%*%freq
    //Q<-Q/scaleQ;
   int loop_k1, loop_k2, loop_i;
    double ** QmatrixTemp = (double **)malloc (nrComb * sizeof (double *));
  
    for (i = 0; i < nrComb; i++)
    {
        QmatrixTemp[i] = (double *)malloc (nrComb * sizeof (double));
    }
     for (loop_k1=0;loop_k1<nrComb;loop_k1++){
        for (loop_k2=0;loop_k2<nrComb;loop_k2++){
            QmatrixTemp[loop_k1][loop_k2]=0.0;
        }
    }
    for (loop_i=0;loop_i<nrComb;loop_i++){ 
        for (loop_j=0;loop_j<nrComb;loop_j++){ 
            QmatrixTemp[loop_i][loop_j]= freqVector[loop_i]*Q[loop_i][loop_j] ;
        } 
    }
    
    // printf("Q matrix\n");
    // for (loop_i=0;loop_i<nrComb;loop_i++){ 
    // for (loop_j=0;loop_j<nrComb;loop_j++){ 
    // printf("%g ",  Q[loop_i][loop_j]);
    // }
    // printf("\n");
    // }
    // printf("\n\n");
    
    
    double  scale=0.0;
    for (loop_i=0;loop_i<nrComb;loop_i++){ 
        scale=scale+QmatrixTemp[loop_i][loop_i];
        
    }
    scale= scale*-1;
    for (loop_i=0;loop_i<nrComb;loop_i++){ 
        for (loop_j=0;loop_j<nrComb;loop_j++){ 
            //printf("Scale %g into ",Q[loop_i][loop_j]);
            Q[loop_i][loop_j]=Q[loop_i][loop_j] /scale;
            //printf("%g using %f\n",Q[loop_i][loop_j],scale);
        } 
    }
    
     //printf("\nQ matrix scaled with %f:\n",scale);
    // for (loop_i=0;loop_i<nrComb;loop_i++){ 
    // for (loop_j=0;loop_j<nrComb;loop_j++){ 
    // printf("%g ",  Q[loop_i][loop_j]);
    // }
    // printf("\n");
    // }
    // printf("\n\n");
    
    
     
    for (i = 0; i < nrComb; i++)
    {
        free (QmatrixTemp[i]);
    }
    free(QmatrixTemp);	
    
    
    transposeMatrix(Q,Qtransposed);
    
    /*******************************************
     ** Do simuations.
     ********************************************/
    
    
    fprintf(stderr,"DEBUG : Simulate according to the Tree.\n");
    simulateTree(Q, Qtransposed,&root, MatrixA, MatrixB, MatrixC,VectCoevComb,freqVector);
    
    fprintf(stderr,"DEBUG : Write in a file.\n");
    
    FILE *fp;
    fp = fopen(outFile, "w+");
    //fprintf (fp,">%s\n%s\n\n",(root.label),(root.data));
    printSimulation(&root,fp);
    fflush (fp);
    fclose (fp);
    
    
    free(AS);
    free(AD);
    free(Ar1);
    free(Ar2);
    
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
    free(freqVector);
    //free memory
    free(VectCoevComb);
    
End_of_Routine:
    /*
     ** Return the status.
     */
    return (status);
}
