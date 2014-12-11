/*
** coev.c: Program for site coevolution analysis.
**
** Linda Dib & Wim Hordijk   Last modified: 11 december 2014
*/

#include "def.h"
#include "tree.h"

/********************************
** Global variables.
********************************/
const char nucleotide[nrNt] = {'A', 'C', 'G', 'T'};
const char aminoacid[nrAa] = {'A','R','N','D','C','E','Q','G','H','I','L','K','M','F','P','S','T','W','Y','V'};

int     IT, sample_freq, print_freq, burnin, col1, col2, method, dataType,
        nrComb, *obsCombs, **profiles, nrProfiles, maxNrProfiles,
        innerLoopIteration, model, nrSimulations, opt;
double  s, d, alpha, beta;
double  r1, r2; 

char    treeFile[100000], alignFile[100000], outFile[100000], *aaComb[nrAaComb];
char   *ntComb[nrNtComb], tree[10*MAX_LINE_LEN];

struct node     root;
struct sequence seqList;


//%%%%%%%%%%%%%%%%%%%%%%%%%%% COMBINATION FUNCTION %%%%%%%%%%%%%%%%%%%%%%%
int isConflict (int index1, int index2){
  int  conflict;
  char c10, c11, c20, c21;

  /*
  ** Get the combinations.
  */
  if (dataType == NT)
  {
    c10 = ntComb[index1][0];
    c11 = ntComb[index1][1];
    c20 = ntComb[index2][0];
    c21 = ntComb[index2][1];
  }
  else
  {
    c10 = aaComb[index1][0];
    c11 = aaComb[index1][1];
    c20 = aaComb[index2][0];
    c21 = aaComb[index2][1];
  }
  
  /*
  ** Compare the combinations.
  */
  if ((c10 == c20) && (c11 == c21))
  {
    // identical
    conflict = 0;
  }
  else if ((c10 == c20) || (c11 == c21))
  {
    // conflict in one position
    conflict = 1;
  }
  else
  {
    // no conflict
    conflict = -1;
  }

  /*
  ** Return the result.
  */
  return (conflict);
} 

/*
** createNtCombs: Create the nucleotide combinations.
*/
void createNtCombs (){
  int i, j, nr;

  /*
  ** Create the labels.
  */
  nr = 0;
  for (i = 0; i < nrNt; i++)
  {
    for (j = 0; j < nrNt; j++)
    {
      ntComb[nr] = (char *)malloc (3*sizeof (char));
      ntComb[nr][0] = nucleotide[i];
      ntComb[nr][1] = nucleotide[j];
      ntComb[nr][2] = '\0';
      nr++;
    }
  }
}

/*
** createAaCombs: Create the amino acid combinations.
*/
void createAaCombs (){
  int i, j, nr;

  /*
  ** Create the labels.
  */
  nr = 0;
  for (i = 0; i < nrAa; i++)
  {
    for (j = 0; j < nrAa; j++)
    {
      aaComb[nr] = (char *)malloc (3*sizeof (char));
      aaComb[nr][0] = aminoacid[i];
      aaComb[nr][1] = aminoacid[j];
      aaComb[nr][2] = '\0';
      nr++;
    }
  }
}

/*
** addProfComb: Recursively add observed combinations to the current profile
**              and see if it creates any conflicts. If not, add it to the
**              list of possible profiles.
**
** Parameters:
**   - curProf: The current profile.
**   - index:   The current index.
*/
void addProfComb (int *curProf, int index){
  int i, j, confl;

  /*
  ** Try each next observed combination.
  */
  for (i = index + 1; i < nrComb; i++)
  {
    if (obsCombs[i] == 1)
    {
      /*
      ** Check for conflicts.
      */
      confl = 0;
      for (j = 0; j <= index; j++)
      {
	if ((curProf[j] == 1) && (isConflict (i, j) != -1))
	{
	  confl = 1;
	  break;
	}
      }
      /*
      ** If no conflict, add the combination and the profile, and recurse.
      */
      if (confl == 0)
      {
	curProf[i] = 1;
	for (j = 0; j < nrComb; j++)
	{
	  profiles[nrProfiles][j] = curProf[j];
	}
	nrProfiles++;
	addProfComb (curProf, i);
	curProf[i] = 0;
      }
    }
  }
}

/*
** getProfiles: Get a list of possible profiles based on the observed
**              combinations.
**
*/
void getProfiles (){
  int i, j, *curProf;

  /*
  ** Allocate and initialize memory for the maximum number of profiles
  ** depending on the data type.
  */
  if (dataType == NT)
  {
    maxNrProfiles = nrNtProf;
  }
  else
  {
    maxNrProfiles = nrAaProf;
  }
  profiles = (int **)malloc (maxNrProfiles * sizeof (int*));
  for (i = 0; i < maxNrProfiles; i++)
  {
    profiles[i] = (int *)malloc (nrComb * sizeof (int));
    for (j = 0; j < nrComb; j++)
    {
      profiles[i][j] = 0;
    }
  }

  /*
  ** Recursively add observed combinations to the profile and check whether it
  ** is possible (i.e., no conflicts).
  */
  curProf = (int *)malloc (nrComb * sizeof (int));
  for (i = 0; i < nrComb; i++)
  {
    curProf[i] = 0;
  }
  nrProfiles = 0;
  for (i = 0; i < nrComb-1; i++)
  {
    if (obsCombs[i] == 1)
    {
      curProf[i] = 1;
      addProfComb (curProf, i);
      curProf[i] = 0;
    }
  }
  free (curProf);
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
/*
** getArguments: Get and check the program arguments.
**
** Parameters:
**   - argc: The number of arguments to the program.
**   - argv: A pointer to the list of arguments.
**
** Returns:
**   If everything is fine:  0.
**   Otherwise:             -1.
*/

int getArguments (int argc, char **argv)
{
  int status, i;

  status = 0;

  /*
  ** Set defaults.
  */
  IT = 10000;          // number of iterations
  sample_freq = 1000;  // sampling frequency (parms are logged to file)
  print_freq = 0;   // print frequency (parms are printed on screen)
  burnin = 0;          // burnin to exclude (keep it 0 for now)
  s = 1;           
  d = 100;             
  r1 = 0.1;  r2 = 0.1;
  strcpy (treeFile, "treeInput.txt");
  strcpy (alignFile, "alignment.txt");
  strcpy (outFile, "output.log");
  col1 = 1;
  col2 = 2;
  method = ML;
  dataType = NT;
  model = JC;
  nrComb = nrNtComb;
  nrSimulations=1;
  opt=1;
  /*
  ** Get and check the arguments.
  */
  i = 1;
  while (i < argc)
  {
       if (strcmp (argv[i], "-method") == 0)
        {
            i++;
            if (strcmp (argv[i], "bayes") == 0)
            {
                method = BAYES;
            }
            else if (strcmp (argv[i], "ml") == 0)
            {
                method = ML;
            }
            else if (strcmp (argv[i], "sm") == 0)
            {
                method = SM;
            }
            else
            {
                status = -1;
                fprintf (stderr, "Invalid method.\n");
                goto End_of_Routine;
            }
            i++;
        }
        else if (strcmp (argv[i], "-data") == 0)
        {
            i++;
            if (strcmp (argv[i], "nt") == 0)
            {
                dataType = NT;
                nrComb = nrNtComb;
            }
            else if (strcmp (argv[i], "aa") == 0)
            {
                dataType = AA;
                nrComb = nrAaComb;
            }
            else
            {
                status = -1;
                fprintf (stderr, "Invalid data type.\n");
                goto End_of_Routine;
            }
            i++;
        }
        
        else if (strcmp (argv[i], "-IT") == 0)
        {
            IT = atoi (argv[++i]);
            if (IT < 1)
            {
                status = -1;
                fprintf (stderr, "Invalid value for argument 'IT'.\n");
                goto End_of_Routine;
            }
            i++;
        }
        else if (strcmp (argv[i], "-sfreq") == 0)
        {
            sample_freq = atoi (argv[++i]);
            if (sample_freq < 0)
            {
                status = -1;
                fprintf (stderr, "Invalid value for argument 'sfreq'.\n");
                goto End_of_Routine;
            }
            i++;
        }
        else if (strcmp (argv[i], "-pfreq") == 0)
        {
            print_freq = atoi (argv[++i]);
            if (print_freq < 0)
            {
                status = -1;
                fprintf (stderr, "Invalid value for argument 'pfreq'.\n");
                goto End_of_Routine;
            }
            i++;
        }
        else if (strcmp (argv[i], "-burnin") == 0)
        {
            burnin = atoi (argv[++i]);
            if (burnin < 0)
            {
                status = -1;
                fprintf (stderr, "Invalid value for argument 'burnin'.\n");
                goto End_of_Routine;
            }
            i++;
        }
        else if (strcmp (argv[i], "-s") == 0)
        {
            s = strtod (argv[++i], NULL);
            if (s < 0)
            {
                status = -1;
                fprintf (stderr, "Invalid value for argument 's'.\n");
                goto End_of_Routine;
            }
            i++;
        }
        else if (strcmp (argv[i], "-d") == 0)
        {
            d = strtod (argv[++i], NULL);
            if (d < 0)
            {
                status = -1;
                fprintf (stderr, "Invalid value for argument 'd'.\n");
                goto End_of_Routine;
            }
            i++;
        }
         else if (strcmp (argv[i], "-r1") == 0)
        {
            r1 = strtod (argv[++i], NULL);
            if (d < 0)
            {
                status = -1;
                fprintf (stderr, "Invalid value for argument 'r1'.\n");
                goto End_of_Routine;
            }
            i++;
        }
         else if (strcmp (argv[i], "-r2") == 0)
        {
            r2 = strtod (argv[++i], NULL);
            if (d < 0)
            {
                status = -1;
                fprintf (stderr, "Invalid value for argument 'r2'.\n");
                goto End_of_Routine;
            }
            i++;
        }
        else if (strcmp (argv[i], "-tree") == 0)
        {
            strcpy (treeFile, argv[++i]);
            i++;
        }
        else if (strcmp (argv[i], "-align") == 0)
        {
            strcpy (alignFile, argv[++i]);
            i++;
        }
        else if (strcmp (argv[i], "-out") == 0)
        {
            strcpy (outFile, argv[++i]);
            i++;
        }
        else if (strcmp (argv[i], "-ns") == 0)
        {
            nrSimulations = atoi (argv[++i]);
            if (nrSimulations >5000 || burnin < 0)
            {
                status = -1;
                fprintf (stderr, "Invalid value for argument 'ns'. It should be between 0 and 5000.\n");
                goto End_of_Routine;
            }
            i++;
        } 
         else if (strcmp (argv[i], "-opt") == 0)
        {
            opt = atoi (argv[++i]);
            i++;
        }   
    else if (strcmp (argv[i], "-cols") == 0)
    {
      col1 = atoi (argv[++i]);
      col2 = atoi (argv[++i]);
      if ((col1 < 0) || (col2 < 0))
      {
	status = -1;
	fprintf (stderr, "Invalid value for argument 'cols'.\n");
	goto End_of_Routine;
      }
      i++;      
    }
    else if (strcmp (argv[i], "-h") == 0)
    {
       printf ("\n%s [OPTIONS]\n\n", argv[0]);
            printf ("  -method s:   Use method s (either 'bayes',  'ml' or 'sm'). Default is 'ml: maximum likelihood'.\n");
            printf ("  -data s:     Use data type s (either 'nt', 'aa' ). Default is 'nt'.\n");
            printf ("  -tree s:     The name of the file containing the input tree in Newick\n");
            printf ("               format. Default is s=treeInput.txt.\n");
            printf ("  -align s:    When method is bayes or ml, the name of the file containing the sequence alignment in\n");
            printf ("               FASTA format. Default is s=alignment.txt.\n");
            printf ("  -cols v1 v2: The columns in the alignment to use for the analysis.\n");
            printf ("               Default is n1=1 and n2=2.\n");
            printf ("  -s v:        ??. Default is v=1.\n");
            printf ("  -d v:        ??. Default is v=1.\n");
            printf ("  -r1 v:      The r1 parameter rates. Default value is 1.\n");
            printf ("  -r2 v:      The r2 parameter rates. Default value is 1.\n");
            printf ("  -out s:      The name of the log file to write the results to. Default\n");
            printf ("               is s=output.log.\n");
            printf ("  -IT n:       When method is bayes, run for n iterations. Default is n=10000.\n");
            printf ("  -sfreq n:    When method is bayes, write every n'th iteration to file. Default is n=1000.\n");
            printf ("  -pfreq n:    When method is bayes, print every n'th iteration on the screen. Default is n=1000.\n");
            printf ("  -burnin n:   Number of burn-in iterations. Default is n=0.\n");
            printf (" -ns: When method is sm, number of simulated pairs. It should be maximum 500.\n");
            printf ("  -h:          Print this help screen and exit.\n");
            
            printf ("\n");
            exit (0);
    }
    else
    {
      status = -1;
      fprintf (stderr, "Unknown argument '%s'.\n", argv[i]);
      goto End_of_Routine;
    }
  }
  End_of_Routine:
  /*
  ** Return the status.
  */
  return (status);
}


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MAIN FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%
int main (int argc, char **argv)
{
  int  status, i, j, num_chars, assign, len;
  char tree[10*MAX_LINE_LEN], dat[3];
  
  /*
  ** Get the arguments.
  */
  status = 0;
    printf( "DEBUG : Get Arguments.\n");
  if (getArguments (argc, argv) == -1)
  {
    status = 1;
    goto End_of_Routine;
  }
   

   printf( "DEBUG : Create combinations.\n");
    /*
     ** Create the nucleotide or amino acid combinations.
     */
    if (dataType == AA)
    {
        createAaCombs();
    }
    else
    {
         createNtCombs ();
    }
    
    if (dataType == AA)
    {
        maxNrProfiles = nrAaProf;
    }
    else
    {
        maxNrProfiles=nrNtProf;
    }
    
    profiles = (int **)malloc (maxNrProfiles * sizeof (int*));
    
    
    obsCombs = (int *)malloc (nrComb * sizeof (int));
    for (i = 0; i < nrComb; i++)
    {
        obsCombs[i] = 0;
    }
    
    ///////TREE READ FROM FILE
    num_chars = readTreeFile(tree);
    
    if(num_chars <0){
        goto End_of_Routine;
    }
    strcpy (root.label, "root");
    strcpy (root.data, "");
    root.parent=NULL;
    root.dist=0.0;
    root.right=NULL;
    root.left=NULL;
    root.probVector = (double *)malloc (nrComb * sizeof (double));
    printf( "DEBUG : Assign Prob vector.\n");
    for (assign=0;assign<nrComb;assign++){
        root.probVector[assign] = 0.0;
    }
    srandom(time(NULL));
    srand ( time(NULL) );
    
      // START//////////////////////
    if (method == BAYES || method == ML){
       printf(  "DEBUG : ReadAlignmentFile\n");
        num_chars =readAlignmentFile();
        if(num_chars <0){
            goto End_of_Routine;
        }
        
        
        len = strlen (seqList.sequence);
       if (col1 > len || col2 > len){
                status = -1;
                fprintf (stderr, "Invalid value for arguments in 'cols'");
                goto End_of_Routine;
            }
            
        printf(  "DEBUG : Create treeStructure with alignmnent data. \n");
        //Create treeStructure with alignmnent data
        if (createTree1 (tree, num_chars) == -1)
        {
            status = 1;
            goto End_of_Routine;
        }
        printf(  "DEBUG : Observed combinations: ");
        for (i = 0; i < nrComb; i++)
        {
            printf(  "%d",obsCombs[i]);
        }
        printf(  "\n");
        printf(  "DEBUG : Get profiles. \n"); 
        getProfiles ();
        
        /*
         ** Call the appropriate function depending on the chosen method.
         */
        if (method == BAYES)
        {
            if (bayes () == -1)
            {
                status = 1;
                goto End_of_Routine;
            }
        }
        else if (method == ML)
        {
            printf(  "DEBUG : call of ML \n");
            if (ml () == -1)
            {
                status = 1;
                goto End_of_Routine;
            }
        }  
  }else {
  if (method == SM) {
        printf( "DEBUG : Simulate Tree:\n");
        nrProfiles=0;
        while(nrProfiles==0){
            
            int index =0;
            while(index <nrComb)
            {   
                obsCombs[index] = 0;
                index++;
            }
            int randObs=rand() % 16;//int randObs=rand() % nrComb; //BECAUSE WE CAN HAVE SEGMENTATION ERROR IF TOO LARGE
            index =0;
            while(index <randObs)
            {
                int valueRand=rand() % nrComb;
                obsCombs[valueRand] = 1;
                index++;
            }
            printf( "DEBUG : Get Profiles:\n");
            
            getProfiles ();
            printf( "DEBUG : Profiles number: %d\n",nrProfiles);
        }
        
        //Create treeStructure without alignmnent data
        if (createTree2 (tree, num_chars) == -1)
        {
            status = 1;
            goto End_of_Routine;
        }
        if (simulate () == -1)
        {
            
            status = 1;
            goto End_of_Routine;
        }
    }
    }
  /*
  ** Free up the allocated memory.
  */
  freeTree(root.left);
    freeTree(root.right);
    free(root.probVector);
    
    if( method != SM){
        freeSequenceList ();
        free (obsCombs);
    }
    
    if (dataType == AA){
        for (i = 0; i < nrAaComb; i++)
        {
            free (aaComb[i]);
        }
    }
    else{
        for (i = 0; i < nrNtComb; i++)
        {
            free (ntComb[i]);
        }
    }
    
    for (i = 0; i < maxNrProfiles; i++)
    {
        free (profiles[i]);
    }
    free (profiles);
  
 End_of_Routine:
  /*
  ** Return the status.
  */
  return (status);
}
