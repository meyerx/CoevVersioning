/*
** coev.c: Program for site coevolution analysis.
**
** Linda Dib & Wim Hordijk   Last modified: 30 April 2013
*/

#include "def.h"


/*
** Global variables.
*/
#define MAX_LINE_LEN 10000

const char nucleotide[nrNt] = {'A', 'C', 'G', 'T'};
const char aminoacid[nrAa] = {'A','R','N','D','C','E','Q','G','H','I','L','K','M','F','P','S','T','W','Y','V'};

int     IT, sample_freq, print_freq, burnin, col1, col2, method, dataType,
        nrComb, *obsCombs, **profiles, nrProfiles, maxNrProfiles,
        innerLoopIteration, model;
double  s, d, alpha, beta;
double  r1, r2; 

char    treeFile[1024], alignFile[1024], outFile[1024], *ntComb[nrNtComb],
       *aaComb[nrAaComb];

struct node     root;
struct sequence seqList;


//%%%%%%%%%%%%%%%%%%%%%%%%%%% COMBINATION FUNCTION %%%%%%%%%%%%%%%%%%%%%%%

/*
** createNtCombs: Create the nucleotide combinations.
*/

void createNtCombs ()
{
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

void createAaCombs ()
{
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
** getComb: Get the index of a given combination.
**
** Parameters:
**   - comb: The nucleotide or amino acid combination.
**
** Returns:
**   The index of the combination (an integer in [0:15] or [0:399]) or
**   -1 if 'comb' is not a valid combination.
*/

int getComb (char *comb)
{
  int index;

  /*
  ** Get the index of the given nucleotide combination.
  */
  index = 0;
  if (dataType == NT)
  {
    while ((index < nrNtComb) && (strcmp (comb, ntComb[index]) != 0))
    {
      index++;
    }
    if (index >= nrNtComb)
    {
      index = -1;
    }
  }
  else
  {
    while ((index < nrAaComb) && (strcmp (comb, aaComb[index]) != 0))
    {
      index++;
    }
    if (index >= nrAaComb)
    {
      index = -1;
    }
  }
  
  /*
  ** Return the result.
  */
  return (index);
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

void addProfComb (int *curProf, int index)
{
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

void getProfiles ()
{
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


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TREE FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%
void insert(struct node *parent, struct node *son){
  if (parent->left == NULL) 
    parent->left = son;
  else if (parent->right == NULL)
    parent->right = son;
  else printf("DEBUG:: ERROR");
}


void printTree(struct node* node) {
  /*
  ** If a leaf node, print the label and data. Otherwise, recurse down the tree.
  */
  if (node->left == NULL)
  {
    printf ("%s %s\n", node->label, node->data);
  }
  else
  {
    printTree (node->left);
    printTree (node->right); 
  }
} 


int readTreeFile(char *tree)
{
  int              num_chars, len;
  char             line[MAX_LINE_LEN];
  FILE            *fp;
  struct sequence *seq, *last;

  /*
  ** Open the tree file for reading.
  */
  num_chars = 0;
  fp = fopen (treeFile, "r");
  if (fp == NULL)
  {
    num_chars = -1;
    fprintf (stderr, "Could not open the tree file '%s'.\n", treeFile);
    goto End_of_Routine;
  }

  /*
  ** Read the Newick tree representation. It can be divided over several
  ** (consecutive) lines.
  */
  tree[0] = '\0';
  while (fgets (line, MAX_LINE_LEN, fp) != NULL)
  {
    /*
    ** Remove the trailing newline character, if present.
    */
    len = strlen (line);
    if (line[len-1] == '\n')
    {
      line[len-1] = '\0';
    }
    /*
    ** Concatenate the next line.
    */
    strcat (tree, line);
  }
  num_chars = strlen (tree);
  fclose (fp);

  /*
  ** Open the alignment file for reading.
  */
  fp = fopen (alignFile, "r");
  if (fp == NULL)
  {
    num_chars = -1;
    fprintf (stderr, "Could not open the alignment file '%s'.\n", alignFile);
    goto End_of_Routine;
  }

  /*
  ** Read the sequences. FASTA file format is assumed.
  */
  seq = NULL;
  while (fgets (line, MAX_LINE_LEN, fp) != NULL)
  {
    /*
    ** Remove the trailing newline character, if present.
    */
    len = strlen (line);
    if (line[len-1] == '\n')
    {
      line[len-1] = '\0';
    }
    /*
    ** If a new sequence is encountered, create a new list element and save
    ** the sequence label.
    */
    if (line[0] == '>')
    {
      last = seq;
      if (seq == NULL)
      {
	seq = &seqList;
      }
      else
      {
	seq = (struct sequence*)malloc (sizeof (struct sequence));
	last->next = seq;
      }
      seq->label = strdup (&(line[1]));
      seq->sequence = NULL;
      seq->next = NULL;
    }
    /*
    ** Otherwise, save the actual sequence.
    */
    else
    {
      if (seq->sequence == NULL)
      {
	seq->sequence = strdup (line);
      }
      else
      {
	len = strlen (seq->sequence) + strlen (line) + 1;
	seq->sequence = (char *)realloc (seq->sequence, len*sizeof (char));
	strcat (seq->sequence, line);
      }
    }
  }
  
 End_of_Routine:
  /*
  ** Close the file and return the length of the Newick tree representation.
  */
  fclose (fp);
  return (num_chars);
}


/*
** getData: Get the nucleotide or amino acid combination for a given species.
**
** Parameters:
**   - label: The species label.
**   - data:  A character array (of at least length 3) to put the observed
**            combination into.
**
** Returns:
**   If the species label exists:  0.
**   Otherwise:                   -1.
*/
int getData (char *label, char *data)
{
  int              status, i;
  struct sequence *seq;

  status = 0;
  
  /*
  ** Find the sequence label.
  */
  seq = &(seqList);
  while (strcmp (seq->label, label) != 0)
  {
    
    if (seq->next == NULL)
    {
      status = -1;
      fprintf (stderr, "Unmatched sequence label in tree: '%s'.\n", label);
      goto End_of_Routine;
    }
    seq = seq->next;
  }

  /*
  ** Get the data from the requested columns.
  */
  data[0] = seq->sequence[col1-1];
  data[1] = seq->sequence[col2-1];
  data[2] = '\0';

  /*
  ** Remember the observed combination.
  */
  i = getComb (data);
  if (i == -1)
  {
    status = -1;
    fprintf (stderr, "Invalid nucleotide or amino acid combination in columns %d and %d of species %s.\n", col1, col2, label);
    goto End_of_Routine;
  }
  obsCombs[i] = 1;
  
 End_of_Routine:
  /*
  ** Return the status.
  */
  return (status);
}


int createTree(char* treeParenth, int num_chars)
{
  int          i, index, assign, indexname, indexdistance, status;
  char         name[128], distance[32];
  struct node *current, *son;

  status = 0;
  i = 0;
  current = NULL;
  while (i< num_chars && treeParenth[i]!=';'){
    if (treeParenth[i]=='(' ){	
      if (current == NULL){
	//printf("HERE1 %d\n",i);			
	current = &root;
      }else{
	//printf("HERE2 %d\n",i);
	son = (struct node*)malloc (sizeof (struct node));
	strcpy(son->data, "internal");
	strcpy(son->label, "internal");
	son->parent=current;
	son->dist=0.0;// will be changed afterwards
	son->right=NULL;
	son->left=NULL;
	son->probVector = (double *)malloc (nrComb * sizeof (double));
	for (assign=0;assign<nrComb;assign++){ 
	  son->probVector[assign]=0.0;
	}
	insert(current,son);
	current=son;
      }
      i++;
    }else if (treeParenth[i]==',' || treeParenth[i]==')'){
      //printf("HERE3 %d\n",i);			
      if (treeParenth[i]==')'){
	while (treeParenth[i]!=':' && treeParenth[i]!=',' && treeParenth[i]!=';'){ 
	  i++;
	}
      }
      else{
	i++;
      }
    }
    else if (treeParenth[i]==':'){
      //printf("HERE4 %d\n",i);			
      indexdistance=0;
      while (treeParenth[i]!=',' && treeParenth[i]!=')' && treeParenth[i]!=';' ){
	if(treeParenth[i]!= ':') {
	  distance[indexdistance]=treeParenth[i];
	  indexdistance++;
	}
	i++;
      }
      distance[indexdistance]='\0';
      current->dist=atof(distance);
      current=current->parent;
			
    }else {
      //printf("HERE5 %d\n",i);
      indexname = 0;
      while (treeParenth[i]!=':'){
	name[indexname]=treeParenth[i];
	indexname++;
	i++;
      }
      name[indexname]='\0';
			
      indexdistance=0;
      while (treeParenth[i]!=',' && treeParenth[i]!=')' && treeParenth[i]!=';' ){
	if(treeParenth[i]!= ':') {
	  distance[indexdistance]=treeParenth[i];
	  indexdistance++;
	}		
	i++;
      }	
      distance[indexdistance]='\0';
      son = (struct node*)malloc(sizeof(struct node));
      strcpy (son->label, name);
      if (getData (name, son->data) == -1)
      {
	status = -1;
	goto End_of_Routine;
      }
      son->parent=current;
      son->dist=atof(distance);
      son->right=NULL;
      son->left=NULL;
      index = getComb (son->data);
      if (index == -1)
      {
	status = -1;
	fprintf (stderr, "Invalid nucleotide combination in given columns.\n");
	goto End_of_Routine;
      }
      son->probVector = (double *)malloc (nrComb * sizeof (double));
      for (assign = 0; assign < nrComb; assign++)
      { 
	son->probVector[assign] = 0.0;
      }
      son->probVector[index] = 1.0;
      insert(current,son);
    }
  }

 End_of_Routine:
  /*
  ** Return the status.
  */
  return (status);
}


void freeTree(struct node* n) { 
  if (n != NULL) { 
    if((n->left == NULL) && (n->right == NULL))// is leaf
    {
      free (n->probVector);
      free(n);
    }else{
      freeTree( n->left);
      freeTree( n->right);
      free (n->probVector);
      free(n);
    }
  }
}


void freeSequenceList ()
{
  struct sequence *seq, *next;

  /*
  ** Free up all list elements and their contents.
  */
  seq = seqList.next;
  free (seqList.label);
  free (seqList.sequence);
  while (seq != NULL)
  {
    next = seq->next;
    free (seq->label);
    free (seq->sequence);
    free (seq);
    seq = next;
  }
}


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
  print_freq = 1000;   // print frequency (parms are printed on screen)
  burnin = 0;          // burnin to exclude (keep it 0 for now)
  s = 1;           
  d = 1;             
  r1 = 1;  r2 = 1;
  strcpy (treeFile, "treeInput.txt");
  strcpy (alignFile, "alignment.txt");
  strcpy (outFile, "output.log");
  col1 = 1;
  col2 = 2;
  method = ML;
  dataType = NT;
  model = JC;
  nrComb = nrNtComb;
  
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
      if (r1 < 0)
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
      if (r2 < 0)
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
      printf ("  -method s:   Use method s (either 'bayes' or 'ml'). Default is 'ml: maximum likelihood'.\n");
      //printf ("  -data s:     Use data type s (either 'nt' or 'aa'). Default is 'nt'.\n");
      printf ("  -IT n:       Run for n iterations. Default is n=10000.\n");
      printf ("  -sfreq n:    Write every n'th iteration to file. Default is n=1000.\n");
      printf ("  -pfreq n:    Print every n'th iteration on the screen. Default is n=1000.\n");
      printf ("  -burnin n:   Number of burn-in iterations. Default is n=0.\n");
      printf ("  -s v:        ??. Default is v=1.\n");
      printf ("  -d v:        ??. Default is v=1.\n");
      printf ("  -r1 v:       ??. Default is v=1.\n");
      printf ("  -r2 v:       ??. Default is v=1.\n");
      printf ("  -tree s:     The name of the file containing the input tree in Newick\n");
      printf ("               format. Default is s=treeInput.txt.\n");
      printf ("  -align s:    The name of the file containing the sequence alignment in\n");
      printf ("               FASTA format. Default is s=alignment.txt.\n");
      printf ("  -out s:      The name of the log file to write the results to. Default\n");
      printf ("               is s=output.log.\n");
      printf ("  -cols n1 n2: The columns in the alignment to use for the analysis.\n");
      printf ("               Default is n1=1 and n2=2.\n");
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
  if (getArguments (argc, argv) == -1)
  {
    status = 1;
    goto End_of_Routine;
  }

  /*
  ** Create the nucleotide or amino acid combinations.
  */
  if (dataType == NT)
  {
    createNtCombs ();
  }
  else
  {
    createAaCombs ();
  }
  obsCombs = (int *)malloc (nrComb * sizeof (int));
  for (i = 0; i < nrComb; i++)
  {
    obsCombs[i] = 0;
  }
    
  ///////TREE READ FROM FILE
  num_chars = readTreeFile(tree);
  len = strlen (seqList.sequence);
  if ((col1 > len) || (col2 > len))
  {
    status = 1;
    fprintf (stderr, "Column value larger than sequence length.\n");
    goto End_of_Routine;
  }

  strcpy (root.data, "root");
  root.parent=NULL;
  root.dist=0.0;
  root.right=NULL;
  root.left=NULL;
  root.probVector = (double *)malloc (nrComb * sizeof (double));
  for (assign=0;assign<nrComb;assign++){ 
    root.probVector[assign] = 0.0;
  }
  if (createTree (tree, num_chars) == -1)
  {
    status = 1;
    goto End_of_Routine;
  }
  getProfiles ();
  srandom(time(NULL));
  printf("Pair %d %d number of profiles: %d\n",col1,col2,nrProfiles);
  for (i = 0; i < nrComb; i++)
  {
    //printf("%d ",obsCombs[i]);
  }
  //printf("\n");
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
    if (ml () == -1)
    {
      status = 1;
      goto End_of_Routine;
    }
  }
  
  /*
  ** Free up the allocated memory.
  */
  freeTree(root.left);
  freeTree(root.right);
  free(root.probVector);
  freeSequenceList ();
  free (obsCombs);
  if (dataType == NT)
  {
    for (i = 0; i < nrNtComb; i++)
    {
      free (ntComb[i]);
    }
  }
  else
  {
    for (i = 0; i < nrAaComb; i++)
    {
      free (aaComb[i]);
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
