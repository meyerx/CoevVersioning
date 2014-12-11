/*
 ** model.h: Header file for the model.c functions.
 **
** Linda Dib & Wim Hordijk   Last modified: 11 december 2014
 */
#include "tree.h"
#include "def.h"
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TREE FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%


/*
 ** Alignement function
 */

int getComb (char *comb){
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
    int index;
    
    /*
     ** Get the index of the given nucleotide combination.
     */
    index = 0;
    if (dataType == AA)
    {while ((index < nrAaComb) && (strcmp (comb, aaComb[index]) != 0))
    {
        index++;
    }
        if (index >= nrAaComb)
        {
            index = -1;
        }
        
    }
    else
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
    
    /*
     ** Return the result.
     */
    return (index);
}

int readAlignmentFile(){ char             line[MAX_LINE_LEN];
    int              num_chars, len;
    struct sequence *seq, *last;
    FILE            *fp;
    /*
     ** Read the sequences. FASTA file format is assumed.
     */
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

int getData (char *label, char *data){
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
        fprintf (stderr, "Invalid nucleotide or amino acid combination of species %s.\n", label);
        goto End_of_Routine;
    }
    obsCombs[i] = 1;
    
End_of_Routine:
    /*
     ** Return the status.
     */
    return (status);
}

void freeSequenceList (){
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

int getNucleotide (char *nucleic){
    int index = 0;
    char s[2];
    sprintf(s,"%c",nucleotide[index]);
    while ((index < nrNt) && (strcmp (nucleic, s) != 0))
    {
        index++;
        s[0]='\0';
        sprintf(s,"%c",nucleotide[index]);
    }
    
    
    if (index >= nrNt)
    {
        index = -1;
    }
    /*
     ** Return the result.
     */
    return (index);
}

int getDataSingle (char *label, char *data, int position){
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
    
    data[0] = seq->sequence[position-1];
    data[1] = '\0';
    /*
     ** Remember the observed combination.
     */
    
    i = getNucleotide (data);
    if (i == -1)
    {
        status = -1;
        goto End_of_Routine;
    }
    obsCombs[i] = 1;
    
End_of_Routine:
    /*
     ** Return the status.
     */
    return (status);
}


/*
 ** Tree function
 */
void insert(struct node *parent, struct node *son){
    if (parent->left == NULL)
        parent->left = son;
    else if (parent->right == NULL)
        parent->right = son;
    else printf("DEBUG:: ERROR in INSERT function\n");
}

int readTreeFile(char *tree){
    int              num_chars, len;
    char             line[MAX_LINE_LEN];
    FILE            *fp;
    
    
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
    return (num_chars);
End_of_Routine:
    /*
     ** Close the file and return the length of the Newick tree representation.
     */
    return (1);
}

int createTree1(char* treeParenth, int num_chars){
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

int createTree2(char* treeParenth, int num_chars){
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
                strcpy(son->data, "");
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
            /* if (getData (name, son->data) == -1)
             {
             status = -1;
             goto End_of_Routine;
             }*/
            son->parent=current;
            son->dist=atof(distance);
            son->right=NULL;
            son->left=NULL;
            /* index = getComb (son->data);
             if (index == -1)
             {
             status = -1;
             fprintf (stderr, "Invalid nucleotide combination in given columns.\n");
             goto End_of_Routine;
             }*/
            son->probVector = (double *)malloc (nrComb * sizeof (double));
            for (assign = 0; assign < nrComb; assign++)
            {
                son->probVector[assign] = 0.0;
            }
            //son->probVector[index] = 1.0;
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

void printTree(struct node* node){
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

