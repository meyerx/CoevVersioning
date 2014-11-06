/*
 ** model.h: Header file for the model.c functions.
 **
 ** Linda Dib    Last modified: 22 July 2014
 */

#ifndef _TREE_H_
#define _TREE_H_


struct node{
    char         data[100000];//DNA or AA sequences
    char         label[128];//sequence name
    double       dist, *probVector;
    struct node *parent, *left, *right;
};

struct sequence
{
    char            *label, *sequence;
    struct sequence *next;
};


/*
 ** Alignement function
 */
int getComb (char *comb);
int readAlignmentFile();
int getData (char *label, char *data);
void freeSequenceList ();

int getNucleotide (char *nucleic);
int getDataSingle (char *label, char *data, int position);


/*
 ** Tree function
 */
void printTree(struct node* node);
int readTreeFile(char *tree);
int createTree1(char* treeParenth, int num_chars);
int createTree2(char* treeParenth, int num_chars);
void freeTree(struct node* n) ;
void insert(struct node *parent, struct node *son);

int createTreeSingle(char* treeParenth, int num_chars, int position);


#endif  /* _TREE_H_ */