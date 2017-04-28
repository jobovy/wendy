/*
  bst.h: A binary-search tree for wendy
 */
#ifndef __BST_H__
#define __BST_H__
#include<stdlib.h>
// Basic node: values and index of those values in another array
struct node
{
  int idx;
  double * val;
  struct node *left, *right;
};
// Functions
struct node * bst_newNode(int,double *);
struct node * bst_forceInsert(struct node *, int, double *);
struct node * bst_build(int,int *, double *);
void bst_destroy(struct node *);
struct node * bst_deleteNode(struct node *,double *);
struct node * bst_minValueNode(struct node *);
void bst_inorder(struct node *);
int bst_nNode(struct node *);
#endif /* bst.h */
