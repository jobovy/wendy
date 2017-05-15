/*
  bst.c: A binary-search tree for wendy
 */
#include <stdio.h>
#include <stdlib.h>
#include <bst.h>
// Adding a node as leaf
struct node * bst_newNode(int idx, double * val){
  struct node *temp =  (struct node *)malloc(sizeof(struct node));
  temp->idx = idx;
  temp->val = val;
  temp->left = temp->right = NULL;
  return temp;
}
// Inserting a new value, always inserts, duplicates to the right
struct node * bst_forceInsert(struct node * node, int idx, double * val){
  // Empty tree: return single new node
  if (node == NULL) return bst_newNode(idx,val);
  // Build tree using recurrence
  if (*val < *(node->val))
    node->left= bst_forceInsert(node->left,idx,val);
  else if (*val >= *(node->val))
    node->right= bst_forceInsert(node->right,idx,val);
  return node;
}
// Building a full tree
struct node * bst_build(int N,int * idxs, double * vals){
  struct node *root = NULL;
  int ii;
  for (ii=0; ii < N; ii++)
    root= bst_forceInsert(root,*(idxs+ii),vals+ii);
  return root;
}
// Free memory of a tree
void bst_destroy(struct node * node){
  // Clean up memory recursively
  if (node->left != NULL)
    bst_destroy(node->left);
  if (node->right != NULL)
    bst_destroy(node->right);
  free(node);
}
// Function to delete a node
struct node * bst_deleteNode(struct node * root,double * val){
  // default case
  if (root == NULL) return root;
  // Find the node to delete using recursion
  if (*val < *(root->val))
    root->left= bst_deleteNode(root->left,val);
  else if (*val > *(root->val))
    root->right= bst_deleteNode(root->right,val);
  // Found the node, now delete
  else {
    // If just one child, move child to node
    if (root->left == NULL) {
      struct node *temp= root->right;
      free(root);
      return temp;
    }
    else if (root->right == NULL) {
      struct node *temp= root->left;
      free(root);
      return temp;
    }
    // Two children, need the inorder successor (the smallest value in 
    // the right subtree)
    struct node * temp = bst_minValueNode(root->right);
    root->idx= temp->idx;
    root->val= temp->val;
    // Then delete the inorder successor
    root->right= bst_deleteNode(root->right,temp->val);
  }
  return root;
}
// Function to find the minimum value in a tree
struct node * bst_minValueNode(struct node * node){
  struct node * current = node;
  // Just move left until the end
  while (current->left != NULL)
    current = current->left;
  return current;
}
// Function to print the tree in order LCOV_EXCL_START
void bst_inorder(struct node *root){
  if (root != NULL){
    bst_inorder(root->left);
    printf("%i,%g\n",root->idx,*root->val);
    fflush(stdout);
    bst_inorder(root->right);
  }
}
// Function to print the number of nodes in the tree
int bst_nNode(struct node *root){
  int out= 0;
  if (root != NULL){
    out+= bst_nNode(root->left);
    out+= bst_nNode(root->right);
    out+= 1;
  }
  return out;
}
// LCOV_EXCL_STOP
