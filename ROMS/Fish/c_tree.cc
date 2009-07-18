/*

Fortran interface for C++ tree code.

*/

#include "c_tree.h"
#include <iostream>
#include <vector>
using namespace std;

void print_traverse(Node * n)
{
  if (n->left != NULL)
    {
      Node * r = n->left;
      print_traverse(r);
    } 

  cout << " Distance: " << n->dist << " : Color: ";
  n->red ? cout << "Red\n" : cout << "Black\n";
  cout << " Egg Sum : " << n->egg_sum <<" Eggs in Node: " << n->eggs << endl;
  cout << "Mother Fish: " << n->momfish << endl;
  cout << "*********************\n";

  if (n->right != NULL)
    {
      Node * t = n->right;
      print_traverse(t);
    }
}

vector<RBTree *> trees;

#ifdef PATH
extern "C" void c_tree_init__(int* species) {
#else
extern "C" void c_tree_init_(int* species) {
#endif
    int sp = *species;
    while (sp >= trees.size()) {
       trees.push_back(new RBTree());
    }
    printf("size of trees %d\n", trees.size());
    printf("species %d\n", sp);
}

#ifdef PATH
extern "C" void c_tree_insert__(int* species, double* dist, double* eggs, int* mom) {
#else
extern "C" void c_tree_insert_(int* species, double* dist, double* eggs, int* mom) {
#endif
    trees[*species]->insert(*dist, *eggs, *mom);
}

#ifdef PATH
extern "C" void c_tree_traverse__(int* species) {
#else
extern "C" void c_tree_traverse_(int* species) {
#endif
    print_traverse(trees[*species]->root);
}

#ifdef PATH
extern "C" void c_tree_collect__(int* species, int* nsuper, int* nfound, double* eggs, int* mom) {
#else
extern "C" void c_tree_collect_(int* species, int* nsuper, int* nfound, double* eggs, int* mom) {
#endif
    vector<RBTree *> arr = trees[*species]->collect(*nsuper);
    *nfound = arr.size();
    int ii;
    for (ii=0; ii<arr.size(); ii++) {
        Node* top = arr[ii]->root;
	eggs[ii] = top->egg_sum;
	mom[ii] = top->momfish;
    }
}

#ifdef PATH
extern "C" void c_tree_destroy__() {
#else
extern "C" void c_tree_destroy_() {
#endif
    int ii;
    for (ii=0; ii<trees.size(); ii++) {
        delete trees[ii];
    }
}
