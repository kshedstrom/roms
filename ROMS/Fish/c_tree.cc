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

extern "C" void c_tree_init_(int* species) {
    int sp = *species;
    while (sp >= trees.size()) {
       trees.push_back(new RBTree());
    }
    printf("size of trees %d\n", trees.size());
    printf("species %d\n", sp);
}

extern "C" void c_tree_insert_(int* species, double* dist, double* eggs, int* mom) {
    trees[*species]->insert(*dist, *eggs, *mom);
}

extern "C" void c_tree_traverse_(int* species) {
    print_traverse(trees[*species]->root);
}

extern "C" void c_tree_collect_(int* species, int* nsuper, int* nfound, double* eggs, int* mom) {
    vector<RBTree *> arr = trees[*species]->collect(*nsuper);
    *nfound = arr.size();
    int ii;
    for (ii=0; ii<arr.size(); ii++) {
        Node* top = arr[ii]->root;
	eggs[ii] = top->egg_sum;
	mom[ii] = top->momfish;
    }
}

extern "C" void c_tree_destroy_() {
    int ii;
    for (ii=0; ii<trees.size(); ii++) {
        delete trees[ii];
    }
}
