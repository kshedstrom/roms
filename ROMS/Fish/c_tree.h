/******************************************************************* 
   Authors: Bob Torgerson & Kate Hedstrom
   File Name: mod_tree.h
   Date: June 16th 2009

   Info: This program is meant to create a Red-Black Tree for use as
         a way to group theoretical fish eggs into super individuals

  *******************************************************************/

#ifndef FILE_MOD_TREE_H_INCLUDED
#define FILE_MOD_TREE_H_INCLUDED

#include <iostream>  // for std::cout
#include <vector>    // for std::vertor
#include <stdlib.h>  // for abort

class Node
{
 public:
    //Node Constructor
    //Preconditions: None
    //Postconditions: Node is created with distance, egg count & left/right nodes
    //Arg.List: Distance is var. dist, egg_data is var. eggs, theLeft is a pointer to 
    //          the left node, theRight is a pointer to the right node 
    Node(double distance, double egg_data, Node * theLeft, Node * theRight)
      :dist(distance), eggs(egg_data), left(theLeft), right(theRight)
    {}
	
    //Node Destructor
    //Preconditions: None.
    //Postconditions: The pointers are deleted for this node, which in turn calls
    //                the destructor on those nodes as well
    ~Node()
    {
	delete left;
        delete right;
    }


    double dist; //Distance along shore from theoretical zero these eggs are 
    double eggs;  //Number of eggs at this node
    double egg_sum;  //Count of eggs to leaves from this node (found using sum_up func)
    bool red; //Bool for checking if this is a red node (red: T, black:F)
    int momfish;  //Integer for tracking number of breeding mothers?
    bool init;  //Check if this is the initial value for a new tree

    Node * parent; //Pointer to parent of this child none, NULL if root
    Node * left; //Pointer to the left node, NULL if none
    Node * right; //Pointer to the right node, NULL if none
  
};

/*
   5 Properties of a Red-Black Tree
   
   1. A node is either red or black.
   2. The root is black.
   3. All leaves are.black.
   4. Both children of every red node are black.
   5. Every simple path from a given node to any descendant leaf, contains 
      the same number of black nodes.
*/

//Function for Finding Grandparent of Given Node
//Preconditions: None.
//Postconditions: Returns pointer to node's grandparent in tree structure, or NULL if none
Node * grandparent(Node * n)
{
  if ((n != NULL) && (n->parent != NULL))
        return n->parent->parent;
  else
        return NULL;
}

//Function for Finding "Uncle" of Given Node
//Precondtions: None.
//Postconditions: Returns pointer to node's uncle in tree structure, or NULL if none
Node * uncle(Node * n)
{
  Node * g = grandparent(n);
  if (g == NULL)
        return NULL;
  if (n->parent == g->left)
        return g->right;
  else
        return g->left;
}

//Function for Finding Sibling of Given Node
//Preconditions: None.
//Postconditions: Returns pointer to node's sibling in tree structure, or NULL if none
Node * sibling(Node * n)
{
  if (n == n->parent->left)
        return n->parent->right;
  else
        return n->parent->left;
}


class RBTree
{

public:

Node * root;  //Binary Tree Root
int NSpawners; //Number of spawning fish

//Initializer for Red-Black Binary Tree
//Preconditions: None.
//Postconditions: Empty Red-Black Tree Created
RBTree()
{
    root = new Node(0, 0, NULL, NULL);
    root->parent = NULL;
    NSpawners = 0;
    root->init = true;
}

//Destructor
//Preconditions: RBTree created
//Postconditions: RBTree destroyed, memory freed
~RBTree()
{
  delete root;
}

/*******************************************
---------[Start Public Functions]---------
*******************************************/

//Insert Function for Binary Tree
//Preconditions: None.
//Postconditions: Inserts new node into appropriate spot in binary tree
//Note: This also uses the functions for balancing the tree using the Red-Black Binary Tree 
//      conditions used in the check# functions below
void insert(double distance, double egg_data, int ifish)
{
    if (egg_data > 0)
    {
      Node * step = root;
      Node * curr = new Node(distance, egg_data, NULL, NULL);
      curr->eggs = egg_data;
      curr->egg_sum = egg_data;
      curr->momfish = ifish;
      curr->red = true;
      curr->init = false;
      NSpawners++;
    
      bool brlp = false; //break loop boolean value

      while(brlp == false)
      {

	//if step is place holder distance -1, replace root with new insert 
	//(only happens once)
	if (step->init == true)
	{
	  root = curr;
	  root->parent = NULL;
	  return;
	}
	//if step->dist is greater than curr->dist, traverse left branch
	else if (curr->dist < step->dist)
	{
	   if (step->left == NULL)
	   {  
	     step->left = curr;
	     curr->parent = step;
	     brlp = true;
	   }    
	   step = step->left;
	}
	//if step->dist is less than curr->dist, traverse right branch
	else
	{
	    if (step->right == NULL)
	    {
	      step->right = curr;
	      curr->parent = step;
	      brlp = true;
	    }
	    step = step->right;
	}
      }

      sum_up(curr);

      check1(curr);
    }
}

//Search for Node
//Preconditions: None.
//Postconditions: Returns the first node that is found to have the given distance value,
//                if no such node exists, it returns an empty node to prevent errors on
//                the client side.
Node * search(double distance)
{
  Node * n = this->root;

  while (n != NULL)
  {
    if (distance < n->dist)
       n = n->left;
    else if (distance > n->dist)
       n = n->right;
    else if (distance == n->dist)
       return n; 
  }
  
  std::cout << "Could Not Find Node Containing Distance " << distance << std::endl;
  n = new Node(0, 0, NULL, NULL);
  return n;
}

//Collect Superindividuals
//Preconditions: None.
//Postconditions: Checks for valid NSuper and uses the private split_up function to
//                collect an asked for number of superindividuals and returns a pointer
//                to a std::vector containing these new smaller red black trees.
//Note: Most of the work for splitting this->tree into smaller RB Trees is done in the split_up
//      function.
std::vector<RBTree *> collect(int NSuper)
{ 
  int count = 0;
  int NFound = 0;
  if (NSuper < 0)
  {
     std::cout << "You have given a negative number of super individuals. <INVALID>" << std::endl;
     abort();
  } 
 else if (NSuper < NSpawners)
  {
     std::vector< RBTree *> individuals(NSuper);
      
     if (NSuper != 0)
     {
     RBTree * wurzul; //German for root (new root for split_up function)

     int NSplit = int(NSpawners / NSuper);  // NSplit is the number of nodes per superindividual 
     int NExtra = NSpawners % NSuper;  //This number will allow for all the nodes to fit into any number of super individuals given by the client
     
     split_up(root, individuals, wurzul, count, NSplit, NFound, NExtra);
     
     delete this->root;
     }
     return individuals;
  }
  else
  {
    std::vector< RBTree *> individuals;
    
    RBTree * wurzul;
    
    split_less(root, individuals, wurzul, this->NSpawners, count);
    delete this->root;
    return individuals;
  }
}

private:

/*******************************************
---------[Start Private Functions]---------
*******************************************/

//Check # 1
//Preconditions: None.
//Postconditions: If this node is the root, paint it black. Otherwise call check2().
void check1(Node * n)
{
  if (n->parent == NULL)
  {
        n->red = false;
  }
  else
        check2(n);
}

//Check # 2
//Preconditions: None.
//Postconditions: Returns if n->parent is black. Otherwise call check3().
//Note: This makes sense as  the only time that we will have trouble regarding a node's color
//      is when a parent node is red, and each of its children must therefore be black.
void check2(Node * n)
{
  if (n->parent->red == false)
  {
        return;
  }
  else
        check3(n);
} 

//Check # 3
//Preconditions: None.
//Postconditions: Checks if parent & uncle nodes are red, colors them black if so, 
//                and changes the grandparent node to red, it then verifies this is
//                what should be done. Otherwise run check4(). 
void check3(Node * n)
{
  Node * u = uncle(n), *g = grandparent(n);
  
  if ((u != NULL) && (u->red == true))
  {
    n->parent->red = false;
    u->red = false;
    g->red = true;
    check1(g);
    return;
  }
  else
    check4(n);
}

//Check # 4
//Preconditions: None.
//Postconditions: If a grandparent node exists, checks the relationship between n and its parent.
//                Rotates the nodes correspondingly, otherwise go to check5().
void check4(Node * n)
{
  Node * g = grandparent(n);
  
  if ((g != NULL) && (n == n->parent->right) && (n->parent == g->left))
  {
    rotate_left(n->parent);
    n = n->left;
  }
  
  else if ((g != NULL) && (n == n->parent->left) && (n->parent == g->right))
  {
    rotate_right(n->parent);
    n = n->right;
  }
  else
    check5(n);
}

//Check # 5
//Preconditions: None.
//Postconditions: Set n->parent's color to black, and rotate if a grandparent node exists.
void check5(Node * n)
{
  Node * g = grandparent(n);
  
  n->parent->red = false;
  if (g != NULL)
  {
  g->red = true;
  if ((n == n->parent->left) && (n->parent == g->left))
   { 
      rotate_right(g);
      return;
   }
  else
   { 
      rotate_left(g);
      return;
   }
  }
}

//Split Up Function for More Superindividuals than Spawners Available
//Preconditions: None.
//Postconditions: The original RB Tree will be cut into an asked for number of super
//                individuals and returned in a std::vector of RB Trees. However, due
//                to the fact that this is meant to be used when more SI are asked for than
//                spawners, this will return the additional number of nodes with initial starting
//                values of distance = -1. 
//Note: This returns the individual vector with empty batches as asked for rather than a message stating that this
//      individual has no eggs & an initial distance of -1.
void split_less(Node * n, std::vector<RBTree *> & ind, RBTree * & wur, int spawn, int & i)
{
  if (spawn != 0)
  {
    if (n->left != NULL)
       split_less(n->left, ind, wur, spawn, i);
   
    wur = new RBTree();
    wur->insert(n->dist, n->eggs, n->momfish);
    ind.push_back(wur);
    i++;

    if (n->right != NULL)
       split_less(n->right, ind, wur, spawn, i);

    if (i == spawn)
    {
      ind.resize(spawn);
    }
  }
}

//Split Up Function
//Preconditions: None.
//Postconditions: The original RB Tree will be cut into a given number of super individuals
//                and returned in a std::vector of RB trees, all of equal size if the number of
//                spawners % the number of super individuals wanted is equal to zero, but if this is not
//                equal to zero, the first RB tree in the vector will contain the spill over of the total number
//                of nodes.
//Note: Function to be used ONLY by the collect function for spliting the RBTree
//      into smaller super individuals
void split_up (Node * n, std::vector<RBTree *> & ind, RBTree * &  wur, int & i, int split, int & found, int & extra)
{
  if (n->left != NULL)
      split_up(n->left, ind, wur, i, split, found, extra);

  if ((extra != 0) && (found == split))
      extra--;
  else if (found >= split)
  {
    found = 0;
    i++;
  }

  if (found == 0)
      wur = new RBTree();

  found++;
  wur->insert(n->dist, n->eggs, n->momfish);
  ind[i] = wur;
  
  if (n->right != NULL)
      split_up(n->right, ind, wur, i, split, found, extra);
 
}

//Sum Up Function
//Preconditions: None.
//Postconditions: Counts the number of eggs from given node to root, adding the updated egg_sum to each
//                node up the tree.
void sum_up(Node * n)
{
  double egg_left = 0;
  double egg_right = 0;
  
  Node * p = n->parent;
  
  while (p != NULL)
  {
    if (p->left != NULL)
    {
      egg_left = p->left->egg_sum;
    }

    if (p->right != NULL)
    {
      egg_right = p->right->egg_sum;
    }
   
    p->egg_sum = egg_left + egg_right + p->eggs;
 
    p = p->parent;

  }
}

//Sum Us Function
//Preconditions: None.
//Postconditions: Egg_sum for all affected nodes are changed to reflect this rotation
//Note: This function is used only after a rotation either left or right has taken place
void sum_us(Node * n)
{
  double egg_left = 0;
  double egg_right = 0;  

  if (n->left != NULL)
  {
      n->left->egg_sum = n->left->eggs;

      if (n->left->left != NULL)
	    n->left->egg_sum += n->left->left->egg_sum;
      if (n->left->right != NULL)
	    n->left->egg_sum += n->left->right->egg_sum;
      
      egg_left = n->left->egg_sum;
  }
  if (n->right != NULL)
  {
      n->right->egg_sum = n->right->eggs;
 
      if (n->right->left != NULL)
	   n->right->egg_sum += n->right->left->egg_sum;
      if (n->right->right != NULL)
	   n->right->egg_sum += n->right->right->egg_sum;   
 
      egg_right = n->right->egg_sum;
  }

     
     n->egg_sum = egg_left + egg_right + n->eggs;
     if (n->parent != NULL)
           sum_us(n->parent);
}

//Rotate Left Function
//Preconditions: (n->right != NULL)
//Postconditions: Nodes n & o are rotated left in the tree. Any nodes that were originally
//                attached to the left pointer of o, are now being pointed to by the right pointer of n.
void rotate_left(Node * n)
{
  Node * o = n->right;
 
  n->right = o->left; 

  if (o->left != NULL)
  {
     n->right = o->left;
     n->right->parent = n;
  }
  
  o->parent = n->parent;
    
  if (n->parent == NULL)
     root = o;
  else if (n == n->parent->left)
     n->parent->left = o;
  else
     n->parent->right = o;
    
  o->left = n;
  n->parent = o;
  
  check1(n);
  sum_us(n);
  return;
}

//Rotate Right Function
//Preconditions: (n->left != NULL)
//Postconditions: Nodes n & o are rotated right in the tree. Any nodes that were originally 
//                attached to the right pointer of o, are now being poinited to by the left pointer of n.
void rotate_right(Node * n)
{
  Node * o = n->left;
  
  n->left = o->right;

  if (o->right != NULL)
  {
    n->left->parent = n;
  }

  o->parent = n->parent;

  if (o->parent == NULL)
    root = o;
  else if (n == n->parent->left)
    n->parent->left = o;
  else
    n->parent->right = o;

  o->right = n;
  n->parent = o;

  check1(n);
  sum_us(n);
  return;
}

};


#endif
