/*****************************************************************
   Author: Bob Torgerson & Kate Hedstrom
   File Name: eggbatch.h
   Date: July 22nd 2009
   
   Info: This program is meant to be a way of storing egg batches

*****************************************************************/

#ifndef FILE_EGG_BATCH_H_INCLUDED
#define FILE_EGG_BATCH_H_INCLUDED

#include <iostream>   // for std::cout
#include <vector>     // for std::vector
#include <algorithm>  // for std::sort
#include <stdlib.h>   // for abort
using std::vector;

class EBV;

class Batch
{

public:

//Batch Constructor
//Preconditions: None.
//Postconditions: New Batch is created of given distance & egg count
Batch(double distance, double egg_data)
    :dist(distance), eggs(egg_data)
{}


//Batch Destructor
//No Dynamic Memory Allocated, None Released
~Batch()
{}

double dist; // Distance along coast
double eggs; // Number of eggs at this distance
int momfish; //Integer for tracking mother fish
double parent;  //Tracks Parent EBV
double previous; //Used for Split-Up

};

bool BatchSortPredicate(Batch * a, Batch * b)
{
  return a->dist < b->dist;
}

class EBV
{
  
public:

vector<Batch *> coast;
int NSpawners;  // Number of spawning fish
double egg_sum; //Total number of eggs along coast
double number; //This is the number of this EBV

//Constructor for EBV
//Precondtions: None.
//Postconditions: NSpawners and egg_sum variables set to zero.
EBV()
  : NSpawners(0), egg_sum(0), number(-1)
{}

//No Dynamic Memory Allocated, None Released
~EBV()
{}

/******************************************
---------[Start Public Functions]----------
******************************************/

//Insert Function
//Preconditions: None.
//Postconditions: Inserts a new Batch at the end of an Egg Batch Vector (EBV)
void insert(double distance, double egg_data, int ifish)
{
  if (egg_data > 0)
  {
    Batch * current = new Batch(distance, egg_data);
    current->momfish = ifish;
    current->parent = number;
    current->previous = 0;
    NSpawners++;
    egg_sum += current->eggs;

    coast.push_back(current);
  }
}

//Sort Function
//Preconditions: None.
//Postconditions: The vector is sorted by distance values in each Batch
void sort()
{
  std::sort(this->coast.begin(), this->coast.end(), BatchSortPredicate);
}

//Print Function
//Preconditions: None.
//Postconditions: Prints coast distance & number of eggs of a batch, and prints total # of spawners and eggs in EBV
//Note: This function could be taken out of the program if unwanted or thought unnecessary
void printthis()
{
  for (int i = 0; i < this->coast.size(); i++)
  {
   std::cout << "Coast Distance: " << this->coast[i]->dist << " : Eggs: " << this->coast[i]->eggs << std::endl;
  }
  std::cout << "Number of spawners : " << this->NSpawners << "  Number of Eggs Total: " << this->egg_sum << std::endl;
  std::cout << "-----------------------\n";
}    

//Collect Function
//Preconditions: None.
//Postconditions: If NSuper is a value capable of being used i.e.(NSuper <= NSpawners), this returns a std::vector
//                of super individuals or EBVs.
//Note: Uses the private function split_up to perform the actual work of splitting this EBV into super individuals
vector< EBV * > collect(int NSuper)
{
  if (NSuper < 0)
  { 
     std::cout << "You have given a negative number of super individuals. <INVALID>";
     abort();
  }
 else if (NSuper < NSpawners)
  {   
     vector< EBV * >individuals;

     if (NSuper != 0)
     {
     EBV * cove;

     int NSplit = int(NSpawners/NSuper);
     int NExtra = NSpawners % NSuper;
     double eggsin = egg_sum/NSuper; //This will be the number of eggs expected in each batch
     
     split_up(this->coast, individuals, cove, eggsin, NSplit,NSuper, NExtra);
     }   

     return individuals;
  }
  else
  {
    vector<EBV *> individuals;
    
    EBV * cove;

    split_less(this->coast, individuals, cove);

    return individuals;
  }
}      

private:

/******************************************
---------[Start Private Functions]----------
******************************************/

//Split Up Function
//Preconditions: None.
//Postconditions: This function performs the initial splitting of the vectors to create the correct number of
//                superindividuals, and passes the individuals to split_batch. 
void split_up(vector<Batch *> cos, vector<EBV * > & ind, EBV * cov, double eggz, int split, int super, int extra)
{
  int found = 0;
  int i = 0;
  double j = 0;
  //This comment below will allow for the minimum number of eggs expected to be printed 
  //std::cout << "Expected Minimum Number of Eggs Per Batch: " << eggz / 10 << std::endl;

  while (super != 0)
  {
    if (found == 0)
    {
       cov = new EBV();
       cov->number = j;
       j++;
    }
    if ((extra != 0) && (found == split))
    {
      cov->insert(cos[i]->dist, cos[i]->eggs, cos[i]->momfish);
      i++;
      extra--;
    }
    else if (found >= split)
    {
      cov->number;
      ind.push_back(cov);
      super--;
      found = 0;
    }
    else
    {
      cov->insert(cos[i]->dist, cos[i]->eggs, cos[i]->momfish);
      found++;
      i++;
    }
  } 

  split_batch(ind, eggz);
}

//Split Up Function for More Superindividuals than Spawners Available
//Preconditions: None.
//Postconditions: Creates a superindividual vector which has only the actual number of
//                spawners rather than adding additional empty eggbatches to fill in the
//                extra superindividuals.
void split_less(vector<Batch *> cos, vector<EBV *> & ind, EBV * cov)
{
  int i = 0;
  
  while (i != NSpawners)
  {
    cov = new EBV();
    cov->insert(cos[i]->dist, cos[i]->eggs, cos[i]->momfish);
    ind.push_back(cov);
    i++;
  }
}

//Split Batch Function
//Preconditions: None.
//Postconditions: This function is called by the split_up function, and does the work of making sure each
//                individual contains at least the average number of eggs off by no more than a factor of ten
//                if possible.
//Note:           In the worst case scenario, this is not always an achievable goal where we want a certain number
//                of superindividuals, but cannot split the eggs to make this goal possible. This should not happen
//                in many general cases.
void split_batch(vector <EBV *> & ind, double eggz)
{
  int sizeint;
  int count;
  double brake;
  double lake;
  double rake;
  while (brake != -1)
  {
    brake = 0;
    rake = 0;
    lake = 0;
    for (int i = 0; i < ind.size(); i++)
    {
      if (i != ind.size()-1)
       {
        count = 0;
	if ((ind[i]->egg_sum < (eggz/2)) || (ind[i+1]->egg_sum >= eggz/5))
	{
	while ((ind[i]->egg_sum < (ind[i+1]->egg_sum *0.9999)) && (ind[i+1]->NSpawners != 1) && ((ind[i+1]->coast[count]->previous != ind[i]->number) || ind[i+1]->coast[count]->previous == 0))
        {
          ind[i]->insert(ind[i+1]->coast[count]->dist, ind[i+1]->coast[count]->eggs, ind[i+1]->coast[count]->momfish);
          sizeint = ind[i]->coast.size()-1;
	  ind[i]->coast[sizeint]->previous = ind[i+1]->number;
	  ind[i]->coast[sizeint]->parent = ind[i]->number;
	  ind[i+1]->NSpawners--;
          ind[i+1]->egg_sum -= ind[i+1]->coast[count]->eggs;
          count++;
	  brake++;
	  lake++;
        }
        if (count != 0)
            ind[i+1]->coast.erase(ind[i+1]->coast.begin(), ind[i+1]->coast.begin()+count);
	}
      }

 if (i != 0)
      {
	sizeint = ind[i-1]->coast.size()-1;
	if ((ind[i]->egg_sum < (eggz/2)) || (ind[i-1]->egg_sum >= eggz/5))
	{
	while ((ind[i]->egg_sum < (ind[i-1]->egg_sum * 0.9999)) && (ind[i-1]->NSpawners != 1) && ((ind[i-1]->coast[sizeint]->previous != ind[i]->number)  || ind[i-1]->coast[sizeint]->previous == 0))
        {
          sizeint = ind[i-1]->coast.size()-1;
          ind[i]->insert(ind[i-1]->coast[sizeint]->dist, ind[i-1]->coast[sizeint]->eggs, ind[i-1]->coast[sizeint]->momfish);
          ind[i]->sort();
	  ind[i]->coast[0]->previous = ind[i-1]->coast[sizeint]->parent;
	  ind[i]->coast[0]->parent = ind[i]->number;
          ind[i-1]->NSpawners--;
          ind[i-1]->egg_sum -= ind[i-1]->coast[sizeint]->eggs;
          ind[i-1]->coast.pop_back();
	  brake++;
	  rake++;
        }
	}
      }
      
    }

    if ((brake == 0) || (lake == rake))
        brake = -1;
  }
}

};

#endif 
