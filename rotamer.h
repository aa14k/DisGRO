// rotamer.h
// class for side-chain rotamer. Not frequently used.

#ifndef _ROTAMER_
#define _ROTAMER_

#include <vector>
#include <map>
#include <string>

using namespace std;

// Data files
#define FILE_ROTPT_DUN "data/rotPt_Dun.dat"
#define FILE_ROTLIB_R "data/rotLib_R.txt"

// class for a rotable bond
class RotBond {
 public:
  double angle;
  double range;
}; // RotBond

// class for one rotamer, store the number of states (r1,r2..) 
// and the ranges of the state if possible

class Rotamer{
 public:
  short aaType;     // the amino acid type of this rotamer
  int intRep;       // integer representation
  double prob;       // the probability of the rotamer
  RotBond chi[5];  // 

  static int numRotBond[20];  // class variable for 20 aa.

  // default constructor
  // n here is the index of the amino acid
  Rotamer(int n){
  }
  Rotamer(){
  }

  Rotamer* operator=(const Rotamer& r){
    intRep=r.intRep;
    prob=r.prob;
    for(int i=0;i<numRotBond[aaType];++i){
      chi[i].angle=r.chi[i].angle;
      chi[i].range=r.chi[i].range;
    }
    return this;
  }
  void operator<<(const Rotamer& r){
  }
}; // Rotamer

typedef vector < Rotamer > ROTVEC;
typedef map <long, ROTVEC> LRVMAP;

class RotLib
{
 public:
	 static int scLib;

	// Richardson Rotamer Lib, scLib equals to 1

	static vector<Rotamer> RLR[20];

	// Dunbrack rotamer lib, scLib equals to 2

  	static LRVMAP RLD[20]; 

  	static void read_rotlib(char* filename, int lib,string dir, int outInfo);
}; // RotLib


#endif

