// dsm.h
//header file for discrete state representation.

#ifndef _REPRST_
#define _REPRST_

#include "atom.h"

// Data files
#define FILE_SCTORSION2 "data/SCT_PF.txt"

// Parameters
#define MAX_NUM_MODEL 5
#define NUM_RES_TP 20
#define MAX_NUM_STATE 30
#define NUM_BB_ANG 1000000  // maximum number of backbone torsion angles 
							// for one amino acid obtained from PDB.
// side chain torsion
class SCT {
 public:
  double torsion[4];
  double bfac[4];
}; // SCT

// side chain representation
class SCR {
 public:
  static vector<SCT> SCANGPDB[NUM_RES_TP];
  static FIMAP SCRMAP[NUM_RES_TP];   // same as BBRMAP
  static void InitSCAng(int, char* scFile, string dir);
  static void testSCRep(int resType, int numStates, char* outFile);
}; // SCR

#endif
