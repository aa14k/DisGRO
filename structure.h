// structure.h
/*
header file for protein structure. It contains definition for Point, Atom, Residue and Structure.

*/


#ifndef _STRUCTURE_
#define _STRUCTURE_

#include "residue.h"
#include "potential.h"

#define CUB_ADJ 500
#define MAX_NUM_SC_ST 50

class Structure{
 public:
  int _numRes;
  // residues of the conformation, start from position 1.
  Residue* _res;
  // number of chains in the structure
  int _numChain;
  // the position of end residues for each chain, maximum number of
  // chains is 20, the first chain is one.
  int _firstResidues[20];
  // store the names of the chain, from the first character of the
  // string, assuming each chain is represented by a single character as
  // defined in PDB format.
  string _chainName;
  // Name of the chain
  string _ProtName;
  // store the priority of residues and their positions
  map <double,int> priMap;
  // weight of the chain
  double _weight;
  // weight of side chain sampling
  double _weight_2;
  // value of individual type of energy
  double _enArr[ENERGY_TYPES];
  // overall energy calculated from the components in _enArr
  double _energy;
  // side chain probability of (selected) rotamers 
  double* _scProb;
  // use the native sidechain conformation as a rotamer
  bool _useNatSC;
  //global RMSD comparing to native structure, backbone only
  double G_rmsd;
  //global RMSD, all atom;
  double allatm_rmsd;
  //Dihedral use the mid point of the loop
  double LoopTorsion;
  //Distance between Midpoint and loop anchors
  double MEdis1;
  double MEdis2;
  //RMSD of each loop
  double _single_rmsd[50];
  // residues to be sampled
  vector<int> _toBeSampled;
  // true if the structure has no backbone gaps
  bool Closed;
  bool Success;
  // predicted secondary structure
  vector<int> _ssPred;
  // probability or reliability of prediction
  vector< vector<double> > _ssPredProb;
  // predicted solvent accessibility
  vector<int> _saPred;
  // the oxygen atom at the end of the chain. It can copied to the last
  // atom of the last residue
  Atom _ATM_OXT;
  Atom _Confcenter;
 // The min and max coordinate value of the whole structure, max-x,min-x, max-y, min-y, max-z, min-z;
  double _Fringevtx[6];
  //
  // -1 means no checking of energy, 
  // a non-negative value means to check a particular energy term
  int _checkEn;  
  // this map stores the detail energy contribution from each
  // atom-pairs or residues
  LLIFMAP _EnProfile[ENERGY_TYPES];		
  LLIFMAP _HB_Profile; 
  
  // energy profile for h-bond, keep both side chain and backbone h-bond
  // information
  // this map stores the count of each Energy term, the detailed term
  // under each energy type.
  // For HALP term, there are 22 atom types as in atomProp2.txt
  IIMAP _AATermCnt;
  IIMAP _HBTermCnt;
  IIMAP _BBTTermCnt;
  IIMAP _SCTTermCnt;

  // temperature
  static double _T;
  // the sequence of the structure
  static string _sequence;
  vector <string> missSeq;
  vector <int> missSeqPos;
  vector <int> missSeqfrom1;
  // Key functions
  Structure(int n = 1);
  Structure(const Structure& S);
  Structure* operator=(const Structure& S);
  ~Structure();

  void Destruct();
  void init(int n);
  void copyStructure(const Structure& C, int start, int end);
  bool readPdb(string, SSET& SelRes, int outInfo);
  bool readPdb(vector <string>&, SSET& SelRes);
  void PredictSS();
  static int writeData(char *ptr, size_t size, size_t nmemb, string* s);
  void calSS();
  void addH();
  void calCenter(int Start, int End);
  void calCenter(int Start, int End, bool Sidechain);
  void writePdb(string filename, int start, int end, int type);
  void writePdb(string filename, int start, int end, int l_start, int l_end,int type);
  void writePdb(string filename, int start, int end, int type, int md_num);
  void copy(Structure&,int start, int end);
  void output(int start,int end,int tp);
  void outputAngles(char* filename, int start, int end, int outType);
  void calBBCo(int resIndex, Residue& res, Residue& res2, double phi, double psi, double omega);
  void calSCCo(int resIndex, double* torAngles, Residue& res);
  void analyticClosure(int start, int Start, int End, int* List, int List_size, bool Ellipsoid);
  void analyticClosure_h(int start, double len_change, double bon_change, double t_change, int Start, int End, int* List, int List_size, bool Ellipsoid);
  void computeEnergy();
  double RMSD(Structure other, int start, int end);
  bool IsClosed(int End);
  bool IsClosed(int End, Residue NewRes);
//  void grow_sc(int start, int end, int type, int numStates, int growType,
//	       vector <int>& resToGrow, Structure& Conf);
  void grow_sc(int start, int end, int type, int numStates, int growType,
	       vector <int>& resToGrow);
  void calSC(int start, int end);
  void readSeq(char*);
  void SinglecalCenter(Residue& _res, int type);
  vector<int> contactNum(int start, int end,double radius);
  void StoreSequence();
  void WriteSequence(string file);
  void checkNAN(int , int);
  void calConfCenter(int start, int end);
};

#endif
