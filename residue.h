/* residue .h
Class Reidue.
*/

#ifndef _RESIDUE_
#define _RESIDUE_

#include "atom.h"

// Data files
#define FILE_ATOMPROP "data/atomProp2.txt"

#define NUM_BB_ATOM 6	// N, CA, C, O, H, CB, Gly has pseudo CB and Pro has a pseudo H
#define CC_DIS_CUT 12	// residue center to center distance cutoff, 
// if smaller than 12 A, the two residues may have atom-atom contacts
#define MAX_NUM_ATOM_PER_RES 18     // 17 + 1, 1 is the frist place to store the residue index

enum {ATM_N, ATM_CA, ATM_C, ATM_O, ATM_H, ATM_CB};
enum {ALA, CYS, ASP, GLU, PHE, GLY, HIS, ILE, LYS, LEU, MET, ASN, PRO, GLN, ARG,
      SER, THR, VAL, TRP, TYR};
      
// enum type of atoms for functional groups of each residue
enum {C_SG=6, D_CG=6, D_OD1=7, D_OD2=8, E_CD=7, E_OE1=8, E_OE2=9, F_CG=6,
      F_CD1=7, F_CD2=8, F_CE1=9, F_CE2=10, F_CZ=11};
enum {H_CG=6,H_ND1=7,H_CD2=8,H_CE1=9,H_NE2=10,I_CG1=6,I_CG2=7,I_CD1=8,K_NZ=9,
      Q_CD=7};
// enum type for secondary structure
enum {SS_H, SS_E, SS_C};
enum {};
enum {};
// there is no OXT atom type in regular residues, The OXT atom can be added to the last residue 
// with proper type and increment its _numAtom by 1.
// For the first N atom, its _type can also be changed. 

// this class is used in side-chain modeling, where the interaction
// energy of rotamers are
// stored in each rotamer, so that when they are grown, the energy of
// their interating rotamer can be updated quickly
class PNT_EN {
 public:
  int rotP;   // position of a rotamer residue *100 + rotamer
  double energy;  // energy
};

typedef vector <PNT_EN> PEVEC;
class Structure;

class Residue {
 public:
  Atom* _atom;     // _atom[0] is ca, _atom[1] is cb, etc.
  short _type;     // residue type represented by int
  short _posn;     // from 1 to NumRes
  short _numAtom;  // number atoms in residue, keep track length of _atom
  // label for secondary structure, 0 for Helix, 1 for Sheet, and 2 for coils.
  short _ss;
  double _phi;      // phi angle
  double _psi;      //  psi angle
  double _omega;    //  amide torsion angle, omega
  double _scChi[5]; // side chain Chi angles
  // used to label side-chain state, initilaized to -1, -2 means
  // side-chain is fixed.
  int _scState;
  int _pdbIndex;   // index in pdb file
  // map for rotamers and their energy, used in hard-shelled model in
  // entropy calculation
  IFMAP _rotE;

  Structure* _parent;   // point to the parent structure
  // the following attributes are used for fast detection of collisions
  Atom _bbc;            // backbone center including Cb
  Atom _scc;            // side-chain center
  Atom _center;         // center of the whole residue
  Atom _rna_center[3];  // center of RNA, 0 is phosphate group, 1 is sugar ring, 2 is the base group
  Atom _SC;             // store the coordinates for seudo side chain atom SC in fragment growth 
  Atom _FG;             // store the center of functional group of side chains
  ISET _res_adj;  // integer set for adjacent residue centers, residue index
  // adjacent centers to bbc, for bbc of other residues, residue
  // index*10, for scc of other residues, index*10+1
  ISET _bb_adj;		// integer set for adjacent backbone centers
  ISET _sc_adj;   // adjacent centers to scc, same rule as above

  // class variable for names
  static string Name1[25];  // one letter name including bases
  static string Name3[25];  // three letter name
  static SIMAP AIMap;   // Amino acid to integer map
  static SSMAP AAMap;  // amino acid three letter name and one letter name map
  static SIMAP SIMap;  // Secondary structure one letter name to integer map
  // map for atoms in each amino acid, for example, N in ASP is DN and
  // its mapped integer value is the index of the atom in atom array _atom
  static SIMAP AtomMap;
  // map for atom index in the atom level distance potential function.
  static SIMAP AtomIndexMap;
  static ISMAP ResMap;
  // parameter between atoms are stored in static variables, the position of each atom is fixed, so that the properties can be located for any atom of a given resiude. 
  static string cType[NUM_RES_TYPE][MAX_NUM_ATOM_RES];   // the string of the atom name
  static int vdwType[NUM_RES_TYPE][MAX_NUM_ATOM_RES];    // the integer type of the atom name
  static int prev_atom[NUM_RES_TYPE][MAX_NUM_ATOM_RES][3];
  static vector < int > FunctionalGroupAtoms[NUM_RES_TYPE];
  // Bond length: length from indexed atom to previous atom
  static double bond_length[NUM_RES_TYPE][MAX_NUM_ATOM_RES];
  // Bond angle: angle formed by indexed atom and previous two atoms
  static double bond_angle[NUM_RES_TYPE][MAX_NUM_ATOM_RES];
  //static double bond_lengthBK[NUM_RES_TYPE][MAX_NUM_ATOM_RES];
  //static double bond_angleBK[NUM_RES_TYPE][MAX_NUM_ATOM_RES];
  static double torsion[NUM_RES_TYPE][MAX_NUM_ATOM_RES];
  static double size[NUM_RES_TYPE];
  static double sc_size[NUM_RES_TYPE];
  static double bb_size;
  static int numAtom[NUM_RES_TYPE];
  // index of functional atoms for each residue
  static IVEC funcAtom[NUM_RES_TYPE];

  Residue(int n = MAX_NUM_ATOM_RES);
  Residue(const Residue& R);
  Residue* operator=(const Residue& r);
  ~Residue();
  void Destruct();
  void init(int n);
  static void InitMap();
  static void InitPar(char* parFile,string&,int outInfo);
  static void RNAPar(char* parFile,string&,int outInfo);  
  void out();
  void cal_scc();
  void cal_fg();
  vector<HBond> AllHBond(bool Given = true, bool Taken = true,
			 int StartAtom = 0);
};

struct ResAtomIdxPair{
   int ResIdx;      // residue index
   int AtmIdx;      // atom index
};
#endif
