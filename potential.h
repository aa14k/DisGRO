// potential.h
// classes for potential function terms.

#ifndef _POTENTIAL_
#define _POTENTIAL_


#include <string>
#include "reprst.h"
#include "residue.h"

using namespace std;

// Parameters
#define MAX_NUM_ROT 81
#define MAX_ATOM_TYPE 30
#define HADIP_ATM_TYPE 210
#define HADIP_ATM_PR 22260	// number of atom pairs in high-resolution atom-level potential (HALP), total number of atoms is set to 210 (the true number is a little less than 210)
// so the total number is 210*210/2+210=22260
#define SIMPL_ATM_TP 99 // number of atom types for simplified atom types, there are 19 side chain types plus C CA N O
#define SIMPL_DIS_BIN 20	// from 0 to 10 A with 0.5 A per bin
//#define HADIP_DIS_BIN 1000	// number of distance bins for HADIP potential, the resolution is 0.01 A \
  // from 1.5 A to 11.5 A, the upper limit is set to 10A currently
#define HADIP_DIS_BIN 65   // from 1.5 to 8 A with 0.1 A as bin length, another version with largest distance bin at 15 A
#define LOODIS_ATM_TP 20 //20 amino acid types

#define LOODIS_DIS_BIN 80 // from 0 - 8 A, 0.1 A per bin 
#define PF_DIS_CUT_SQUARE 64	   
#define PF_DIS_CUT 8	    // atom-atom distance cutoff when the distance is greater than PF_DIS_CUT, interaction energy is ignored
//#define LOODIS_DIS_BIN 150 // from 0 - 8 A, 0.1 A per bin 
#define START_DIS 1.5		// starting distance of HADIP is set to 1.5 A
//#define H_INTV 0.01		// the interval of HADIP is set to 0.01A
#define H_INTV 0.1          // another version
#define H_INLO 0.1          // the interval of Loopdistance is set to 0.1
//#define H_INTV 0.5
#define MAX_ENERGY 10		// maximum energy for a contact pair
#define MAX_SCT_ENERGY 5    // maximum energy for SCT term
#define B_T_INT 4           // backbone torsion interval in degrees
#define SC_T_INT 4          // side chain torsion interval in degrees
#define LOODIS_PF_DIS_CUT 15	    // atom-atom distance cutoff for LOODIS
#define SIMPL_DIS_CUTOFF 10 // distance cutoff for SIMPL potential
#define SIMPL_INTV 0.5      // distance interval is set to 0.5 A from 0 to 10 A for SIMPL potential
// if the centers of two residues smaller than this cutoff then the
// two residues may have atom-atom contact.
// this is also used to update _res_adj
#define RES_CENT_DIS_CUTOFF 12

// this is also used in determine whether two centers are close
#define CUB_SIZE 5.5
#define VDW_CLASH_CUTOFF 0.60
// enum type for energy types
// E_VDWA: atractive part of van der Waals energy
// E_VDWR: repulsive part of van der Waals energy
// E_SOL: solvation
// E_ADP: atom-atom distance potential (statistical, high resolution)
// E_RP: residue pair energy using functional groups FGP.txt
// E_BBT: backbone torsion energy
// E_SCA: side chain atom-atom distance energy, this is used in side
//        chain sampling to store the side chain energy of SIMPL
// E_SCT: side chain torsion energy
// E_HBB: backbone hydrogen bond energy
// E_HBS: sidechain hydrogen bond energy
// E_RAMA: Rama term of Rosetta
// E_SS: secondary structure term
// E_CNT: a fully parameterized VDW
// E_VAA: VDW within the same amino acid, this term affects side chain
// conformations the most.
// E_HAV: h-bond adjusted VDW. For donor and acceptor atom, their VDW
// energy is adjusted
// E_CON:  contact numbers
// E_EP:  loop entropy
// E_SRC:  short-range correlation

enum{E_VDWA, E_VDWR, E_VAA, // VDW (0,1,2)
     E_SOL,                 // SOL (3)
     E_HALP,                // HALP (4)
     E_RP,                  // RP (5)
     E_HBB, E_HBS,          // HB (6,7)
     E_BBT,                 // BBT (8)
     E_SCT,                 // SCT (9)
     E_RAMA,                // RAMA H (10)
     E_SIMPL,               // SIMPL (11)
     E_SS,                  // SS (12)
     E_CNT,                 // CNT (13)
     E_ROT,                 // ROT (14)
     E_HAV,		    // HAV (15)
     E_RAMAE,               // RAMA E (16)
     E_RAMAC,		    // RAMA C (17)
     E_CON,		    // CON (18)
     E_EP,                  // ENTROPY (19)
     E_SRC,                 // SRC (20)
     E_MODEL,               // MODEL (21)
     E_TEMP,                // Template-based residue-residue distance terms (22)
     E_HB,                  // Hydrogen bond (23)
     E_LOODIS               // LOODIS (24)
};
#define ENERGY_TYPES 25

// Energy modes
enum{EM_VDW, EM_SOL, EM_HALP, EM_RP, EM_HB, EM_BBT, EM_SCT, EM_RAMA, EM_SIMPL,
     EM_SS, EM_CNT, EM_ROT, EM_CON, EM_EP,EM_SRC, EM_MODEL, EM_TEMP, EM_LOODIS};
#define ENERGY_MODES 18

class Point;
class Structure;
class Rotamer;
class RotLib;
class SCT;

// all parameters are stored in class PF
class PF {
 public:

  // whether B-factor will be used or not
  static int _useBfactor;
  // whether weight for contacts per protein is used in contact potential
  static int _contWt;	    
  static double T;			
  static bool COUNT;
  // initial model
  static Structure ModelStart;

  // Whether to calculate each energy mode
  static bool   cal[ENERGY_MODES];

  // Intercept
  static double intercept;

  // Other energy types
  static double Rama[20][3][36][36];
  static double ResPair[20][20][20];
  static double HB_dis[5][32];
  static double HB_angle_bb[9][72];
  static double HB_angle_sc[10][18];
  static double SRC_below4[20][20];
  static double SRC_above4[20][20];
  static double VDWA[22][22][81];
  static double VDWR[22][22][81];
  static double  LOODIS[20][20][LOODIS_DIS_BIN];
  static vector<double> CNT_bucket;
  static string SSDir;
  //static double EnergyWeight[25];

  // use distance bin to deal with vdw and solvation allow further
  // optimization by statistical method.

  // the first subscript is equal to number of possible atom pairs,
  // the second is distance bin, from 1 A to 6 A with every 0.01 A per bin.
  static double VDW[(MAX_ATOM_TYPE*MAX_ATOM_TYPE-MAX_ATOM_TYPE)/2 +
		   MAX_ATOM_TYPE][500];
  static double SOL[(MAX_ATOM_TYPE*MAX_ATOM_TYPE-MAX_ATOM_TYPE)/2 +
		   MAX_ATOM_TYPE][500];
  static double ROT[NUM_RES_TP][MAX_NUM_ROT];

  // VDW_AA store for each amino acids the atom pairs that need to be
  // calculated for VDW energy
  static vector < vector <int> > VDW_AA[20];

  // high-resolution atom-level distance dependent statistical potential
  // Atom Pair Index Locator to save about half of space for HADIP
  // array, store the index of atom pairs in this array,
  static vector < vector <int> > APIL;
  // for example, at APIL[10][11] the atom pair index of atom 10 and
  // atom 11 is stored.
  static vector < vector < double > > HADIP;
  // simplified atom types distance dependent potential
  static vector < vector < vector < double > > > SIMPL;
  // loop distance dependent potential
 // static vector < vector < vector < double > > > LOODIS;
  // map atom types to integers for SIMPL potential
  static SIMAP simplMap;
  // map of integers to atom types for SIMPL potentil
  static ISMAP simplMapr;
  static vector < vector < vector < vector < double > > > > angCount;
  // side chain torsion potential energy parameters
  static IFMAP SCTP[25];
  
  // back bone torsion potential energy parameters
  static IFMAP BBTP[25];
  // CNT potential terms
  static map <string, IDMAP> Parameter;



  static IDMAP TEMP;
  static string TempFile;
  static IDMAP TEMP_CONT;

  // index of atom pairs for HADIP, the first two dimensions are residues types
  // the second two dimensions are atom types.
  static void InitPar(string, string&, int outInfo);

  // rot_e calculate the rotamer term of the potential function
  // It can also be used to assign rotameric state of each side-chain
  // for a conformation
  double cal_rot_e(Structure& Conf, RotLib& R);
  
  // calculate the residue pair term
  double res_pair_e(Structure& Conf);
  
  void assign_ss(Structure& Conf);
  
  // rama_e calculate the rama term of the potential function
  // It can also be used to calculate the phi and psi angles for a conformation
  double rama_e(Structure& Conf);
  
  // calculate the h-bond energy
  double hb_e(Structure& Conf,double*,double*);
  
  // dfire type potential
  double dfire(Structure&);

  // initialize parameters in SIMPL potential
  static void initSIMPL(string& dir);

  // initialize parameters in LOODIS potential
  static void initLOODIS(string filename); 
  
  // initialize parameters in CNT potential
  static void initCNT(string filename);
  
  // initialize short-range correlation
   static void initSRC(string filename);
 
  //initialize vdw energy
  static void initVDWtable(string filename);

   // initialize temperature schedule
   static void initTemp(string filename);
   static int TempSched[100];
   static float TempValues[100];
   static int curTempSched;
};


int AAmap(string aname);

double hb_energy(Point*, int type, int tp_ss, int tp_sp, int tp_en, double dis);

double vdw_sol(Point& p1, Point& p2, double dis, double r_ij, int tp1, int tp2);

double pnorm(double x, double mean, double sd);

#endif
