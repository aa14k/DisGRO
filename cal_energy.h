// cal_energy.h 
// header file for cal_energy.cpp
// Dec 20, 2005

#ifndef _CALE_
#define _CALE_

#include <vector>

using namespace std;

#include "structure.h"
#include "util.h"
#include "potential.h"
#include "rotamer.h"
//#include "smc.h"

typedef vector <double> FVEC;

// set the weight for different energy terms (from Baker Protein 2005,59:15-29)
//#define W_RAMA 0.1
//#define W_ROT 0.1          // weight for rotamer
//#define W_ATR 0.8          // weight for attraction part of vdw
//#define W_REP 0.65       // weight for repulsion part of vdw
//#define W_REP 0.1
//#define W_SOL 0.65         // weight for solvation
//#define W_PAIR 0.5        // weight for residue pair
//#define W_HB_BB 0.8        // weight for hydrogen bond
//#define W_HB_SC 0.5

class Structure;
void calE(Structure& Conf, int Start, int End, bool type);
void calE_list(Structure& Conf, int* List, int Start, int End, int List_size, bool type);

IDMAP calEraw(Structure& Conf);
void computeEnergy(Structure& Conf, int start, int end);

void vdw_sol(Structure& Conf, double& e_atr, double& e_rep, double& e_sol,
	     int Start, int End, int StartAtom, double& e_halp, double& e_vaa,
	     double& e_cnt, bool debug);
	     
void vdw_sol_list(Structure& Conf, double& e_atr, double& e_rep, double& e_sol,
             int Start, int End, int StartAtom, double& e_halp, double& e_vaa,
	                  double& e_cnt, int* List, int List_size, bool debug);

// calculate the residue pair term
double res_pair_e(Structure& Conf,int,int,int StartAtom);
// assign secondary structure types
void assign_ss(Structure& Conf,int Start,int End, int type);

// ss_e computes energy based on predicted vs. actual secondary structures
double ss_e(Structure& Conf, int Start, int End);

void one_res_en(Structure& conf, Residue& res, int start, int end, int Start, int End,
		double* energy, int type = 0);

void one_res_en_list(Structure& conf, Residue& res, int Start, int End,  double* energy, int type,  int* List, int List_size);

void aa_int_en(Residue& res1, Residue& res2, Atom& atom1, Atom& atom2, double dis,
	       double* energy, int type = 0);
void frag_en(Structure& conf, Residue* res, int start, int end, int length,
	     double* energy);
void two_res_en(Structure& conf, Residue& res1, Residue& res2, int start,
		int end, double* energy);
void one_res_en_SIMPL(Residue& res1, Residue& res2, double* energy, int type);
void one_res_en_sc(Structure& conf, Residue& res, int position, int Start,  int End, int type,
		   double* energy);
void one_res_en_sc_list(Structure& conf, Residue& res, int position, int Start, int End, int type, int* Reslist, int List_size,
                   double* energy); 
void vdw_within_aa(Structure& Conf, Residue& res, double& energy);
double contact_num_e(Structure& Conf, int Start, int End);
void printEnergy(ostream& out, Structure& Conf);
void printRecord(ostream& out, Structure& Conf, string other);
void printRecord(ostream& out, Structure& Conf, string other, bool weighted, bool counts);
double src_e(Structure& Conf, int Start, int End) ;
double model_e(Structure& Conf, int Start, int End);
double temp_e(Structure& Conf, int Start, int End);
double loodis_e(Structure& Conf, int Start, int End, bool type);
double loodis_e_list(Structure& Conf, vector<int> start, vector<int> end, int* List, int List_size, bool type);
void aa_vdw_enum(short type1, short type2, double dis, ofstream& outfile);
void Clash_detection_list(Structure& conf, int Start, int End, vector<int>& ResIdx, vector<int>& ClashNum, int* List, int List_size );
void BBClash_detection_list(Structure& conf, int Start, int End, vector<int>& ResIdx, vector<int>& ClashNum, int* List, int List_size );
void BBClash_detection(Structure& conf, int Start, int End, vector<int>& ResIdx, vector<int>& ClashNum);
int Res_clash_detection_list(Structure& conf,  Residue& res, int Start, int End,  int* List, int List_size);



// Abstract class for an energy mode
class EnergyMode {
 public:
  // Name of this mode
  string Name;
  // Dense section of energy values -- these should all be populated and given
  // descriptive string keys
  SDMAP Dense;
  // Sparse section of energy values -- these need not all be populated and
  // are given integer keys to save space
  IDMAP Sparse;
  // Initializes the energy mode -- sets any parameters, loads any files.
  virtual void Init() = 0;
  // Computes the energy value on a fragment
  virtual void Calc(Structure Conf, int Start, int End) = 0;
};

class ETerm {
 public:
  int index;
  int bin;
  int sign;
};
#endif
