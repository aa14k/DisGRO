/* smc.h
header file for smc.cpp, sequential Monte Carlo functions
*/
#ifndef _SMC_
#define _SMC_

#include "structure.h"
#include "rotamer.h"
#include "potential.h"
#include "reprst.h"
#include "util.h"
#include "cal_energy.h"
#include "sample_states.h"

// Empirical distance files
#define FILE_FRAG_C_CA "data/frag.C.CA_pdf_32_19.txt"
#define FILE_FRAG_N_C "data/frag.N.C_pdf_32_19.txt"
#define FILE_LOOPGEO "data/LoopGeo_37_pdf_21.txt"


using namespace std;


// Distribution value: frequency and escore
struct DistVal {
  int    Freq; // Frequency
  double MI;   // Mutual information
};

struct LoopLocaInfo {
  int Start;
  int End;
  int CurrentPos; 
};

// Distribution: map of integers to DistVal's
typedef map<int, DistVal> Dist;

class Folding;

class SMC {
 public:

  // Totals
  int NumConf; // number of conformations in SMC simulation
   
  //bool sortfun (combiEIndex x, combiEIndex y);
  // The current conformation and the array of all sampled conformations
  Structure Conf;
  Structure _conf;
  Structure * ConfArr;
  vector <Structure> LoopStore;
  // States for SMC sampling
  vector<int> NumAngleStates;    // # angles for SMC sampling at each residue
  vector<int> NumDistanceStates; // # distances for SMC sampling
  vector<int> prev_NumDistanceStates;
  vector<int> prev_NumAngleStates;    // # angles for SMC sampling at each residue
  // Some parameters for SMC sampling
  int RandSeed;      // random seed
  string ProtName;   // name of the protein
  string Dir;        // directory storing the files  
  int IsResamp;      // IsResamp = 1: do resampling
  double disConsWt;  // weight for distance constraint in loop modeling
  int AngType;       // torsion angle type
 // int DistanceStyle; // Style of distance select
  double DistanceBy; // Increment in distance distribution
  int DistanceFallback; // Distance sampling when no empiricals fall in range
  bool Close;        // Will analytic closure by performed?
  bool prev_Close;        // Will analytic closure by performed?
  //bool Verbose;      // Should it print out heartbeat information?
  int _isSampSC;    // whether sampling side chain conformations or not
  int numSCStates;	// number of states for side chain sampling. if numSCStates=5, then 5 sidechain conformations will be sampled for each residue
  // the starting position of the protein sequence for SMC sampling
  // for loop modeling, the C atom of the start residue will be
  // sampled. When sampling start residue, the atoms N and CA of start+1
  // are also sampled.
   int Start;
  // The end position of the protein sequence for SMC sampling 
  // for loop modeling, the N and CA atoms of end residue are sampled
  // when sampling residue end-1
   int End;
   bool Backward;
   bool UseBackward;

  bool kmcluster;
  bool Ellipsoid;  // Apply Ellipsoid Algorithm 
  bool Surface;
  bool Eval;
  int* Reslist;
  int List_size;
  int Joint_Angle[20][TORBIN][TORBIN];
  double etedCon[2][20][32][32];
  double minDistcon[2][20];
  double minDistdel[2][20];
  double DistconBy[2][20];
  double DistdelBy[2][20];
  
  double MtorsionEEdis[17][37][37];    // The dihedral angle between the plane of middle point of the loop and the plane extended protein body, start from length =4 
  double minEEdis;
  double minMtorsion;
  double EEdisBy[17];
  double MtorsionBy;
  

  bool foldinglabel;
  bool Refine;
  bool noScore;
  int NumSConf;
  int NDStates_BK;
  int confkeep;    //the number of kept conformations sorted by energy
  int NumClosedconf;
  int outputconf;
  int CA_Constraint;
  Structure InitConf;
  int maxDepth;
  double minRMSD;
  double minEnergy;
  double minERMSD;
  double RMSDsum;
  int min_idx;
  bool TorsionStyle;
  // Constructors
  SMC() {
    SMC(1, 1);
  }
  
  SMC(int r)
  {
    IsResamp=0;
    disConsWt=1;
    End=-1;
    Start=-1;
    Conf.init(r);
    _conf.init(r);
    InitConf.init(r);
    kmcluster = false;
    foldinglabel = false;
    Ellipsoid = false;
    Surface   = false;
    Refine    = false;
    NDStates_BK = 0;
    confkeep =1;
    NumClosedconf = 0;
    outputconf =  0;   
    CA_Constraint   = 1000;
    noScore = false;
    maxDepth = 0;
    minRMSD = 10000;
    minERMSD = 10000;
    minEnergy = 10000;
    RMSDsum = 0;
    TorsionStyle = false;
    min_idx = -1;
    minEEdis = 3.8;
    minMtorsion = 0;
    MtorsionBy = 10;
	Eval = false;
  }

  SMC(int n, int r) {
    IsResamp=0;
    disConsWt=1;
    End=-1;
    Start=-1;
    kmcluster = false;
    foldinglabel = false;
    Ellipsoid = false;
    Surface   = false;
    Refine    = false;
    reshape(n, r);
    NDStates_BK = 0;
    confkeep =1;
    NumClosedconf = 0;
    outputconf = 0;
    CA_Constraint   = 1000;
    noScore = false;
    maxDepth = 0;
    minERMSD = 10000;
    minRMSD = 10000;
    minEnergy = 10000;
    RMSDsum = 0;
    TorsionStyle = false;
    min_idx = -1;
    minEEdis = 3.8;
    minMtorsion = 0;
    MtorsionBy = 10;
	Eval = false;
  }

  // Reshape -- regenerate structures but keep all other params intact
  void reshape(int n, int r) {
    NumConf = n;
    Conf.init(r);
    ConfArr = new Structure[NumConf];
    LoopStore.clear();
   
    for(int nNumStruct = 0; nNumStruct < NumConf; nNumStruct++)
      ConfArr[nNumStruct].init(r);
    if (ConfArr == NULL) {
      cout<<"memory allocation error for ConfArr!"<<endl;
      exit(0);
    }
  }
  

  SMC& operator = (const SMC & other)
  {
      if (this != &other)
      {
        this->IsResamp = other.IsResamp;
        this->disConsWt = other.disConsWt;
        this->End = other.End;
        this->Start = other.Start;
        this->Conf = other.Conf;
        NumConf = other.NumConf;
	if(other.ConfArr)
	{
          ConfArr = new Structure[NumConf];
          for(int i=0;i<NumConf;i++)
          ConfArr[i] = other.ConfArr[i];
	}
	this->_conf = other._conf;
	this->NumAngleStates = other.NumAngleStates;
	this->ProtName = other.ProtName;
	this->Dir = other.Dir;
	this->NumDistanceStates = other.NumDistanceStates;
	this->RandSeed = other.RandSeed;
	this->Close = other.Close;
	AngType = other.AngType;
	DistanceBy = other.DistanceBy;
	Eval = other.Eval;
	InitConf = other.InitConf;
	 numSCStates= other. numSCStates;
	_isSampSC = other._isSampSC;
	DistanceFallback = other.DistanceFallback;
	LoopStore = other.LoopStore;
	kmcluster = other.kmcluster;
	NDStates_BK = other.NDStates_BK;
    confkeep    = other.confkeep;
    NumClosedconf = other.NumClosedconf;
	outputconf  = other.outputconf;
	Ellipsoid = other.Ellipsoid;
	Surface = other.Surface;
	Reslist = other.Reslist;
	Refine	   = other.Refine;
    CA_Constraint   = other.CA_Constraint;
    noScore = other.noScore;
	maxDepth = other.maxDepth;
	List_size = other.List_size;
    minERMSD = other.minERMSD;
	minRMSD = other.minRMSD;
    minEnergy = other.minEnergy;
    RMSDsum = other.RMSDsum;
	TorsionStyle = other.TorsionStyle;
	min_idx = other.min_idx;
    }
     return *this;
  }

  
  // Destructor.
  ~SMC() {
    delete [] ConfArr;
    delete [] Reslist;
  }

  void smc();                 // overall Sequential Monte Carlo algorithm
  bool grow_one(Structure& Conf, int Position, int tmpEnd, Atom& EndPt);   // grow one residue
  void select_one_conf(Structure& conf, Residue* Res, double* energy, int start, int end, int& selected, int, double&, int position);
  void output();
  double eteDis(Residue&, Residue&);
  // for loop modeling
  void fragdis(string fname, int label);
  void geometryinfo(string fname);
  void loop_one(Structure& conf, Residue* Res, double* energy, int start, int end, int& selected, int numStates, double& p_selected, int position);
  void output_loop(string);
  void SinglecalCenter(Residue&);
  void PreProcess();
  void sample_distance( Atom& b, Atom& c, Atom& B,
		       double theta, double lcd,
		       Atom& p, double lcon, int label, int rem, bool verbose);
  bool geometryProb(int length, double EEdis, double Mtorsion);
  int DisGroup(double dis);
  double GroupDis(int group);
  string convertDouble(double value);
  void Wholeproc( vector<Structure>& Topconflist);
  void simpBBT_Init(string fname);

};

#endif
