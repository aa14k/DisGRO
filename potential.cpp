// potential.cpp

#include "potential.h"
#include "atom.h"
#include "structure.h"
#include "rotamer.h"
#include "util.h"
#include "cal_energy.h"
#include "reprst.h"
#include <sstream>

int PF::_useBfactor = 0;
int PF::_contWt = 0;
bool PF::COUNT;
double PF::T=1;
bool PF::cal[ENERGY_MODES];
string PF::SSDir;
Structure PF::ModelStart;

double PF::Rama[20][3][36][36];
double PF::ResPair[20][20][20];
double PF::HB_dis[5][32];
double PF::HB_angle_bb[9][72];
double PF::HB_angle_sc[10][18];
double PF::VDW[(MAX_ATOM_TYPE*MAX_ATOM_TYPE-MAX_ATOM_TYPE)/2 +
	      MAX_ATOM_TYPE][500];
double PF::SOL[(MAX_ATOM_TYPE*MAX_ATOM_TYPE-MAX_ATOM_TYPE)/2 +
	      MAX_ATOM_TYPE][500];
double PF::ROT[20][MAX_NUM_ROT];
double PF::SRC_below4[20][20];
double PF::SRC_above4[20][20];
double PF::VDWA[22][22][81];
double PF::VDWR[22][22][81];

double PF::LOODIS[20][20][LOODIS_DIS_BIN];
int PF::TempSched[100];
float PF::TempValues[100];
int PF::curTempSched;
vector < double > PF::CNT_bucket;

vector < vector <int> > PF::APIL;
vector < vector <double> > PF::HADIP;
vector < vector < vector < double > > > PF::SIMPL;
SIMAP PF::simplMap;
ISMAP PF::simplMapr;
vector< vector < vector < vector < double > > > > PF::angCount;
IFMAP PF::SCTP[25];   // side chain torsion potential energy parameters
IFMAP PF::BBTP[25];
vector < vector <int> > PF::VDW_AA[20];
map <string, IDMAP> PF::Parameter;
IDMAP PF::TEMP;
string PF::TempFile;
IDMAP PF::TEMP_CONT;
void PF::InitPar(string ParFile, string& dir, int outInfo) {
  int j,k,l;
  char tmpLine[500], filename[500], longLine[10000];
  string tmpStr;
  vector <string> strVec; // used in reading files
  ifstream inFile;
    // Set up maps to fold parameters
    IDMAP EMap, AMap, BMap, HMap, SMap;
    Parameter["E"] = EMap;
    Parameter["A"] = AMap;
    Parameter["B"] = BMap;
    Parameter["H"] = HMap;
    Parameter["S"] = SMap;
    for (int i = 0; i < ENERGY_TYPES; i++)
      PF::Parameter["E"][i] = 0.0;
    // Read parameters and load into maps                                           
    if (ParFile != "none") {                                                        
      vector <string> ParFileV = FileLines(ParFile);                                
      for (int i = 0; i < ParFileV.size(); i++) {
		if(ParFileV[i][0] == '#') continue;
        SVEC bob; 
        split(ParFileV[i], ":", bob);
        string id = ParFileV[i].substr(0, 1);                                       
        bob[0].erase(0, 1);                                                         
        int    key   = atoi(bob[0].c_str());
        double value = atof(bob[1].c_str());
        if (value != 0)
          Parameter[id][key] = value;
      }
    }

  if (cal[EM_VDW]) {
    // initialize VDW_AA for VDW energy within an amino acids
    vector <int> tmpV;
    // set number of atoms for each AA that need to calculate VDW
    int VDW_CNT[20] = {
      0, 1, 3, 4, 3,
      0, 3, 3, 5, 3,
      3, 3, 0, 4, 6,
      1, 1, 1, 4, 3
    }; 
    for (int i = 0; i < 20; ++i) {
      // push 7 vector<int> to each VDW_AA[i], most of them will be empty
      // VDW_AA[i][ATM_N] will store the atoms that need to calculate
      // VDW energy with ATM_N of residue i
      for (j = 0; j<7; ++j)
      {
	VDW_AA[i].push_back(tmpV);
      }
      if (VDW_CNT[i] >= 1) {
	for (k = NUM_BB_ATOM; k<Residue::numAtom[i]; k++)
	  VDW_AA[i][ATM_O].push_back(k);
      }
      if(VDW_CNT[i] >= 3){
	if (i == 7) {
	  // isoleucine
	  for(k = (NUM_BB_ATOM+2); k<Residue::numAtom[i]; ++k){
	    VDW_AA[i][ATM_N].push_back(k);
	    VDW_AA[i][ATM_C].push_back(k);
	  }
	}
	else {
	  for(k = (NUM_BB_ATOM+1); k<Residue::numAtom[i]; ++k){
	    VDW_AA[i][ATM_N].push_back(k);
	    VDW_AA[i][ATM_C].push_back(k);
	  }
	}
      }
      if(VDW_CNT[i] >= 4){
	for(k = (NUM_BB_ATOM+2); k<Residue::numAtom[i]; ++k){
	  VDW_AA[i][ATM_CA].push_back(k);
	}
      }
      else if(VDW_CNT[i] >= 5){
	for(k = (NUM_BB_ATOM+3); k<Residue::numAtom[i]; ++k){
	  VDW_AA[i][ATM_CB].push_back(k);
	}
      }
      else {
	for (k = (NUM_BB_ATOM+4); k < Residue::numAtom[i]; k++)
	  VDW_AA[i][6].push_back(k);
      }
    }
  }
}


// initialize temperature schedule
void PF::initTemp(string filename){
  vector <string> inFile= FileLines(filename);
  for (int i = 0; i < inFile.size(); i++) {
    SVEC holder;
    split(inFile[i]," ",holder);
    TempSched[i]=atoi(holder[0].c_str());
    TempValues[i]=atof(holder[1].c_str());
  }
  TempSched[inFile.size()]=1000000000;

  // set initial temperature
  T = TempValues[0];
  curTempSched=0;

}

// initialize short-range correlation term look-up table
void PF::initSRC(string filename){
  int indic=0;

  int index1;
  int index2;
	
  vector <string> inFile= FileLines(filename);
  for (int i = 0; i < inFile.size(); i++) {
    if (indic==0) {
      SVEC holder;
      split(inFile[i]," ",holder);
      index1= AAmap(holder[0].c_str());
      index2=AAmap(holder[1].c_str());
      indic=1;
    }
    else if (indic==1) {
      SVEC holder;
      split(inFile[i],"  ",holder);
      SRC_below4[index1][index2] = atof(holder[0].c_str());
      SRC_above4[index1][index2] = atof(holder[1].c_str());
      indic=0;
    }
  }	
}

int AAmap(string aname) {
  if (aname=="ALA") return ALA;
  if (aname=="CYS") return CYS;
  if (aname=="ASP") return ASP;
  if (aname=="GLU") return GLU;
  if (aname=="PHE") return PHE;
  if (aname=="GLY") return GLY;
  if (aname=="HIS") return HIS;
  if (aname=="ILE") return ILE;
  if (aname=="LYS") return LYS;
  if (aname=="LEU") return LEU;
  if (aname=="MET") return MET;
  if (aname=="ASN") return ASN;
  if (aname=="PRO") return PRO;
  if (aname=="GLN") return GLN;
  if (aname=="ARG") return ARG;
  if (aname=="SER") return SER;
  if (aname=="THR") return THR;
  if (aname=="VAL") return VAL;
  if (aname=="TRP") return TRP;
  if (aname=="TYR") return TYR;
}	

//initialize potential function LooDis
void PF::initLOODIS(string filename) {
  int i,j,k;
  double pt;
  char tmpLine[500];
  string tmpStr;
  vector <string> strVec;       // used in reading files
  string tmp;
  ifstream inFile;

  inFile.open(filename.c_str(), ios::in);
  
  if(!inFile.is_open()) {
      cout<<"Error!!! Cannot open LOODIS potential file! "<<filename<<endl;
      exit(0);
  }
					
  while(!inFile.eof()){
    inFile.getline(tmpLine,500);
    tmpStr=tmpLine;
	if(tmpStr[0] == '#') continue;
    if(tmpStr.length()<3) continue;
    strVec.clear();
    split(tmpStr," ",strVec);
	if(strVec.size() == 6)
	{
     i = atoi(strVec[1].c_str());
     j = atoi(strVec[2].c_str());
     k = atoi(strVec[0].c_str());
     pt = atof(strVec[5].c_str());
     LOODIS[i][j][k] = pt;
     LOODIS[j][i][k] = pt;
    }

  }
   inFile.close();
}

//init vdw lookup table
void PF::initVDWtable(string filename) {
  int i,j,k;
  double pt1, pt2;
  char tmpLine[500];
  string tmpStr;
  vector <string> strVec;       // used in reading files
  string tmp;
  ifstream inFile;

  //strcpy(filename,FILE_SIMPL);

  inFile.open(filename.c_str(), ios::in);
  
  if(!inFile.is_open()) {
      cout<<"cannot open VDW potential file! "<<filename<<endl;
      exit(0);
  }
  cout << "intialize   VDW	" << filename <<endl;
					
  while(!inFile.eof()){
    inFile.getline(tmpLine,500);
    tmpStr=tmpLine;
	if(tmpStr[0] == '#') continue;
    if(tmpStr.length()<3) continue;
    strVec.clear();
    split(tmpStr," ",strVec);
    i = atoi(strVec[0].c_str());
    j = atoi(strVec[1].c_str());
    k = atoi(strVec[2].c_str());
    pt1 = atof(strVec[3].c_str());
    pt2 = atof(strVec[4].c_str());
   
    VDWA[i][j][k] = pt1;
    VDWA[j][i][k] = pt1;
    VDWR[i][j][k] = pt2;
    VDWR[j][i][k] = pt2;
  }
  inFile.close();
}


// calculate hbond energy from all atoms involved.
// The atoms in the pntArr are Donor, H, Acceptor, Acceptor Base and the atom preceeding Acceptor Base.
double hb_energy(Point* pntArr, int type, int tp_ss, int type_sp, int tp_en, double dis)
{
  int type_dis=0;
  double e_hb=0,chi,theta,psi;
  Point p_hbondH=pntArr[1],p_donor=pntArr[0],p_acc=pntArr[2],p_ab=pntArr[3],p_ab_prev=pntArr[4];

  if(dis<2.1)
    type_dis=0; // short
  else type_dis=1; // long
 	
  // calculate the angles
  chi=torsion(p_ab_prev,p_ab,p_acc,p_hbondH);
  theta=angle(p_acc,p_hbondH,p_donor);
  psi=angle(p_ab,p_acc,p_hbondH);
 
  int dis_label;
  
  dis_label=3+type_sp;
  e_hb+=PF::HB_dis[dis_label][int((dis-1.4)/0.05)];
  
  // angle terms
  // prevent array subscripts from overflow
  if(chi==180) chi-=0.0001;
  if(psi==180) psi-=0.0001;
  if(theta==180) theta-=0.0001;

  if(type_sp==0){ // sp2
    if(type_dis==0) { // short
      // chi
      e_hb+=PF::HB_angle_sc[0][int(abs(chi)/10)];
      // psi
      e_hb+=PF::HB_angle_sc[1][int(psi/10)];
      // theta
      e_hb+=PF::HB_angle_sc[2][int(theta/10)];	    
    }
    else { // long
      // chi
      e_hb+=PF::HB_angle_sc[3][int(abs(chi)/10)];
      // psi
      e_hb+=PF::HB_angle_sc[4][int(psi/10)];
      // theta
      e_hb+=PF::HB_angle_sc[5][int(theta/10)];	    
    }
  }
  else {  // sp3, no chi
    if(type_dis==0) { // short
      // psi
      e_hb+=PF::HB_angle_sc[6][int(psi/10)];
      // theta
      e_hb+=PF::HB_angle_sc[7][int(theta/10)];	  
    }
    else { // long
      // psi
      e_hb+=PF::HB_angle_sc[8][int(psi/10)];
      // theta
      e_hb+=PF::HB_angle_sc[9][int(theta/10)];	  
    }
  }
  return e_hb;
}

// calculate vdw and solvation force given coordinates of twp atoms and their type, if r_ij == 2.95 then it is the case of hbond
double vdw_sol(Point& p1, Point& p2, double dis, double r_ij,int tp1,int tp2) {
  double e_atr=0,e_rep=0,e_sol=0,e_tmp,e_ij,slope,y_intercept;

  // VDW

  if(r_ij!=2.95){
    r_ij=Atom::radius[tp1]+Atom::radius[tp2];
  }
  e_ij=sqrt(Atom::welldepth[tp1]*Atom::welldepth[tp2]);
  if((r_ij/dis)<1.12){ // attraction part
    e_tmp=(pow(r_ij/dis,12)-2*pow(r_ij/dis,6))*e_ij;
    e_atr+=e_tmp;
  }

  // repulsion part
  else if((r_ij/dis)<1.33){
    e_tmp=(pow(r_ij/dis,12)-2*pow(r_ij/dis,6))*e_ij;
    e_rep+=e_tmp;
  }
  else {
    slope=-12*e_ij*(33.383)*(1/r_ij);  // 33.383=(1.33^13-1.33^7)
    y_intercept=-1*slope*(r_ij/1.33)+e_ij*19.565; // 19.565=(1.33^12-2(1.33^6))
    e_tmp=(y_intercept + dis*slope);
    e_rep+=e_tmp;
  }

  // Solvation

  r_ij=Atom::radius[tp1]+Atom::radius[tp2];
  e_tmp=(-1*(Atom::s_dgfree[tp1]*pow(EXPO,-(double)dis*dis)*Atom::s_volume[tp2]/(2*pow((double)PI,(double)1.5)*Atom::s_lambda[tp1]*r_ij*r_ij)+Atom::s_dgfree[tp2]*pow(EXPO,(double)-dis*dis)*Atom::s_volume[tp1]/(2*pow(PI,1.5)*Atom::s_lambda[tp2]*r_ij*r_ij)));

  e_sol+=e_tmp;

  return PF::Parameter["E"][E_SOL] * e_sol + PF::Parameter["E"][E_VDWA] * e_atr + PF::Parameter["E"][E_VDWR] * e_rep;


}

double pnorm(double x)
{
  double result;
  result = 0.5*(1+sqrt(1-pow(EXPO,(-2*x*x/PI))));
  return result;
}

double pnorm(double x, double mean, double sd)
{
  double v;

  if(x > mean) {
    v = (x - mean)/sd;
    return 0.5*(1+sqrt(1-pow(EXPO,(-2*v*v/PI))));
  }
  else { 
    v = (mean - x)/sd;
    return 1 - 0.5*(1+sqrt(1-pow(EXPO,(-2*v*v/PI))));
  }
}

void mapAtoms(SIMAP& aiMap) {

  // insert N, CA, C, O, H into aiMap
  aiMap.insert(SIMAP::value_type("N",0));
  aiMap.insert(SIMAP::value_type("CA",1));
  aiMap.insert(SIMAP::value_type("C",2));
  aiMap.insert(SIMAP::value_type("O",3));
  aiMap.insert(SIMAP::value_type("H",4));
	    
  aiMap.insert(SIMAP::value_type("ACB",5));
  aiMap.insert(SIMAP::value_type("CCB",6));
  aiMap.insert(SIMAP::value_type("CSG",7));
  aiMap.insert(SIMAP::value_type("DCB",8));
  aiMap.insert(SIMAP::value_type("DCG",9));
  aiMap.insert(SIMAP::value_type("DOD1",10));
  aiMap.insert(SIMAP::value_type("DOD2",11));
  aiMap.insert(SIMAP::value_type("ECB",12));
  aiMap.insert(SIMAP::value_type("ECG",13));
  aiMap.insert(SIMAP::value_type("ECD",14));
  aiMap.insert(SIMAP::value_type("EOE1",15));
  aiMap.insert(SIMAP::value_type("EOE2",16));
  aiMap.insert(SIMAP::value_type("FCB",17));
  aiMap.insert(SIMAP::value_type("FCG",18));
  aiMap.insert(SIMAP::value_type("FCD1",19));
  aiMap.insert(SIMAP::value_type("FCD2",20));
  aiMap.insert(SIMAP::value_type("FCE1",21));
  aiMap.insert(SIMAP::value_type("FCE2",22));
  aiMap.insert(SIMAP::value_type("FCZ",23));
  aiMap.insert(SIMAP::value_type("HCB",24));
  aiMap.insert(SIMAP::value_type("HCG",25));
  aiMap.insert(SIMAP::value_type("HND1",26));
  aiMap.insert(SIMAP::value_type("HCD2",27));
  aiMap.insert(SIMAP::value_type("HCE1",28));
  aiMap.insert(SIMAP::value_type("HNE2",29));
  aiMap.insert(SIMAP::value_type("HHD1",30));
  aiMap.insert(SIMAP::value_type("HHE2",31));
  aiMap.insert(SIMAP::value_type("ICB",32));
  aiMap.insert(SIMAP::value_type("ICG1",33));
  aiMap.insert(SIMAP::value_type("ICG2",34));
  aiMap.insert(SIMAP::value_type("ICD1",35));
  aiMap.insert(SIMAP::value_type("KCB",36));
  aiMap.insert(SIMAP::value_type("KCG",37));
  aiMap.insert(SIMAP::value_type("KCD",38));
  aiMap.insert(SIMAP::value_type("KCE",39));
  aiMap.insert(SIMAP::value_type("KNZ",40));
  aiMap.insert(SIMAP::value_type("LCB",41));
  aiMap.insert(SIMAP::value_type("LCG",42));
  aiMap.insert(SIMAP::value_type("LCD1",43));
  aiMap.insert(SIMAP::value_type("LCD2",44));
  aiMap.insert(SIMAP::value_type("MCB",45));
  aiMap.insert(SIMAP::value_type("MCG",46));
  aiMap.insert(SIMAP::value_type("MSD",47));
  aiMap.insert(SIMAP::value_type("MCE",48));
  aiMap.insert(SIMAP::value_type("NCB",49));
  aiMap.insert(SIMAP::value_type("NCG",50));
  aiMap.insert(SIMAP::value_type("NOD1",51));
  aiMap.insert(SIMAP::value_type("NND2",52));
  aiMap.insert(SIMAP::value_type("N1HD2",53));
  aiMap.insert(SIMAP::value_type("N2HD2",54));
  aiMap.insert(SIMAP::value_type("PCB",55));
  aiMap.insert(SIMAP::value_type("PCG",56));
  aiMap.insert(SIMAP::value_type("PCD",57));
  aiMap.insert(SIMAP::value_type("QCB",58));
  aiMap.insert(SIMAP::value_type("QCG",59));
  aiMap.insert(SIMAP::value_type("QCD",60));
  aiMap.insert(SIMAP::value_type("QOE1",61));
  aiMap.insert(SIMAP::value_type("QNE2",62));
  aiMap.insert(SIMAP::value_type("Q1HE2",63));
  aiMap.insert(SIMAP::value_type("Q2HE2",64));
  aiMap.insert(SIMAP::value_type("RCB",65));
  aiMap.insert(SIMAP::value_type("RCG",66));
  aiMap.insert(SIMAP::value_type("RCD",67));
  aiMap.insert(SIMAP::value_type("RNE",68));
  aiMap.insert(SIMAP::value_type("RCZ",69));
  aiMap.insert(SIMAP::value_type("RNH1",70));
  aiMap.insert(SIMAP::value_type("RNH2",71));
  aiMap.insert(SIMAP::value_type("RHE",72));
  aiMap.insert(SIMAP::value_type("R1HH1",73));
  aiMap.insert(SIMAP::value_type("R2HH1",74));
  aiMap.insert(SIMAP::value_type("R1HH2",75));
  aiMap.insert(SIMAP::value_type("R2HH2",76));
  aiMap.insert(SIMAP::value_type("SCB",77));
  aiMap.insert(SIMAP::value_type("SOG",78));
  aiMap.insert(SIMAP::value_type("TCB",79));
  aiMap.insert(SIMAP::value_type("TOG1",80));
  aiMap.insert(SIMAP::value_type("TCG2",81));
  aiMap.insert(SIMAP::value_type("VCB",82));
  aiMap.insert(SIMAP::value_type("VCG1",83));
  aiMap.insert(SIMAP::value_type("VCG2",84));
  aiMap.insert(SIMAP::value_type("WCB",85));
  aiMap.insert(SIMAP::value_type("WCG",86));
  aiMap.insert(SIMAP::value_type("WCD1",87));
  aiMap.insert(SIMAP::value_type("WCD2",88));
  aiMap.insert(SIMAP::value_type("WNE1",89));
  aiMap.insert(SIMAP::value_type("WCE2",90));
  aiMap.insert(SIMAP::value_type("WCE3",91));
  aiMap.insert(SIMAP::value_type("WCZ2",92));
  aiMap.insert(SIMAP::value_type("WCZ3",93));
  aiMap.insert(SIMAP::value_type("WCH2",94));
  aiMap.insert(SIMAP::value_type("WHH2",95));
  aiMap.insert(SIMAP::value_type("YCB",96));
  aiMap.insert(SIMAP::value_type("YCG",97));
  aiMap.insert(SIMAP::value_type("YCD1",98));
  aiMap.insert(SIMAP::value_type("YCD2",99));
  aiMap.insert(SIMAP::value_type("YCE1",100));
  aiMap.insert(SIMAP::value_type("YCE2",101));
  aiMap.insert(SIMAP::value_type("YCZ",102));
  aiMap.insert(SIMAP::value_type("YOH",103));

}
