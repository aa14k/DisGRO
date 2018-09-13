// util.h
// header file for util.cpp

#ifndef _UTIL_
#define _UTIL_

#include <sstream>
#include "atom.h"
#include "structure.h"
#include "cal_energy.h"

struct combiEIndex{
  double energy;
  int index;
};

double angle(Point,Point,Point);
double torsion(Point a, Point b, Point c, Point d);
int split(const string& input, const string& delimiter, vector<string>& results);
void insert_p_lvmap(Point& tmpP,LVMAP& mcMap, int key, double cub_size);
string itoa(int num);
string ftoa(double num);
void calCo(Point* prev_atoms, double length, double bAngle, double tAngle, Point& n_res);
void calCo(Atom* prev_atoms, double, double, double, Atom&);
// calculate backbone torsion angles
void cal_dha(Structure& Conf,int,int);
// calculate the side chain torsion angles
void cal_sc_angles(Structure& Conf, int, int);
// calculate atom-atom distances
void cal_aadis(Structure& Conf, int start, int end, char* filename, double disCutoff, int cType);
template<class T> void quicksort(T a[], const int& leftarg, const int& rightarg);
double rmsd(Structure& S1, Structure& S2, int type, int normal);
double rmsd(Structure& S1, Structure& S2, int Start1, int End1, int Start2,
	   int End2, int type,int normal);
double rmsd(Structure& S1,Structure& S2, int,int,int,int,int type,int normal, double**,double*,double**);
double Rg(Structure& Conf);
vector<int> ExpandNumStates(vector<int> NumStates, int FragLength, string type);
int SampleOne(vector<double> prob);
double VectorSum(vector<double> vec);
vector<string> FileLines(string inFile);
int BucketInd(double element, vector<double> bucket);
string File2ProtName(string fname);
bool FileExists(string fname);
double Root_MSD(Structure& A,Structure& B, int,int,int,int, int, int);
double Normal_RMS(vector<Atom> A, vector<Atom>B);
double Root_MSD(vector<Residue>& A, vector<Residue>& B,int);
double MSD(Structure& A, Structure& B, int Start1, int End1, int Start2,
           int End2, int& size,  int type, int normal);
double entropy(double* Prob, int size);
void box_MullerNsample(double normal_Svalue[],double mean[], int size, double sigma);
void box_MullerNsample_single(double& normal_Svalue,double mean, double sigma);
double distance_to_line(Point begin, Point end, Point x );
void MeanSQ(Residue& A,Residue& B, double& sq_Sum, int& size, int type);
void Loop_DecoyPDB_read(string PathparPDB, vector < vector <string> >& TMPCC, int conf_Upbound);
void SCE_Minimization_list(Structure& conf, int start, int end, vector<int>& ResIdx, vector<int>& ClashNum,  int* List, int List_size, double change_rot);
double Ellipsoid_Detect (Structure& Conf, int Start, int End, bool type);
void read_Surface(string surfaceFile, vector< vector<int> >&SurfaceList,Structure& Conf, int outInfo);

  bool sortfun_E (combiEIndex x, combiEIndex y);
  bool sortfun_R (combiEIndex x, combiEIndex y);
  bool sortfun_StrE(Structure x, Structure y);
//  bool sortfun_LoostE(Structure x, Structure y);

template <typename T>
std::string to_string(T const& value) {
  stringstream sstr;
  sstr << value;
  return sstr.str();
}

#endif
