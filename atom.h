/* atom.h
The basic header file for all programs.
Class Atom inherited from Point.

*/

#ifndef _ATOM_
#define _ATOM_

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <string>
#include <sstream>
#include <cctype>
#include <set>
#include <map>
#include <list>
#include <vector>
#include <fstream>
#include <iomanip>
#include <cstdio>
#include <ctime>
#include <new>
#include <cstring>
#include <math.h>

using namespace std;
class Residue;

#define MAX_ATOM_TYPE 30  // maximum number of atom types allowed
#define MAX_NUM_ATOM_RES 30
#define NUM_RES_TYPE 30
#define MAX_NUM_RES 1500
#define MAX_HBOND_PER_DONOR 1
#define MAX_HBOND_PER_ACCEPTOR 2
#define MAX_DISTANCE_HBOND 3.0
#define PI 3.14159265358979
#define EXPO 2.718281828
#define UNDEF -12345   // undefined value
#define DEFINED 1   // defined value
#define MAX_E_CUTOFF 10 // maximum energy value for atom-atom clash
#define intrand(a,b) ((int)(a+(b-a)*((double)rand())/RAND_MAX)) // [a, b-1]
#define frand(a,b) ((rand()+.5)/RAND_MAX*(b-a)+a)
#define drand(a,b) (((double)rand()+.5)/RAND_MAX*(b-a)+a)
#define MAX(a,b) (a>=b? a : b)
#define MIN(a,b) (a<=b? a : b)
#define BBTbinSize 5
#define TORBIN 360/BBTbinSize

// defined types
typedef set <int> ISET;                  // set of integers
typedef set <string> SSET;               // set of strings
typedef list <int> ILIST;                // list of integers
typedef vector <double> FVEC;             // vector of double
typedef vector <double> DVEC;            // vector of double
typedef vector <int> IVEC;               // vector of int
typedef vector <string> SVEC;            // vector of string
typedef map <int, int> IIMAP;            // integer to integer map
typedef map <long, vector<int> > LVMAP;  // map of long int and vector of int
typedef map <int, double> IFMAP;          // map of int and double
typedef map <int, double> IDMAP;         // map of int and double
typedef map <double, int> FIMAP;          // double to integer map
typedef map <double, int> DIMAP;         // double to integer map
typedef map <long, string> ISMAP;        // long to string map
typedef map <string, string> SSMAP;      // string to string map
typedef map <string, int> SIMAP;         // string to int map
typedef map <string, double> SFMAP;       // string to double map
typedef map <string, double> SDMAP;      // string to double map
typedef map <int, vector<double> > IVMAP; // int to vector of double map 
typedef map <int, ISET> ISTMAP;          // map of int and integer set, ISET
typedef map <long long int, double> LLIFMAP; // long long int and double map

class Point {
 public:
  double x,y,z;                 // coordinates on 3D
  
 Point():x(0),y(0),z(0){}
  Point(double X,double Y,double Z) {x=X;y=Y;z=Z;}
  Point* operator=(const Point& r){ x=r.x;y=r.y;z=r.z; return this;}
  // operation can be carried on Point
  Point operator+(Point p){  return Point(x+p.x,y+p.y,z+p.z);}
  Point operator-(Point p){  return Point(x-p.x,y-p.y,z-p.z);}
  // cross product
  Point operator*(Point v){  return Point(y*v.z-z*v.y,z*v.x-x*v.z,x*v.y-y*v.x);}
  Point operator/(double d){return Point(x/d,y/d,z/d);}
  Point operator*(double d) { return Point(x*d,y*d,z*d);}
  void minus(double d);
  void plus(double d);
  void multiply(double d);
  void divide(double d);
  double square() { return x*x+y*y+z*z;}
  Point squarev() const{ return Point(x*x,y*y,z*z);}
  double sum() { return x+y+z;}
  // absolute value of a point(or vector)
  double pabs();
  // ^ is used as dot operator
  double operator^(const Point& v);
  // distance to another point
  double dis(const Point& v);
  double disquare(const Point& v);
  // output coordinates x,y,z
  void out();
  // randomize the coordinates between 0 and 1
  void randomizeCo(){x = drand(0,1); y = drand(0,1); z = drand(0,1); }
  // functions to calculate bond angle and torsion angles
  // the current point is at the end of the three or four points
  double angle(Point& , Point& );
  double torsion(Point& , Point& , Point&);
};

class Atom: public Point {
 public:
  short    _type;    // for van der wall radii
  string   _name;    // name of the atom
  double    _Bfactor; // store the temperature factor in pdb files
  int      _posn;    // position in the parent residue
  Residue* _parent;  // pointer to the parent residue
  short int   _state;         //This attribute can be used for multiple purposes
                              //0: default
			      //1: this atom is missing from native

  // HBond-related terms
  vector<Atom*> HGiven; // Atoms donated H to (max. 1 right now)
  vector<Atom*> HTaken; // Atoms taken H from (max. 2 right now)

  // parameters for different atoms are declared as static
  static double vdw_adj;        // adjust the vdw radius
  static double radius[MAX_ATOM_TYPE];
  static double welldepth[MAX_ATOM_TYPE];
  static double s_volume[MAX_ATOM_TYPE];
  static double s_lambda[MAX_ATOM_TYPE];
  static double s_dgfree[MAX_ATOM_TYPE];
  static short acceptor[MAX_ATOM_TYPE];
  static short donor[MAX_ATOM_TYPE];
  static short hbondH[MAX_ATOM_TYPE];
  static double R_SC[NUM_RES_TYPE];

  Atom(double X = 0., double Y = 0., double Z = 0., int I = UNDEF);
  Atom(const Point& P);
  Atom(const Atom& A);
  void init(double X = 0., double Y = 0., double Z = 0., int I = UNDEF);
  Atom* operator=(const Atom& A);
  Atom* operator=(const Point& P);
  void CopyPos(const Atom& A);
  void reset();
  static void InitPar(char*,string&);
  void out();
  void SetParent(Residue& R);
  bool AnyHBond();
  static bool ShareHBond(const Atom& A, const Atom& B);
  static bool DonatedHBond(const Atom& A, const Atom& B);
};

// This class represents an HBond as a pair of atoms
class HBond {
 public:
  Atom* DonorAtom;
  Atom* AcceptorAtom;
};

typedef vector <Point> PTVEC;
typedef vector <Atom> ATMVEC;

#endif
