/* atom.cpp
*/

#include "atom.h"

//// Stuff from the Point class ////
// absolute value of a point(or vector)
double Point::pabs() {  return sqrt(x*x+y*y+z*z);}
// ^ is used as dot operator
double Point::operator^(const Point& v) {  return x*v.x+y*v.y+z*v.z;}
// distance to another point
double Point::dis(const Point& v){  return sqrt((x-v.x)*(x-v.x)+(y-v.y)*(y-v.y)+(z-v.z)*(z-v.z));}
double Point::disquare(const Point& v){  return ((x-v.x)*(x-v.x)+(y-v.y)*(y-v.y)+(z-v.z)*(z-v.z));}
void Point::minus(double d) {  x=x-d;  y=y-d;  z=z-d;}
void Point::plus(double d) {  x=x+d;    y=y+d;    z=z+d;}
void Point::multiply(double d) {    x*=d;    y*=d;    z*=d;}
void Point::divide(double d) {  x/=d; y/=d; z/=d;}
 // output coordinates x,y,z
void Point::out() {
  cout << "(" << x << ", " << y << ", " << z << ") ";
}

double Point::angle(Point& p1, Point& p2) {
  double d12 = this->dis(p1);
  double d23 = this->dis(p2);
  double dot, angle;
  dot = (p1-(*this))^((*this)-p2);
  angle = 180 * acos(dot/(d12*d23))/PI;
  return (180-angle);
}

double Point::torsion(Point& p1, Point& p2, Point& p3) {
  double xij,yij,zij,
    xkj,ykj,zkj,
    xkl,ykl,zkl,
    dxi,dyi,dzi,
    gxi,gyi,gzi,
    bi,bk,ct,
    boi2,boj2,
    z1,z2,ap,s,
    bioj,bjoi;
  double xi=p1.x,yi=p1.y,zi=p1.z;
  double xj=p2.x,yj=p2.y,zj=p2.z;
  double xk=p3.x,yk=p3.y,zk=p3.z;
  double xl=x,yl=y,zl=z;

  /* Calculate the vectors C,B,C                                       */
  xij = xi - xj;
  yij = yi - yj;
  zij = zi - zj;
  xkj = xk - xj;
  ykj = yk - yj;
  zkj = zk - zj;
  xkl = xk - xl;
  ykl = yk - yl;
  zkl = zk - zl;

  dxi = yij * zkj - zij * ykj;     /* Normal to plane 1                */
  dyi = zij * xkj - xij * zkj;
  dzi = xij * ykj - yij * xkj;
  gxi = zkj * ykl - ykj * zkl;     /* Mormal to plane 2                */
  gyi = xkj * zkl - zkj * xkl;
  gzi = ykj * xkl - xkj * ykl;

   /* Calculate the length of the two normals                           */
  bi = dxi * dxi + dyi * dyi + dzi * dzi;
  bk = gxi * gxi + gyi * gyi + gzi * gzi;
  ct = dxi * gxi + dyi * gyi + dzi * gzi;

  boi2 = 1./bi;
  boj2 = 1./bk;
  bi   = (double)sqrt((double)bi);
  bk   = (double)sqrt((double)bk);

  z1   = 1./bi;
  z2   = 1./bk;
  bioj = bi * z2;
  bjoi = bk * z1;
  ct   = ct * z1 * z2;
  if (ct >  1.0)   ct = 1.0;
  if (ct < (-1.0)) ct = -1.0;
  ap   = acos(ct);

  s = xkj * (dzi * gyi - dyi * gzi)
    + ykj * (dxi * gzi - dzi * gxi)
    + zkj * (dyi * gxi - dxi * gyi);
  
  if (s < 0.0) ap = -ap;

  ap = (ap > 0.0) ? PI-ap : -(PI+ap);
  
  return ap / PI * 180;
}

// Stuff for the Atom class//
double Atom::vdw_adj=1.0;
double Atom::radius[MAX_ATOM_TYPE];
double Atom::welldepth[MAX_ATOM_TYPE];
double Atom::s_volume[MAX_ATOM_TYPE];
double Atom::s_lambda[MAX_ATOM_TYPE];
double Atom::s_dgfree[MAX_ATOM_TYPE];
short Atom::acceptor[MAX_ATOM_TYPE];
short Atom::donor[MAX_ATOM_TYPE];
short Atom::hbondH[MAX_ATOM_TYPE];
double Atom::R_SC[NUM_RES_TYPE] = {1.53, 2.14, 2.49, 3.16, 3.42, 0.0, 3.17, 2.35, 3.57, 2.65, 3.00, 2.52, 1.88, 3.15, 4.17, 1.95, 1.95, 1.97, 3.89, 3.78, 0,0,0,0,0,0,0,0,0,0};  // new assignment calculated recently

// Simple initialization
Atom::Atom(double X, double Y, double Z, int I) {
  init(X, Y, Z, I);
}
Atom::Atom(const Point& P) {
  init(P.x, P.y, P.z);
}
Atom::Atom(const Atom& A) {
  init();
  *this = A;
}
void Atom::init(double X, double Y, double Z, int I) {
  x = X;
  y = Y;
  z = Z;
  _type = I;
  _Bfactor = -1;
  _state = 0;
}

// This is the full assignment operator, copying all members
Atom* Atom::operator=(const Atom& A) {
  x = A.x;
  y = A.y;
  z = A.z;
  _type = A._type;
  _name = A._name;
  _Bfactor = A._Bfactor;
  _posn = A._posn;
  _state = A._state;
  return this;
}
Atom* Atom::operator=(const Point& P) {
  x = P.x;
  y = P.y;
  z = P.z;
  return this;
}

// This is a limited assignment operator, copying only the spatial position
void Atom::CopyPos(const Atom& A) {
  x = A.x;
  y = A.y;
  z = A.z;
}

// Set the parent of an Atom to the input residue
void Atom::SetParent(Residue& R) {
  _parent = &R;
}

void Atom::reset() {
  init();
}

void Atom::out() {
  cout << "(" << x << ", " << y << ", " << z << ", " << _type << ") ";
}

bool Atom::AnyHBond() {
  // Tests whether this atom has any sort of HBond
  return (HGiven.size() > 0 || HTaken.size() > 0);
}

bool Atom::ShareHBond(const Atom& A, const Atom& B) {
  // Tests whether the two input atoms share any HBonds
  return (DonatedHBond(A, B) || DonatedHBond(B, A));
}

bool Atom::DonatedHBond(const Atom& A, const Atom& B) {
  // Tests whether Atom A has donated an H to Atom B in an HBond
  if (A.HGiven.size() > 0 && B.HTaken.size() > 0)
    // This atom has donated to an HBond, and A has accepted an HBond
    for (int i = 0; i < B.HTaken.size(); i++)
      if (A.HGiven[0] == B.HTaken[i])
	return true;
}

void Atom::InitPar(char* parFile, string& dir) {
  int i,numLine;
  ifstream inFile;
  char tempCh[200],tmpLine[500];
  cout<<"Reading Atom and Residue Parameters in "<<parFile<<endl;
  inFile.open(parFile,ios::in);
  if(!inFile.is_open()){
    strcpy(tempCh,dir.c_str());
    strcat(tempCh,parFile);
    inFile.open(tempCh,ios::in);    
    if(!inFile.is_open()) {  
      cout<<"cannot open parameter file "<<tempCh<<endl;
      exit(0);
    }
  }

  // read parameters on atoms
  inFile.clear();
  inFile.seekg(0,ios::beg);
  while(!inFile.eof()){
    inFile.getline(tmpLine,1000);
    sscanf(tmpLine,"%*s %s %d",tempCh,&numLine);
    if(strcmp(tempCh,"atom_parameters")==0){ // read vdw parameters
      for(i=0;i<numLine;++i){
        inFile.getline(tmpLine,1000);
		sscanf(tmpLine, "%lf, %lf, %lf, %lf, %lf, %d %d %d",
	       &Atom::radius[i+1],
	       &Atom::welldepth[i+1],
	       &Atom::s_volume[i+1],
	       &Atom::s_lambda[i+1],
	       &Atom::s_dgfree[i+1],
	       (int*)&Atom::acceptor[i+1],
	       (int*)&Atom::donor[i+1],
	       (int*)&Atom::hbondH[i+1]);
	Atom::radius[i+1]=Atom::vdw_adj*Atom::radius[i+1];
      }      
    }
    else continue;
  }
  inFile.close();
}
