// sc_util.cpp
// useful functions for side-chain modeling program

#include "util.h"
#include "structure.h"
#include "matrix.h"

// angle formed by ab and bc
double angle(Point a, Point b, Point c) {
  //  return(acos(((a - b)^(c - b)) / (a.dis(b) * c.dis(b))) * 180. / PI);
  double d12 = b.dis(a);
  double d23 = b.dis(c);
  double dot,angle;
  dot = (a-b)^(b-c);
  angle = 180 * acos(dot/(d12*d23))/PI;
  return (180-angle);
}

// torsion angle formed by four points a,b,c,and d.
double torsion(Point a, Point b, Point c, Point d) {
  double xij,yij,zij,
    xkj,ykj,zkj,
    xkl,ykl,zkl,
    dxi,dyi,dzi,
    gxi,gyi,gzi,
    bi,bk,ct,
    boi2,boj2,
    z1,z2,ap,s,
    bioj,bjoi;
  double xi=a.x,yi=a.y,zi=a.z;
  double xj=b.x,yj=b.y,zj=b.z;
  double xk=c.x,yk=c.y,zk=c.z;
  double xl=d.x,yl=d.y,zl=d.z;


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

  /* Calculate the normals to the two planes n1 and n2
     this is given as the cross products:
     AB x BC
     --------- = n1
     |AB x BC|

     BC x CD
     --------- = n2
     |BC x CD|
  */
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

 // return(ap);
  return((ap/PI)*180);
}

// this split function take a string and split it into pieces
// separated by delimiter and the results is stored in results. It
// returns the number of strings found which is different from
// delimiter.
int split(const string& Input, const string& delimiter,
	  vector<string>& results) {

  int i,isDel=0;
  int sizeS2 = delimiter.size();
  string s,del;
  vector<int> positions;
  string input=Input;
  del=delimiter;
  if(input.size()==0)
    return 0;

  if(sizeS2==1 || (sizeS2>1&&input.substr(0,1) == "1")) {
    // only one character in delimiter, simple case
    if(sizeS2>1){
      del=delimiter.substr(1,delimiter.length()-1);
      sizeS2--;
    }
    i=0;
    input=input+del;
    if(input.substr(i,sizeS2)==del){
      isDel=1;
      i=i+sizeS2;
    }
    else {
      positions.push_back(i);
    }
    while(i<input.size()){
      if(isDel==0){
	if(input.substr(i,sizeS2)==del){
	  positions.push_back(i);
	  isDel=1;
	}
      }
      else if(input.substr(i,sizeS2)!=del){
	positions.push_back(i);
	isDel=0;
      }
      i=i+sizeS2;
    }
    
    if(positions.size()==0)
      return 0;
    
    if(positions.size()==1){
      if(isDel==0){
	results.push_back(input);
	return 1;
      }
      else {
	return 0;
      }
    }

    //if(positions.size()%1==1){// odd number
    //  positions.push_back(input.size());
    //}
    //cout<<"size "<<positions.size()<<" ";
    for( i=0; i < positions.size(); i=i+2 ){
      s=input.substr(positions[i],positions[i+1]-positions[i]);
      results.push_back(s);
    }
    
    return results.size();
  }
  else {  // more than one character in delimiter
    del=delimiter.substr(1,delimiter.length()-1);
    input=input+del;
    sizeS2=1;
    i=0;
    if(del.find(input.substr(i,sizeS2),0)!=string::npos){
      isDel=1;
      i=i+sizeS2;
    }
    else {
      positions.push_back(i);
    }
    while(i<input.size()){
      if(isDel==0){
	if(del.find(input.substr(i,sizeS2),0)!=string::npos){
	  positions.push_back(i);
	  isDel=1;
	}
      }
      else if(del.find(input.substr(i,sizeS2),0)==string::npos){
	positions.push_back(i);
	isDel=0;
      }
      i=i+sizeS2;
    }
    
    if(positions.size()==0)
      return 0;
    
    if(positions.size()==1){
      if(isDel==0){
	results.push_back(input);
	return 1;
      }
      else {
	return 0;
      }
    }
       
    for( i=0; i < positions.size(); i=i+2 ){
      s=input.substr(positions[i],positions[i+1]-positions[i]);
      results.push_back(s);
    }
    
    return results.size();
  }
}

// insert a point (27 integers) and its associated key to a LVMAP, <long, vector<int> > 
void insert_p_lvmap(Point& tmpP,LVMAP& mcMap, int key, double cub_size=5.5)
{
  int l,m,n,index;
  LVMAP::iterator mcItr;
  vector<int> *tmpVec;
  // insert key into mcMap
  for(l=-1;l<=1;++l){
    for(m=-1;m<=1;++m){
      for(n=-1;n<=1;++n){
        index=(int)(l+(tmpP.x/cub_size+CUB_ADJ))*1000000+(int)(m+(tmpP.y/cub_size+CUB_ADJ))*1000+(int)(n+(tmpP.z/cub_size+CUB_ADJ));
        if((mcItr=mcMap.find(index))!=mcMap.end()){
          mcItr->second.push_back(key);
        }
        else {
          tmpVec = new vector<int>;
          tmpVec->push_back(key);
          mcMap.insert(LVMAP::value_type(index,(*tmpVec)));
        } // else
      } // for(int n=-1;
    } // for(int m=-1;
  }  // for(int l=-1;
}

string itoa(int num) {
  stringstream converter;
  converter << num;
  return converter.str(); 
}

string ftoa(double num) {
  stringstream converter;
  converter << num;
  return converter.str(); 
}

// calculate position of n_res given the previous positions of three atoms and the bond length, bond angle and torsion angle
void calCo(Point* prev_atoms, double length, double bAngle, double tAngle,
	   Point& n_res) {
  Point su,u2,u3,_SvdV;
  Point _u,_v,_w,cur_res,last_res,bfl_res;
  double d,dis2;
  cur_res  = prev_atoms[2];
  last_res = prev_atoms[1];
  bfl_res  = prev_atoms[0];
  _SvdV    = bfl_res-last_res;
  su       = cur_res-last_res;
  d        = sqrt(su.x*su.x+su.y*su.y+su.z*su.z);
  _u       = su/d;
  d        = _SvdV^_u;
  u3       = _u*d;
  dis2     = _SvdV.dis(u3);
  _v       = (_SvdV-u3)/dis2;
  _w       = _u*_v;
  n_res    = cur_res + _u*length*cos(PI-bAngle)
    + _v*length*sin(PI-bAngle)*cos(tAngle)
    + _w*length*sin(PI-bAngle)*sin(tAngle);

  //double tor = torsion(prev_atoms[0],prev_atoms[1],prev_atoms[2],n_res);
  //cout<<"torsion "<<tor<<endl;
  /*
    if(isnan(n_res.x)){	
    cout<<"calCo ";
    cur_res.out();
    last_res.out();
    bfl_res.out();
    cout<<" "<<length<<" "<<bAngle<<" "<<tAngle<<" "<<dis2<<" "<<d<<endl;
    }
  */
}

// calCo for atoms
void calCo(Point* prev_atoms, double length, double bAngle, double tAngle,
	   Atom& n_res) {
  Point su,u2,u3,_SvdV;
  Point _u,_v,_w,cur_res,last_res,bfl_res;
  double d,dis2;
  cur_res  = prev_atoms[2];
  last_res = prev_atoms[1];
  bfl_res  = prev_atoms[0];
  _SvdV    = bfl_res-last_res;
  su       = cur_res-last_res;
  d        = sqrt(su.x*su.x+su.y*su.y+su.z*su.z);
  _u       = su/d;
  d        = _SvdV^_u;
  u3       = _u*d;
  dis2     = _SvdV.dis(u3);
  _v       = (_SvdV-u3)/dis2;
  _w       = _u*_v;
  n_res    = cur_res + _u*length*cos(PI-bAngle)
    + _v*length*sin(PI-bAngle)*cos(tAngle)
    + _w*length*sin(PI-bAngle)*sin(tAngle);
}

// calculate backbone dihedral angles, phi, psi and omega
void cal_dha(Structure& Conf,int Start,int End)
{
  int i,j;
  double phi,psi;

  // calculate phi and psi angles without considering the number of chains
  if(Start==1){
    Conf._res[1]._phi=999.99;
    Conf._res[1]._psi=torsion(Conf._res[1]._atom[0],Conf._res[1]._atom[1],Conf._res[1]._atom[2],Conf._res[2]._atom[0]);
    Conf._res[1]._omega=torsion(Conf._res[1]._atom[1],Conf._res[1]._atom[2],Conf._res[2]._atom[0],Conf._res[2]._atom[1]);
    Start++;
  }
  for(i=Start;i<=End;++i){
    if(i==Conf._numRes) continue;
    else Conf._res[i]._phi = torsion(Conf._res[i-1]._atom[2],Conf._res[i]._atom[0],Conf._res[i]._atom[1],Conf._res[i]._atom[2]);
    if(Conf._res[i]._pdbIndex != (Conf._res[i+1]._pdbIndex-1)) {  // not continuous at this position
      Conf._res[i]._psi = 999.99;
      Conf._res[i]._omega = 999.99;
    }
    else {
      Conf._res[i]._psi = torsion(Conf._res[i]._atom[0],Conf._res[i]._atom[1],Conf._res[i]._atom[2],Conf._res[i+1]._atom[0]);
      Conf._res[i]._omega = torsion(Conf._res[i]._atom[1],Conf._res[i]._atom[2],Conf._res[i+1]._atom[0],Conf._res[i+1]._atom[1]);
    }
  if(End==Conf._numRes){
    j=Conf._numRes;
    Conf._res[j]._phi=torsion(Conf._res[j-1]._atom[2],Conf._res[j]._atom[0],Conf._res[j]._atom[1],Conf._res[j]._atom[2]);
    Conf._res[j]._psi=999.99;
    Conf._res[j]._omega=999.99;
  }
 }
}
// function for calculating side chain torsion angles
void cal_sc_angles(Structure& Conf, int start, int end)
{
  int i,j,p3,p2,p1,count;
  double tor,bondAngle;
  for(i=start;i<=end;++i)
    {
      count=0;
      for(j=NUM_BB_ATOM;j<Conf._res[i]._numAtom;++j)  // caution needs to be taken here since H atom may or may not be counted. If H is counted then this is fine, otherwise j<=Conf._res[i]._numAtom should be used instead
	{
	  //cout<<j<<" "<<Residue::torsion[Conf._res[i]._type][j]<<endl;
	  if(Residue::torsion[Conf._res[i]._type][j]==-1234)
	    { // calculate torsion angle
	      p1 = Residue::prev_atom[Conf._res[i]._type][j][0];
	      p2 = Residue::prev_atom[Conf._res[i]._type][j][1];
	      p3 = Residue::prev_atom[Conf._res[i]._type][j][2];
	      bondAngle = angle(Conf._res[i]._atom[j],Conf._res[i]._atom[p1],Conf._res[i]._atom[p2]);
	      if(Conf._res[i]._atom[j]._type!=UNDEF&&Conf._res[i]._atom[p3]._type!=UNDEF&&Conf._res[i]._atom[p2]._type!=UNDEF&&Conf._res[i]._atom[p1]._type!=UNDEF){
		tor = torsion(Conf._res[i]._atom[j],Conf._res[i]._atom[p3],Conf._res[i]._atom[p2],Conf._res[i]._atom[p1]);
		Conf._res[i]._scChi[count++] = tor;
	      }
	      else {
		Conf._res[i]._scChi[count++] = 999.99;
	      }
	    }	
	}  	
    }
}

// calculating atom-atom distance and consider each observation as
// normal distribution
void cal_aadis(Structure& Conf, int start, int end, char* filename,
	       double disCutoff, int cType=1) {
  int i, j, k, l, count=0;
  double dis, sumB=0, sumBsquare=0,Bfactor, sdB, avgB;
  ofstream fout;
  fout.open(filename);
  
  // normalize B-factors
  if(cType==0) {
    for(i=start;i<=end-2;++i) {
      for(k=0;Conf._res[i]._atom[k]._type!=UNDEF;++k) {
	if(k==4||Conf._res[i]._atom[k]._type>=22) continue;  // ignore H atoms
	Bfactor = Conf._res[i]._atom[k]._Bfactor;
	sumB += Bfactor;
	sumBsquare += (Bfactor*Bfactor);
	count++;
      }
    }
    
    avgB = sumB/count;
    // standard deviation of Bfactors
    sdB = sqrt((sumBsquare - count*avgB*avgB)/(count-1));
    
    for(i=start;i<=end-2;++i) {
      for(k=0;Conf._res[i]._atom[k]._type!=UNDEF;++k) {
	if(k==4||Conf._res[i]._atom[k]._type>=22) continue;  // ignore H atoms
	Bfactor = Conf._res[i]._atom[k]._Bfactor;
	Conf._res[i]._atom[k]._Bfactor = (Bfactor - avgB)/sdB;
      }
    }
  } // end of if(cType==0)
  
  // calculate distances
  if(cType == 1){
    for(i=start;i<=end-2;++i) {
      for(k=0;k<=Conf._res[i]._numAtom && k<MAX_NUM_ATOM_RES;++k) {
	if(Conf._res[i]._atom[k]._type==UNDEF) {continue;}
	if(k==4||Conf._res[i]._atom[k]._type>=22) continue; // ignore H atoms
	for(j=i+2;j<=end;++j) {
	  for(l=0;l<=Conf._res[j]._numAtom && l <  MAX_NUM_ATOM_RES;++l) {
	    if(Conf._res[j]._atom[l]._type==UNDEF) {continue;}
	    //cout<<i<<" "<<j<<" "<<k<<" "<<l<<" "<<Conf._res[i]._atom[k]._name<<" "<<Conf._res[j]._atom[l]._name<<endl;
	    if(l==4||Conf._res[j]._atom[l]._type>=22) continue; // ignore H atoms
	    dis=Conf._res[i]._atom[k].dis(Conf._res[j]._atom[l]);
	    if(dis<disCutoff) {
	      fout<<setw(6)<<Conf._res[i]._pdbIndex<<setw(4)<<Residue::Name3[Conf._res[i]._type];
	      fout<<setw(5)<<Conf._res[i]._atom[k]._name;
	      fout<<setw(6)<<Conf._res[j]._pdbIndex<<setw(4)<<Residue::Name3[Conf._res[j]._type];
	      fout<<setw(5)<<Conf._res[j]._atom[l]._name;
	      fout<<setw(8)<<setprecision(4)<<dis;
	      fout<<setw(8)<<Conf._res[i]._atom[k]._Bfactor;
	      fout<<setw(8)<<Conf._res[j]._atom[l]._Bfactor<<endl;
	    }
	  }
	}
      }
    }
  }
  else if(cType == 2) { // calculate distance for SC atoms (seudo-side chain atom)
    // calculate SC atom from all side chain atoms including CB and assign it to CB

    Point prev_atom[3]; 
    for(i=start;i<=end;++i) {
      if(Conf._res[i]._type == GLY || Conf._res[i]._type == ALA) continue;
      prev_atom[0] = Conf._res[i]._atom[ATM_N];
      prev_atom[1] = Conf._res[i]._atom[ATM_C];
      prev_atom[2] = Conf._res[i]._atom[ATM_CA];
      calCo(prev_atom,
	    Atom::R_SC[Conf._res[i]._type],
	    Residue::bond_angle[Conf._res[i]._type][ATM_CB],
	    PI * 122.55 / 180,
	    Conf._res[i]._SC);
      Conf._res[i]._atom[ATM_CB].CopyPos(Conf._res[i]._SC);
    }

    // calculate the distance of SC to all others and output to file
    for(i=start;i<=end-2;++i) {
      for(k=0;k<=ATM_CB;++k) {
	if(Conf._res[i]._atom[k]._type==UNDEF) {continue;}
	if(k==4||Conf._res[i]._atom[k]._type>=22) continue; // ignore H atoms
	for(j=i+2;j<=end;++j) {
	  for(l=0;l<=ATM_CB;++l) {
	    if(Conf._res[j]._atom[l]._type==UNDEF) {continue;}
	    //cout<<i<<" "<<j<<" "<<k<<" "<<l<<" "<<Conf._res[i]._atom[k]._name<<" "<<Conf._res[j]._atom[l]._name<<endl;
	    if(l==4||Conf._res[j]._atom[l]._type>=22) continue; // ignore H atoms
	    dis = Conf._res[i]._atom[k].dis(Conf._res[j]._atom[l]);
	    if(dis<disCutoff) {
	      fout<<setw(6)<<Conf._res[i]._pdbIndex<<setw(4)<<Residue::Name3[Conf._res[i]._type];
	      if(k==ATM_CB) fout<<setw(5)<<"SC";
	      else	fout<<setw(5)<<Conf._res[i]._atom[k]._name;
	      fout<<setw(6)<<Conf._res[j]._pdbIndex<<setw(4)<<Residue::Name3[Conf._res[j]._type];
	      if(l==ATM_CB) fout<<setw(5)<<"SC";
	      else	fout<<setw(5)<<Conf._res[j]._atom[l]._name;
	      fout<<setw(8)<<setprecision(4)<<dis;
	      fout<<setw(8)<<Conf._res[i]._atom[k]._Bfactor;
	      fout<<setw(8)<<Conf._res[j]._atom[l]._Bfactor<<endl;
	    }
	  }
	}
      }
    }
  }
}

template<class T> void quicksort(T a[], const int& leftarg,
				 const int& rightarg) {
  if (leftarg < rightarg) {
    
    T pivotvalue = a[leftarg];
    int left = leftarg - 1;
    int right = rightarg + 1;
    
    for(;;) {
      
      while (a[--right] > pivotvalue);
      while (a[++left] < pivotvalue);
      
      if (left >= right) break;
      
      T temp = a[right];
      a[right] = a[left];
      a[left] = temp;
    }
    
    int pivot = right;
    quicksort(a, leftarg, pivot);
    quicksort(a, pivot + 1, rightarg);
  }
}

void P_PT(Point p0, Point p1,double matrix[4][4]);

double rmsd(Structure& S1, Structure& S2, int type,int normal) {
 double RMSD =  rmsd(S1, S2, 1, S1._numRes, 1, S2._numRes, type, normal);
 return RMSD;
}

double rmsd(Structure& S1, Structure& S2, int Start1, int End1, int Start2,
	   int End2, int type, int normal) {
  double** u;
  double** v;
  double* w;
  u = new double*[4];
  v = new double*[4];
  for (int ii = 0; ii < 4; ++ii) {
    u[ii] = new double[4];
    v[ii] = new double[4];
  }
  w = new double[4];
  return rmsd(S1, S2, Start1, End1, Start2, End2, type, normal, u, w, v);
}

double rmsd(Structure& S1, Structure& S2, int Start1, int End1, int Start2,
	   int End2, int type, int normal, double** u, double* w, double** v) {
  Point X_sum,Y_sum;
  Point X_square_sum,Y_square_sum;
  double X_norm,Y_norm;
  double rmsd=0; 
  double Y_X_T_sum[4][4];
  double u1[4][4];

  double a1,a2;
  int i,j,k,numAtom;
  Point X_mean,Y_mean;

  if((End1-Start1)!=(End2-Start2)){
    cout << "rmsd calculation error! Atom numbers not equal " << Start1 << " "
	 << End1 << " " << Start2 << " " << End2 << endl;
    return 100;
  }

  if(abs(End1-Start1)<1)
    return 0;
  j=0;
  X_sum=Point(0,0,0);
  X_square_sum=Point(0,0,0);
  Y_sum=Point(0,0,0);
  Y_square_sum=Point(0,0,0);
  
  for(i=1;i<4;i++)
    for(j=1;j<4;j++)
      Y_X_T_sum[i][j]=0;

  if(type==1){ // Ca RMSD (cRMSD)
    numAtom=End1-Start1+1;
    for(i=Start1;i<=End1;i++)   {
      X_sum=X_sum+S1._res[i]._atom[ATM_CA];
      X_square_sum=X_square_sum+S1._res[i]._atom[ATM_CA].squarev();
    }
    
    for(i=Start2;i<=End2;i++)   {
      Y_sum=Y_sum+S2._res[i]._atom[ATM_CA];
      Y_square_sum=Y_square_sum+S2._res[i]._atom[ATM_CA].squarev();
    }
    for(i=0;i<numAtom;i++)   {
      P_PT(S2._res[i+Start2]._atom[ATM_CA],S1._res[i+Start1]._atom[ATM_CA],u1);
      for(j=1;j<4;j++)
	for(k=1;k<4;k++)
	  Y_X_T_sum[j][k]+=u1[j][k];   
    }

  } else if(type==0){ // all backbone atoms
    numAtom=(End1-Start1+1);
    for(i=Start1;i<=End1;i++)   {
      X_sum=X_sum+S1._res[i]._atom[ATM_N];
      X_sum=X_sum+S1._res[i]._atom[ATM_CA];
      X_sum=X_sum+S1._res[i]._atom[ATM_C];
      X_square_sum=X_square_sum+S1._res[i]._atom[ATM_N].squarev();
      X_square_sum=X_square_sum+S1._res[i]._atom[ATM_CA].squarev();
      X_square_sum=X_square_sum+S1._res[i]._atom[ATM_C].squarev();
    }
    
    for(i=Start2;i<=End2;i++)   {
      Y_sum=Y_sum+S2._res[i]._atom[ATM_N];
      Y_sum=Y_sum+S2._res[i]._atom[ATM_CA];
      Y_sum=Y_sum+S2._res[i]._atom[ATM_C];
      Y_square_sum=Y_square_sum+S2._res[i]._atom[ATM_N].squarev();
      Y_square_sum=Y_square_sum+S2._res[i]._atom[ATM_CA].squarev();
      Y_square_sum=Y_square_sum+S2._res[i]._atom[ATM_C].squarev();
    }
    
    for(i=0;i<numAtom;i++)   {
      P_PT(S2._res[i+Start2]._atom[ATM_N],S1._res[i+Start1]._atom[ATM_N],u1);
      for(j=1;j<4;j++)
	for(k=1;k<4;k++)
	  Y_X_T_sum[j][k]+=u1[j][k];
      P_PT(S2._res[i+Start2]._atom[ATM_CA],S1._res[i+Start1]._atom[ATM_CA],u1);
      for(j=1;j<4;j++)
	for(k=1;k<4;k++)
	  Y_X_T_sum[j][k]+=u1[j][k];
      P_PT(S2._res[i+Start2]._atom[ATM_C],S1._res[i+Start1]._atom[ATM_C],u1);
      for(j=1;j<4;j++)
	for(k=1;k<4;k++)
	  Y_X_T_sum[j][k]+=u1[j][k];
      
    }
      numAtom*=3;   
  }
  X_norm=X_square_sum.sum()-X_sum.square()/(double)(numAtom);
  Y_norm=Y_square_sum.sum()-Y_sum.square()/(double)(numAtom);

  P_PT(Y_sum, X_sum, u1);
  for (i = 1; i < 4; i++)
    for (j = 1; j < 4; j++) {
      u[i][j] = Y_X_T_sum[i][j] - u1[i][j] / (double)(numAtom);
    }
  svdcmp(u,3,3,w,v);
  a1=X_norm+Y_norm;
  a2=2*traceVector(w,3);
	
  rmsd = (X_norm + Y_norm - 2*traceVector(w,3))/(double)(numAtom-1);
  
  if(rmsd<0)
    rmsd=0;
  
  rmsd=sqrt(rmsd);
 int N = End1 - Start1 +1;
 if ( normal == 1)
 {
    rmsd = rmsd/(1+log(sqrt((double)N/6))) ;
           return (double)rmsd;
        }
 else
  return (double)rmsd;
  
}

// calculate CA radius of gyration of a protein structure
double Rg(Structure& Conf)
{	
  int i,NumRes;
  double cx=0,cy=0,cz=0,r=0,rg;
	
  NumRes = Conf._numRes;

  for(i=1;i<=NumRes;++i){
    cx+=Conf._res[i]._atom[1].x;
    cy+=Conf._res[i]._atom[1].y;
    cz+=Conf._res[i]._atom[1].z;
  }
  cx=cx/NumRes;
  cy=cy/NumRes;
  cz=cz/NumRes;

  for(i=1;i<=NumRes;++i){
    r=r+pow(Conf._res[i]._atom[1].x-cx,2)+pow(Conf._res[i]._atom[1].y-cy,2)+pow(Conf._res[i]._atom[1].z-cz,2);
  }

  rg=sqrt(r/NumRes);
  return rg;
}

vector<int> ExpandNumStates(vector<int> NumStates, int FragLength, string type) {
  // Take the inputted NumStates; expand to fragment length if necessary
  vector<int> Rst;
  if (NumStates.size() == 1) {
    Rst.resize(FragLength);
    for (int j = 0; j < FragLength; j++)
      Rst[j] = NumStates[0];
  } else if (NumStates.size() == FragLength) {
    Rst = NumStates;
  } else {
    cout << "The fragment length, " << FragLength << ", does not match the "
	 << "number, " << NumStates.size() << ", of inputted numbers of "
	 << type << " states." << endl;
    exit(0);
  }
  return(Rst);
}

int SampleOne(vector<double> prob) {
  double total = 0.0;
  for (int i = 0; i < prob.size(); i++)
    total += prob[i];
  double chosentotal = drand(0, total);
  total = 0.0;
  for (int i = 0; i < prob.size(); i++) {
    total += prob[i];
    if (total > chosentotal)
    {
      return(i);
      }
  }
  cout << chosentotal << " " << total << "BADBADBAD, prob size: " << prob.size()
       << " probabilities: " << endl;
  for (int i = 0; i < prob.size(); i++)
    cout << prob[i] << " ";
  cout << endl;
  exit(0);
}

double VectorSum(vector<double> vec) {
  double sum = 0.0;
  for (int i = 0; i < vec.size(); i++)
    sum += vec[i];
  return sum;
}

vector<string> FileLines(string inFile) {
  // Returns a vector with strings representing the lines of a file
  vector<string> inFileV;
  ifstream inFileS;
  inFileS.open(inFile.c_str(), ios::in);
  inFileS.clear();
  inFileS.seekg(0);
  while (!inFileS.eof()) {
    string line;
    getline(inFileS, line);
    if (line.size() > 0)
      inFileV.push_back(line);
  }
  return inFileV;
}

int BucketInd(double element, vector<double> bucket) {
  int ind = UNDEF;
  for (int i = 0; i < bucket.size() - 1; i++)
    if (element > bucket[i] & element < bucket[i + 1]) {
      ind = i;
      break;
    }
  return(ind);
}

string File2ProtName(string fname) {
  vector<string> fnameTok;
  split(fname, "/", fnameTok);
  string lastTok = fnameTok[fnameTok.size() - 1];
  int end = 0;
  for (end = 0; end <= 6; end++)
    if (lastTok.substr(end, 1) == ".")
      break;
  return (lastTok.substr(0, end));
}

bool FileExists(string fname) {
  bool flag = false;
  fstream fin;
  fin.open(fname.c_str(), ios::in);
  if (fin.is_open())
    flag = true;
  fin.close();
  return flag;
}

//Calculate rmsd for two structures with out superpose
//type =1, N,CA,C,O Model; type =2, CB added; type =3, all heavy atom
double Root_MSD(Structure& A, Structure& B, int Start1, int End1, int Start2,
	   int End2, int type, int normal) {
       double sq_Sum=0;
       double RMS =0;
       int size = 0;
       
       //cout << Start1 << " " << End1 << " " << Start2 << " " << End2 << endl;
       for(int i = Start1; i <= End1;i++)
       {
          int j = Start2 + i - Start1;       
          MeanSQ(A._res[i], B._res[j], sq_Sum, size, type);
        }
        sq_Sum = sq_Sum/size;
        RMS = sqrt(sq_Sum);
        int N = End1 - Start1 +1;
        if ( normal == 1)
        {
           RMS = RMS/(1+log(sqrt((double)N/6))) ;
           return RMS;
        }
        else
         return RMS;

}

double MSD(Structure& A, Structure& B, int Start1, int End1, int Start2,
	   int End2, int& size,  int type, int normal) {
       double sq_Sum=0;
       for(int i = Start1; i <= End1;i++)
       {
          int j = Start2 + i - Start1;       
          MeanSQ(A._res[i], B._res[j], sq_Sum, size, type);
       }
       return sq_Sum;
}


void MeanSQ(Residue& A, Residue& B, double& sq_Sum, int& size,int type) 
{
          
          double x,y,z;
          //cout << A._type << " and " << B._type << endl;
          //cout << Residue::Name3[A._type] << " at position " << A._posn << " and " << Residue::Name3[B._type] << " at position " << B._posn << endl; 
          if (A._type  != B._type)
          {
             cout << "Not the same residues when calculate RMSD!!!!!!"<<endl;
             cout << Residue::Name3[A._type] << " at position " << A._posn << " and " << Residue::Name3[B._type] << " at position " << B._posn << endl; 
             exit(0);
          }
          for(int k=0; k < 4; k++)
          {
              x = A._atom[k].x - B._atom[k].x;
              y = A._atom[k].y - B._atom[k].y;
              z = A._atom[k].z - B._atom[k].z;
              sq_Sum += x*x+y*y+z*z;
              size ++;
          }
          if (type == 2)
          {
              if (A._type != 5)
              {
                 x = A._atom[5].x - B._atom[5].x;
                 y = A._atom[5].y - B._atom[5].y;
                 z = A._atom[5].z - B._atom[5].z;
                 sq_Sum += x*x+y*y+z*z;
                 size ++;
               }
           }
           if (type == 3)
           {
              if(A._numAtom != B._numAtom)
              {
                cout << " The num of atoms are unequal!!" <<endl;
                exit(0);
              }
              else
              {
                   for(int k=5; k < A._numAtom; k++)
                   {
                      if( A._atom[k]._type < 21 && B._atom[k]._type < 21 && A._atom[k]._type != UNDEF &&  B._atom[k]._type != UNDEF)
                      {   
                           x = A._atom[k].x - B._atom[k].x;
                   	       y = A._atom[k].y - B._atom[k].y;
              		       z= A._atom[k].z - B._atom[k].z;
                           sq_Sum += x*x+y*y+z*z;
                           size ++;
                       }
                       else
                          continue;
                   }
              }
            }


 }





//simple RMS without superimpose in Atom level, flexible 
double Normal_RMS(vector<Atom> A, vector<Atom>B) {
  int A_size,B_size;
  double RMS=0;
  double x,y,z;
  double sq_Sum=0;
  A_size = A.size();
  B_size = B.size();
  if(A_size !=B_size)
    cout<<"RMS size error"<<endl;
  else {
    for(int i=0;i<A_size;++i) {
      x=A[i].x-B[i].x;
      y=A[i].y-B[i].y;
      z=A[i].z-B[i].z;
      sq_Sum += x*x+y*y+z*z;
    }
    sq_Sum = sq_Sum/A_size;
    RMS = sqrt(sq_Sum);
  }
  return RMS;
}

//simple RMS without superimpose in Residue level, flexible 
double Root_MSD(vector<Residue>& A, vector<Residue>& B, int type) {
      
       double rms =0;
       double sq_Sum = 0;
       int size = 0;
       if(A.size() == B.size())
       { 

	 for(int i = 0; i < A.size() ;i++)
          MeanSQ(A[i], B[i], sq_Sum, size, type);
         rms = sqrt(sq_Sum/size);
        // rms = sqrt(sq_Sum);
         return rms;
       }
      else
      {
        cout <<"sequence size unequal!!!"<<endl;
        exit(0);
      }
}


  void box_MullerNsample (double normal_Svalue[],double mean[], int size, double sigma)
  {
     double U, V;
     sigma = sigma *PI/180;
     for (int i=0;i < (size/2)+1;i++)
     {
       U = frand(0,1);
       V = frand(0,1);
       if(i < (size/2))
       {
       normal_Svalue[2*i] = sqrt((-2)*log(U))*cos(2*PI*V)*sigma + mean[2*i];
       normal_Svalue[2*i+1] = sqrt((-2)*log(U))*sin(2*PI*V)*sigma + mean[2*i+1];
       }
       else
       {
         if(size%2 == 1)
	 {	   
           normal_Svalue[2*i] = sqrt((-2)*log(U))*cos(2*PI*V)*sigma + mean[2*i];
	 }
       }
      }
     
   }


  void box_MullerNsample_single(double& normal_Svalue,double mean, double sigma)
  {
     double U, V;
     sigma = sigma *PI/180;
     U = frand(0,1);
     V = frand(0,1);
     normal_Svalue = sqrt((-2)*log(U))*cos(2*PI*V)*sigma + mean;
     if(normal_Svalue > PI)
     {
      normal_Svalue = (normal_Svalue - mean) - mean;
     }
  }
     
  bool sortfun_E (combiEIndex x, combiEIndex y) {return (x.energy<y.energy);}

  bool sortfun_StrE(Structure x, Structure y) {return (x._energy < y._energy);}


//for get the atomname in each line in PDB file
  string PDB_atomname_Extraction(string tmpStr)
  {
  	string atomName;
	if(tmpStr.substr(12,1) != " ")
	  atomName += tmpStr.substr(12,4);
        else if(tmpStr.substr(14,1) == " ")
          atomName += tmpStr.substr(13,1);
	else if(tmpStr.substr(15,1) == " ")			
	  atomName += tmpStr.substr(13,2);                       
	else atomName += tmpStr.substr(13,3);          
	return atomName;
 }

//Read the PDB file generated by myself only store loop region
 void Loop_DecoyPDB_read(string PathparPDB, vector < vector <string> >& TMPCC, int conf_Upbound)
 {
      vector <string> TMPC;
      string tmpStr;
      int confCount = 0;
      ifstream INPUTFILE;
      INPUTFILE.open(PathparPDB.c_str(), ios::in);
      if(!INPUTFILE.is_open()) 
      {
                 cout<<"cannot open pdb file "<<PathparPDB<<endl;
                 exit(0);
      }
      
      while (!INPUTFILE.eof()){
         getline(INPUTFILE,tmpStr);
         if(tmpStr.substr(0,6) == "MODEL ")
	     { 
           if(confCount < conf_Upbound)
		   {
                        confCount++;
		   }
		   else
			break;
        }
	    else if(tmpStr.substr(0,6) == "ATOM  ")
	    {
		   TMPC.push_back(tmpStr);
	    }
	    else if(tmpStr.substr(0,6) == "ENDMDL")
	    {
		    TMPCC.push_back(TMPC);
                    TMPC.clear();

	    }
            else 
                 continue;
      }

           int Size = TMPCC.size();
           INPUTFILE.close();
  }

// Side chain Energy Minimization and clash solver

void SCE_Minimization_list(Structure& conf, int start, int end, vector<int>& ResIdx, vector<int>& ClashNum,  int* List, int List_size, double basic_rot)
{
   int i,j,k,direction,posn,resType, collnum, minIdx, rot_per_unit;
   int minCollnum = 10000;
   int size =  ResIdx.size();
   double pri_x, pri_y, pri_z, R, change_rot = 0;
   Atom rot_axis, tmp;
   Residue tmpRes, tmpRes2;
   double **trans_matrix;
   trans_matrix = matrix(1,3,1,3);

   for(i = 0; i< size; i++)
   {
     minCollnum = 10000;
     if(ClashNum[i] == 0)
     continue;
     posn = ResIdx[i];
     resType = conf._res[posn]._type;
     tmpRes = conf._res[posn];
     if(resType == 12)   // Proline has special treatment
     {
       rot_axis.x = conf._res[posn]._atom[ATM_CA].x-conf._res[posn]._atom[ATM_N].x;
       rot_axis.y = conf._res[posn]._atom[ATM_CA].y-conf._res[posn]._atom[ATM_N].y;
       rot_axis.z = conf._res[posn]._atom[ATM_CA].z-conf._res[posn]._atom[ATM_N].z;
       R= sqrt(rot_axis.x*rot_axis.x+rot_axis.y*rot_axis.y+rot_axis.z*rot_axis.z);
       rot_axis = rot_axis/R;
       for(k = 1 ; k <= 35; k++)
       {
          for(direction = -1; direction <= 1; direction++)
	  {
	   if(direction == 0)
	   continue;
	   rot_per_unit = k*direction;
           change_rot = k*basic_rot*direction;
           RotMatrix(trans_matrix, (double)rot_axis.x, (double)rot_axis.y, (double)rot_axis.z, (double)change_rot);      
           for(j = NUM_BB_ATOM - 1;j < conf._res[posn]._numAtom; ++j)
           {
            if(conf._res[posn]._atom[j]._type==UNDEF || conf._res[posn]._atom[j]._type>=22) continue;
	    else{
		tmp = conf._res[posn]._atom[j]-conf._res[posn]._atom[ATM_CA];
		pri_x = tmp.x;
		pri_y = tmp.y;
		pri_z = tmp.z;
		tmp.x = (double)(trans_matrix[1][1]*pri_x+trans_matrix[1][2]*pri_y+trans_matrix[1][3]*pri_z);
		tmp.y = (double)(trans_matrix[2][1]*pri_x+trans_matrix[2][2]*pri_y+trans_matrix[2][3]*pri_z);
		tmp.z = (double)(trans_matrix[3][1]*pri_x+trans_matrix[3][2]*pri_y+trans_matrix[3][3]*pri_z);
		tmpRes._atom[j] = conf._res[posn]._atom[ATM_CA] + tmp;
	    }
           }
           collnum = Res_clash_detection_list(conf, tmpRes, start, end,  List, List_size);
           if(collnum == 0)
           {
	         minCollnum = 0;
	         tmpRes2 = tmpRes;
	         minIdx = rot_per_unit;
             break;
           }
	   else
	   {
	      if(collnum < minCollnum)
	      {
	        minCollnum = collnum;
	        tmpRes2 = tmpRes;
	        minIdx = rot_per_unit;
	      }
	   }
	 }
	 if(minCollnum == 0)
	 break;
       }
       conf._res[posn] = tmpRes2;
       ClashNum[i] = minCollnum;
     }
     else
     {
       rot_axis.x = conf._res[posn]._atom[ATM_CB].x-conf._res[posn]._atom[ATM_CA].x;
       rot_axis.y = conf._res[posn]._atom[ATM_CB].y-conf._res[posn]._atom[ATM_CA].y;
       rot_axis.z = conf._res[posn]._atom[ATM_CB].z-conf._res[posn]._atom[ATM_CA].z;
       R= sqrt(rot_axis.x*rot_axis.x+rot_axis.y*rot_axis.y+rot_axis.z*rot_axis.z);
       rot_axis = rot_axis/R;
       for(k = 1 ; k <= 35; k++)
       {
          for(direction = -1; direction <= 1; direction++)
	  {
	   if(direction == 0)
	   continue;
	   rot_per_unit = k*direction;
           change_rot = k*basic_rot*direction;
           RotMatrix(trans_matrix, (double)rot_axis.x, (double)rot_axis.y, (double)rot_axis.z, (double)change_rot);      
           for(j = NUM_BB_ATOM;j < conf._res[posn]._numAtom; ++j)
           {
            if(conf._res[posn]._atom[j]._type==UNDEF || conf._res[posn]._atom[j]._type>=22) continue;
	    else{
		tmp = conf._res[posn]._atom[j]-conf._res[posn]._atom[ATM_CB];
		pri_x = tmp.x;
		pri_y = tmp.y;
		pri_z = tmp.z;
		tmp.x = (double)(trans_matrix[1][1]*pri_x+trans_matrix[1][2]*pri_y+trans_matrix[1][3]*pri_z);
		tmp.y = (double)(trans_matrix[2][1]*pri_x+trans_matrix[2][2]*pri_y+trans_matrix[2][3]*pri_z);
		tmp.z = (double)(trans_matrix[3][1]*pri_x+trans_matrix[3][2]*pri_y+trans_matrix[3][3]*pri_z);
		tmpRes._atom[j] = conf._res[posn]._atom[ATM_CB] + tmp;
	    }
           }
           collnum = Res_clash_detection_list(conf, tmpRes, start, end,  List, List_size);
	       if(collnum == 0)
           {
	         minCollnum = 0;
	         tmpRes2 = tmpRes;
	         minIdx = rot_per_unit;
             break;
           }
	       else{
	        if(collnum < minCollnum)
	        {
	            minCollnum = collnum;
	            tmpRes2 = tmpRes;
	            minIdx = rot_per_unit;
	        }
	   }
	 }
	 if(minCollnum == 0)
	 break;
       }
       conf._res[posn] = tmpRes2;
       ClashNum[i] = minCollnum;
     }
   }
  
  free_matrix(trans_matrix,1,3,1,3);
}   
   

double Ellipsoid_Detect (Structure& Conf, int Start, int End, bool type) {
	
 	//determine 2a	
	double ConstL = Residue::bond_length[Conf._res[Start]._type][ATM_C] + 3* (End-Start);  // 3A is the C-C distance
	double ConstLtotal = 0;
	if(type == false)
	 ConstLtotal = ConstL;	
  	else
	{	
	      int MAX_SC_SIZE = 7;
 	      double dist = Conf._res[Start]._atom[ATM_CA].dis(Conf._res[End]._atom[ATM_C]);

	      ConstLtotal = sqrt((sqrt((ConstL/2)*(ConstL/2) - (dist/2)*(dist/2)) +  MAX_SC_SIZE) * (sqrt((ConstL/2)*(ConstL/2) - (dist/2)*(dist/2)) +  MAX_SC_SIZE) + (dist/2)*(dist/2)) * 2;
	}
	return ConstLtotal;
}


void read_Surface(string surfaceFile,vector < vector<int> > &SurfaceList, Structure& Conf, int outInfo) {
  int tmpInt,resType, atomCount;
  char tmpLine[1000];
  string atomName, resName;
  ifstream inFile;
  SSMAP::iterator aaItr;
  SIMAP::iterator aiItr;
  if(outInfo==1)
    cout<<"Reading protein's coordinates..."<<endl;
  inFile.open(surfaceFile.c_str(), ios::in);
  if (!inFile.is_open()) {
    cout << "cannot open surface file " << surfaceFile << endl;
    exit(0);
  }
   // Read the PDB file line by line
  while (!inFile.eof()) {
    inFile.getline(tmpLine, 1000);
    string tmpStr = tmpLine;
    if (tmpStr.substr(0,6) == "ATOM  ") {
      // This is an atom line, so we execute the main routine to add this atom
      // to the structure
      
      // Check if this atom line is malformed
      for (int i = 0; i < tmpStr.size(); i++)
	if (!isgraph(tmpStr[i]) & !isspace(tmpStr[i])) {
	  cout << surfaceFile << " has a malformed atom line in residue "
	       << " and it is " << tmpStr[i] << " at character " << i << endl;
	  return;
        }
      if(tmpStr.substr(13, 1) == "D")
	  continue;
      
      // ignore alternative atom positions
      if (tmpStr.substr(16, 1) != "A" && tmpStr.substr(16, 1) != " ")
	  continue;

      // ignore alternative residues
      if (tmpStr.substr(26, 1) != " ")
	  continue;
      
      // Get the residue sequence number
      tmpInt = atoi(tmpStr.substr(22, 4).c_str());
      resName = tmpStr.substr(17,3);
	if ((aaItr = Residue::AAMap.find(resName)) == Residue::AAMap.end()) {
	  // if the residue is unknown$
	  continue;
	}
      
      resType=Residue::AIMap.find(aaItr->second)->second; // value of aaItr was obtained earlier
      atomName=Residue::Name1[resType];
      
      if(tmpStr.substr(12,1) != " ")
	atomName += tmpStr.substr(12,4);
      else if(tmpStr.substr(14,1)==" ")
	atomName+=tmpStr.substr(13,1);
      else if(tmpStr.substr(15,1)==" ")
	atomName+=tmpStr.substr(13,2);
      else atomName+=tmpStr.substr(13,3);
      if((aiItr=Residue::AtomMap.find(atomName))!=Residue::AtomMap.end())
	// backbone atoms
	atomCount = aiItr->second;
      for(int i=1; i <= Conf._numRes; i++) 
      {
      	   if( tmpInt == Conf._res[i]._pdbIndex && resType == Conf._res[i]._type)
	   {
		int surf_res_index = i-1;
		SurfaceList[surf_res_index].push_back(atomCount);
	   }
      }     
   }
  }
}

double distance_to_line(Point begin, Point end, Point x ){
   //translate the begin to the origin
      Point axis = end - begin;
      Point D1 = x - begin;
      double area = (D1*axis).pabs();
      return area / axis.pabs();
}


double entropy(double* Prob, int size){
      double partitionF = 0;
      double entropy = 0;
      for(int i=0; i < size; i++)
      partitionF += Prob[i];
      
      for(int i=0; i < size; i++)
      {
        Prob[i] = Prob[i]/partitionF;
        entropy += Prob[i]*(-1)*log(Prob[i]);
      }        
   return entropy;
}
