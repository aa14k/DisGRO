// structure.cpp

#include "structure.h"
#include "util.h"
#include "reprst.h"
#include "rotamer.h"
#include "cal_energy.h"
#include "loop_closure_c/sturm.h"
#include "loop_closure_c/tripep_closure.h"
#include "sample_states.h"
#include <math.h>

#define CUB_SIZE 5.5

double Structure::_T = 1;
string Structure::_sequence = "";

Structure::Structure(int n) {
  init(n);
}
Structure::Structure(const Structure& S) {
  init(1);
  *this = S;
}
Structure* Structure::operator=(const Structure& S) {
  copyStructure(S, 1, S._numRes);
  return this;
}

void Structure::Destruct() {
//  delete [] _res->_parent;
  delete [] _res;
  delete [] _scProb;
}
Structure::~Structure() {
//  cout << "DESTRUCTi Structure!!!!!!"<<endl;
  Destruct();
}

void Structure::init(int n) {
  _numRes    = n; 
  _weight    = 0;
  _weight_2  = 0;
  _energy    = 0;
  _chainName = "";
  G_rmsd     = 100000;
  allatm_rmsd = 100000;
  _checkEn	 = -1;
  for(int i=0; i<6;i++)
  {
   if(i%2 == 0)
   _Fringevtx[i] = -30000;
   else
   _Fringevtx[i] = 30000;
  }
  _res       = new Residue[n+1];
  _scProb    = new  double[n+1];
  Atom p(0, 0, 0, UNDEF);
  for (int i = 0; i <= n; i++) {
    _res[i]._parent = this;
    _res[i]._bbc  = p;
    _res[i]._scc  = p;
    _res[i]._posn = i;
  }
  _Confcenter = p;
  // initialize _res[0] and _res[1], _res[0] is just a place holder,
  // the structure starts from _res[1]
  _res[0]._atom[ATM_N].CopyPos(Atom(-1,     1,    1));
  _res[0]._atom[ATM_CA].CopyPos(Atom( 0,     1,    1));
  _res[0]._atom[ATM_C].CopyPos(Atom( 0,     1,    0));
  _res[1]._atom[ATM_N].CopyPos(Atom( 0,     0,    0));
  _res[1]._atom[ATM_CA].CopyPos(Atom( 1.458, 0,    0));
  _res[1]._atom[ATM_H].CopyPos(Atom( 0,    -1.01, 0));
  _ATM_OXT._type = 20;
  for (int i = 0; i < ENERGY_TYPES; i++)
    _enArr[i] = 0;
}

void Structure::copyStructure(const Structure& S, int start, int end) {
  // Copy a source structure C to the target structure, from start to end
 // int n = end -start +1;
  if (_numRes != S._numRes) {
    // The function is being used to assign a source structure with a different
    // number of residues, indicating the intent to clear out all members of
    // the target structure.
    if (start != 1 | end != S._numRes) {
     cout << "Cannot call copyStructure on fragments when _numRes incompatible" << endl;
      exit(0);
    }
    Destruct();
    init(S._numRes);
  }
  // Now copy over all the members from the source structure
  _numRes   = S._numRes;
  _weight   = S._weight;
  _weight_2 = S._weight_2;
  _energy   = S._energy;
  for (int i = start; i <= end; i++) {
    _res[i]    = S._res[i];
    _scProb[i] = S._scProb[i];
  }
  _useNatSC    = S._useNatSC;
  _toBeSampled = S._toBeSampled;
  priMap       = S.priMap;
  _chainName   = S._chainName;
  _ProtName    = S._ProtName;
  G_rmsd       = S.G_rmsd;
  allatm_rmsd  = S.allatm_rmsd;
  LoopTorsion  = S.LoopTorsion;
  MEdis1       = S.MEdis1;
  MEdis2       = S.MEdis2;
  _numChain    = S._numChain;
  missSeq      = S.missSeq;
  missSeqPos   = S.missSeqPos;
  missSeqfrom1 = S.missSeqfrom1;
  for (int i = 0; i < 20; i++)
    _firstResidues[i] = S._firstResidues[i];
  for (int i = 0; i < ENERGY_TYPES; i++)
    _enArr[i] = S._enArr[i];
  // copy _HB_Profile
  _HB_Profile = S._HB_Profile;
  // Update adj_list
  if (start > 1 || end < _numRes)
    calCenter(start, end);
  Closed     = S.Closed;
  Success    = S.Success;
  _Confcenter = S._Confcenter;
  for(int i=0; i< 6; i++)
  _Fringevtx[i] = S._Fringevtx[i];
}


bool Structure::readPdb(vector <string>& stringlist,SSET& SelRes)
{
  int tmpInt,resType;
  string tmpStr;
  string atomName;
  ifstream inFile;
  SSMAP::iterator aaItr;
  SIMAP::iterator aiItr;

  // assign all atoms as missing at the beginning
 for (int i = 0; i < _numRes; ++i)
   for (int j = 0; j < _res[i]._numAtom; ++j)
     _res[i]._atom[j]._state = 1; 


  int resNum = 0;
  int prev_num=-1000,atomCount,numChain=0;
  
  string resName,chainName,tmpResName,prev_chain="**";
  
  int skipped = 0;
  for(int i=0;i<stringlist.size();i++)
    {
      string tmpStr = stringlist[i];
      if (tmpStr.size() == 0) {
	skipped++;
	if (skipped > 1000) {
	  cout << "The input structure has a bad chunk in it somewhere, "
	       << "which has caused over 1,000 blank lines in the file." << endl;
	  return false;
	}
	continue;
      }
      
      if (tmpStr.substr(0,6) == "ATOM  ") {
	// This is an atom line, so we execute the main routine to add this atom
	// to the structure
	
	// Check if this atom line is malformed
	for (int i = 0; i < tmpStr.size(); i++)
	  if (!isgraph(tmpStr[i]) & !isspace(tmpStr[i])) {
	    return false;
	  }
	if (tmpStr.size() < 54 || tmpStr.size() > 90) {
	  cout << tmpStr << " has a malformed atom line in residue " << resNum
	       << " that, with " << tmpStr.size() << " characters, falls outside"
	       << " the line width region" << endl;
	  return false;
	}
	// ignore Hydrogen atoms
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
	
	if(SelRes.size()!=0) { // only read those residues in SelRes if it is not empty
	  tmpResName=itoa(tmpInt);
	  chainName=tmpStr.substr(21,1);
	  tmpResName =chainName + tmpResName;
	  if(SelRes.find(tmpResName)==SelRes.end()) continue; // read only selected residues, others are omitted.
	}
	chainName=tmpStr.substr(21,1);
	
	if(prev_num==-1000 || prev_num != tmpInt || chainName != prev_chain) {
	  // check for the next residue
	  resName = tmpStr.substr(17,3);
	  if ((aaItr = Residue::AAMap.find(resName)) == Residue::AAMap.end()) {
	    // if the residue is unknown$
	    continue;
	  }
	  if (prev_num != -1000 && (tmpInt - prev_num) != 1) {
	    // broken chain
	    cout << "readPdb found missing residue between " << prev_num
		 << " and " << tmpInt << endl;
//		 return false;
	  }
	  prev_num=tmpInt;
	  resNum++;
	  _res[resNum]._numAtom = 0;
	  _res[resNum]._posn = resNum;
	  // this is set to 1 because H atom occupies the fifth position
	  // of the atom array.  in most X-ray structures, there is no H
	  // atoms. But sometimes, a pdb file does have H atoms In those
	  // cases, when H atom is read, _numAtom does not increase.
	  
	  _res[resNum]._pdbIndex = tmpInt;
	  if(chainName != prev_chain){
	    _chainName += chainName;
	    _firstResidues[numChain++] = resNum;
	  }
	  prev_chain = chainName;
	  resType = Residue::AIMap.find(aaItr->second)->second; // value of aaItr was obtained earlier
	  _res[resNum]._type = resType;
	}
	atomName = Residue::Name1[resType];
	
	if(tmpStr.substr(12,1) != " ")
	  atomName += tmpStr.substr(12,4);
	else if(tmpStr.substr(14,1) == " ")
	  atomName += tmpStr.substr(13,1);
	else if(tmpStr.substr(15,1) == " ")
	  atomName += tmpStr.substr(13,2);
	else atomName += tmpStr.substr(13,3);
	if((aiItr = Residue::AtomMap.find(atomName)) != Residue::AtomMap.end()) {
	  // backbone atoms
	  atomCount = aiItr->second;
	  _res[resNum]._atom[atomCount]._type = Residue::vdwType[resType][atomCount];
	  _res[resNum]._atom[atomCount]._posn = atomCount;
	  if (tmpStr.size() >= 67)
	    _res[resNum]._atom[atomCount]._Bfactor =
	      atof(tmpStr.substr(60, 6).c_str());
	  // assign the atom name back to their name in PDB files like CA, CB...
	  
	  _res[resNum]._atom[atomCount]._name = atomName.substr(1);
	  _res[resNum]._atom[atomCount].x = atof(tmpStr.substr(30,8).c_str());
	  _res[resNum]._atom[atomCount].y = atof(tmpStr.substr(38,8).c_str());
	  _res[resNum]._atom[atomCount].z = atof(tmpStr.substr(46,8).c_str());
	  _res[resNum]._atom[atomCount]._state = 0;
	  _res[resNum]._numAtom++;
	}
	else if (atomName.substr(1,3) == "OXT")
	  continue;
	else {
	  continue;
	}
      _res[resNum]._atom[4]._type = UNDEF;
      }
      
    }
  // Collection of sanity checks for this structure
  if (resNum >= MAX_NUM_RES) {
    cout << " has " << resNum << " residues, but the most allowed "
  	 << "is " << MAX_NUM_RES << endl;
    return false;
  }
  if (resNum == 0) {
    cout << " has no atoms" << endl;
    return false;
  }
  _numRes = resNum;
  calCenter(1,_numRes);
  calSC(1,_numRes);
  assign_ss(*this, 1, this->_numRes,0);
  for (int i = 1; i <= _numRes; i++){
    // check missing atoms and if any initialize them by setting random coordinates
    for(int j = 0; j < Residue::numAtom[_res[i]._type]; j++){
      if(_res[i]._atom[j]._type == UNDEF){
	_res[i]._atom[j]._type = Residue::vdwType[_res[i]._type][j];
	_res[i]._atom[j].randomizeCo();
      }
    }
    _res[i]._numAtom = Residue::numAtom[_res[i]._type];
  } 
  // calculate center of functional group
  for(int i = 1; i <= _numRes; i++)
    _res[i].cal_fg();
  return true;
}


bool Structure::readPdb(string protFile, SSET& SelRes, int outInfo) {
  int tmpInt,resType;
  char tmpLine[1000];
  string atomName;
  ifstream inFile;
  SSMAP::iterator aaItr;
  SIMAP::iterator aiItr;


  // assign all atoms as missing at the beginning
 for (int i = 0; i < _numRes; ++i)
   for (int j = 0; j < _res[i]._numAtom; ++j)
     _res[i]._atom[j]._state = 1; 

  // residue number count starts from 1
  int resNum = 0;
 
  // Read the protein's coordinates 
  if(outInfo==1)
    cout<<"Reading protein's coordinates..."<<endl;
  
  int prev_num=-1000,atomCount,numChain=0;
  string resName,chainName,tmpResName,prev_chain="**";
  
  // Open the PDB file
  inFile.open(protFile.c_str(), ios::in);
  if (!inFile.is_open()) {
    cout << "cannot open pdb file " << protFile << endl;
    exit(0);
  }
  
  // Read the PDB file line by line
  int skipped = 0;
  while (!inFile.eof()) {
    inFile.getline(tmpLine, 1000);
    string tmpStr = tmpLine;
    
    // Is this line empty? If so, skip it, but stop after 1,000 skips (error)
    if (tmpStr.size() == 0) {
      skipped++;
      if (skipped > 1000) {
	cout << "The input structure has a bad chunk in it somewhere, "
	     << "which has caused over 1,000 blank lines in the file." << endl;
	return false;
      }
      continue;
    }
    
    
    if (tmpStr.substr(0,6) == "ATOM  ") {
      // This is an atom line, so we execute the main routine to add this atom
      // to the structure
      // Check if this atom line is malformed
      for (int i = 0; i < tmpStr.size(); i++)
	if (!isgraph(tmpStr[i]) & !isspace(tmpStr[i])) {
	  cout << protFile << " has a malformed atom line in residue " << resNum
	       << " and it is " << tmpStr[i] << " at character " << i << endl;
	  return false;
	}
      if (tmpStr.size() < 54 || tmpStr.size() > 90) {
	cout << protFile << " has a malformed atom line in residue " << resNum
 	     << " that, with " << tmpStr.size() << " characters, falls outside"
	     << " the line width region" << endl;
 	return false;
      }
      // ignore Hydrogen atoms
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
      
      if(SelRes.size()!=0) {
	tmpResName=itoa(tmpInt);
	chainName=tmpStr.substr(21,1);
	tmpResName =chainName + tmpResName;
	if(SelRes.find(tmpResName)==SelRes.end()) continue; // read only selected residues, others are omitted.
      }
      chainName=tmpStr.substr(21,1);
      
      if(prev_num==-1000||prev_num!=tmpInt||chainName!=prev_chain) {
	// check for the next residue
	resName = tmpStr.substr(17,3);
	if ((aaItr = Residue::AAMap.find(resName)) == Residue::AAMap.end()) {
	  // if the residue is unknown$
	  continue;
	}
	if(prev_num != -1000 && (tmpInt - prev_num) != 1) { // broken chain
	  cout<<"readPdb found missing residue between "<<prev_num<<" and "<<tmpInt<<endl;
	  if(!missSeq.empty() && !missSeqPos.empty())
	  {		
		for(int t = 0; t < missSeqPos.size();t++)
		{	     
		    if(tmpInt - missSeqPos[t] == 3 || tmpInt -missSeqPos[t] ==5)		   
		    {
			missSeqfrom1.push_back(resNum);
			SIMAP::iterator siItr;
			for(int i = 1; i<= missSeq[t].length() ; ++i){
   			   if((siItr=Residue::AIMap.find(missSeq[t].substr(i-1,1))) !=	Residue::AIMap.end())
			   {
				resNum++;     				
				_res[resNum]._type = siItr->second;
				_res[resNum]._numAtom = 0;
				_res[resNum]._posn = resNum;
				_res[resNum]._pdbIndex=1000+resNum;
			    }
				
   			   else {
      			      cout<<"Residue type error! "<<_sequence.substr(i-1,1)<<endl;
      				exit(0);
    			    }
			}
			missSeqfrom1.push_back(resNum+1);
		   }
		}
	   }
	}
	prev_num=tmpInt;
	resNum++;
	
//	_res[resNum]._numAtom = 0;
	_res[resNum]._posn = resNum;
	// this is set to 1 because H atom occupies the fifth position
	// of the atom array.  in most X-ray structures, there is no H
	// atoms. But sometimes, a pdb file does have H atoms In those
	// cases, when H atom is read, _numAtom does not increase.
	
	_res[resNum]._pdbIndex=tmpInt;
	if(chainName!=prev_chain){
	  _chainName+=chainName;
	  _firstResidues[numChain++]=resNum;
	}
	prev_chain=chainName;
	resType=Residue::AIMap.find(aaItr->second)->second; // value of aaItr was obtained earlier
	_res[resNum]._type=resType;
        _res[resNum]._numAtom = Residue::numAtom[_res[resNum]._type];
      }
      atomName=Residue::Name1[resType];
      
      if(tmpStr.substr(12,1) != " ")
	atomName += tmpStr.substr(12,4);
      else if(tmpStr.substr(14,1)==" ")
	atomName+=tmpStr.substr(13,1);
      else if(tmpStr.substr(15,1)==" ")
	atomName+=tmpStr.substr(13,2);
      else atomName+=tmpStr.substr(13,3);
      if((aiItr=Residue::AtomMap.find(atomName))!=Residue::AtomMap.end()) {
	// backbone atoms
	atomCount = aiItr->second;
	_res[resNum]._atom[atomCount]._type=Residue::vdwType[resType][atomCount];
	_res[resNum]._atom[atomCount]._posn=atomCount;
	if (tmpStr.size() >= 67)
	  _res[resNum]._atom[atomCount]._Bfactor =
	    atof(tmpStr.substr(60, 6).c_str());
	// assign the atom name back to their name in PDB files like CA, CB...
	_res[resNum]._atom[atomCount]._name = atomName.substr(1);
	_res[resNum]._atom[atomCount].x = atof(tmpStr.substr(30,8).c_str());
	_res[resNum]._atom[atomCount].y = atof(tmpStr.substr(38,8).c_str());
	_res[resNum]._atom[atomCount].z = atof(tmpStr.substr(46,8).c_str());
	
	_res[resNum]._atom[atomCount]._state = 0;
     }
      else if (atomName.substr(1,3) == "OXT")
	continue;
      else {
	  continue;
      }
    }
    
  }
  // Collection of sanity checks for this structure
  if (resNum >= MAX_NUM_RES) {
    cout << protFile << " has " << resNum << " residues, but the most allowed "
	 << "is " << MAX_NUM_RES << endl;
    return false;
  }
  if (resNum == 0) {
    cout << protFile << " has no atoms" << endl;
    return false;
  }
  
  if (outInfo == 1)
    cout << "Number of residues is " << resNum << "." << endl;
  _numRes=resNum;
  
  inFile.close();
  
  calCenter(1,_numRes);
  calSC(1,_numRes);
  assign_ss(*this, 1, this->_numRes,0);
  
  // the number of atoms for each residue is reset here using the
  // fixed values in Residue::numAtom[] When adding H atoms, addition
  // of the backbone H atom should not increase _numAtom, but
  // sidechain H atoms should 
  
  for (int i = 1; i <= _numRes; i++){
    _res[i]._numAtom = Residue::numAtom[_res[i]._type];
    for(int j = 0; j < _res[i]._numAtom; j++)
    {
       _res[i]._atom[j]._type=Residue::vdwType[_res[i]._type][j];
       _res[i]._atom[j]._posn=j;
    }
  } 
  // calculate center of functional group
  for(int i = 1; i <= _numRes; i++)
    _res[i].cal_fg();
  
  // Finally, store the sequence string and protein name
  StoreSequence();
  _ProtName = File2ProtName(protFile);
  
  // If we're here, everything went correctly
  return true;
}

// read a sequence from a file
// initialize all the residues and atoms at the same time, 
// atom coordinates are all set to zero except res[0] and res[1], which are initialized when Structure is allocated
void Structure::readSeq(char* protFile) {
  char tmpLine[1000];
  ifstream inFile;
  inFile.open(protFile,ios::in);
  if(!inFile.is_open()) {
    cout<<"cannot open pdb file "<<protFile<<endl;
    exit(0);
  }

  inFile.getline(tmpLine,1000);
  _sequence = tmpLine;

  SIMAP::iterator siItr;
  _numRes = _sequence.length();
  cout<<"Number of residues "<<_numRes<<endl;
  // set the number of atoms for each residue
  for(int i = 1; i<=_numRes; ++i){
    if((siItr=Residue::AIMap.find(_sequence.substr(i-1,1))) !=
       Residue::AIMap.end())
      _res[i]._type = siItr->second;
    else {
      cout<<"Residue type error! "<<_sequence.substr(i-1,1)<<endl;
      exit(0);
    }
    _res[i]._posn = i;
    _res[i]._numAtom = Residue::numAtom[_res[i]._type];
    for(int j = 0; j<_res[i]._numAtom; ++j){
      _res[i]._atom[j]._type = Residue::vdwType[_res[i]._type][j];
      _res[i]._atom[j]._name = Residue::cType[_res[i]._type][j];
    }
  }
}

// Helper function for secondary structure prediction
int Structure::writeData(char *data, size_t size, size_t nmemb, string *s) {
  if (s == NULL)
    return 0;
  s->append(data, size*nmemb);
  return size * nmemb;
}

void Structure::StoreSequence() {
  _sequence = "";
  for (int i = 1; i <= _numRes; i++)
    _sequence = _sequence + Residue::Name1[_res[i]._type];
}

void Structure::WriteSequence(string file) {
  ofstream fout;
  fout.open(file.c_str());
  fout << _sequence << endl;
  fout.close();
}

// fetch the predicted secondary structure of the sequence
void Structure::PredictSS() {
  // First check if this structure is in the SS cache
  string SableName   = PF::SSDir + "/" + _ProtName + ".sable.txt";
  string PsipredName = PF::SSDir + "/" + _ProtName + ".psipred.txt";
  if (!FileExists(SableName) | !FileExists(PsipredName)) {
    cout << "Couldn't find secondary structure file. "
	 << "Need to unset SS and CON energies or produce SS file." << endl;
    exit(0);
  }
  vector<string> rst1 = FileLines(SableName);
  vector<string> rst2 = FileLines(PsipredName);
  vector<string> ss, ssph, sspe, sspc, sa;
  split(rst2[0], "\t", ss);
  split(rst2[1], "\t", ssph);
  split(rst2[2], "\t", sspe);
  split(rst2[3], "\t", sspc);
  split(rst1[4], "\t", sa);
  _ssPred     = vector<int>(_numRes + 1);
  _ssPredProb = vector< vector<double> >(_numRes + 1);
  _saPred     = vector<int>(_numRes + 1);
  for (int i = 0; i < _numRes; i++) {
    _ssPred[i+1]           = Residue::SIMap.find(ss[i])->second;
    _ssPredProb[i+1]       = vector<double>(3);
    _ssPredProb[i+1][SS_H] = atof(ssph[i].c_str());
    _ssPredProb[i+1][SS_E] = atof(sspe[i].c_str());
    _ssPredProb[i+1][SS_C] = atof(sspc[i].c_str());
    _saPred[i+1]           = atoi(sa[i].c_str());
  }
}

// find disulfide bonds in a structure, the residue is labelled by
// changing _scState to -2 (or other negative value)
void Structure::calSS() {
  int i,j,tmpInt;
  string tmpStr;

  for(i=1;i<=_numRes;++i){
    if(_res[i]._type==(Residue::AIMap.find("C")->second)) {
      for(j=i+3;j<=_numRes;++j) {
        if(_res[j]._type==(Residue::AIMap.find("C")->second)) {
          tmpStr="CSG";
          tmpInt=Residue::AtomMap.find(tmpStr)->second;
          if(_res[i]._atom[tmpInt].dis(_res[j]._atom[tmpInt])<2.15){
            _res[i]._scState=-2;
            _res[j]._scState=-2;
            break;
          }
        }
      }
    }
  }
}

void Structure::addH()
{
  int i,tmpInt;
  string tmpStr;
  SIMAP::iterator aiItr;
  Point tmpArr[3];

  // checked by eye examination of helix.pdb with 14 residues
  for(i=1;i<=_numRes;++i){
    if(_res[i]._type>=20) continue;  // only add H to amino acids
    tmpStr=Residue::Name1[_res[i]._type];
    // add backbone H atom
    if(_res[i]._atom[4]._type == UNDEF){
      if(i==1){// there may be three H atoms on the terminal N whose index is i, Here only one hydrogens are added,
	tmpStr+="H";
	if((aiItr=Residue::AtomMap.find(tmpStr))!=Residue::AtomMap.end())
	  tmpInt=aiItr->second;
	else {
	  cout<<"error at first residue!"<<endl;
	}
	tmpArr[0]=_res[i]._atom[2];
	tmpArr[1]=_res[i]._atom[1];
	tmpArr[2]=_res[i]._atom[0];
	// add one H among three possible ones, one H may have h-bond with backbone oxygen, one H may be close to the side-chain preventing it forming h-bond with other atoms, the third one is added here. The -75 is got from sheet.pdb. 
	calCo(tmpArr,1.01,PI*(-75)/180,0,_res[i]._atom[tmpInt]);
	_res[i]._atom[tmpInt]._type=25;
	// add the second H after the last atom
      }
      else if(tmpStr=="P"||_res[i]._type==(Residue::AIMap.find("P")->second)){  // put proline H 
	tmpStr+="H"; 
	if((aiItr=Residue::AtomMap.find(tmpStr))!=Residue::AtomMap.end()){
	  tmpInt=aiItr->second;
	  _res[i]._atom[tmpInt]=_res[i]._atom[0];
	  _res[i]._atom[tmpInt]._type=25;
	}
	else {
	  _res[i]._atom[4]=_res[i]._atom[0];
	  _res[i]._atom[4]._type=25;
	}
   }
   else {
	tmpStr+="H";
	if((aiItr=Residue::AtomMap.find(tmpStr))!=Residue::AtomMap.end())
	  tmpInt=aiItr->second;
	else {
	  //cout<<"Pro found at position "<<i<<endl;
	}
	if(_res[i]._atom[tmpInt]._type==UNDEF){
	  tmpArr[0]=_res[i-1]._atom[1];
	  tmpArr[1]=_res[i-1]._atom[2];
	  tmpArr[2]=_res[i]._atom[0];
	  calCo(tmpArr,1.01,PI*123/180,0,_res[i]._atom[tmpInt]);
	  _res[i]._atom[tmpInt]._type=25;
	}
	else {
	  cout<<"defined H "<<_res[i]._atom[tmpInt]._type<<endl;
	}
      }
    }
    // add the rest of H atoms
    // comment this out because there is a bug here
    for(int j=5;j<Residue::numAtom[_res[i]._type];++j){
      // bug is in the next line, 
      // _res[i]._atom[j]._type can be UNDEF, which makes Atom::hbondH[_res[i]._atom[j]._type] out of array index range
      if(Atom::hbondH[Residue::vdwType[_res[i]._type][j]]==1 && _res[i]._atom[j]._type==UNDEF) { // add hbondH
	tmpArr[0]=_res[i]._atom[Residue::prev_atom[_res[i]._type][j][0]];
	tmpArr[1]=_res[i]._atom[Residue::prev_atom[_res[i]._type][j][1]];
	tmpArr[2]=_res[i]._atom[Residue::prev_atom[_res[i]._type][j][2]];
	calCo(tmpArr,1.01,PI*Residue::bond_angle[_res[i]._type][j]/180,Residue::torsion[_res[i]._type][j],_res[i]._atom[j]);
	_res[i]._atom[j]._type=22;
      }
		if(isnan(_res[i]._atom[j].x))
		{
			_res[i]._atom[j].x =0;
			_res[i]._atom[j].y =0;
			_res[i]._atom[j].z =0;
		}
    }
   	
  } 
}

void Structure::calCenter(int Start, int End, bool Sidechain) {
    // for proteins

    // calculate backbone (N CA C O H Cb) and side-chain center
    for (int i = Start; i <= End; i++) {
      _res[i]._bbc    = Point(0, 0, 0);
      _res[i]._scc    = Point(0, 0, 0);
      _res[i]._center = Point(0, 0, 0);
      int numBBAtom = 0;
      int numSCAtom = 0;

    for (int j = 0; j < NUM_BB_ATOM; j++) {
	if (_res[i]._atom[j]._type > 20)
	  continue;
	if(_res[i]._atom[j].x== 0 && _res[i]._atom[j].y== 0 && _res[i]._atom[j].z== 0) continue;
	if (_res[i]._atom[j]._type != UNDEF) {
	  numBBAtom++;
	  _res[i]._bbc    = _res[i]._bbc    + _res[i]._atom[j];
	  _res[i]._center = _res[i]._center + _res[i]._atom[j];
	}
      }
      // because we count H in, therefore j <= _res[i]._numAtom
    if(Sidechain)
    for (int j = NUM_BB_ATOM; j <= _res[i]._numAtom; j++) {
	if (_res[i]._atom[j]._type > 20)
	  continue;
	if(_res[i]._atom[j].x== 0 && _res[i]._atom[j].y== 0 && _res[i]._atom[j].z== 0) continue;
	if (_res[i]._atom[j]._type != UNDEF) {
	  numSCAtom++;
	  _res[i]._scc    = _res[i]._scc    + _res[i]._atom[j];
	  _res[i]._center = _res[i]._center + _res[i]._atom[j];
	}
      }

      if (numBBAtom != 0)
	_res[i]._bbc = _res[i]._bbc / numBBAtom;
      else
	_res[i]._bbc._type = UNDEF;

      if(Sidechain)
      {
       if (numSCAtom != 0)
	_res[i]._scc = _res[i]._scc / numSCAtom;
       else if (_res[i]._type == 0)
	// ALA
	_res[i]._scc = _res[i]._atom[5];
       else if (_res[i]._type == 5)
	// GLY
	_res[i]._scc = _res[i]._atom[1];
       else
	_res[i]._scc._type = UNDEF;
      }

      if ((numBBAtom + numSCAtom) != 0) {
        if(Sidechain)
	_res[i]._center = _res[i]._center / (numBBAtom+numSCAtom);
	else
        _res[i]._center = _res[i]._bbc;
	_res[i]._center._type = 0;
      }
      else
	_res[i]._center._type=UNDEF;
    }
}

void Structure::calCenter(int Start, int End) {
  if (_res[1]._type >= 20 && _res[1]._type <= 24) {
    // for RNA or DNA
    int i,j,numBaseAtom;
    int numPAtom = 3;
    int numSugarAtom = 9;
    string tMpatom_name;
    for (i = 1; i <= _numRes; i++) {
      for (int k = 0; k < 3; k++)
	_res[i]._rna_center[k]=Point(0,0,0);
      numBaseAtom = 0;
      for(j=0;j<_res[i]._numAtom;++j) {
	if(_res[i]._atom[j]._type == UNDEF) continue;
	tMpatom_name = _res[i]._atom[j]._name.substr(0,3);
	if(tMpatom_name == "OP3")                                      //op3 is ignored
	  {	continue;}
	else if(tMpatom_name == "P" || tMpatom_name == "OP1" ||tMpatom_name == "OP2")  //p group
	  {
	    _res[i]._rna_center[0] = _res[i]._rna_center[0]+_res[i]._atom[j];
	  }
	else if(tMpatom_name.substr(2) == "'") {
	  //sugar ring
	  _res[i]._rna_center[1] = _res[i]._rna_center[1]+_res[i]._atom[j];
	} else {
	  //base group
	  _res[i]._rna_center[2] = _res[i]._rna_center[2]+ _res[i]._atom[j];
	  numBaseAtom++;
	}
      }
      _res[i]._rna_center[0] = _res[i]._rna_center[0]/numPAtom;
      _res[i]._rna_center[1] = _res[i]._rna_center[1]/numSugarAtom;
      _res[i]._rna_center[2] = _res[i]._rna_center[2]/numBaseAtom;
       }	

  } else {
    // for proteins
    // calculate backbone (N CA C O H Cb) and side-chain center
    for (int i = Start; i <= End; i++) {
      _res[i]._bbc    = Point(0, 0, 0);
      _res[i]._scc    = Point(0, 0, 0);
      _res[i]._center = Point(0, 0, 0);
      int numBBAtom = 0;
      int numSCAtom = 0;

    for (int j = 0; j < NUM_BB_ATOM; j++) {
	if (_res[i]._atom[j]._type > 20)
	  continue;
	if(_res[i]._atom[j].x== 0 && _res[i]._atom[j].y== 0 && _res[i]._atom[j].z== 0) continue;
	if (_res[i]._atom[j]._type != UNDEF) {
	  numBBAtom++;
	  _res[i]._bbc    = _res[i]._bbc    + _res[i]._atom[j];
	  _res[i]._center = _res[i]._center + _res[i]._atom[j];
	}
      }
    // because we count H in, therefore j <= _res[i]._numAtom
    for (int j = NUM_BB_ATOM; j <= _res[i]._numAtom; j++) {
	if (_res[i]._atom[j]._type > 20)
	  continue;
	if(_res[i]._atom[j].x== 0 && _res[i]._atom[j].y== 0 && _res[i]._atom[j].z== 0) continue;
	if (_res[i]._atom[j]._type != UNDEF) {
	  numSCAtom++;
	  _res[i]._scc    = _res[i]._scc    + _res[i]._atom[j];
	  _res[i]._center = _res[i]._center + _res[i]._atom[j];
	}
      }

      if (numBBAtom != 0)
	_res[i]._bbc = _res[i]._bbc / numBBAtom;
      else
	_res[i]._bbc._type = UNDEF;
      if (numSCAtom != 0)
	_res[i]._scc = _res[i]._scc / numSCAtom;
      else if (_res[i]._type == 0)
	// ALA
	_res[i]._scc = _res[i]._atom[5];
      else if (_res[i]._type == 5)
	// GLY
	_res[i]._scc = _res[i]._atom[1];
      else
	_res[i]._scc._type = UNDEF;

      if ((numBBAtom + numSCAtom) != 0) {
	_res[i]._center = _res[i]._center / (numBBAtom+numSCAtom);
	_res[i]._center._type = 0;
      }
      else
	_res[i]._center._type=UNDEF;
    }
    
    // calculate the distance of residue centers
    // _res_adj contain all the contacting residues for easier
    // searching hbond pairs

    for (int i = Start; i <= End; i++) {
      if (_res[i]._center._type != UNDEF) {
	for (int j = 1; j <= _numRes; j++) {
	  if (j == i)
	    continue;
	  if (_res[i]._center.dis(_res[j]._center) < 
	      (Residue::size[_res[i]._type] +
	       Residue::size[_res[j]._type] + CUB_SIZE)) {
		_res[i]._res_adj.insert(j);
		_res[j]._res_adj.insert(i);	
	  } else {
		  _res[i]._res_adj.erase(j);
		  _res[j]._res_adj.erase(i);
	  }
	}
      }
    }
  }
}

// write pdb file for the whole structure or a part of it,tp=0, all
// atom, tp=1 backbone only
void Structure::writePdb(string filename, int start, int end, int tp) {
  int i,j,atomNum=1,type;
  string Str;
  ofstream fout;
  fout.open(filename.c_str());

  if(start<1){
    cout << "starting position cannot smaller than 1 and reset to 1! "
	 << start << endl;
    start=1;
  }
  if(end>_numRes){
    cout<<"end position ("<<end<<") cannot exceed the total number of residues of the structure ("<<_numRes<<"). It is reset to the total number of residues."<<endl;
    end=_numRes;
  }
  for(i=start;i<=end;i++){
    type=_res[i]._type;
    Str=Residue::Name3[type];
    if(_res[i]._atom[1].x == 0 && _res[i]._atom[1].y == 0 && _res[i]._atom[1].z == 0)
    {
      fout<<"ATOM  "<<setw(5)<<atomNum++<<" ";
      fout<<" ";
      fout.setf(ios::left);
      fout<<setw(3)<<"H"<<" ";
      fout.unsetf(ios::left);
      fout<<Str<<"  "<<setw(4)<<i<<"    ";
      fout<<setiosflags(ios::fixed)<<setprecision(3)<<setw(8);
      fout<<_res[i]._atom[1].x<<setw(8)<<_res[i]._atom[1].y<<setw(8)<<_res[i]._atom[1].z;
      fout<<setiosflags(ios::fixed)<<setprecision(2)<<setw(6)<<1.0<<setw(6)<<1.0<<endl;
      
    }
    else{
     for(j=0;j<_res[i]._numAtom;++j){
      if(_res[i]._atom[j]._type==UNDEF ||_res[i]._atom[j]._type >=22) continue;
      if(tp==1&&j==NUM_BB_ATOM) break;
      if(type == GLY && Residue::cType[type][j] == "CB") continue;
      
      fout<<"ATOM  "<<setw(5)<<atomNum++<<" ";
      if(Residue::cType[type][j].length()<4) 
        fout<<" ";
      fout.setf(ios::left);
      fout<<setw(3)<<Residue::cType[type][j]<<" ";
      fout.unsetf(ios::left);
      fout<<Str<<"  "<<setw(4)<<i<<"    ";
      fout<<setiosflags(ios::fixed)<<setprecision(3)<<setw(8);
      fout<<_res[i]._atom[j].x<<setw(8)<<_res[i]._atom[j].y<<setw(8)<<_res[i]._atom[j].z;
      fout<<setiosflags(ios::fixed)<<setprecision(2)<<setw(6)<<1.0<<setw(6)<<1.0<<endl;
    }
   }
  }
  fout.close();
}


void Structure::writePdb(string filename, int start, int end, int l_start, int l_end, int tp) {
  int i,j,atomNum=1,type;
  string Str;
  ofstream fout;
  fout.open(filename.c_str());

  if(start<1){
    cout << "starting position cannot smaller than 1 and reset to 1! "
	 << start << endl;
    start=1;
  }
  if(end>_numRes){
    cout<<"end position ("<<end<<") cannot exceed the total number of residues of the structure ("<<_numRes<<"). It is reset to the total number of residues."<<endl;
    end=_numRes;
  }
  for(i=start;i<=end;i++){
    type=_res[i]._type;
    Str=Residue::Name3[type];
    if(i >= l_start && i <= l_end)
    {
     for(j=0;j< NUM_BB_ATOM;++j){
      if(_res[i]._atom[j]._type==UNDEF ||_res[i]._atom[j]._type >=22) continue;
      if(tp==1&&j==NUM_BB_ATOM) break;
      if(type == GLY && Residue::cType[type][j] == "CB") continue;
      
      fout<<"ATOM  "<<setw(5)<<atomNum++<<" ";
      if(Residue::cType[type][j].length()<4) 
        fout<<" ";
      fout.setf(ios::left);
      fout<<setw(3)<<Residue::cType[type][j]<<" ";
      fout.unsetf(ios::left);
      fout<<Str<<"  "<<setw(4)<<i<<"    ";
      fout<<setiosflags(ios::fixed)<<setprecision(3)<<setw(8);
      fout<<_res[i]._atom[j].x<<setw(8)<<_res[i]._atom[j].y<<setw(8)<<_res[i]._atom[j].z;
      fout<<setiosflags(ios::fixed)<<setprecision(2)<<setw(6)<<1.0<<setw(6)<<1.0<<endl;
    }
   }
   else
   {
     for(j=0;j<_res[i]._numAtom;++j){
      if(_res[i]._atom[j]._type==UNDEF ||_res[i]._atom[j]._type >=22) continue;
      if(tp==1&&j==NUM_BB_ATOM) break;
      if(type == GLY && Residue::cType[type][j] == "CB") continue;
      
      fout<<"ATOM  "<<setw(5)<<atomNum++<<" ";
      if(Residue::cType[type][j].length()<4) 
        fout<<" ";
      fout.setf(ios::left);
      fout<<setw(3)<<Residue::cType[type][j]<<" ";
      fout.unsetf(ios::left);
      fout<<Str<<"  "<<setw(4)<<i<<"    ";
      fout<<setiosflags(ios::fixed)<<setprecision(3)<<setw(8);
      fout<<_res[i]._atom[j].x<<setw(8)<<_res[i]._atom[j].y<<setw(8)<<_res[i]._atom[j].z;
      fout<<setiosflags(ios::fixed)<<setprecision(2)<<setw(6)<<1.0<<setw(6)<<1.0<<endl;

   }
  }
  }
  fout.close();
}



void Structure::writePdb(string filename, int start, int end, int tp, int md_num) {
  int i,j,atomNum=1,type;
  string Str;
  ofstream fout;
  fout.open(filename.c_str(),ios::app);

  if(start<1){
    cout<<"starting position cannot smaller than 1 and reset to 1! " <<start<<endl;
    start=1;
  }
  if(end>_numRes){
    cout<<"end position ("<<end<<") cannot exceed the total number of residues of the structure ("<<_numRes<<"). It is reset to the total number of residues."<<endl;
    end=_numRes;
  }
  
  fout << "MODEL"<<setw(9)<< md_num<<setw(6)<<"RMS"<<setw(6)<< setprecision(3) << G_rmsd<<setw(9)<<setprecision(4)<<_energy<<endl;
  for(i=start;i<=end;i++){
    type=_res[i]._type;
    Str=Residue::Name3[type];
    for(j=0;j<_res[i]._numAtom;++j){
      if(_res[i]._atom[j]._type==UNDEF ||_res[i]._atom[j]._type >= 22) continue;
      if(tp==1&&j==NUM_BB_ATOM) break;
      if(type == GLY && Residue::cType[type][j] == "CB") continue;
      fout<<"ATOM  "<<setw(5)<<atomNum++<<" ";
      if(Residue::cType[type][j].length()<4) 
        fout<<" ";
      fout.setf(ios::left);
      fout<<setw(3)<<Residue::cType[type][j]<<" ";
      fout.unsetf(ios::left);
      fout<<Str<<" A"<<setw(4)<<i<<"    ";
      fout<<setiosflags(ios::fixed)<<setprecision(3)<<setw(8);
      fout<<_res[i]._atom[j].x<<setw(8)<<_res[i]._atom[j].y<<setw(8)<<_res[i]._atom[j].z<<endl;
    }
  }
  fout << "ENDMDL" <<endl; 
  fout.close();
}


// Calculate backbone atom coordinates of a residue given the torsion angles
// The calculated coordinates are stored in atomArr
void Structure::calBBCo(int resInd, Residue& res, Residue& res2, double phi, double psi,
			double omega) {
  Point prev_atom[3];
  // These are defined because of the "cutting" form of growth
  int resType_prev = _res[resInd]._type;     // For C, O, CB, SC
  int resType_next = _res[resInd + 1]._type; // For N, CA, H

  // ATM_C, use phi angle, torAngles[0]
  prev_atom[0] = _res[resInd-1]._atom[ATM_C];
  prev_atom[1] = _res[resInd]._atom[ATM_N];
  prev_atom[2] = _res[resInd]._atom[ATM_CA];
  calCo(prev_atom,
	Residue::bond_length[resType_prev][ATM_C],
	Residue::bond_angle[resType_prev][ATM_C],
	phi,
	res._atom[ATM_C]);

  // ATM_O, use psi+PI
  prev_atom[0] = _res[resInd]._atom[ATM_N];
  prev_atom[1] = _res[resInd]._atom[ATM_CA];
  prev_atom[2] = res._atom[ATM_C];
  calCo(prev_atom,
	Residue::bond_length[resType_prev][ATM_O],
	Residue::bond_angle[resType_prev][ATM_O],
	psi + PI,
	res._atom[ATM_O]);

  // ATM_N of next residue, use psi angle, here the ATM_N of next
  // residue is stored at res2
  // this atom WAS actually stored at the same residue, res, before
  prev_atom[0] = _res[resInd]._atom[ATM_N];
  prev_atom[1] = _res[resInd]._atom[ATM_CA];
  prev_atom[2] = res._atom[ATM_C];
  calCo(prev_atom,
	Residue::bond_length[resType_next][ATM_N],
	Residue::bond_angle[resType_next][ATM_N],
	psi, res._atom[ATM_N]);
	
  // ATM_CA of next residue use omega angle, set to PI, stored as res2 
  // again it was stored at res
  prev_atom[0] = _res[resInd]._atom[ATM_CA];
  prev_atom[1] = res._atom[ATM_C];
  prev_atom[2] = res._atom[ATM_N];
  calCo(prev_atom,
	Residue::bond_length[resType_next][ATM_CA],
	Residue::bond_angle[resType_next][ATM_CA],
	omega, res._atom[ATM_CA]);
  // ATM_CB if not GLY, previous atoms are N, C, CA of the same residue, torsion angle is set to 122.55
  if (_res[resInd]._type != GLY) {
    prev_atom[0] = _res[resInd]._atom[ATM_N];
    prev_atom[1] = res._atom[ATM_C];
    prev_atom[2] = _res[resInd]._atom[ATM_CA];
    calCo(prev_atom,
	  Residue::bond_length[resType_prev][ATM_CB],
	  Residue::bond_angle[resType_prev][ATM_CB],
	  PI * 122.55 / 180,
	  res._atom[ATM_CB]);
    res._numAtom = 6;

    // calculate the seudo side chain atom SC and place it at res._SC
    prev_atom[0] = _res[resInd]._atom[ATM_N];
    prev_atom[1] = res._atom[ATM_C];
    prev_atom[2] = _res[resInd]._atom[ATM_CA];
    calCo(prev_atom,
	  Atom::R_SC[resType_prev],
	  Residue::bond_angle[resType_prev][ATM_CB],
	  PI * 122.55 / 180,
	  res._SC);
  } else {
    res._numAtom = 5;
  }
}

// calculate side chain atom coordinates of a residue, the type of the residue res._type should be given
void Structure::calSCCo(int resIndex, double* torAngles, Residue& res)
{
   
  int j,rotAngleCnt=0, numSCHeavyAtom = 0;
  Point tmpArr[3];
  double bond_length,bond_angle,torsion_angle;
  res._scc.reset();
  if(res._type == 0)
   res._scc = res._atom[5];
  else if(res._type == 5)
    res._scc = res._atom[1];
  else
   for(j=NUM_BB_ATOM;j<Residue::numAtom[res._type];++j){  // loop through each atom of the side-chain
    tmpArr[0] = res._atom[Residue::prev_atom[res._type][j][0]];
    tmpArr[1] = res._atom[Residue::prev_atom[res._type][j][1]];
    tmpArr[2] = res._atom[Residue::prev_atom[res._type][j][2]];
    bond_length=Residue::bond_length[res._type][j];
    bond_angle=Residue::bond_angle[res._type][j];
    if(Residue::torsion[res._type][j]==-1234){
      if(rotAngleCnt>=Rotamer::numRotBond[res._type]){
	cout<<"rotatable angle count error! "<<res._type<<" "<<Rotamer::numRotBond[res._type]<<" "<<rotAngleCnt<<endl;
	exit(0);
      }
      torsion_angle = torAngles[rotAngleCnt++];
    }
    else torsion_angle=Residue::torsion[res._type][j];	 
    calCo(tmpArr,bond_length,bond_angle,torsion_angle,res._atom[j]);
    if(isnan(res._atom[j].x) || isnan(res._atom[j].y) || isnan(res._atom[j].z)){
      cout<<"calSCCo error at "<<res._posn<<" "<<j<<endl; 
      exit(0);
    }
	if(res._atom[j]._type < 22){
		res._scc = (res._scc + res._atom[j]);
		numSCHeavyAtom++;
	}
  }
  res._numAtom = Residue::numAtom[res._type];
  if(res._type == 0 || res._type ==5)
    res._center = res._bbc; 
  else
  {
     res._scc = res._scc/numSCHeavyAtom;
     res._center = (res._bbc * NUM_BB_ATOM + res._scc * numSCHeavyAtom) / (NUM_BB_ATOM + numSCHeavyAtom);
  }
}


// grow side chains given backbone structure
// resToGrow contain the indices of sidechains to grow
// type == 0: grow sidechains from start to end including start and end
// type == 1: grow sidechains whose indices are in resToGrow
void Structure::grow_sc(int start, int end, int type, int numStates,
			int growType, vector <int>& resToGrow) {
  int i,j,size,n,m,selected;
  double p_selected;
  //Residue tmpRes[MAX_NUM_ST];	
  double enArr[MAX_NUM_SC_ST][ENERGY_TYPES];
  double energy[MAX_NUM_SC_ST];
  vector<double> prob(MAX_NUM_SC_ST,0.0);
  double sum = 0;
  double angles[MAX_NUM_SC_ST][6];
  Residue tmpRes[MAX_NUM_SC_ST];
  _weight_2 = 0;   // this is used for weight side chain sampling
  _toBeSampled.clear();	
  if(type == 0){  // grow from start to end
    for(i=start;i<=end;++i){
      if(_res[i]._type != GLY && _res[i]._type != ALA)
      {
	_toBeSampled.push_back(i);
      }
    }
  }
  else if(type == 1){  // use resToGrow which may not be a consecutive sequence of residues
    _toBeSampled = resToGrow;
  }
  else {
    cout<<"wrong type in grow_sc! "<<type<<endl;
    return;
  }
	
  for(i = 0; i<numStates;++i)
    for(j = 0; j<ENERGY_TYPES;++j)
      enArr[i][j] = 0;
  for(m=0;m < _toBeSampled.size();m++){
      n = _toBeSampled[m];
    sample_sc_angles((*this)._res[n], angles, numStates, growType);
    // calculate the coordinates of the side chains
	 // set this so that all the atom types of _res[n] will be copied to tmpRes[i]
    _res[n]._numAtom = Residue::numAtom[_res[n]._type];
    for(i=0;i<numStates;++i){
      tmpRes[i] = _res[n];
      calSCCo(n, angles[i], tmpRes[i]);
    }
    // calculate energy of those sampled side chain conformations

    for(i=0;i<numStates;++i){
      one_res_en_sc(*this, tmpRes[i], n, start, end, 1, enArr[i]);
      energy[i] = 0;
      for(j=0;j<ENERGY_TYPES;++j) {
	energy[i] += enArr[i][j];
      }
      //energy[i] = one_res_en_sc(*this, tmpRes[i], n, 1);		
    }

    // select one of them
    double minE = 10000;
    for (i = 0; i < numStates; i++)
      if (energy[i] < minE)
	minE = energy[i];
    sum = 0;
    for (i = 0; i < numStates; i++) {
      prob[i] = pow(EXPO, double(0.5*(minE - energy[i])/_T));
      sum += prob[i];
      if(isinf(prob[i])){
	  cout << "grow_sc probability error! " << energy[i] << " " << prob[i]
	     << endl;
	exit(0);
      }
    }
    selected = SampleOne(prob);
    p_selected = prob[selected]/sum;
    // update the conformation
    _res[n] = tmpRes[selected];
    if (tmpRes[selected]._numAtom == 8 &&
	abs(tmpRes[selected]._atom[8].x) > 0.01) {
      cout << "This is VERY bad at " << n << endl;
      exit(0);
    }
    _weight_2 += log(1/(p_selected*(numStates/(2*PI)))); 
    _energy += energy[selected];
    for(int i = 0; i<ENERGY_TYPES; ++i)
      _enArr[i] += enArr[selected][i];
    if(PF::cal[EM_RP])
      _res[n].cal_fg();  // update functional group center
    _res[n].cal_scc();  // update sidechain center 
  }
  // grow funtional group, _FG for ALA residues
  for(int i = start ; i <= end; ++i){
	  if(_res[i]._type == ALA)
		_res[i].cal_fg();
  }
}

// output conformation, if type==0 all atom, type==1 backbone
void Structure::output(int start, int end, int tp)
{
  int i,j,type,atomNum=1;
  string Str;

  for(i=start;i<=end;++i){
    cout<<_res[i]._numAtom<<endl;
    type=_res[i]._type;
    Str=Residue::Name3[type];
    for(j=0;j<_res[i]._numAtom;++j){
      if(tp==1&&i==NUM_BB_ATOM) break;
      cout<<"ATOM  "<<setw(5)<<atomNum++<<" ";
      if(Residue::cType[type][j].length()<4) 
        cout<<" ";
      cout.setf(ios::left);
      cout<<setw(3)<<Residue::cType[type][j]<<" ";
      cout.unsetf(ios::left);
      cout<<Str<<" A"<<setw(4)<<i<<"    ";
      cout<<setiosflags(ios::fixed)<<setprecision(3)<<setw(8);
      cout<<_res[i]._atom[j].x<<setw(8)<<_res[i]._atom[j].y<<setw(8)<<_res[i]._atom[j].z<<endl;
    }    
  }

}

// output backbone and side chain torsion angles  
void Structure::outputAngles(char* filename, int start, int end, int outType=0)
{
  int i,j,k,count,p1,p2,p3;
  double Bfactors[4][4];
  double tor;
  ofstream fout;
  fout.open(filename);
  if(outType==0)   // output torsion angles
    {
      for(i=start;i<=end;++i)
	{
	  // output backbone torsion angles
	  fout<<setw(5)<<_res[i]._pdbIndex<<" "<<setw(4)<<Residue::Name3[_res[i]._type]<<" ";
	  fout<<setw(8)<<setprecision(5)<<_res[i]._phi<<" "<<setw(8)<<_res[i]._psi<<" "<<setw(8)<<_res[i]._omega<<" ";
	  // output B-factors
	  fout<<_res[i-1]._atom[2]._Bfactor<<" ";
	  fout<<_res[i]._atom[0]._Bfactor<<" ";
	  fout<<_res[i]._atom[1]._Bfactor<<" ";
	  fout<<_res[i]._atom[2]._Bfactor<<" ";
	  fout<<_res[i+1]._atom[0]._Bfactor<<" ";
	  fout<<_res[i+1]._atom[1]._Bfactor<<endl;
	  // output Chi angles and B-factors
	  count=0;		
	  for(j=NUM_BB_ATOM;j<_res[i]._numAtom||_res[i]._atom[j]._type!=UNDEF;++j) { // here the condition is to accomodate the case that H atom is not counted sometimes, but it still occupy a place on the _atom array.
	    if(Residue::torsion[_res[i]._type][j]==-1234) {
	      p1 = Residue::prev_atom[_res[i]._type][j][0];
	      p2 = Residue::prev_atom[_res[i]._type][j][1];
	      p3 = Residue::prev_atom[_res[i]._type][j][2];
	      if(_res[i]._atom[j]._type!=UNDEF&&_res[i]._atom[p3]._type!=UNDEF&&_res[i]._atom[p2]._type!=UNDEF&&_res[i]._atom[p1]._type!=UNDEF)
		tor = torsion(_res[i]._atom[j],_res[i]._atom[p3],_res[i]._atom[p2],_res[i]._atom[p1]);
	      else tor = 999.99;    // 999.99 means undefined torsion angles
	      if(count==0)
		 fout<<setw(5)<<_res[i]._pdbIndex<<" "<<setw(4)<<Residue::Name3[_res[i]._type]<<" Chi ";
	      // fout<<" Chi "<<++count<<"  "<<setw(8)<<setprecision(5)<<tor<<" ";
	      fout<<setw(8)<<setprecision(5)<<tor<<" ";
	      Bfactors[count][0] = _res[i]._atom[j]._Bfactor;
	      Bfactors[count][1] = _res[i]._atom[p3]._Bfactor;
	      Bfactors[count][2] = _res[i]._atom[p2]._Bfactor;
	      Bfactors[count][3] = _res[i]._atom[p1]._Bfactor;
	      count++;
	    }
	  }
	  for(j=0;j<count;++j){
	    for(k=0;k<4;++k){
	      fout<<Bfactors[j][k]<<" ";
	    }
	  }
	  if(count>0)
	    fout<<endl;	
	}
    }
}

void Structure::analyticClosure(int start, int Start, int End, int* List, int List_size, bool Ellipsoid) {
  // Finds an analytic closure using the method of Seok et al
  
  //!-----------------------------------------------------------------------
  //! This is a sample driver routine to reconstruct tripeptide loops
  //! from the coordinates in a pdb file using a canonical bond lengths and
  //! angles.
  //!-----------------------------------------------------------------------
  //  use in_out
  //  use tripep_closure
  //!-----------------------------------------------------------------------
  //  implicit none
  //  real(dp) :: time1, time2
  //  character(len=100) :: in_pdb, out_prefix, out_pdb
  char out_prefix[100], out_pdb[100];
  //  integer :: n0, n_soln, unit = 17, i, j, k, n1, n2, order(max_soln)
  int n_soln;
  //  real(dp) :: r_n(3,5), r_a(3,5), r_c(3,5)
  double r_n[5][3], r_a[5][3], r_c[5][3];
  //  real(dp) :: r_soln_n(3,3,max_soln), r_soln_a(3,3,max_soln), r_soln_c(3,3,max_soln)
  double r_soln_n[max_soln][3][3], r_soln_a[max_soln][3][3], r_soln_c[max_soln][3][3], r_soln_o[max_soln][3][3], r_soln_cb[max_soln][3][3];
  //  real(dp) :: rmsd, sum, dr(3)
  double rmsd, sum, dr[3];
  //  real(dp) :: r0_n(3,3), r0_a(3,3), r0_c(3,3)
  double r0_n[3][3], r0_a[3][3], r0_c[3][3];
  
  //  logical :: calc_rmsd, write_out_pdb
  int  calc_rmsd, write_out_pdb, chosen;
  //  real(dp) :: b_len(6), b_ang(7), t_ang(2)
  double b_len[6], b_ang[7], t_ang[2];

  double MinE = 10000;   
  double Sum = 0; 
  //!-----------------------------------------------------------------------

  //  call my_timer(time1)

  //  !! initialize: set up input parameters
  //  ! bond lengths
  //  b_len(1:6) = (/ b_ac, b_cn, b_na, b_ac, b_cn, b_na /)
  b_len[0] = Residue::bond_length[_res[start]._type][ATM_C];
  b_len[1] = Residue::bond_length[_res[start+1]._type][ATM_N];
  b_len[2] = Residue::bond_length[_res[start+1]._type][ATM_CA];
  b_len[3] = Residue::bond_length[_res[start+1]._type][ATM_C];
  b_len[4] = Residue::bond_length[_res[start+2]._type][ATM_N];
  b_len[5] = Residue::bond_length[_res[start+2]._type][ATM_CA];
  
  //  ! bond angles
  //  b_ang(1:7) = (/ ang_nac, ang_acn, ang_cna, ang_nac, ang_acn, ang_cna, ang_nac /)
  b_ang[0] = Residue::bond_angle[_res[start]._type][ATM_C];
  b_ang[1] = Residue::bond_angle[_res[start+1]._type][ATM_N];
  b_ang[2] = Residue::bond_angle[_res[start+1]._type][ATM_CA];
  b_ang[3] = Residue::bond_angle[_res[start+1]._type][ATM_C];
  b_ang[4] = Residue::bond_angle[_res[start+2]._type][ATM_N];
  b_ang[5] = Residue::bond_angle[_res[start+2]._type][ATM_CA];
  b_ang[6] = Residue::bond_angle[_res[start+2]._type][ATM_C];
   
  //  ! peptide torsion angles
  //  t_ang(1:2) = pi
  t_ang[0] = PI;
  t_ang[1] = PI;
  //  call initialize_loop_closure(b_len, b_ang, t_ang)
  initialize_loop_closure(b_len, b_ang, t_ang);

  //  ! output options
  //  calc_rmsd = .true.
  calc_rmsd = 0;
  //  write_out_pdb = .false.
  write_out_pdb = 0;

  // Populate r_n, r_a, r_c from structure
  for (int i = 0; i < 5; i++) {
    r_n[i][0] = _res[i + start - 1]._atom[ATM_N].x;
    r_n[i][1] = _res[i + start - 1]._atom[ATM_N].y;
    r_n[i][2] = _res[i + start - 1]._atom[ATM_N].z;
    r_a[i][0] = _res[i + start - 1]._atom[ATM_CA].x;
    r_a[i][1] = _res[i + start - 1]._atom[ATM_CA].y;
    r_a[i][2] = _res[i + start - 1]._atom[ATM_CA].z;
    r_c[i][0] = _res[i + start - 1]._atom[ATM_C].x;
    r_c[i][1] = _res[i + start - 1]._atom[ATM_C].y;
    r_c[i][2] = _res[i + start - 1]._atom[ATM_C].z;
  }

  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) {
      r0_n[i][j]=r_n[i+1][j];
      r0_a[i][j]=r_a[i+1][j];
      r0_c[i][j]=r_c[i+1][j];
    }
  solve_3pep_poly(r_n[1], r_a[1], r_a[3], r_c[3], r_soln_n, r_soln_a, r_soln_c,
		  &n_soln);
   

  int k = 0;
  if (n_soln > 0) {

    double dis = 10000, disquare = 10000;
    double* Energy = new double [n_soln];
    vector <double>  prob;

    double** enArr = new double* [n_soln];
    for(k=0; k < n_soln; k++)
     enArr[k] = new double [ENERGY_TYPES];
	
    Point prev_atom[3];
    double psi;

    for(k=0; k < n_soln; k++)
    {
       Energy[k] = 0;
       for( int j=0; j < ENERGY_TYPES; j++)
        enArr[k][j] = 0;
    }
   
   // Choose a random solution to populate back to the structure
   for(k=0; k < n_soln; k++)
   {
    for (int i = 0; i < 3; i++) {
      _res[i + start]._atom[ATM_N].x = r_soln_n[k][i][0];
      _res[i + start]._atom[ATM_N].y = r_soln_n[k][i][1];
      _res[i + start]._atom[ATM_N].z = r_soln_n[k][i][2];
      _res[i + start]._atom[ATM_CA].x = r_soln_a[k][i][0];
      _res[i + start]._atom[ATM_CA].y = r_soln_a[k][i][1];
      _res[i + start]._atom[ATM_CA].z = r_soln_a[k][i][2];
      _res[i + start]._atom[ATM_C].x = r_soln_c[k][i][0];
      _res[i + start]._atom[ATM_C].y = r_soln_c[k][i][1];
      _res[i + start]._atom[ATM_C].z = r_soln_c[k][i][2];

    }

	for(int i = 0; i< 3; ++i) {
		// add other atoms

		// ATM_O, use psi+PI
	
		psi = torsion(_res[i+start]._atom[ATM_N], _res[i+start]._atom[ATM_CA], _res[i+start]._atom[ATM_C], _res[i+start+1]._atom[ATM_N])*PI/180;
		prev_atom[0] = _res[i+start]._atom[ATM_N];
		prev_atom[1] = _res[i+start]._atom[ATM_CA];
		prev_atom[2] = _res[i+start]._atom[ATM_C];
		calCo(prev_atom,
			Residue::bond_length[_res[i+start]._type][ATM_O],
			Residue::bond_angle[_res[i+start]._type][ATM_O],
			psi + PI,
			_res[i+start]._atom[ATM_O]);
               r_soln_o[k][i][0] = _res[i + start]._atom[ATM_O].x;
               r_soln_o[k][i][1] = _res[i + start]._atom[ATM_O].y;
               r_soln_o[k][i][2] = _res[i + start]._atom[ATM_O].z;

		// ATM_CB
  
		if (_res[i+start]._type != GLY) {
			prev_atom[0] = _res[i+start]._atom[ATM_N];
			prev_atom[1] = _res[i+start]._atom[ATM_C];
			prev_atom[2] = _res[i+start]._atom[ATM_CA];
			calCo(prev_atom,
				Residue::bond_length[_res[i+start]._type][ATM_CB],  
				Residue::bond_angle[_res[i+start]._type][ATM_CB],
				PI * 122.55 / 180,
				_res[i+start]._atom[ATM_CB]); 
               r_soln_cb[k][i][0] = _res[i + start]._atom[ATM_CB].x;
               r_soln_cb[k][i][1] = _res[i + start]._atom[ATM_CB].y;
               r_soln_cb[k][i][2] = _res[i + start]._atom[ATM_CB].z;
		} 
		else {
			_res[i+start]._numAtom = 5;
		}
		SinglecalCenter(_res[i+start], 1);
	}
       for(int i = 0; i< 3; ++i)
       {
	if(Ellipsoid)
	{
          one_res_en_list(*this, _res[i+start], Start, End, enArr[k],0, List, List_size);
    }
	else
	{
	    one_res_en(*this, _res[i+start], 1, i+start - 2, Start, End,  enArr[k], 1);
	    if (End != _numRes)
	    one_res_en(*this, _res[i+start], End +1  , _numRes, Start, End,  enArr[k], 2);
          for( int j = 0; j < _res[i+start]._numAtom; j++) 
    	  {
	      if (_res[i+start]._atom[j]._type == UNDEF || _res[i+start]._atom[j]._type >= 22)
   	      continue;
	      for(int q=1; q<=2; q++)
	      {
	          disquare = (_res[i+start]._atom[j].x - _res[End]._atom[q].x)*(_res[i+start]._atom[j].x - _res[End]._atom[q].x)
			 +(_res[i+start]._atom[j].y - _res[End]._atom[q].y)*(_res[i+start]._atom[j].y - _res[End]._atom[q].y) + (_res[i+start]._atom[j].z - _res[End]._atom[q].z)*(_res[i+start]._atom[j].z - _res[End]._atom[q].z);		                       if (disquare <= PF_DIS_CUT_SQUARE)
	       {
	          dis = sqrt(disquare);
	          int disInd = (int)(dis/H_INLO);
	          enArr[k][E_LOODIS] += PF::LOODIS[_res[i+start]._atom[j]._type-1][_res[End]._atom[q]._type-1][disInd];
               }
	      }
            }
	}
       }
      for (int j = 0; j < ENERGY_TYPES; j++)
      {
	   Energy[k] += enArr[k][j];
      }
      if(Energy[k] < MinE)
       MinE = Energy[k];
       
    }
    for(k =0; k <n_soln; k++)
    {  
       prob.push_back(pow(EXPO, (double(MinE - Energy[k])*0.5)));
	 if (isinf(prob[k])) {
	    cout << "grow_one probability error in state " << k << ": "
		 << " " << Energy[k] << " " << prob[k] << endl;
         exit(0);
         }
	 Sum += prob[k];
    }
    if (Sum == 0)
      chosen = intrand(0, n_soln);
    else
      chosen = SampleOne(prob);

    delete [] Energy;
    for(k =0; k <n_soln; k++)
     delete [] enArr[k];
    delete [] enArr;

    for (int i = 0; i < 3; i++) {
      _res[i + start]._atom[ATM_N].x = r_soln_n[chosen][i][0];
      _res[i + start]._atom[ATM_N].y = r_soln_n[chosen][i][1];
      _res[i + start]._atom[ATM_N].z = r_soln_n[chosen][i][2];
      _res[i + start]._atom[ATM_CA].x = r_soln_a[chosen][i][0];
      _res[i + start]._atom[ATM_CA].y = r_soln_a[chosen][i][1];
      _res[i + start]._atom[ATM_CA].z = r_soln_a[chosen][i][2];
      _res[i + start]._atom[ATM_C].x = r_soln_c[chosen][i][0];
      _res[i + start]._atom[ATM_C].y = r_soln_c[chosen][i][1];
      _res[i + start]._atom[ATM_C].z = r_soln_c[chosen][i][2];
      _res[i + start]._atom[ATM_O].x = r_soln_o[chosen][i][0];
      _res[i + start]._atom[ATM_O].y = r_soln_o[chosen][i][1];
      _res[i + start]._atom[ATM_O].z = r_soln_o[chosen][i][2];
      // ATM_CB
      if (_res[i+start]._type != GLY) {
        _res[i + start]._atom[ATM_CB].x = r_soln_cb[chosen][i][0];
        _res[i + start]._atom[ATM_CB].y = r_soln_cb[chosen][i][1];
        _res[i + start]._atom[ATM_CB].z = r_soln_cb[chosen][i][2];
	//		_res[i+start]._numAtom = 6;
      } 
      else {
        _res[i+start]._numAtom = 5;
      }
       SinglecalCenter(_res[i+start], 1);
    }
  }
  if (calc_rmsd) {
    for (int k = 0; k < n_soln; k++) {
      for (int i = 0; i < 3; i++)
	for (int j = 0; j < 3; j++) {
	  r_n[i+1][j] = r_soln_n[k][i][j];
	  r_a[i+1][j] = r_soln_a[k][i][j];
	  r_c[i+1][j] = r_soln_c[k][i][j];
	} 
      sum = 0.0e0;
      for(int i = 0; i < 3; i++) {
	for (int j = 0; j < 3; j++)
	  dr[j] = r_soln_n[k][i][j] - r0_n[i][j];
	sum += dot_product(dr,dr);
	for (int j = 0; j < 3; j++)
	  dr[j] = r_soln_a[k][i][j] - r0_a[i][j];
	sum += dot_product(dr,dr);
	for (int j = 0; j < 3; j++)
	  dr[j] = r_soln_c[k][i][j] - r0_c[i][j];
	sum += dot_product(dr,dr);
      }
      //           rmsd = sqrt(sum/9.0d0)
      rmsd = sqrt(sum/9.0e0);
      if (write_out_pdb) {
      }

      //        end do
    }
    //     end if
  }
  //  end do
}


void Structure::analyticClosure_h(int start, double len_change, double bon_change, double t_change, int Start, int End, int* List, int List_size, bool Ellipsoid) {
  // Finds an analytic closure using the method of Seok et al
  
  //!-----------------------------------------------------------------------
  //! This is a sample driver routine to reconstruct tripeptide loops
  //! from the coordinates in a pdb file using a canonical bond lengths and
  //! angles.
  //!-----------------------------------------------------------------------
  //  use in_out
  //  use tripep_closure
  //!-----------------------------------------------------------------------
  //  implicit none
  //  real(dp) :: time1, time2
  //  character(len=100) :: in_pdb, out_prefix, out_pdb
  char out_prefix[100], out_pdb[100];
  //  integer :: n0, n_soln, unit = 17, i, j, k, n1, n2, order(max_soln)
  int n_soln;
  //  real(dp) :: r_n(3,5), r_a(3,5), r_c(3,5)
  double r_n[5][3], r_a[5][3], r_c[5][3];
  //  real(dp) :: r_soln_n(3,3,max_soln), r_soln_a(3,3,max_soln), r_soln_c(3,3,max_soln)
  double r_soln_n[max_soln][3][3], r_soln_a[max_soln][3][3], r_soln_c[max_soln][3][3],
  r_soln_o[max_soln][3][3], r_soln_cb[max_soln][3][3];
  //  real(dp) :: rmsd, sum, dr(3)
  double rmsd, sum, dr[3];
  //  real(dp) :: r0_n(3,3), r0_a(3,3), r0_c(3,3)
  double r0_n[3][3], r0_a[3][3], r0_c[3][3];
  
  //  logical :: calc_rmsd, write_out_pdb
  int  calc_rmsd, write_out_pdb, chosen;
  //  real(dp) :: b_len(6), b_ang(7), t_ang(2)
  double b_len[6], b_ang[7], t_ang[2];
  
  double MinE = 10000;
  double Sum = 0;

  //!-----------------------------------------------------------------------

  //  call my_timer(time1)

  //  !! initialize: set up input parameters
  //  ! bond lengths
  //  b_len(1:6) = (/ b_ac, b_cn, b_na, b_ac, b_cn, b_na /)
  b_len[0] = frand((Residue::bond_length[_res[start]._type][ATM_C] - len_change), (Residue::bond_length[_res[start]._type][ATM_C] + len_change) );
  b_len[1] = frand((Residue::bond_length[_res[start+1]._type][ATM_N] - len_change), (Residue::bond_length[_res[start+1]._type][ATM_N] + len_change) );
  b_len[2] = frand((Residue::bond_length[_res[start+1]._type][ATM_CA] - len_change), (Residue::bond_length[_res[start+1]._type][ATM_CA] + len_change) );
  b_len[3] = frand((Residue::bond_length[_res[start+1]._type][ATM_C] - len_change), (Residue::bond_length[_res[start+1]._type][ATM_C] + len_change));
  b_len[4] = frand((Residue::bond_length[_res[start+2]._type][ATM_N] - len_change), (Residue::bond_length[_res[start+1]._type][ATM_N] + len_change));
  b_len[5] = frand((Residue::bond_length[_res[start+2]._type][ATM_CA]- len_change), (Residue::bond_length[_res[start+2]._type][ATM_CA]+ len_change)) ;
  
  //  ! bond angles
  //  b_ang(1:7) = (/ ang_nac, ang_acn, ang_cna, ang_nac, ang_acn, ang_cna, ang_nac /)
  b_ang[0] = frand( (Residue::bond_angle[_res[start]._type][ATM_C] - bon_change), (Residue::bond_angle[_res[start]._type][ATM_C] + bon_change));
  b_ang[1] = frand( (Residue::bond_angle[_res[start+1]._type][ATM_N] - bon_change), (Residue::bond_angle[_res[start+1]._type][ATM_N] + bon_change) );
  b_ang[2] = frand( (Residue::bond_angle[_res[start+1]._type][ATM_CA] - bon_change), (Residue::bond_angle[_res[start+1]._type][ATM_CA] + bon_change));
  b_ang[3] = frand( (Residue::bond_angle[_res[start+1]._type][ATM_C] - bon_change), (Residue::bond_angle[_res[start+1]._type][ATM_C] + bon_change));
  b_ang[4] = frand( (Residue::bond_angle[_res[start+2]._type][ATM_N] - bon_change), (Residue::bond_angle[_res[start+2]._type][ATM_N] + bon_change));
  b_ang[5] = frand( (Residue::bond_angle[_res[start+2]._type][ATM_CA] -bon_change), (Residue::bond_angle[_res[start+2]._type][ATM_CA] + bon_change));
  b_ang[6] = frand( (Residue::bond_angle[_res[start+2]._type][ATM_C] -bon_change), (Residue::bond_angle[_res[start+2]._type][ATM_C] + bon_change));
  //  ! peptide torsion angles
  //  t_ang(1:2) = pi
  t_ang[0] = frand((PI - t_change), (PI + t_change));
  t_ang[0] = t_ang[0] > PI ? (t_ang[0] - 2*PI):t_ang[0];
  t_ang[1] = frand((PI - t_change), (PI + t_change));
  t_ang[1] = t_ang[1] > PI ? (t_ang[1] - 2*PI):t_ang[1];
  //  call initialize_loop_closure(b_len, b_ang, t_ang)
  initialize_loop_closure(b_len, b_ang, t_ang);

  //  ! output options
  //  calc_rmsd = .true.
  calc_rmsd = 0;
  //  write_out_pdb = .false.
  write_out_pdb = 0;

  // Populate r_n, r_a, r_c from structure
  for (int i = 0; i < 5; i++) {
    r_n[i][0] = _res[i + start - 1]._atom[ATM_N].x;
    r_n[i][1] = _res[i + start - 1]._atom[ATM_N].y;
    r_n[i][2] = _res[i + start - 1]._atom[ATM_N].z;
    r_a[i][0] = _res[i + start - 1]._atom[ATM_CA].x;
    r_a[i][1] = _res[i + start - 1]._atom[ATM_CA].y;
    r_a[i][2] = _res[i + start - 1]._atom[ATM_CA].z;
    r_c[i][0] = _res[i + start - 1]._atom[ATM_C].x;
    r_c[i][1] = _res[i + start - 1]._atom[ATM_C].y;
    r_c[i][2] = _res[i + start - 1]._atom[ATM_C].z;
  }

  //     ! save initial coord for rmsd 
  //     r0_n(:,:) = r_n(:,2:4)
  //     r0_a(:,:) = r_a(:,2:4)
  //     r0_c(:,:) = r_c(:,2:4)
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) {
      r0_n[i][j]=r_n[i+1][j];
      r0_a[i][j]=r_a[i+1][j];
      r0_c[i][j]=r_c[i+1][j];
    }
  //     ! call tripeptide loop closure routine
  solve_3pep_poly(r_n[1], r_a[1], r_a[3], r_c[3], r_soln_n, r_soln_a, r_soln_c,
		  &n_soln);
  int k = 0;
  if (n_soln > 0) {
    // Choose a random solution to populate back to the structure
    double dis = 10000, disquare = 10000;
    double* Energy = new double [n_soln];
    vector <double>  prob;
    double** enArr = new double* [n_soln];
    for(k=0; k < n_soln; k++)
      enArr[k] = new double [ENERGY_TYPES];
    Point prev_atom[3];
    double psi;
    
    for(k=0; k < n_soln; k++)
    {
	 Energy[k] = 0;
	 for( int j=0; j < ENERGY_TYPES; j++)
	 enArr[k][j] = 0;
    }

   for(k=0; k < n_soln; k++)
   {
      for (int i = 0; i < 3; i++) {
           _res[i + start]._atom[ATM_N].x = r_soln_n[k][i][0];
	   _res[i + start]._atom[ATM_N].y = r_soln_n[k][i][1];
	   _res[i + start]._atom[ATM_N].z = r_soln_n[k][i][2];
	   _res[i + start]._atom[ATM_CA].x = r_soln_a[k][i][0];		        
	   _res[i + start]._atom[ATM_CA].y = r_soln_a[k][i][1];			 
           _res[i + start]._atom[ATM_CA].z = r_soln_a[k][i][2];				        
	   _res[i + start]._atom[ATM_C].x = r_soln_c[k][i][0];
           _res[i + start]._atom[ATM_C].y = r_soln_c[k][i][1];					     
	   _res[i + start]._atom[ATM_C].z = r_soln_c[k][i][2];
      }

	for(int i = 0; i< 3; ++i) {
		// add other atoms

		// ATM_O, use psi+PI
	
		psi = torsion(_res[i+start]._atom[ATM_N], _res[i+start]._atom[ATM_CA], _res[i+start]._atom[ATM_C], _res[i+start+1]._atom[ATM_N])*PI/180;
		prev_atom[0] = _res[i+start]._atom[ATM_N];
		prev_atom[1] = _res[i+start]._atom[ATM_CA];
		prev_atom[2] = _res[i+start]._atom[ATM_C];
		calCo(prev_atom,
			Residue::bond_length[_res[i+start]._type][ATM_O],
			Residue::bond_angle[_res[i+start]._type][ATM_O],
			psi + PI,
			_res[i+start]._atom[ATM_O]);
                r_soln_o[k][i][0] = _res[i + start]._atom[ATM_O].x;
		r_soln_o[k][i][1] = _res[i + start]._atom[ATM_O].y;	         
		r_soln_o[k][i][2] = _res[i + start]._atom[ATM_O].z;
		// ATM_CB
  
		if (_res[i+start]._type != GLY) {
			prev_atom[0] = _res[i+start]._atom[ATM_N];
			prev_atom[1] = _res[i+start]._atom[ATM_C];
			prev_atom[2] = _res[i+start]._atom[ATM_CA];
			calCo(prev_atom,
				Residue::bond_length[_res[i+start]._type][ATM_CB],  
				Residue::bond_angle[_res[i+start]._type][ATM_CB],
				PI * 122.55 / 180,
				_res[i+start]._atom[ATM_CB]); 
                 r_soln_cb[k][i][0] = _res[i + start]._atom[ATM_CB].x;
                 r_soln_cb[k][i][1] = _res[i + start]._atom[ATM_CB].y;
                 r_soln_cb[k][i][2] = _res[i + start]._atom[ATM_CB].z;
		} 
		else {
			_res[i+start]._numAtom = 5;
		}
		SinglecalCenter(_res[i+start], 1);
	}
    for(int i = 0; i< 3; ++i)
	{
        if(Ellipsoid)
	    {
	     one_res_en_list(*this, _res[i+start], Start, End, enArr[k],0, List, List_size);
        }
	  else
	  {                         
	     one_res_en(*this, _res[i+start], 1, i+start - 1, Start, End,  enArr[k], 1);
	     if (End != _numRes)
	     one_res_en(*this, _res[i+start], End +1  , _numRes, Start, End,  enArr[k], 2);
             for( int j = 0; j < _res[i+start]._numAtom; j++) 
	     {
	       if (_res[i+start]._atom[j]._type == UNDEF || _res[i+start]._atom[j]._type >= 22)
   	        continue;
	       for(int q=1; q<=2; q++)
	       {
	            disquare = (_res[i+start]._atom[j].x - _res[End]._atom[q].x)*(_res[i+start]._atom[j].x - _res[End]._atom[q].x)
			 +(_res[i+start]._atom[j].y - _res[End]._atom[q].y)*(_res[i+start]._atom[j].y - _res[End]._atom[q].y) 
			 + (_res[i+start]._atom[j].z - _res[End]._atom[q].z)*(_res[i+start]._atom[j].z - _res[End]._atom[q].z);		                  
               if (disquare <= PF_DIS_CUT_SQUARE)
	       {
	          dis = sqrt(disquare);
	          int disInd = (int)(dis/H_INLO);
	          enArr[k][E_LOODIS] += PF::LOODIS[_res[i+start]._atom[j]._type-1][_res[End]._atom[q]._type-1][disInd];
               }
	       }
              }
           }
          }
          for (int j = 0; j < ENERGY_TYPES; j++)
	  {
	     Energy[k] += enArr[k][j];
	  }
	  if(Energy[k] < MinE)
	        MinE = Energy[k];
    }
    for(k =0; k <n_soln; k++)
    {  
       prob.push_back(pow(EXPO, (double(MinE - Energy[k])*0.5)));
	 if (isinf(prob[k])) {
	    cout << "grow_one probability error in state " << k << ": "
		 << " " << Energy[k] << " " << prob[k] << endl;
         exit(0);
         }
	 Sum += prob[k];
    }
    if (Sum == 0)
      chosen = intrand(0, n_soln);
    else
      chosen = SampleOne(prob);
    delete [] Energy;
    for(k =0; k <n_soln; k++)
     delete [] enArr[k];
    delete [] enArr;

    for (int i = 0; i < 3; i++) {
      _res[i + start]._atom[ATM_N].x = r_soln_n[chosen][i][0];
      _res[i + start]._atom[ATM_N].y = r_soln_n[chosen][i][1];
      _res[i + start]._atom[ATM_N].z = r_soln_n[chosen][i][2];
      _res[i + start]._atom[ATM_CA].x = r_soln_a[chosen][i][0];
      _res[i + start]._atom[ATM_CA].y = r_soln_a[chosen][i][1];
      _res[i + start]._atom[ATM_CA].z = r_soln_a[chosen][i][2];
      _res[i + start]._atom[ATM_C].x = r_soln_c[chosen][i][0];
      _res[i + start]._atom[ATM_C].y = r_soln_c[chosen][i][1];
      _res[i + start]._atom[ATM_C].z = r_soln_c[chosen][i][2];
      _res[i + start]._atom[ATM_O].x = r_soln_o[chosen][i][0];
      _res[i + start]._atom[ATM_O].y = r_soln_o[chosen][i][1];
      _res[i + start]._atom[ATM_O].z = r_soln_o[chosen][i][2];
      // ATM_CB
      if (_res[i+start]._type != GLY) {
        _res[i + start]._atom[ATM_CB].x = r_soln_cb[chosen][i][0];
        _res[i + start]._atom[ATM_CB].y = r_soln_cb[chosen][i][1];
        _res[i + start]._atom[ATM_CB].z = r_soln_cb[chosen][i][2];
	//		_res[i+start]._numAtom = 6;
      } 
      else {
        _res[i+start]._numAtom = 5;
      }
       SinglecalCenter(_res[i+start], 1);
    }
    // _weight += log((double)n_soln);
  }
  //     if (calc_rmsd) then
  //        ! backbone coord in pdb 
  //        do k = 1, n_soln
  if (calc_rmsd) {
    for (int k = 0; k < n_soln; k++) {
      //           ! copy coord
      //           do i = 2, 4
      //              r_n(:,i) = r_soln_n(:,i-1,k)
      //              r_a(:,i) = r_soln_a(:,i-1,k)
      //              r_c(:,i) = r_soln_c(:,i-1,k)
      //           end do
      for (int i = 0; i < 3; i++)
	for (int j = 0; j < 3; j++) {
	  r_n[i+1][j] = r_soln_n[k][i][j];
	  r_a[i+1][j] = r_soln_a[k][i][j];
	  r_c[i+1][j] = r_soln_c[k][i][j];
	} 
      //           ! rmsd
      //           sum = 0.0d0
      sum = 0.0e0;
      //           do i = 2, 4
      //              dr(:) = r_n(:,i) - r0_n(:,i-1)
      //              sum = sum + dot_product(dr,dr)
      //              dr(:) = r_a(:,i) - r0_a(:,i-1)
      //              sum = sum + dot_product(dr,dr)
      //              dr(:) = r_c(:,i) - r0_c(:,i-1)
      //              sum = sum + dot_product(dr,dr)
      //           end do
      for(int i = 0; i < 3; i++) {
	for (int j = 0; j < 3; j++)
	  dr[j] = r_soln_n[k][i][j] - r0_n[i][j];
	sum += dot_product(dr,dr);
	for (int j = 0; j < 3; j++)
	  dr[j] = r_soln_a[k][i][j] - r0_a[i][j];
	sum += dot_product(dr,dr);
	for (int j = 0; j < 3; j++)
	  dr[j] = r_soln_c[k][i][j] - r0_c[i][j];
	sum += dot_product(dr,dr);
      }
      rmsd = sqrt(sum/9.0e0);
      //           write(*,"(f15.3)") rmsd
      //           if (write_out_pdb) then
      //              out_pdb = trim(out_prefix)??'_'??char(ichar('0')+k)??'.pdb'
      //              write(*,*) ' Writing soln ', k, ' in file ', trim(out_pdb)
      //              call write_pdb_backbone(out_pdb, res_name, r_n, r_a, r_c, n0-1, n0+3)
      //           end if
      if (write_out_pdb) {
      }
      //        end do
    }
    //     end if
  }
  //  end do
}
double Structure::RMSD(Structure other, int start, int end) {
  // Compute the backbone RMSD of this structure against "other"
  double sum = 0;
  for (int i = start; i <= end; i++)
    for (int j = 0; j < 3; j++)
      sum += (_res[i]._atom[j] - other._res[i]._atom[j]).square();
  sum /= end - start + 1;
  return(sqrt(sum));
}

bool Structure::IsClosed(int End) {
  // Tests whether this structure is closed, using a given endpoint
  // Also tests the N-CA distance, since this may change as the part of the
  // closure
  // procedure where the CA is 
  double CurDistance1 = _res[End]._atom[ATM_N].dis(_res[End]._atom[ATM_CA]);
  double CurDistance2 = _res[End]._atom[ATM_CA].dis(_res[End]._atom[ATM_C]);
  double CurAngle1 = angle(_res[End]._atom[ATM_N],
			   _res[End]._atom[ATM_CA],
			   _res[End]._atom[ATM_C]);
  double CurAngle2 = angle(_res[End]._atom[ATM_CA],
			   _res[End]._atom[ATM_C],
			   _res[End+1]._atom[ATM_N]);
  double CurTorsion1 = torsion(_res[End]._atom[ATM_CA],
			       _res[End]._atom[ATM_C],
			       _res[End+1]._atom[ATM_N],
			       _res[End+1]._atom[ATM_CA]);
  return (CurDistance1 > 1.358 && CurDistance1 < 1.558 &&
	  CurDistance2 > 1.425 && CurDistance2 < 1.625 &&
	  CurAngle1 > 86.1 && CurAngle1 < 136.1 &&
	  CurAngle2 > 95 && CurAngle2 < 135 &&
	  (CurTorsion1 > 160 | CurTorsion1 < -160));
}

bool Structure::IsClosed(int End, Residue NewRes) {
  // Tests whether this structure would be closed if it had NewRes at
  // position End
  Residue oldRes[2];
  oldRes[0] = _res[End];
  oldRes[1] = _res[End + 1];
  _res[End]._atom[ATM_C] = NewRes._atom[ATM_C];
  _res[End + 1]._atom[ATM_N] = NewRes._atom[ATM_N];
  _res[End + 1]._atom[ATM_CA] = NewRes._atom[ATM_CA];
  bool Result = IsClosed(End);
  _res[End]     = oldRes[0];
  _res[End + 1] = oldRes[1];
  return(Result);
}

// calculate the coordinates of SC atom, the squdo-sidechain atom
void Structure::calSC(int start, int end) {
  Point prev_atom[3]; 
  for(int i = start; i <= end; ++i) {
    if(_res[i]._type == GLY) continue;
    else if(_res[i]._type == ALA) {
      _res[i]._SC = _res[i]._atom[ATM_CB];
      _res[i]._SC._name = "SC";
    } else {
      prev_atom[0] = _res[i]._atom[ATM_N];
      prev_atom[1] = _res[i]._atom[ATM_C];
      prev_atom[2] = _res[i]._atom[ATM_CA];
      calCo(prev_atom,
	    Atom::R_SC[_res[i]._type],
	    Residue::bond_angle[_res[i]._type][ATM_CB],
	    PI*122.55/180,_res[i]._SC);
      // set the type of SC atom the same as the CB of ALA atom
      _res[i]._SC._type = 5;
    }
  }
}


// calculate centers for a single residue
// type = 1 calculate only backbone center, _bbc, and residue center _center
// type = 0 calculate backbone center, side chain center and residue center 
void Structure::SinglecalCenter(Residue& _res, int type) {
  int numBBAtom = 0;
  int numSCAtom = 0;
  int j;
  _res._bbc = Point(0,0,0);
  _res._scc = Point(0,0,0);
  _res._center = Point(0,0,0);
  for(j=0;j<NUM_BB_ATOM;++j){
    if(_res._atom[j]._type>20) continue;
    if(_res._atom[j].x== 0 && _res._atom[j].y== 0 && _res._atom[j].z== 0) continue;
    if(_res._atom[j]._type!=UNDEF){
      numBBAtom++;
      _res._bbc=_res._bbc+_res._atom[j];
      _res._center=_res._center+_res._atom[j];
    }
  }
  if(type == 0) {
    for(j = NUM_BB_ATOM; j < _res._numAtom; ++j){
      if(_res._atom[j]._type>20) continue;
      if(_res._atom[j].x== 0 && _res._atom[j].y== 0 && _res._atom[j].z== 0) continue;
      if(_res._atom[j]._type!=UNDEF){
        numSCAtom++;
        _res._scc=_res._scc+_res._atom[j];
        _res._center=_res._center+_res._atom[j];
      }
    }
  }
  if(numBBAtom!=0) {
    _res._bbc=_res._bbc/numBBAtom;
  }
  else _res._bbc._type=UNDEF;
  if(numSCAtom!=0){
    _res._scc=_res._scc/numSCAtom;
  }
  else if(_res._type==0){  // ALA
    _res._scc=_res._atom[5];
  }
  else if(_res._type==5){ // GLY
    _res._scc=_res._atom[1];
  }
  else _res._scc._type=UNDEF;

  if((numBBAtom+numSCAtom)!=0){
    _res._center=_res._center/(numBBAtom+numSCAtom);
    _res._center._type=0;
  }
  else _res._center._type=UNDEF;
}

// calculate the contact numbers
vector<int> Structure::contactNum(int start, int end, double radius)
{
	if (end > _numRes)  end=_numRes;
	vector<int> contactNumList(_numRes+1);
	
	
	for(int i=start; i<=end;i++) {
		for (int j=0; j<=_numRes;j++){
			double distance = _res[i]._center.dis(_res[j]._center);
			if   ( (distance < radius)  && (i != j) ) {
				contactNumList[i]++;
			}
		}
	}
	return contactNumList;
}

// check whether any atom coordinates are nan
void Structure::checkNAN(int start, int end)
{
  for(int i = 1; i<_numRes; ++i){
    for(int j = 0; j<_res[i]._numAtom; ++j){
      if(isnan(_res[i]._atom[j].x) || isnan(_res[i]._atom[j].y) || isnan(_res[i]._atom[j].z)) {
	cout<<"Found NAN at residue "<<i<<" atom "<<j<<" "<<start<<" "<<end<<endl;
	exit(0);
      }
      if(_res[i]._atom[j]._type != Residue::vdwType[_res[i]._type][j]){
	cout<<"Type error at residue "<<i<<" atom "<<j<<" "<<_res[i]._atom[j]._type<<endl;
	exit(0);
      }
    }
  }
}

// roughly calculate the center of the whole protein except the loop region
void Structure::calConfCenter(int start, int end)
{
  int totalatom = 0;
  for(int i = 1; i < start; i++)
  {
     for (int j = 0; j < _res[i]._numAtom; j++)
     {
       if (_res[i]._atom[j]._type > 20)
       continue;
       if (_res[i]._atom[j]._type != UNDEF) {
           _Confcenter = _Confcenter +  _res[i]._atom[j];
          totalatom++;
         if(_res[i]._atom[j].x > _Fringevtx[0])
	   _Fringevtx[0] = _res[i]._atom[j].x;
         if(_res[i]._atom[j].x < _Fringevtx[1])
	   _Fringevtx[1] = _res[i]._atom[j].x;
         if(_res[i]._atom[j].y > _Fringevtx[2])
	   _Fringevtx[2] = _res[i]._atom[j].y;
         if(_res[i]._atom[j].y < _Fringevtx[3])
	   _Fringevtx[3] = _res[i]._atom[j].y;
         if(_res[i]._atom[j].z > _Fringevtx[4])
	   _Fringevtx[4] = _res[i]._atom[j].z;
         if(_res[i]._atom[j].z < _Fringevtx[5])
	   _Fringevtx[5] = _res[i]._atom[j].z;
       }
     }
  }

  for(int i = end+1; i <= _numRes; i++)
  {
     for (int j = 0; j < _res[i]._numAtom; j++)
     {
       if (_res[i]._atom[j]._type > 20)
       continue;
       if (_res[i]._atom[j]._type != UNDEF) {
           _Confcenter = _Confcenter + _res[i]._atom[j];
          totalatom++;
         if(_res[i]._atom[j].x > _Fringevtx[0])
	   _Fringevtx[0] = _res[i]._atom[j].x;
         if(_res[i]._atom[j].x < _Fringevtx[1])
	   _Fringevtx[1] = _res[i]._atom[j].x;
         if(_res[i]._atom[j].y > _Fringevtx[2])
	   _Fringevtx[2] = _res[i]._atom[j].y;
         if(_res[i]._atom[j].y < _Fringevtx[3])
	   _Fringevtx[3] = _res[i]._atom[j].y;
         if(_res[i]._atom[j].z > _Fringevtx[4])
	   _Fringevtx[4] = _res[i]._atom[j].z;
         if(_res[i]._atom[j].z < _Fringevtx[5])
	   _Fringevtx[5] = _res[i]._atom[j].z;
       }
     }

  }
  _Confcenter = _Confcenter/totalatom;
}

