// reprst.cpp 

#include "reprst.h"
#include "util.h"
#include "residue.h"
#include "rotamer.h"
#include "potential.h"
      
vector<SCT> SCR::SCANGPDB[NUM_RES_TP];
FIMAP SCR::SCRMAP[NUM_RES_TP];

// initialize side chain angle distribution, which is used in sampling
// side chain torsion angles
void SCR::InitSCAng(int scType, char* scFile, string dir) {
  int i,j,aaIndex,numRot;
  char filename[100], tmpLine[500];
  ifstream inFile;
  string tmpStr,aaName;
  vector <string> strVec;
  SIMAP::iterator SIitr;
  double torsion[4],bFacAve[4];
  SCT tmpSCT;
  // read side chain torsion angles from scFile
  inFile.open(scFile,ios::in);
  if(!inFile.is_open()){
    strcpy(filename,dir.c_str());
    strcat(filename,scFile);
    inFile.open(filename,ios::in);
    inFile.clear();
    inFile.seekg(0,ios::beg);
    
    if(!inFile.is_open()){
      cout<<"cannot open side chain torsion angle file at "<<filename<<endl;
      return;
    }
  }
    int key;
    double prob, cumProb, count, totalCnt[NUM_RES_TP];
    vector<int> scDataKey[NUM_RES_TP];
    vector<double> scDataCnt[NUM_RES_TP];

    for (i = 0; i < NUM_RES_TP; i++)
      totalCnt[i] = 0;
    while (!inFile.eof()) {
      inFile.getline(tmpLine,500);
      tmpStr = tmpLine;
      if(tmpStr[0] == '#') continue;
      strVec.clear();
      split(tmpStr, " ", strVec);
      aaIndex = atoi(strVec[0].c_str());
      numRot = Rotamer::numRotBond[aaIndex];
      if(strVec.size() < (numRot + 3)) continue;
      key=0;
      for (j=1; j <= numRot; ++j)
	key = key*100 + atoi(strVec[j].c_str());
      count = atof(strVec[numRot+2].c_str());
      scDataKey[aaIndex].push_back(key);
      scDataCnt[aaIndex].push_back(count);
      totalCnt[aaIndex] += count;
    }
    for(i=0;i<NUM_RES_TP;++i){
      if(scDataKey[i].size()==0) continue;
      cumProb=0;
      for(j=0;j<scDataKey[i].size();++j){
	prob = scDataCnt[i][j]/totalCnt[i];
	cumProb += prob;
	SCRMAP[i].insert(FIMAP::value_type(cumProb,scDataKey[i][j]));
      }

      // make sure the last one has cumulative probability equals to 1
      if (cumProb < 1)
	SCRMAP[i].insert(FIMAP::value_type(1,
					   scDataKey[i][scDataKey[i].size()-1]));
    }
}

// test continuous side chain representation by sampling a number of angles and check
// whether they have the same distribution as observed from PDB
void SCR::testSCRep(int resType, int numStates, char* outFile)
{
  int i,j,numRot,key;
  double r,angles[6];
  FIMAP::iterator fiItr;
  ofstream fout;
  fout.open(outFile);
  numRot = Rotamer::numRotBond[resType];
  for(i=0;i<numStates;++i){
    r = frand(0,1);
    if((fiItr = SCR::SCRMAP[resType].upper_bound(r))!=SCR::SCRMAP[resType].end()){
      key = fiItr->second;
      for(j=(numRot-1);j>=0;--j){
	angles[j] = (key%100*SC_T_INT + frand(0,SC_T_INT) - 180);
	key = key/100;
	fout<<angles[j]<<" ";
      }
      fout<<endl;
    }
  }
}
