// rotamer.cpp

#include "rotamer.h"
#include "residue.h"
#include <new>

int RotLib::scLib;
vector<Rotamer> RotLib::RLR[20];
LRVMAP RotLib::RLD[20];

int Rotamer::numRotBond[20]={
    0,1,2,3,2,    // A C D E F 
    0,2,2,4,2,    // G H I K L
    3,2,2,3,4,    // M N P Q R
    1,1,1,2,2     // S T V W Y
};

void RotLib::read_rotlib(char* scFile, int scLib, const string dir,int outInfo)
{
  int i,j,numRot,resInt,tmpInt;
  char tmpLine[1000],tempCh[200];
  string tmpStr,rotName;
  ifstream inFile;
  double tmpF;
  Rotamer* tmpRot;
  ROTVEC* tmpRotVec;    
  int rotCount=1,aaIndex;
  double valueHolder[10];
  char ChHolder[3][50];
  string strHolder[5];
  long key;   // index key for SRL_D
  LRVMAP::iterator rvItr;
  // read side chain rotamer library of Dunbrack
  if(scLib==2){
    strcpy(scFile, FILE_ROTPT_DUN);
    if(outInfo==1)
    cout<<"Reading side chain library of Dunbrack's ..."<<endl;

    // Read the side chain rotamer library //
    strcpy(tempCh,dir.c_str());
    strcat(tempCh,scFile);
    inFile.open(tempCh,ios::in);
    if(!inFile.is_open()){
      cout<<"cannot find rotamer library "<<scFile<< " at "<<tmpLine<<endl;
      exit(0);
    }

    while(!inFile.eof()){
      inFile.getline(tmpLine,1000);      
      tmpStr=tmpLine;
	  if(tmpStr[0] == '#') continue;
      if(tmpStr.length()<10)
	break;
      // read one line, 
      // strHolder[0] store the name of the amino acid
      // strHolder[1] and strHolder[2] store the phi and psi angles of the backbone
      // valueHolder[0] store the probability of the rotamer
      // valueHolder[1-4] store the X angles
      // valueHolder[5-8] store the range of X angles
      // some values in valueHolder[] may be 0.0 if that angle does not exist.
      sscanf(tmpLine,"%s %s %s %lf %f %f %f %f %f %f %f %f",ChHolder[0],ChHolder[1],ChHolder[2],&valueHolder[0],&valueHolder[1],&valueHolder[2],&valueHolder[3],&valueHolder[4],&valueHolder[5],&valueHolder[6],&valueHolder[7],&valueHolder[8]);  
      strHolder[0]=ChHolder[0];strHolder[1]=ChHolder[1];strHolder[2]=ChHolder[2];
      aaIndex=Residue::AIMap.find(strHolder[0])->second;
      tmpRot= new Rotamer(Rotamer::numRotBond[aaIndex]);
      tmpRot->prob=valueHolder[0];
      for(i=0;i<Rotamer::numRotBond[aaIndex];++i){
        tmpRot->chi[i].angle=(PI/180)*valueHolder[i+1];
        tmpRot->chi[i].range=(PI/180)*valueHolder[i+4+1];
      }
      key=(atoi(strHolder[1].c_str())+180)*1000+(atoi(strHolder[2].c_str())+180);
      rvItr=RLD[aaIndex].find(key);
      if(rvItr!=RLD[aaIndex].end()){
	rvItr->second.push_back(*tmpRot);
      }
      else {
        tmpRotVec = new ROTVEC;
	if(tmpRotVec==NULL){
	  cout<<"no memory!"<<endl;
	  break;
	}
	tmpRotVec->push_back(*tmpRot);
        RLD[aaIndex].insert(LRVMAP::value_type(key,*tmpRotVec));
      }      
    }
  } // finish reading Dunbrack rotamer library
  
    // Read the side chain rotamer library of Richardson
  
  else if(scLib==1){
    strcpy(scFile, FILE_ROTLIB_R);
    if(outInfo==1)
    cout<<"Reading side chain library of Richardson's ..."<<endl;

    strcpy(tempCh,dir.c_str());
    strcat(tempCh,scFile);
    inFile.open(tempCh,ios::in);
    if(!inFile.is_open()){
      cout<<"cannot find rotamer library "<<scFile<< " at "<<tmpLine<<endl;
      exit(0);
    }
    while(!inFile.eof()){
      inFile.getline(tmpLine,1000);      
      tmpStr=tmpLine;
	  if(tmpStr[0] == '#') continue;
      if(tmpStr.length()<10)
	break;
      sscanf(tmpLine,"%s %lf %f %f %f %f",ChHolder[0],&valueHolder[0],&valueHolder[1],&valueHolder[2],&valueHolder[3],&valueHolder[4]); 
      strHolder[0]=ChHolder[0];
      aaIndex=Residue::AIMap.find(strHolder[0])->second;
      tmpRot= new Rotamer(Rotamer::numRotBond[aaIndex]);
      tmpRot->prob=-1*log(valueHolder[0]/100);
      for(i=0;i<Rotamer::numRotBond[aaIndex];++i){
        tmpRot->chi[i].angle=(PI/180)*valueHolder[i+1];
        tmpRot->chi[i].range=20;
      }
      // insert the rotamer to library
      RLR[aaIndex].push_back(*tmpRot);
    }
  }
}
