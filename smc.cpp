// SMC code

#include <algorithm>
#include "cal_energy.h"
#include "smc.h"
#include "sample_states.h"
#include "util.h"
#include <iostream>
#include <fstream>
#include <cstring>
#include <string>
#include <sstream>

using namespace std;


void k_mean_clustering(vector<vector<double> >& matrix,int K,vector<int>& centerMin,vector<int>& idclusterMin);

void SMC::Wholeproc(vector<Structure>& Topconflist)
{
  int NClosed =0, NEnergy = 0;
  vector<combiEIndex> Erank;
  vector<combiEIndex> Erank_small;
  vector<int> ResIdx;
  vector<int> ClashNum; 
  combiEIndex temp;
  double MinEnergy = 1000000, minE = 10000, minEidx = 0;
  int _model_num = 0, t=0, Meet_Num = 0, size_keep = 0, index = 0;
  double change_rot = PI/36;  //5 degree
  bool BBclash_flag = false;

  
  if(!noScore)
  {
     if(Surface || Ellipsoid)
       calE_list(Conf, Reslist, Start, End, List_size, true);
     else 
      calE(Conf, Start, End, true);
  }

  cout<<"DiSGro in progress ..."<<endl;
  PreProcess();
  time_t t1 = time(NULL);
  int Counting = 0;
  /********Do SMC*************/
  for(t=0; t< NumConf;t++)
  {
    smc();
    if(_conf.Success)
    {
     Counting ++;
     temp.energy = _conf._energy;
     temp.index = NumClosedconf;
     if(!noScore)
     Erank.push_back(temp);
   }
    if(NumClosedconf == 5000)
    break;
  }
  time_t t2 = time(NULL);
  cout << "Conformational sampling done, time cost  "<< (t2-t1)<<" s"<<endl;
  cout << "Counting: " <<Counting<< endl;
  Structure tmpconf;
  Structure copyconf;
  tmpconf = Conf;
  copyconf = LoopStore[0];
	
  // rank by energy
  time_t t3 = time(NULL);
   if (Erank.size()) {
       for(int j=0;j<Erank.size();j++)
         Erank[j].index = j;
       size_keep = (Erank.size() > confkeep)? confkeep : Erank.size();
       sort(Erank.begin(), Erank.end(), sortfun_E);
       for(int j=0;j<size_keep;j++)
	   {
	    index = Erank[j].index;
	    for(int k=Start;k <= End;k++)
		{
	      tmpconf._res[k] = LoopStore[index]._res[k-Start+1];
        }
	   tmpconf._energy  = Erank[j].energy;
       BBclash_flag = false;
	   if(Ellipsoid)
	   BBClash_detection_list(tmpconf, Start, End, ResIdx, ClashNum, Reslist, List_size);
	   else
	   BBClash_detection(tmpconf, Start, End, ResIdx, ClashNum);
 
       for(int t=0 ; t< ResIdx.size(); t++)
       {
	      if(ClashNum[t] > 10)
	      {
	         BBclash_flag = true;
	         break;
	      }
       }
	   ResIdx.clear();
	   ClashNum.clear();
	   if(BBclash_flag)
	   continue;
	  
       if(_isSampSC)
       {
	    tmpconf.grow_sc(Start, End, 0, numSCStates, AngType, tmpconf._toBeSampled);
	    if(Ellipsoid)
        {
	      Clash_detection_list(tmpconf, Start, End, ResIdx, ClashNum, Reslist, List_size);
	      SCE_Minimization_list(tmpconf, Start, End, ResIdx,ClashNum, Reslist, List_size, change_rot);
	      ClashNum.clear();
	      ResIdx.clear();
	      Clash_detection_list(tmpconf, Start, End, ResIdx, ClashNum, Reslist, List_size);
	    }
	    
		if(Eval)
		{
          tmpconf._energy = 0;
	      tmpconf.calCenter(Start,End,true);	   
	      if(Surface || Ellipsoid)
	      calE_list(tmpconf, Reslist, Start, End, List_size, true);
	      else
	      calE(tmpconf, Start, End, true);
		}
		else
        tmpconf._energy  = Erank[j].energy;
      }
	  else
        tmpconf._energy  = Erank[j].energy;
	    if(tmpconf._energy < minE)
	    {
	       minE = tmpconf._energy;
	       minEidx = j;	   
		   if(outputconf == 1)
		   {
	         Topconflist.clear();		   
	         Topconflist.push_back(tmpconf);		   
		   }
	    } 
	 
	    Erank[j].energy = tmpconf._energy;
	    if(outputconf > 1)
		{
	     for(int k=1; k <= End - Start + 1 ;k++)
	      copyconf._res[k] = tmpconf._res[k + Start - 1];
          copyconf._energy = tmpconf._energy;
	      Topconflist.push_back(copyconf);
		}
	    Erank_small.push_back(Erank[j]);
      }
      LoopStore.clear();
      
	  if(Eval)
	  {
        sort(Erank_small.begin(), Erank_small.end(), sortfun_E);
        if(outputconf > 1)
        sort(Topconflist.begin(),Topconflist.end(),sortfun_StrE);
	  }
      
      int final_size = Erank_small.size(); 
      if(!kmcluster)
	    cout << "Calculation Done!!"<<endl;
	    
	  //RMSD Calculations
	  double rms = 0;
	  double rms_norm = 0;
	  double Rms = 0;
	  int Size;
	  vector<double> onedim;	
	  vector<vector<double> > Matrix;
	  
	  
	  ofstream myfile;
	  ofstream myfile1;
	  myfile.open("RSMD.pdb");
	  myfile1.open("RSMD_norm.xlsx");
	  if(Topconflist.size() < 30000)
         Size = Topconflist.size();
       else
         Size = 30000;  
	  for(int i=0;i<Size;i++)
       {   
		//rms = Root_MSD(Topconflist[i], Conf, 1, End - Start + 1, 1, End - Start + 1,1,0);
		//rms_norm = Root_MSD(Topconflist[i], Conf, 1, End - Start + 1, 1, End - Start + 1,1,1);
		Rms = Root_MSD(Conf, Topconflist[i], Start, End, 1, End - Start + 1,1,0);
		myfile << i << ", " << Rms <<endl;
		//myfile1 << i << ", " << rms_norm << endl;        
	   }
  }

    //K-mean clustering
    if(kmcluster == true)
    {
	   cout << "Do k-mean clustering ..."<<endl;
       vector<vector<double> > Matrix;
       vector<int> centermin;
       vector<int> clustermin;
       vector<double> onedim;
       vector<Structure> tmpConflist;
       double rms = 0;
       double minRMS = 10000;
       double aveRMS = 0;
       double sumRMS = 0;
       int Size;
       int ksize = 50;

       if(Topconflist.size() < 30000)
         Size = Topconflist.size();
       else
         Size = 30000;
       if(Size < 2*ksize)
       {
         cout << "Not Enough Samples!!!!!!!!!!!"<<endl;
		 exit(0);
       }
       for(int i=0;i<Size;i++)
       {
         for(int j=0;j<Size;j++)
	     {
             if(i<j)
	         {
                rms = Root_MSD(Topconflist[i], Topconflist[j], 1, End - Start + 1, 1, End - Start + 1,1,0);
                 
		        onedim.push_back(rms);
             }
	         else if (i == j)
		       onedim.push_back(0);
	         else
		       onedim.push_back(Matrix[j][i]);
	     }
	     Matrix.push_back(onedim);
	     onedim.clear();
       }

     k_mean_clustering(Matrix, ksize, centermin,clustermin);
	 for(int i=0; i < centermin.size(); i++)
     {
	  tmpConflist.push_back(Topconflist[centermin[i]]);
	 }
     Topconflist.clear();
     for(int i=0; i < ksize; i++)
     Topconflist.push_back(tmpConflist[i]);
	 cout <<"Calculation Done !!"<<endl;
   }
}

// sequential monte carlo
void SMC::smc() {
  Structure tmploop(End-Start+1);
  _conf = Conf;
  NumSConf = 100;
  _conf._energy = 0;
  vector<int> prev_NumDistanceStates = NumDistanceStates;
  vector<int> prev_NumAngleStates = NumAngleStates;
  bool prev_Close = Close;
  bool grow_Success = false, grow_preSuccess = false;  // label the residue successfully generated or not
  double eedis = _conf._res[Start]._atom[ATM_CA].dis(_conf._res[End]._atom[ATM_CA]);    //End to End distance between start and end residues
  double Ecutoff = 350;
  int confnum_ub = 10000;
  _conf.Closed = false;
  _conf.Success = false;
  for (int j = 0; j < ENERGY_TYPES; j++)
     _conf._enArr[j] = 0;
  if(End == _conf._numRes)
  {
      int i;
      for(i=Start;i<End;i++)
      {
   	    grow_Success = grow_one(_conf, i, End, _conf._res[End]._atom[ATM_C]);
        
	    if(!grow_Success)
        {
	        _conf.Success = false; 
			return;
        }
      }
      
  }
  else
  {
   for(int i=Start;i<End-2;i++)
   {
      grow_Success = grow_one(_conf, ((Backward) ? (Start + End - i) : i), End, _conf._res[End]._atom[ATM_C]);
	  if(!grow_Success)
      {
	    _conf.Success = false ;
	    return;
      }
   }
  }
   /********** Set closure indicators *************/
   Point next_atom[3];
   if(End-Start < 4)
   {
     grow_preSuccess = grow_one(_conf, End-2, End, _conf._res[End]._atom[ATM_C]);
     if(grow_preSuccess)
     {
     grow_Success = grow_one(_conf, End-1, End, _conf._res[End]._atom[ATM_C]);
      
      next_atom[0] = _conf._res[End+1]._atom[ATM_CA];
      next_atom[1] = _conf._res[End+1]._atom[ATM_N];
      next_atom[2] = _conf._res[End]._atom[ATM_C];
      calCo(next_atom,
      Residue::bond_length[_conf._res[End]._type][ATM_C],
      Residue::bond_angle[_conf._res[End+1]._type][ATM_N],PI, _conf._res[End]._atom[ATM_CA]);
    }

    if(grow_Success && grow_preSuccess)
     _conf.Closed = _conf.IsClosed(End);
   }
     if (Close && (_conf.Closed == false)) {
	// The algorithm assumes the final CA is placed in accordance
	// with the following atoms, which grow_one does not guarantee
	// currently. Here we place it using torsion angle omega = PI.
    _conf.analyticClosure(End - 2, Start, End, Reslist, List_size, Ellipsoid);
     }

      /********** Set closure indicators *************/
      if (Close && ((Start != 1 && UseBackward) || End != Conf._numRes))
      {
       if(grow_Success && grow_preSuccess)
       {
	     _conf.Closed = _conf.IsClosed(End);
       }  
       else if(!grow_Success && grow_preSuccess)
       {	
        _conf.Closed = _conf.IsClosed(End-1);
       }
       else
       {
        _conf.Closed = _conf.IsClosed(End-2);
		  
	    if(!_conf.Closed)
	    {
          for( int i =0; i < 300; i++)
	      {
	         _conf.analyticClosure_h(End - 2, 0.06, PI/32, PI/32, Start, End, Reslist, List_size, Ellipsoid);
             _conf.Closed = _conf.IsClosed(End-2);
	         if(_conf.Closed)
	         break;
	      }
	    }
       }
      }

      if ((Start == 1 && UseBackward) || End == Conf._numRes) {
      // Replace original values of overwritten members
      NumDistanceStates = prev_NumDistanceStates;
      NumAngleStates    = prev_NumAngleStates;
      Close = prev_Close;
      }
     
      /************** If End = final residue, do some special processing  ***********************/
      if (End == Conf._numRes) {
	// Replace original values of overwritten members
	double rPhiAng = frand((-PI), PI);
	double rPsiAng = frand((-PI), PI);
	Point prev_atom[3];				
	Residue& ttRes = _conf._res[End];  // the last residue		
	// place final ATM_C of chain 
	prev_atom[0] = _conf._res[End-1]._atom[ATM_C];
	prev_atom[1] = _conf._res[End]._atom[ATM_N];
	prev_atom[2] = _conf._res[End]._atom[ATM_CA];
	calCo(prev_atom,
	      Residue::bond_length[_conf._res[End]._type][ATM_C],
	      Residue::bond_angle[_conf._res[End]._type][ATM_C],
	      rPhiAng,
	      _conf._res[End]._atom[ATM_C]);
	// place ATM_O atom on ATM_C
	prev_atom[0] = _conf._res[End]._atom[ATM_N];
	prev_atom[1] = _conf._res[End]._atom[ATM_CA];
	prev_atom[2] = _conf._res[End]._atom[ATM_C];
	calCo(prev_atom,
	      Residue::bond_length[_conf._res[End]._type][ATM_O],
	      Residue::bond_angle[_conf._res[End]._type][ATM_O],
	      rPsiAng,
	      _conf._res[End]._atom[ATM_O]);
		
	// place ATM_OXT on ATM_C
	prev_atom[0] = _conf._res[End]._atom[ATM_N];
	prev_atom[1] = _conf._res[End]._atom[ATM_CA];
	prev_atom[2] = _conf._res[End]._atom[ATM_C];
	calCo(prev_atom,
	      Residue::bond_length[_conf._res[End]._type][ATM_O],
	      Residue::bond_angle[_conf._res[End]._type][ATM_O],
	      rPsiAng+PI,
	      _conf._ATM_OXT);
     
	// sample side chains
	if(ttRes._type == GLY) {
	  ttRes._numAtom = 6;   // 5 backbone atoms plus OXT	
	}
	else{
	  // place ATM_CB on ATM_CA
	  prev_atom[0] = _conf._res[End]._atom[ATM_N];
	  prev_atom[1] = _conf._res[End]._atom[ATM_C];
	  prev_atom[2] = _conf._res[End]._atom[ATM_CA];
	  calCo(prev_atom,
		Residue::bond_length[_conf._res[End]._type][ATM_CB],
		Residue::bond_angle[_conf._res[End]._type][ATM_CB],
		PI*122.55/180,
		_conf._res[End]._atom[ATM_CB]);
	}
	// copy ATM_OXT to the last residue. This copy needs to be done here
	ttRes._atom[ttRes._numAtom] = _conf._ATM_OXT;
   }
   if(_conf.Closed || (Start == 1 && UseBackward) || End == Conf._numRes)
       NumClosedconf++;
    if(!noScore)
	if (_conf.Closed) {
    // save energy
	   _conf._energy = 0;
	   _conf.calCenter(Start, End, false);
  	  if(Surface || Ellipsoid)
        calE_list(_conf, Reslist, Start, End, List_size, false);
      else
        calE(_conf, Start, End, false);
	}
   //store conformations
   if(_conf.Closed || (Start == 1 && UseBackward) || End == Conf._numRes)
   {
      if(_conf._energy < minEnergy + Ecutoff && LoopStore.size() < confnum_ub )   //control the ensemble size
      {
        if(minEnergy > _conf._energy)
		{
		    minEnergy = _conf._energy;
	    }
        tmploop._energy = _conf._energy;
        
	    if(!noScore)                      //Don't do k-mean clustering, keep top low energy conformations
        {
	      for (int i=Start; i < End +1; i++)
	      {
	    	 tmploop._res[i-Start+1] = _conf._res[i];

	     }
           LoopStore.push_back(tmploop);
        }
        _conf.Success = true; // do not change
      }
	  else
        _conf.Success = false;

  }

 }
  





// A number of candidate conformations for a residue are sampled
// one of them is chosen and its weight is updated.
bool SMC::grow_one(Structure& CurConf, int Position, int tmpEnd, Atom& EndPt) {
  // Initalize positional values
  int RelPosition = Position - Start;
  int RemLength   = tmpEnd - Position;
  int LargeNumDistanceStates = 160;
  int phi_int, psi_int;
  
  // Initialize numbers of states
  int CurNumAngleStates    =    NumAngleStates[RelPosition];
  int CurNumDistanceStates = NumDistanceStates[RelPosition];
  int NumStates = CurNumAngleStates + CurNumDistanceStates;
  if(!TorsionStyle)
  LargeNumDistanceStates = NumStates;
    if (Backward) {
      RelPosition = tmpEnd - Position;
      RemLength = Position - Start;
    }

  // The current residue
  int CurResType = CurConf._res[Position]._type;
  // Placeholders for sampling
  double bb_angles[NumStates][3];
  double larbb_angles[LargeNumDistanceStates][3];

  double energy[NumStates];
  double axisdist[NumStates];  // the distance between current atom and the protein axis, specific used for beta barrel membrane protein
  double enArr[NumStates][ENERGY_TYPES];
  int Phi[LargeNumDistanceStates];
  int Psi[LargeNumDistanceStates];
  int TorProb[LargeNumDistanceStates];
  int deadlabel[LargeNumDistanceStates];
  int binNum = 180/BBTbinSize;
  Residue tmpRes[NumStates];		// position n
  Residue tmpRes2[NumStates];		// position n+1
//  Residue tmpRes3[TotalStates];		// for side chain position n
  vector<double> prob(NumStates, 0.0);
  double EEdisC = CurConf._res[Position]._atom[ATM_CA].dis(EndPt);
  double EEdisN = 0;
  double minEEDisC = minDistcon[1][RemLength];
  double maxEEDisC = minDistcon[1][RemLength]+ 31*DistconBy[1][RemLength];
  double minEEDisN  = minDistcon[0][RemLength-1];
  double maxEEDisN  = minDistcon[0][RemLength-1] + 31*DistconBy[0][RemLength-1];

  // Sample angles
  if (CurNumAngleStates > 0)
  {
       if(AngType == 3)
       {
          bb_angles[0][0] = torsion(Conf._res[Position-1]._atom[2],
			      Conf._res[Position]._atom[0],
			      Conf._res[Position]._atom[1],
			      Conf._res[Position]._atom[2]) * PI/180;
          bb_angles[0][1] = torsion(Conf._res[Position]._atom[0],
			      Conf._res[Position]._atom[1],
			      Conf._res[Position]._atom[2],
			      Conf._res[Position+1]._atom[0]) * PI/180;
          bb_angles[0][2] = torsion(Conf._res[Position]._atom[1],
			      Conf._res[Position]._atom[2],
			      Conf._res[Position+1]._atom[0],
			      Conf._res[Position+1]._atom[1]) * PI/180;
           sample_bb_angles(CurConf, (Backward) ? (Position - 1) : Position, bb_angles, CurNumAngleStates, AngType);
       }
   else
       sample_bb_angles(CurConf, (Backward) ? (Position - 1) : Position, bb_angles, CurNumAngleStates, AngType);
  }


  // Sample end-to-end distances and convert to positions
  if (CurNumDistanceStates > 0) {
    // Loop through states
    for (int i = 0; i < LargeNumDistanceStates; i++) {
      // Initialize four residues for this sample
        Residue sam;
     // Sample new points for C atom based on distance to final C
        TorProb[i] = -1;
    if( EEdisC >= minEEDisC  && EEdisC <= maxEEDisC)
	{
         sample_distance(CurConf._res[Position]._atom[ATM_N],
		      CurConf._res[Position]._atom[ATM_CA],
		      EndPt,
		      Residue::bond_angle[CurResType][ATM_C],
		      Residue::bond_length[CurResType][ATM_C],
		      sam._atom[ATM_C], EEdisC,
		      1, RemLength,
		      false);
	}
	else
	{
	return false;
	}
    // Sample new points for N atom based on distance to final C
	if(isnan(sam._atom[ATM_C].x))
	{
	 TorProb[i] = 0;
     deadlabel[i] = 1;
	 continue;
	}
    EEdisN = sam._atom[ATM_C].dis(EndPt);
    if(EEdisN >= minEEDisN && EEdisN <= maxEEDisN)
	{
	sample_distance(CurConf._res[Position]._atom[ATM_CA],
			sam._atom[ATM_C],
			EndPt,
			Residue::bond_angle[CurResType][ATM_N],
			Residue::bond_length[CurResType][ATM_N],
			sam._atom[ATM_N], EEdisN,
			0, RemLength - 1, // -1 since it's NEXT N
			false);
	}
	else
	{
	TorProb[i] = 0;
	deadlabel[i] = 1;
	continue;
	}

	if(isnan(sam._atom[ATM_N].x))
	{
	 TorProb[i] = 0;
	 deadlabel[i] = 1;
	 continue;
	}
      // Populate remainder of bb_angles with the torsion angles
      // implied by the positional selections for C and N. This is for
      // unity with the existing code base. Omega = pi by default.
	larbb_angles[i][0] =
	  torsion(CurConf._res[Position-1]._atom[ATM_C],
		  CurConf._res[Position]._atom[ATM_N],
		  CurConf._res[Position]._atom[ATM_CA],
		  sam._atom[ATM_C]) * PI/180;
	larbb_angles[i][1] =
	  torsion(CurConf._res[Position]._atom[ATM_N],
		  CurConf._res[Position]._atom[ATM_CA],
		  sam._atom[ATM_C],
		  sam._atom[ATM_N]) * PI/180;
        box_MullerNsample_single(larbb_angles[i][2],PI,4);
	if(!TorsionStyle)
	{
	 bb_angles[i][0] = larbb_angles[i][0];
	 bb_angles[i][1] = larbb_angles[i][1];
	 bb_angles[i][2] = larbb_angles[i][2];

	}
     }
    }

  if(TorsionStyle)
  {
    int Total = 0;
    for(int i=0; i < LargeNumDistanceStates; i++) {

    if(TorProb[i] == -1)
    {
     phi_int = (int)(larbb_angles[i][0]*180/PI/BBTbinSize)+binNum;  
     psi_int = (int)(larbb_angles[i][1]*180/PI/BBTbinSize)+binNum; 

    if(larbb_angles[i][0]*180/PI/BBTbinSize+binNum - phi_int < 0)
    phi_int--;
    if(larbb_angles[i][1]*180/PI/BBTbinSize+binNum - psi_int < 0)
    psi_int--;
    Phi[i] = phi_int;
    Psi[i] = psi_int;

    TorProb[i] = Joint_Angle[CurResType][phi_int][psi_int];
    Total += TorProb[i]; 
    }
   }
  
   for(int i=0; i < NumStates; i++) {
   if(Total == 0)
   {
      int r= rand() % LargeNumDistanceStates;
      while(deadlabel[r] == 1 && r < LargeNumDistanceStates)
      r++;
      
      bb_angles[i][0] = frand((double)((Phi[r]-binNum)*PI/binNum),(double)((Phi[r]-binNum+1)*PI/binNum));
      bb_angles[i][1] = frand((double)((Psi[r]-binNum)*PI/binNum),(double)((Psi[r]-binNum+1)*PI/binNum));
      bb_angles[i][2] = larbb_angles[r][2];
      continue;
   }
   else
   {
    int r = rand() % Total;
    int cur = 0;
    for(int j=0; j< LargeNumDistanceStates; j++)
    {
      cur += TorProb[j]; 
      if( r < cur)
      {
         bb_angles[i][0] = frand((double)((Phi[j]-binNum)*PI/binNum), (double)((Phi[j]-binNum+1)*PI/binNum));
         bb_angles[i][1] = frand((double)((Psi[j]-binNum)*PI/binNum), (double)((Psi[j]-binNum+1)*PI/binNum));
         bb_angles[i][2] = larbb_angles[j][2];
      	 break;
      } 
     }
    }
   }
  }

  double dis = 10000, disquare = 10000;
 
 for(int i = 0; i < NumStates; i++)
  {
   energy[i] = 0;
   axisdist[i] = 0;
   for (int j = 0; j < ENERGY_TYPES; j++)
       enArr[i][j] = 0;
  }

  // Loop through states, converting sampled angles to positions
  for(int i = 0; i < NumStates; i++)
  {
    tmpRes[i] = CurConf._res[Position];

    if(Backward)
       tmpRes2[i] = CurConf._res[Position-1];
    else
       tmpRes2[i] = CurConf._res[Position+1];
    tmpRes[i]._phi   = bb_angles[i][0];
    tmpRes[i]._psi   = bb_angles[i][1];
    tmpRes[i]._omega = bb_angles[i][2];
    CurConf.calBBCo(Position, tmpRes[i], tmpRes2[i],
      bb_angles[i][0], bb_angles[i][1], bb_angles[i][2]);
    // Calculate residue center
    CurConf.SinglecalCenter(tmpRes[i], 1);
     
    // calculate energy
    if(Backward) {
    one_res_en(CurConf, tmpRes[i], Position + 1, CurConf._numRes, Start, End,  enArr[i],     
                   2);
    }
    else{ 
       if(Ellipsoid)
       {
        // Protein Body
         one_res_en_list(CurConf, tmpRes[i], Start, End, enArr[i],1, Reslist, List_size);
       }
       else
       {
        one_res_en(CurConf, tmpRes[i], 1, Position - 1, Start, End,  enArr[i], 1);
        if (End != CurConf._numRes)
        one_res_en(CurConf, tmpRes[i], End +1  , CurConf._numRes, Start, End,  enArr[i], 2);
       
        for( int j = 0; j < tmpRes[i]._numAtom; j++) 
        {
         if (tmpRes[i]._atom[j]._type == UNDEF || tmpRes[i]._atom[j]._type >= 22)
	 continue;
         if (tmpRes[i]._atom[j].x == 0 && tmpRes[i]._atom[j].y == 0 && tmpRes[i]._atom[j].z == 0) continue;
	 for(int k = 1; k<= 2;k++)
	 {
            disquare = (tmpRes[i]._atom[j].x-CurConf._res[End]._atom[k].x)*(tmpRes[i]._atom[j].x-CurConf._res[End]._atom[k].x)+
	 	    (tmpRes[i]._atom[j].y-CurConf._res[End]._atom[k].y)*(tmpRes[i]._atom[j].y-CurConf._res[End]._atom[k].y)+
                    (tmpRes[i]._atom[j].z-CurConf._res[End]._atom[k].z)*(tmpRes[i]._atom[j].z-CurConf._res[End]._atom[k].z);
            if (disquare <= PF_DIS_CUT_SQUARE)
	    {
	       dis = sqrt(disquare);
	       int disInd = (int)(dis/H_INLO);
	       enArr[i][E_LOODIS] += PF::LOODIS[tmpRes[i]._atom[j]._type-1][CurConf._res[End]._atom[k]._type-1][disInd];
	    }
	 }
        }
       }
	
       disquare = (tmpRes[i]._atom[1].x-CurConf._res[Position]._atom[0].x)*(tmpRes[i]._atom[1].x-CurConf._res[Position]._atom[0].x)+
	 	  (tmpRes[i]._atom[1].y-CurConf._res[Position]._atom[0].y)*(tmpRes[i]._atom[1].y-CurConf._res[Position]._atom[0].y)+
                  (tmpRes[i]._atom[1].z-CurConf._res[Position]._atom[0].z)*(tmpRes[i]._atom[1].z-CurConf._res[Position]._atom[0].z);
       if (disquare <= PF_DIS_CUT_SQUARE)
       {
	  dis = sqrt(disquare);
	  int disInd = (int)(dis/H_INLO);
	  enArr[i][E_LOODIS] += PF::LOODIS[tmpRes[i]._atom[1]._type-1][CurConf._res[Position]._atom[0]._type-1][disInd];
       }
    }
	for(int j=0; j <ENERGY_TYPES; j++)
	{
		energy[i] += enArr[i][j];
	}
  }  

  double minE = 100000;
  double minDist = 100000;
  for(int i = 0; i < NumStates; ++i)
  {
    if(energy[i] < minE)
      minE = energy[i];
    if(axisdist[i] < minDist)
      minDist = axisdist[i]; 
  }

  // Now find the probabilities of states and select one
  double sum=0, sumE = 0, sumD = 0;
 
    for(int i =0; i <NumStates; i++)
    {
        prob[i] = pow(EXPO, double( (minE - energy[i])*0.5));
	    if (isinf(prob[i])) {
      		cout << "grow_one probability error in state " << i << ": "
	        << " " << energy[i] << " " << prob[i] << endl;
     	 exit(0);
   	 }
         sum += prob[i];
    }
   

   // in case that all energies are very high, prob[i] are all zero
    // assign all prob[i] to 1
    if (sum == 0)
    {
      if(Refine)
       return false;
      else
      {
       for (int i = 0; i < NumStates; i++)
       {
          prob[i] = 1; 
       }
      }
    }
    int chosen = SampleOne(prob);

  // Update the weight of the chain
  CurConf._weight += log(VectorSum(prob) / prob[chosen]);
  // keep track of the individual energy types
  for (int i = 0; i < ENERGY_TYPES; ++i) {
    CurConf._enArr[i] += enArr[chosen][i];
    CurConf._energy += enArr[chosen][i];
  }

  CurConf._res[Position]._atom[ATM_C].CopyPos(tmpRes[chosen]._atom[ATM_C]);
  CurConf._res[Position]._atom[ATM_O].CopyPos(tmpRes[chosen]._atom[ATM_O]);
  CurConf._res[Position+1]._atom[ATM_N].CopyPos(tmpRes[chosen]._atom[ATM_N]);
  CurConf._res[Position+1]._atom[ATM_CA].CopyPos(tmpRes[chosen]._atom[ATM_CA]);

  if (CurConf._res[Position]._type != GLY) {
    CurConf._res[Position]._atom[ATM_CB].CopyPos(tmpRes[chosen]._atom[ATM_CB]);
    CurConf._res[Position]._SC = tmpRes[chosen]._SC;
  }
  
  CurConf.SinglecalCenter(CurConf._res[Position], 1);
  if (CurConf._res[Position]._type == -1) {
    cout << "grow_one atom type error " << Position << endl;
    exit(0);
  }

  return true;
}





// end to end distance in loop modeling
double SMC::eteDis(Residue& tmpRes, Residue& endRes) {
  return tmpRes._atom[ATM_CA].dis(endRes._atom[ATM_C]);
}

// Compute group for distance
int SMC::DisGroup(double dis) {
  return int(dis*(1/DistanceBy));          //solve the rounding problem
}

// Compute distance for group
double SMC::GroupDis(int group) {
  return group * DistanceBy;
}



void SMC::fragdis(string fname, int label) {
  // Populate etedP collection of distance lengths
  // Conditional distribution
  ifstream input (fname.c_str());
  int Countset = 0, FraglenCount = -1;
  while(!input.eof()) {
    // Read line
    string line;
    getline(input, line);
    if(line[0] == '#') continue;
    // Tokenize by tabs
    vector<string> token;
    split(line, " ", token);
    if (token.size() == 4)
    {
       FraglenCount++;
       Countset = 0;
       minDistcon[label][FraglenCount] = atof(token[0].c_str());
       DistconBy[label][FraglenCount] = atof(token[1].c_str());
       minDistdel[label][FraglenCount] = atof(token[2].c_str());
       DistdelBy[label][FraglenCount] = atof(token[3].c_str());
    }
    else if (token.size() == 32) 
    {      
      for(int i=0; i < 32; i++)
      {
        etedCon[label][FraglenCount][Countset][i] = atof(token[i].c_str());
      }
      Countset++;

    }
  }
}


void SMC::geometryinfo(string fname) {
  // Populate etedP collection of distance lengths
  // Conditional distribution
  ifstream input (fname.c_str());
  int Countset = 0, FraglenCount = -1;
  while(!input.eof()) {
    // Read line
    string line;
    getline(input, line);
    if(line[0] == '#') continue;

    // Tokenize by tabs
    vector<string> token;
    split(line, " ", token);
    if (token.size() == 1)
    {
       FraglenCount++;
       Countset = 0;
       EEdisBy[FraglenCount] = atof(token[0].c_str());
    }
    else if (token.size() == 37) 
    {      
      for(int i=0; i < 37; i++)
      {
        MtorsionEEdis[FraglenCount][Countset][i] = atof(token[i].c_str());
      }
      Countset++;

    }
  }
}

bool SMC::geometryProb(int length, double EEdis, double Mtorsion) {
     int FragCount = length - 4;
     double distby = EEdisBy[FragCount];
     double x = EEdis - minEEdis;
     int lower = int(x/distby + 0.00001);
     double x1 = lower*distby;
     double x2 = x1 + distby;
     int lowerY = int(Mtorsion/MtorsionBy +0.00001);
     double y1 = lowerY*MtorsionBy;
     double y2 = y1 + MtorsionBy;
     double pdf[37];
     double maxpdf = 0;
     double mypdf = 0;
     for(int i=0; i < 37; i++)
     {
       pdf[i] = ((x2-x)*MtorsionEEdis[FragCount][lower][i] + (x-x1)*MtorsionEEdis[FragCount][lower+1][i])*(1/distby);
       if(pdf[i] > maxpdf)
       maxpdf = pdf[i];
     }
     mypdf = ((y2 - Mtorsion)*pdf[lowerY] + (Mtorsion - y1)*pdf[lowerY+1])*(1/MtorsionBy);
     double myprob = mypdf/maxpdf;
     if(myprob > 1) myprob = 1;
     if(myprob < 0.5)
     myprob = myprob/10;
     double randr = frand(0,1);
     
     if(randr < myprob)
     return true;
     else
     return false;

}

void SMC::simpBBT_Init(string fname)
{
     ifstream input (fname.c_str());
     if(!input.is_open()){
		 cout <<"Error!! Cannot open BBT sampling file!!   "<<fname<<endl;
		 exit(0);
	 }
     int resType, phi_angle,psi_angle;
     while(!input.eof()) {
         // Read line
        string line;
        getline(input, line);
        if(line[0] == '#') continue;
	    if(strlen((char*)line.c_str()) == 0) break;
	 // Tokenize by tabs
        vector<string> token;
        split(line, "\t", token);
		if (token.size() == 4)
		{
	      resType = atoi(token[0].c_str());
          phi_angle = atoi(token[1].c_str())+180;
          psi_angle = atoi(token[2].c_str())+180;
	      Joint_Angle[resType][phi_angle/BBTbinSize][psi_angle/BBTbinSize] = atoi(token[3].c_str());
		}
    }
}
   
string SMC::convertDouble (double value) {
  std::ostringstream o;
  if(!(o<<value))
    return "";
  return o.str();
}

// Pre decide some of the useful parameters, e.g. Close
void SMC::PreProcess()
{
    // Edit end fragment
    if ((Start == 1 && UseBackward) || End == Conf._numRes) {
       cout << "Editing end fragment! " <<endl;
      // If we're editing a start or end fragment, we switch distance states
      // to angle states
      for (int i = 0; i < NumAngleStates.size(); i++)
        NumAngleStates[i] = NumDistanceStates[i] + NumAngleStates[i];
      for (int i = 0; i < NumDistanceStates.size(); i++)
        NumDistanceStates[i] = 0;
      Close = false;
    }
   // Whether to grow backwards
    Backward = (Start == 1 && UseBackward) ? 1 : 0;

}   


