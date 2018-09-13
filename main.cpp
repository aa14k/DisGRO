// main program for folding programs

#include "residue.h"
#include "cal_energy.h"
#include "smc.h"
#include "reprst.h"
#include "potential.h"
#include "util.h"
#include <algorithm>

using namespace std;

void Usage();


int main(int argc, char* argv[]) {

  // Initializations of default parameter values
    int cMODE = 1, NumConf = 100000, scLib = 1, NumRes = 0,
    OutputPDB = 0,  outInfo = 0, start = -1, end =-1, RandSeed = time(NULL), DistanceFallback = 2,
	isContinuous = 1, angType = 2, numSCStates = 5, MAX_NUM_CONF = 5000000, 
	isSampSC = 0, Confkeep =1, CA_Constraint = 1000;
    bool Close = true, UseSingleGrow = true, Kcluster = false, Eval = false, 
	UseBackward = false, OptimizePT = false, noScore = false, Ellipsoid = false, Surface = false,
	Refine = false, TorsionStyle = true;
    string ParameterFile = "data/parameter.txt", DumpSeq = "",
    InitialPdb = "", CommandLog = "command.log", OutDir = ".", ProtName = "", surface_file = "" ,
	OutFile = "", TempSched = "";
    double Radius = 10, dis = 0;
  
  vector <int> sTart;
  vector <int> eNd;
  vector <double> EllipConst2A; 

  int NumState = 5;
  // Initialize number of angle and distance states to 5 and 0
  vector<int> NumAngleStates, NumDistanceStates;
  NumAngleStates.push_back(0);
  NumDistanceStates.push_back(8);

  // protFile stores protein's coordinates;
  // OutFile stores name of output file;
  char resFile, parFile[50], tmpLine[200], selResFile[200];
 

  string tmpStr, dir = "./", protFile = "";
  SSET tmpsSet;
  double disConsWt = 0;
  strcpy(parFile, FILE_ATOMPROP);
  strcpy(selResFile, "no");

  scLib = 1;
  Atom::vdw_adj = 1;
  // Initialize energy modes and weights: we start off weighting energy types
  // equally on their native scales, and including none of them, and an
  // intercept of zero
  for (int i = 0; i < ENERGY_MODES; i++)
    PF::cal[i] = false;
  PF::cal[EM_LOODIS] = true;

  // There are no arguments, so something went awry
  if (argc == 1)
    Usage();

  // Get the mode
  int ii = 1;
  while (ii < argc) {
    if (strcmp(argv[ii], "-mode") == 0) {
      // mode
      // sequential Monte Carlo for both backbone and side chains
      if(strcmp(argv[ii+1], "smc") == 0) cMODE = 1;
      // calculating backbone and side chain torsion angles
      else if(strcmp(argv[ii+1], "torsion") == 0) cMODE = 2;
      else {
	  cout << "wrong mode type!\n";
	  Usage();
      }
      ii += 2;
    }
    else if (strcmp(argv[ii], "-n") == 0) {
      // sample size
      NumConf = atoi(argv[ii+1]);
      if (NumConf > MAX_NUM_CONF) {
        cout << "The maximum sample size is " << MAX_NUM_CONF << "." << endl;
        cout << "Please decrease the input sample size (number of confs)"
	     << "or modify the program." << endl;
        exit(0);
      }
      ii += 2;
    }
    else if (strcmp(argv[ii], "-f") == 0) {
      // protein coordinate file
      protFile = argv[ii+1];
      if (ProtName == "")
	ProtName = File2ProtName(protFile);
      ii += 2;
    }
    else if (strcmp(argv[ii], "-OutFile") == 0) {
      // output file name
      OutFile = argv[ii+1];
      ii += 2;
    }
    else if (strcmp(argv[ii], "-Surface") == 0) {
      // surface list file, ended by ".surface"
      surface_file = argv[ii+1];
      Surface = true;
      ii += 2;
    }
    else if (strcmp(argv[ii], "-nas") == 0) {
      // number of angle states

      vector<string> NumAngleStatesStr;
      int ntoken = split(argv[ii+1], ",", NumAngleStatesStr);
      NumAngleStates.resize(ntoken);
      for (int j = 0; j < ntoken; j++)
	  NumAngleStates[j] = atoi(NumAngleStatesStr[j].c_str());
      ii += 2;

    } else if (strcmp(argv[ii], "-nds") == 0) {
      // number of distance states

      vector<string> NumDistanceStatesStr;
      int ntoken = split(argv[ii+1], ",", NumDistanceStatesStr);
      NumDistanceStates.resize(ntoken);
      for (int j = 0; j < ntoken; j++)
	NumDistanceStates[j] = atoi(NumDistanceStatesStr[j].c_str());
      ii += 2;

    }
    else if (strcmp(argv[ii],"-t") == 0) {
      // temperature
      PF::T = atof(argv[ii+1]);
      ii += 2;
    }
    else if (strcmp(argv[ii], "-va") == 0) { // vdw radius adjustment
      Atom::vdw_adj=atof(argv[ii+1]);
      ii += 2;
    }
    else if(strcmp(argv[ii],"-useSinGr") == 0){
      //use single chain to grow, don't copy the primary conformations
      UseSingleGrow = true;
      ii++;
    }
    else if(strcmp(argv[ii], "-dumpseq") == 0) {
      DumpSeq = argv[ii+1];
      ii += 2;
    }
    else if (strcmp(argv[ii], "-close") == 0) {
      // use analytic closure in SMC fragment regrowth
      Close = bool(atoi(argv[ii + 1]));
      ii += 2;
    }
    else if (strcmp(argv[ii], "-ProtName") == 0) {
      // use analytic closure in SMC fragment regrowth
      ProtName = argv[ii + 1];
      ii += 2;
    }
    else if (strcmp(argv[ii], "-nscc") == 0) {
      // number of side chain states in side chain sampling
      numSCStates = atoi(argv[ii+1]);
      isSampSC = 1;
      ii += 2;
    }
    else if (strcmp(argv[ii], "-e") == 0) {
       if (strcmp(argv[ii+1], "loodis") == 0) 
      	PF::cal[EM_LOODIS] = true;
       else {
        cout << "Unrecognized energy term: " << argv[ii+1] << " !" << endl;
      }
      ii += 2;
    }
    else if (strcmp(argv[ii], "-Count") ==0) {
      PF::COUNT = true;
      ii++;
    }
    else if (strcmp(argv[ii], "-ellip") ==0) {
      Ellipsoid = true;
      ii++;
    }
    else if (strcmp(argv[ii], "-refine") == 0) {
      Refine = true;
      CA_Constraint = atoi(argv[ii+1]);
      ii += 2;
    }
    else if (strcmp(argv[ii], "-dir") == 0) {
      dir = argv[ii+1];
      ii += 2;
    }
    else if (strcmp(argv[ii], "-OutDir") == 0) {
      OutDir = argv[ii+1];
      ii += 2;
    }
    else if (strcmp(argv[ii], "-selres") == 0) {
      // files for selected residues, other residues will b
      strcpy(selResFile, argv[ii+1]);
      ii += 2;
    }
    else if (strcmp(argv[ii], "-outinfo") == 0) {
      outInfo = atoi(argv[ii+1]);
      ii += 2;
    }
    else if (strcmp(argv[ii], "-angt") == 0) {
      // torsion angle representation type
      angType = atoi(argv[ii+1]);
      ii += 2;
    }
    else if (strcmp(argv[ii], "-noscore") == 0) {
	  // Do not store conformations
      noScore = true;
      ii++;
     }	     
    else if (strcmp(argv[ii], "-kcluster") == 0) {
	  // K-mean clustering
      Kcluster = true;
      ii++;
    }
    else if(strcmp(argv[ii], "-ts") == 0) {
      TorsionStyle = true;
      ii++;
    }
    else if (strcmp(argv[ii], "-eval") == 0) {
      Eval = true;
      ii++;
    }
    else if (strcmp(argv[ii], "-confkeep") == 0) {
      Confkeep = atoi(argv[ii+1]);
      ii += 2;
    }
    else if (strcmp(argv[ii], "-start") == 0) {
      start = atoi(argv[ii+1]);
      sTart.push_back(start);
      ii += 2;
    }
    else if(strcmp(argv[ii],"-end")==0) {
      end = atoi(argv[ii+1]);
      eNd.push_back(end);
      ii += 2;
    } 
    else if (strcmp(argv[ii], "-pdbout") == 0) {
	  //Number of conformations need to be outputted
      OutputPDB = atoi(argv[ii+1]);
      ii += 2;
    } else {
      string tmpStr;
      tmpStr = argv[ii];
      if (tmpStr.substr(tmpStr.length() - 4, 4) == ".pdb") {
	  protFile = tmpStr;
  	  ii++;
      }
      else {
	cout << "Unrecognized command: " << argv[ii] << endl;
        Usage();
      }
    }
  }

  // output command line
  if (CommandLog != "") {
    ofstream comOut(CommandLog.c_str(), ios::out | ios::app);
    for (int i = 0; i < argc; i++)
      comOut << argv[i] << " ";
    comOut << endl;
    comOut.close();
  }
  // seed random numbers
  srand(RandSeed);

  // initialize class variables
  Atom::InitPar(parFile, dir);
  Residue::InitMap();
  Residue::InitPar(parFile, dir, outInfo);
  PF::InitPar(ParameterFile, dir, outInfo);
  if (PF::cal[EM_LOODIS] == true) {
    string LOOPTfile = "data/LOODIS_ed4_8_V3.txt";
    PF::initLOODIS(LOOPTfile);
  }
  // If all we want to do is dump the sequence, do it now and exit
  if (DumpSeq != "") {
    Structure tmpStruct(MAX_NUM_RES);
    tmpStruct.readPdb((char*)protFile.c_str(), tmpsSet, 0);
    tmpStruct.WriteSequence(DumpSeq);
    exit(0);
  }
  
  
  // This is the core code, launching the appropriate algorithm
  if (cMODE == 1) {
    // We're in either fold or grow modes, so we first do a common set of
    // things necessary to initialize the SMC object
    // initialize side chain angles
    if (isSampSC) {
      string file_sct = FILE_SCTORSION2;
      SCR::InitSCAng(angType, (char*)file_sct.c_str(), dir);
    }

    Structure tmpStruct(MAX_NUM_RES);
    int looplen =0;
    int total_looplen = 0;
    // Read the initial structure
	tmpStruct.readPdb((char*)protFile.c_str(), tmpsSet, outInfo);  // if the protFile has been regulated to start from 1, use regulated start number, otherwise, can use the index stored in the original PDB, e.g. 1ARP.pdb we can use the number in Sali_loop_9.txt, 127-135, while we use 1ARP_1_1.pdb, we should use 119-127
    
    // Set up the SMC object
    SMC smc(tmpStruct._numRes);
    if (!UseSingleGrow)
      smc = SMC(NumConf, tmpStruct._numRes);

// Set start and ending positions for fragment to be regrown (only for smc)
   if (start == -1)
   {
      cout<<"The start position is required for running DiSGro!!"<<endl;
      cout<<"Example:      ./disgro -mode smc -f 1ctqa.pdb -n 5000 -nds 32 -start 26 -end 37 -eval -confkeep 1000 -ellip -nscc 20 -pdbout 100"<<endl;
      exit(0);
   }
   else
    smc.Start = start;

   if (end == -1)
   {
      cout<<"The end position is required for running DiSGro!!"<<endl;                                                         cout<<"Example:      ./disgro -mode smc -f 1ctqa.pdb -n 5000 -nds 32 -start 26 -end 37 -eval -confkeep 1000 -ellip -nscc 20 -pdbout 100"<<endl;
      exit(0);
   }
   else
    smc.End = end;
  
    // Override default protein name (from file) if given explcitly
    if (ProtName != "")
      tmpStruct._ProtName = ProtName;

    smc.ProtName          = ProtName;
    smc.RandSeed          = RandSeed;
    smc.NumConf           = NumConf;
    smc.disConsWt         = disConsWt;
    smc.Dir               = OutDir;
    smc.AngType           = angType;
    smc.DistanceBy        = 0.1;
    smc.DistanceFallback  = DistanceFallback;
    smc.Close             = Close;
    smc.NumAngleStates    = NumAngleStates;
    smc.NumDistanceStates = NumDistanceStates;
    smc.Conf              = tmpStruct;
    smc.numSCStates	  = numSCStates;
    smc._isSampSC	  = isSampSC;
    smc.UseBackward       = UseBackward;
    smc.confkeep	  = Confkeep;
    smc.Surface		  = Surface;
    smc.Refine		  = Refine;
	smc.Eval          = Eval;
    smc.CA_Constraint	  = CA_Constraint;
    smc.noScore		  = noScore;
    smc.TorsionStyle	  = TorsionStyle;
	smc.outputconf       = OutputPDB;
   
    if(smc.confkeep == 1)
       cout << "Warning!!!!  only keep one conformation at last!!!" <<endl;
    // calculate Ellipsoid constant 2a
    if(Ellipsoid)
    {
        for(int j=0;j<sTart.size();j++)
    	 EllipConst2A.push_back(Ellipsoid_Detect (smc.Conf, sTart[j], eNd[j], true));
	    smc.Ellipsoid = Ellipsoid;	
    }
    vector <int> Tmplist;
    if(Ellipsoid || Surface)
    {
       int resCount = 0;
	//test the space after using ellipsoid and surface
       for(int i=1; i<= tmpStruct._numRes; i++)
       {
           if(smc.Conf._res[i]._center.x== 0 && smc.Conf._res[i]._center.y== 0  && smc.Conf._res[i]._center.z== 0 )
           {
               Tmplist.push_back(i);
	           resCount++;
	           continue;
           }
          
           for(int j = 0; j < sTart.size(); j++)
	       {
	          dis = smc.Conf._res[i]._center.dis(smc.Conf._res[sTart[j]]._atom[ATM_CA]) + smc.Conf._res[i]._center.dis(smc.Conf._res[eNd[j]]._atom[ATM_C]);
	          if(dis < EllipConst2A[j] +  Residue::size[smc.Conf._res[i]._type]+ CUB_SIZE)
	          {
   	                Tmplist.push_back(i);
	                resCount++;
		            break;
	          }
	       }
	   
       }

      smc.Reslist = new int [resCount];
      for(int i=0; i < resCount; i++)
	  smc.Reslist[i] = Tmplist[i];

      smc.List_size = resCount; 
	  cout << "Using REDCELL!! Residue number for calculation is reduced from  "<<tmpStruct._numRes <<" to "<<resCount <<endl;
    }


   // K-mean clustering
    if(Kcluster == true)
      smc.kmcluster = Kcluster;
    // Copy conformation to each element of ConfArr
    if (UseSingleGrow)
      smc._conf = smc.Conf;
    else
      for (int j = 0; j < NumConf; j++)
	smc.ConfArr[j] = smc.Conf;
					     
     //bbt torsion anlge for joint prob sampling 
    smc.simpBBT_Init("data/BBT_phi_psi_pair_NEW.txt");
    if (NumDistanceStates[0] > 0) {
        smc.NDStates_BK = NumDistanceStates[0];
      // Load fragment end-to-end distance distributions
	// Initialize conditional distributions for C and N selection
	smc.fragdis(FILE_FRAG_N_C,  0);
	smc.fragdis(FILE_FRAG_C_CA, 1);
	smc.geometryinfo(FILE_LOOPGEO);
    }
    vector<Structure>Topconflist;
      int FragLength = 0;
      // We're in grow mode
      string topconfpdb = "";
      vector<string> tmpVec;
      split(ParameterFile,"/",tmpVec);
      string parafile_num = tmpVec[1];
      cout <<"Protein Name:        "<< smc.ProtName.substr(0,4) <<endl;
        // Set number of states, populating out to residue length if necessary
      FragLength = smc.End - smc.Start;
      cout << "Start Residue  " << smc.Start <<"   :   End Residue  "<<smc.End <<endl;
      smc.NumAngleStates =
		ExpandNumStates(NumAngleStates, FragLength, "angle");
      smc.NumDistanceStates =
		ExpandNumStates(NumDistanceStates, FragLength, "distance");
      // Launch SMC
      smc.Wholeproc(Topconflist);
      if(OutputPDB > 0)
	  {
		cout << "Output generated conformations to pdb_output."<<endl;
		int ConfNum = (OutputPDB > Topconflist.size())?Topconflist.size():OutputPDB;
        Structure copyconf;
		if(OutputPDB == 1)
		 copyconf = Topconflist[0];
		else
		 copyconf = smc.Conf;
		for(int i=0;i< ConfNum;i++)
        {
           topconfpdb = "topconf.pdb";
	       topconfpdb = "pdb_output/"+ProtName.substr(0,4) + "_" + itoa(smc.Start) +"_" + itoa(smc.End) + "_" + itoa(i+1) + topconfpdb;	    
           if(OutputPDB > 1)
		   { 
			  for(int k = smc.Start; k < smc.End + 1; k++)
				copyconf._res[k] = Topconflist[i]._res[k - smc.Start + 1];
		   }
		   if(isSampSC)
		   copyconf.writePdb(topconfpdb.c_str(),1,copyconf._numRes,0);
		   else
		   copyconf.writePdb(topconfpdb.c_str(),1,copyconf._numRes, smc.Start, smc.End, 0);
        }
	 }
  }
  else if (cMODE == 2) {
    // calculating backbone and side chain torsion angles

    Structure Conf(MAX_NUM_RES);
    Conf.readPdb((char*)protFile.c_str(),tmpsSet,outInfo);
    OutFile = Conf._ProtName+".out";
    cal_dha(Conf,1,Conf._numRes);
    cout << "OUTFILE		"<<OutFile <<endl;
    // side chain torsion angles are calculated in outputAngles function, 
    // however, they are not stored in _scChi[] array in Conf._res[]._atom[].
   Conf.outputAngles((char*)OutFile.c_str(), 1, Conf._numRes,0);
  }
}

void Usage() {
  cout << "disgro: program for simulation of protein structures using multi-dimension Monte Carlo method." << endl;
  cout << "options:\n";
  cout << "\t-mode:\t mode of the computation, currently modes: smc, torsion" << endl;
  cout << "\t\tsmc: loop modeling using sequential chain-groth Monte Carlo." << endl;
  cout << "\t-f:\t input coordinate file of the protein" << endl;
  cout << "\t-n:\t number of conformations" << endl;
  cout << "\t-nds:\t number of sampling states" << endl;
  cout << "\t-start:\t The starting position for growth mode. Required for running Disgro." << endl;
  cout << "\t-end:\t The end position for growth mode. Required for running Disgro." << endl;
  cout << "\t-pdbout: \t the number of the output conformations." << endl;
  cout << "\t-ellip:\t turn on or turn off REsidue-residue Distance Cutoff and ELLipsoid criterion (Redcell)." << endl;
  cout << "\t-nscc:\t number of side chain states" << endl;
  cout << "\t" << endl;
  cout << "example:\t ./disgro -mode smc -f 1ctqa.pdb -n 5000 -nds 32 -start 26 -end 37 -eval -confkeep 1000 -ellip -nscc 20 -pdbout 100" << endl;
  cout << "\t \n";
  cout << "\t \n";
  exit(0);
}

void readSelRes(char* selResFile, SSET& SelRes) {
  ifstream inFile;
  char tmpLine[1000];
  string tmpStr;

  if(strcmp(selResFile,"no")!=0){
    cout << "Reading selected residue file "<<selResFile<<endl;
    inFile.open(selResFile,ios::in); 
    if(!inFile.is_open()){
      cout << "resFile cannot be opened!" << endl;
    }
    else {
      while(!inFile.eof()){
        inFile.getline(tmpLine,1000);
        tmpStr=tmpLine;
        if(tmpStr.length()<2) break;
        SelRes.insert(tmpStr);
      }    
      inFile.close();
    }
  }
}

