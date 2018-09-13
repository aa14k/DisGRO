//// cal_energy.cpp
// calculate the potential energy of a given structure

#include "cal_energy.h"

// this function calculate energy for a fixed structure
void calE_list(Structure& Conf, int* List, int Start, int End, int List_size, bool type) {
  double dis;
  double enArr[ENERGY_TYPES];


  for (int i = 0; i < ENERGY_TYPES; i++)
    enArr[i] = 0;

  double e_hb_bb = 0;
  double e_hb_sc = 0;
  double e_hav = 0;
  double e_hb = 0;
  // the weight for different environment of hbond, same for bb and sc hbond
  double env_wt[3] = {0.18, 0.28, 0.91};
  // Reset energy terms to zero
  Conf._energy = 0;
  for (int i = 0; i < ENERGY_TYPES; i++)
    Conf._enArr[i] = 0;
  // Some terms need the secondary structure info, so do it first
  // Note that we have to step one outside the fragment range to catch all
  // residues whose secondary structures might have changed
  int Start2 = (Start == 1) ? 1 : (Start - 1);
  int   End2 = (End == Conf._numRes) ? End : (End + 1);
  assign_ss(Conf, 1, Conf._numRes, 0);
  // Check if we need to populate the secondary structure
  if (PF::cal[EM_SS] || PF::cal[EM_CON])
    if (Conf._ssPred.size() == 0)
      Conf.PredictSS();
  // All atom-atom interaction terms
  if(type == true)
  if (PF::cal[EM_VDW] || PF::cal[EM_SOL] || PF::cal[EM_HALP] ||
      PF::cal[EM_CNT])
  {  
    vdw_sol_list(Conf, enArr[E_VDWA], enArr[E_VDWR], enArr[E_SOL],
	    Start, End, 0, enArr[E_HALP], enArr[E_VAA], enArr[E_CNT], List, List_size, false); 
    enArr[E_VDWA] = enArr[E_VDWA];
    enArr[E_VDWR] = enArr[E_VDWR];
    enArr[E_VAA] = enArr[E_VAA];
  }
   // loop distance energy
   if(PF::cal[EM_LOODIS])
   {
    enArr[E_LOODIS] = loodis_e(Conf, Start, End, type);
   }
  // Join the energy types with linear weights
  for (int i = 0; i < ENERGY_TYPES; i++) {
    Conf._enArr[i] =  enArr[i]; // first dont consider the weight
    Conf._energy += Conf._enArr[i];
  }
}


void calE(Structure& Conf, int Start, int End, bool type) {
  double dis;
  double enArr[ENERGY_TYPES];

  for (int i = 0; i < ENERGY_TYPES; i++)
    enArr[i] = 0;

  double e_hb_bb = 0;
  double e_hb_sc = 0;
  double e_hav = 0;
  double e_hb = 0;
  // the weight for different environment of hbond, same for bb and sc hbond
  double env_wt[3] = {0.18, 0.28, 0.91};

  // Reset energy terms to zero
  Conf._energy = 0;
  for (int i = 0; i < ENERGY_TYPES; i++)
    Conf._enArr[i] = 0;

  // Some terms need the secondary structure info, so do it first
  // Note that we have to step one outside the fragment range to catch all
  // residues whose secondary structures might have changed
  int Start2 = (Start == 1) ? 1 : (Start - 1);
  int   End2 = (End == Conf._numRes) ? End : (End + 1);
  assign_ss(Conf, 1, Conf._numRes, 0);


  // All atom-atom interaction terms
  if (PF::cal[EM_VDW] || PF::cal[EM_SOL] || PF::cal[EM_HALP] ||
      PF::cal[EM_CNT])
  {  
    vdw_sol(Conf, enArr[E_VDWA], enArr[E_VDWR], enArr[E_SOL],
	    Start, End, 0, enArr[E_HALP], enArr[E_VAA], enArr[E_CNT], false);
    
    enArr[E_VDWA] = enArr[E_VDWA];
    enArr[E_VDWR] = enArr[E_VDWR];
    enArr[E_VAA] = enArr[E_VAA];
  }
   // loop distance energy
   if(PF::cal[EM_LOODIS])
    enArr[E_LOODIS] = loodis_e(Conf, Start, End, type);

  // Join the energy types with linear weights
  for (int i = 0; i < ENERGY_TYPES; i++) {
    Conf._enArr[i] =  enArr[i]; // first dont consider the weight
    Conf._energy += Conf._enArr[i];
  }
}



// template-based energy
// use the distances of CB-CB atoms in template to guide the 
// refinement of target sequence where those CA atoms of aligned
// residues are used
double temp_e(Structure& Conf, int Start, int End){
  IDMAP::iterator idItr;
  IDMAP::iterator idItr2;
  int res1, res2, disbin, key;
  double dis, ttDis, energy=0;
  Atom atom1;
  Atom atom2;
  // not working currently, so just return 0
  return 0;  
  for(idItr=PF::TEMP_CONT.begin();idItr!=PF::TEMP_CONT.end();++idItr){
    res1 = (int)(idItr->first/1000);
    res2 = (int)(idItr->first%1000);
    dis = idItr->second;
    disbin = (int)((dis-2)/2);
    if(disbin<0) {
      cout<<"temp_e disbin error "<<disbin<<endl;
      disbin = 0;
    }
    else if(disbin>2) {
      cout<<"temp_e disbin error "<<disbin<<endl;
      disbin =2;
    }
    if((res1 >= Start && res1 <= End) || (res2 >= Start && res2 <= End)){
      if(Conf._res[res1]._type == GLY)
	atom1 = Conf._res[res1]._atom[ATM_CA];
      else atom1 = Conf._res[res1]._atom[ATM_CB];
      if(Conf._res[res2]._type == GLY)
	atom2 = Conf._res[res2]._atom[ATM_CA];
      else atom2 = Conf._res[res2]._atom[ATM_CB];
      ttDis = abs(dis - atom1.dis(atom2));
      if(Conf._res[res1]._type <= Conf._res[res2]._type)
	key = disbin*10000 + Conf._res[res1]._type * 100 + Conf._res[res2]._type;
      else key = disbin*10000 + Conf._res[res2]._type * 100 + Conf._res[res1]._type;
      if((idItr2 = PF::TEMP.find(key)) != PF::TEMP.end()){
	energy += idItr2->second;
      }
      else {
	cout<<"TEMP index error "<<key<<endl;
      }
    }
  }
  return energy;
}

// function for calculating energy of one residue with the rest of the protein
// when type==0, it calculate energy of this residue with all atoms of
// all residues between start and end
// when type==1, it is for fragment closure of fress, calculating the
// energy of this residue with residues before the fragment
// when type==2, it is for fragment closure of fress, calculating the
// energy of this residue with residues after the fragment
// this function cannot be used to grow the last two residues of the
// fragment, which is done in two_res_en()
// -------------------------
// Note: not all energy terms are used in growing backbone atoms of fragments.
// It may not be that useful to use all the energy terms. A fragment
// needs to be re-evaluated after backbones have been fully grown and
// side chains placed with all the energy terms.  Here VDW terms are
// important since sterics/collision has to be considered during
// backbone growth.
// -------------------------
void one_res_en(Structure& conf, Residue& res, int start, int end, int Start, int End,
		double* energy, int type) {
  int i,j,k,position, atomNum, disInd;
  double dis,e_tmp, phi, psi, disquare;
  // backbone torsion potential

    // get the position of res to check for adjacent residues
    position = res._posn;
    for (i = start; i <= end; i++) {
      if(conf._res[i]._center.x == 0 && conf._res[i]._center.y == 0 && conf._res[i]._center.z == 0) continue;
      if((type == 0 && res._center.dis(conf._res[i]._center) < Residue::size[res._type] + Residue::size[conf._res[i]._type] + CUB_SIZE) ||
	 (type != 0 && res._bbc.dis(conf._res[i]._center) < Residue::bb_size + Residue::size[conf._res[i]._type] + CUB_SIZE)) {
	  
          if(i >= Start && i <= End)
          atomNum =6;
	  else
	  atomNum = conf._res[i]._numAtom;
	for (j = 0; j < NUM_BB_ATOM; j++) {
	  // skip undefined and H atoms
	  if (res._atom[j]._type == UNDEF || res._atom[j]._type >= 22)
	    continue;
          if (res._atom[j].x== 0 && res._atom[j].y== 0 && res._atom[j].z== 0) continue;
	  //if(type != 0 && j >= NUM_BB_ATOM ) continue;  // when sampling only backbone atoms, skip side chain atoms if there are any
	  for(k = 0; k < atomNum; ++k) {
	    if(conf._res[i]._atom[k]._type==UNDEF || conf._res[i]._atom[k]._type >= 22) continue;
            if(conf._res[i]._atom[k].x== 0 && conf._res[i]._atom[k].y== 0 && conf._res[i]._atom[k].z== 0) {
	    continue;
	    }
	    // when calculating energy for the chain before the fragment
	    // do not calculate the energy of C and O of this residue with C and O of residue at end position
	   
	    if (type == 1) {
	      // adjacent residues, do not calculate VDW for atoms
	      // separated by three bonds or less
	      if (position - i == 1) {
	       if (k == ATM_C && (j == ATM_C || j == ATM_CB))
	        continue;
	      }
	    }
	   
            if(type == 3)
	    if (abs(position - i) == 1) {

	      // adjacent residues, do not calculate VDW for atoms
	      // separated by three bonds or less
	      if (i - position == 1) {
		// most of the time it is j - 1 == 1
		if (j == ATM_N && k == ATM_N)
		  continue;
		if (j == ATM_CA && (k == ATM_N || k == ATM_CA))
		  continue;
		if (j == ATM_C &&
		    (k == ATM_N || k == ATM_CA || k == ATM_C || k == ATM_CB))
		  continue;
		if (j == ATM_O && (k == ATM_N || k == ATM_CA))
		  continue;
		if (j == ATM_CB && k == ATM_N)
		  continue;
	      } else if (position - i == 1) {
		// only when i == Start and j == Start - 1
		if (k == ATM_N && j == ATM_N)
		  continue;
		if (k == ATM_CA && (j == ATM_N || j == ATM_CA))
		  continue;
		if (k == ATM_C &&  (j == ATM_N || j == ATM_CA || j == ATM_C || j==ATM_CB))
		  continue;
		if (k == ATM_O && (j == ATM_N || j == ATM_CA))
		  continue;
		if (k == ATM_CB && j == ATM_N)
		  continue;
	      }
	    }
	   if(position > Start)
	   if (type == 1 && i == end && (k == ATM_C || k == ATM_O) &&
		(j == ATM_C || j == ATM_O))
	      continue;
	      disquare = (res._atom[j].x-conf._res[i]._atom[k].x)*(res._atom[j].x-conf._res[i]._atom[k].x)+(res._atom[j].y-conf._res[i]._atom[k].y)*(res._atom[j].y-conf._res[i]._atom[k].y)+
	      (res._atom[j].z-conf._res[i]._atom[k].z)*(res._atom[j].z-conf._res[i]._atom[k].z);

	    if (disquare <= PF_DIS_CUT_SQUARE)
	    { 
	          dis = sqrt(disquare);
              disInd = (int)(dis/H_INLO);
	 	    energy[E_LOODIS] += PF::LOODIS[res._atom[j]._type-1][conf._res[i]._atom[k]._type-1][disInd];
	    }
	  }
	}
      }
    }

}

//type 0 is used for loop closure calculation , type 1 is used for normal chain-growth
void one_res_en_list(Structure& conf, Residue& res, int Start, int End, double* energy, int type, int* List, int List_size) {
  int p,i,j,k,position, atomNum;
  int disInd;
  double dis,e_tmp, phi, psi, disquare;
  // atom-atom distance potential
    // get the position of res to check for adjacent residues
    position = res._posn;
    for (p = 0; p < List_size; p++) {
      i = List[p];
      if(conf._res[i]._center.x == 0)
      {
	  if( i != position)
          continue; 
      }
      if(position != End)
      {	
       if( i == position)
       {
        if(type)
	    {
	     if (res._atom[1]._type == UNDEF || conf._res[i]._atom[0]._type == UNDEF || res._atom[1].x== 0.000 || conf._res[i]._atom[0].x== 0.000)
	     continue;
	     else
	     {
	       disquare = (res._atom[1].x-conf._res[i]._atom[0].x)*(res._atom[1].x-conf._res[i]._atom[0].x)+
			 (res._atom[1].y-conf._res[i]._atom[0].y)*(res._atom[1].y-conf._res[i]._atom[0].y)+
	                 (res._atom[1].z-conf._res[i]._atom[0].z)*(res._atom[1].z-conf._res[i]._atom[0].z);
	      if (disquare <= PF_DIS_CUT_SQUARE)
	      { 
	          dis = sqrt(disquare);
             	  disInd = (int)(dis/H_INLO);
	 	  energy[E_LOODIS] += PF::LOODIS[res._atom[1]._type-1][conf._res[i]._atom[0]._type-1][disInd];
	      }
         }
	    }
       }
       else if( i > position && i < End)  continue;
       else
       {
     // if((type == 0 && res._center.dis(conf._res[i]._center) < Residue::size[res._type] + Residue::size[conf._res[i]._type] + CUB_SIZE) ||
	    if (res._bbc.dis(conf._res[i]._center) < Residue::bb_size + Residue::size[conf._res[i]._type] + CUB_SIZE) {
	    for (j = 0; j < NUM_BB_ATOM; j++) {
	  // skip undefined and H atoms
	    if (res._atom[j]._type == UNDEF || res._atom[j]._type >= 22)
	    continue;
        if (res._atom[j].x== 0.000) continue;
	    for(k = 0; k < conf._res[i]._numAtom; ++k) {
	    
	    if(i == End)
           if(type)
	    if(k != ATM_C && k != ATM_CA)
	    continue;    
	    if(conf._res[i]._atom[k]._type==UNDEF || conf._res[i]._atom[k]._type >= 22) continue;
            if(conf._res[i]._atom[k].x== 0.000) continue;
	    if (position - i == 1) {
	      if(k == ATM_C && (j == ATM_C || j == ATM_CB))
	       continue;
	      if(type)
              {
	        if(position > Start) 
                if( (k == ATM_C || k == ATM_O) && (j == ATM_C || j == ATM_O || j == ATM_CB) )
		   continue;
              }
	      else
	      {
                if( (k <= 5 && j == ATM_N) || ( (k == ATM_CA || k== ATM_C || k == ATM_O) && j == ATM_CA) )
		 continue;
	      }
	    }
	      disquare = (res._atom[j].x-conf._res[i]._atom[k].x)*(res._atom[j].x-conf._res[i]._atom[k].x)+
			 (res._atom[j].y-conf._res[i]._atom[k].y)*(res._atom[j].y-conf._res[i]._atom[k].y)+
	                 (res._atom[j].z-conf._res[i]._atom[k].z)*(res._atom[j].z-conf._res[i]._atom[k].z);
	    if (disquare <= PF_DIS_CUT_SQUARE)
	    { 
	          dis = sqrt(disquare);
             	  disInd = (int)(dis/H_INLO);
	 	  energy[E_LOODIS] += PF::LOODIS[res._atom[j]._type-1][conf._res[i]._atom[k]._type-1][disInd];
	    }
	  }
	 }
    }
    }
   }
   else
   {
      if( i <= position && i >= position - 2) continue;
      if(res._bbc.dis(conf._res[i]._center) < Residue::bb_size + Residue::size[conf._res[i]._type] + CUB_SIZE) {
	for (j = 0; j < NUM_BB_ATOM; j++) {
	  // skip undefined and H atoms
	  if ( res._atom[j]._type == UNDEF || res._atom[j]._type >= 22 )
	    continue;
          if(i >= Start && i < End)    //CA and C have been calculated previously
          {
            if(j == 1 || j == 2)
	    continue;
	  }
        if (res._atom[j].x== 0.000) continue;
  	    for(k = 0; k < conf._res[i]._numAtom; ++k) {
	    if(conf._res[i]._atom[k]._type==UNDEF || conf._res[i]._atom[k]._type >= 22) continue;
            if(conf._res[i]._atom[k].x== 0.000) continue;
	    if (position - i == 1) {
	      if( (k == ATM_C && j == ATM_CB) ||  (k <= 5 && j == ATM_N) )
	       continue;
	    }
	    if(i - position == 1) {
	      if( (k == ATM_N && j <= 5) || (k == ATM_CA && (j == ATM_CA ||  j == ATM_C || j == ATM_O) ) || (j == ATM_C && (k == ATM_C || k == ATM_CB)) )
              continue;
	    }
	    disquare = (res._atom[j].x-conf._res[i]._atom[k].x)*(res._atom[j].x-conf._res[i]._atom[k].x)+
		       (res._atom[j].y-conf._res[i]._atom[k].y)*(res._atom[j].y-conf._res[i]._atom[k].y)+
	               (res._atom[j].z-conf._res[i]._atom[k].z)*(res._atom[j].z-conf._res[i]._atom[k].z);
	    if (disquare <= PF_DIS_CUT_SQUARE)
	    { 
	          dis = sqrt(disquare);
             	  disInd = (int)(dis/H_INLO);
	 	  energy[E_LOODIS] += PF::LOODIS[res._atom[j]._type-1][conf._res[i]._atom[k]._type-1][disInd];

	    }
        }
       }
     }
   }
  }
}


// side chain interaction energy for side chain modeling. 
// only side chain atoms are considered
// Note: side chain center for the residue _scc should be calculated using Residue::cal_scc() 
// before calling this function
void one_res_en_sc(Structure& conf, Residue& res, int position, int Start, int End,  int type,
		   double* energy) {
  int i,j,k,key, disInd;
  double dis,e_tmp, disquare;
  if(res._type == 0 || res._type == 5)
   return;
  else
  {
  // atom-atom distance potential
    for(i=1;i< position;++i){
      if( (res._scc.dis(conf._res[i]._center) < (Residue::sc_size[res._type] + Residue::size[conf._res[i]._type] + CUB_SIZE)) || i >= Start) {
	// Here we can use another distance cutoff instead of CC_DIS_CUT which is 12 for whole residues
	// the distance cutoff here for sidechains can be considerablly smaller and residue dependent
	  for(k = 0; k < conf._res[i]._numAtom; ++k){
	    if(k==4 || conf._res[i]._atom[k]._type==UNDEF || conf._res[i]._atom[k]._type >= 22) continue;
            if(conf._res[i]._atom[k]._name == "") continue;	   
	    for(j = NUM_BB_ATOM;j < res._numAtom;++j){  
	    if(res._atom[j]._type==UNDEF || res._atom[j]._type>=22) continue;
	    
	    disquare = (res._atom[j].x-conf._res[i]._atom[k].x)*(res._atom[j].x-conf._res[i]._atom[k].x)+(res._atom[j].y-conf._res[i]._atom[k].y)*(res._atom[j].y-conf._res[i]._atom[k].y)+
	                  (res._atom[j].z-conf._res[i]._atom[k].z)*(res._atom[j].z-conf._res[i]._atom[k].z);

	    if(disquare <= PF_DIS_CUT_SQUARE)
	    {
	     dis = sqrt(disquare);
	     disInd = (int)(dis/H_INLO);
	     energy[E_LOODIS] += PF::LOODIS[res._atom[j]._type-1][conf._res[i]._atom[k]._type-1][disInd];
	    }
	  }
	}
   }
 } 

    for(i=End+1;i<=conf._numRes;++i){
      if(res._scc.dis(conf._res[i]._center) < (Residue::sc_size[res._type] + Residue::size[conf._res[i]._type] + CUB_SIZE)) {
	// Here we can use another distance cutoff instead of CC_DIS_CUT which is 12 for whole residues
	// the distance cutoff here for sidechains can be considerablly smaller and residue dependent
	  for(k = 0; k < conf._res[i]._numAtom; ++k){
	    if(k==4 || conf._res[i]._atom[k]._type==UNDEF || conf._res[i]._atom[k]._type >= 22) continue;
            if(conf._res[i]._atom[k]._name == "") continue;
	    for(j = NUM_BB_ATOM;j < res._numAtom;++j){  
	    if(res._atom[j]._type==UNDEF || res._atom[j]._type>=22) continue;
	    disquare = (res._atom[j].x-conf._res[i]._atom[k].x)*(res._atom[j].x-conf._res[i]._atom[k].x)+(res._atom[j].y-conf._res[i]._atom[k].y)*(res._atom[j].y-conf._res[i]._atom[k].y)+
	                              (res._atom[j].z-conf._res[i]._atom[k].z)*(res._atom[j].z-conf._res[i]._atom[k].z);
	    if(disquare <= PF_DIS_CUT_SQUARE)
	    {
	     dis = sqrt(disquare);
	     disInd = (int)(dis/H_INLO);
	     energy[E_LOODIS] += PF::LOODIS[res._atom[j]._type-1][conf._res[i]._atom[k]._type-1][disInd];
	    }
	  }
	}
   }
  }

    //C atom in End 
    if(position < End)
    {
      for(i= position +1 ;i<= End;++i){
      if(res._scc.dis(conf._res[i]._center) < (Residue::sc_size[res._type] + Residue::size[conf._res[i]._type] + CUB_SIZE)) {
	  for(k = 0; k < NUM_BB_ATOM; ++k){
	    if(k==4 || conf._res[i]._atom[k]._type==UNDEF || conf._res[i]._atom[k]._type >= 22) continue;
            if(conf._res[i]._atom[k]._name == "") continue;
	    for(j = NUM_BB_ATOM;j < res._numAtom;++j){  
	    if(res._atom[j]._type==UNDEF || res._atom[j]._type>=22) continue;
	    disquare = (res._atom[j].x-conf._res[i]._atom[k].x)*(res._atom[j].x-conf._res[i]._atom[k].x)+(res._atom[j].y-conf._res[i]._atom[k].y)*(res._atom[j].y-conf._res[i]._atom[k].y)+
	                              (res._atom[j].z-conf._res[i]._atom[k].z)*(res._atom[j].z-conf._res[i]._atom[k].z);
	    if(disquare <= PF_DIS_CUT_SQUARE)
	    {
	     dis = sqrt(disquare);
	     disInd = (int)(dis/H_INLO);
	     energy[E_LOODIS] += PF::LOODIS[res._atom[j]._type-1][conf._res[i]._atom[k]._type-1][disInd];
	    }
	   }
	  }
     }
    }
   }
  }
}


// surface side chain interaction energy for side chain modeling. 
// only side chain atoms are considered
// Note: side chain center for the residue _scc should be calculated using Residue::cal_scc() 
// before calling this function
void one_res_en_sc_list(Structure& conf, Residue& res, int position, int Start, int End,  int type, int* List, int List_size,  double* energy) {
  int i,j,k,key,p,q;
  double dis,e_tmp, disquare;
  if(res._type == 0 || res._type == 5)
   return;
  else
  {
  if(PF::cal[EM_VDW] == 1){
    for(p=0;p < List_size;++p){
      i = List[p];
      if(res._scc.dis(conf._res[i]._center) < (Residue::sc_size[res._type] + Residue::size[conf._res[i]._type] + CUB_SIZE)) {
        for(k = 0; k < conf._res[i]._numAtom; ++k){
	    if(k==4 || conf._res[i]._atom[k]._type==UNDEF) continue;
            if(conf._res[i]._atom[k]._name == "") continue;
	    for(j = NUM_BB_ATOM;j < res._numAtom;++j){  
	      if(res._atom[j]._type==UNDEF || res._atom[j]._type>=22) continue;
	    disquare = (res._atom[j].x-conf._res[i]._atom[k].x)*(res._atom[j].x-conf._res[i]._atom[k].x)+(res._atom[j].y-conf._res[i]._atom[k].y)*(res._atom[j].y-conf._res[i]._atom[k].y)+
	                             (res._atom[j].z-conf._res[i]._atom[k].z)*(res._atom[j].z-conf._res[i]._atom[k].z);
	    if(disquare>PF_DIS_CUT_SQUARE) continue;
	    else
	    dis = sqrt(disquare);
	    aa_int_en(res, conf._res[i], res._atom[j], conf._res[i]._atom[k], dis, energy, 1);
	  }
	}
      }
    }

    if(position > Start)
     for(i= Start; i < position; ++i){
	    if(res._scc.dis(conf._res[i]._center) < (Residue::sc_size[res._type] + Residue::size[conf._res[i]._type] + CUB_SIZE)) {
	      for (k = 0; k < conf._res[i]._numAtom; ++k){
	       if(k==4 || conf._res[i]._atom[k]._type==UNDEF) continue;
	       if(conf._res[i]._atom[k]._name == "") continue;
	       for(j = NUM_BB_ATOM;j < res._numAtom;++j){
	     	 if(res._atom[j]._type==UNDEF || res._atom[j]._type>=22) continue;
		 disquare = (res._atom[j].x-conf._res[i]._atom[k].x)*(res._atom[j].x-conf._res[i]._atom[k].x)+
		 (res._atom[j].y-conf._res[i]._atom[k].y)*(res._atom[j].y-conf._res[i]._atom[k].y)+
		 (res._atom[j].z-conf._res[i]._atom[k].z)*(res._atom[j].z-conf._res[i]._atom[k].z);
		 if(disquare>PF_DIS_CUT_SQUARE) continue;
		 else
		 dis = sqrt(disquare);	
                 aa_int_en(res, conf._res[i], res._atom[j], conf._res[i]._atom[k], dis, energy, 1);
              }
	     }
	  }
     }
     //C atom in End 
     if(position < End)
     for(j = NUM_BB_ATOM;j < res._numAtom;++j){  
       if(res._atom[j]._type==UNDEF || res._atom[j]._type>=22) continue;
       disquare = (res._atom[j].x-conf._res[End]._atom[2].x)*(res._atom[j].x-conf._res[End]._atom[2].x)+
                  (res._atom[j].y-conf._res[End]._atom[2].y)*(res._atom[j].y-conf._res[End]._atom[2].y)+
		  (res._atom[j].z-conf._res[End]._atom[2].z)*(res._atom[j].z-conf._res[End]._atom[2].z);
      if(disquare <= PF_DIS_CUT_SQUARE)
      {  
           dis = sqrt(disquare);
           aa_int_en(res, conf._res[End], res._atom[j], conf._res[End]._atom[2], dis, energy, 1);
      }
    }
				      // VDW term within a residue E_VAA
    vdw_within_aa(conf, res, energy[E_VAA]);
  }
 } 
}



// this function is writen for fragment closure using Dill's method
// to add the last two residues, here p1 is position of res1 in the sequence 
// and p2 is position of res2 in the sequence
void two_res_en(Structure& conf, Residue& res1, Residue& res2, int p1, int p2, double* energy)
{
  int i,j,k;
  double dis,e_tmp;
  // energy of res1

  for(i=1;i<=conf._numRes;++i){
    if(i==p1) continue;  
    if(res1._center.dis(conf._res[i]._center)<CC_DIS_CUT){
      for(j=0;j<res1._numAtom;++j){
	if(res1._atom[j]._type==UNDEF || res1._atom[j]._type>=22) continue;
	for(k=0;k<conf._res[i]._numAtom;++k){
	  if(conf._res[i]._atom[k]._type==UNDEF || conf._res[i]._atom[k]._type>=22 ) continue;
	  // for residue p1-1, do not calculate the energy of C and O of res1 with C and O of conf._res[p1-1]
	  if(i==(p1-1) && (k==ATM_C||k==ATM_O) && (j==ATM_C||j==ATM_O)) continue;
	  // for energy between res1 and res2, only calculate the C and O of p1 with the rest of p2 except C and O
	  if(i==p2 && ((j!=ATM_C && j!=ATM_O)||(k==ATM_C||k==ATM_O))) continue;
	  // for energy between res1 and p2+1, only calculate the C and O of res1 with C and O of p2+1
	  if(i==(p2+1) && ((j!=ATM_C && j!=ATM_O) || (k!=ATM_C && k!=ATM_O))) continue;
	  dis = res1._atom[j].dis(conf._res[i]._atom[k]);
	  if(dis>PF_DIS_CUT) continue;
	  aa_int_en(res1, conf._res[i], res1._atom[j],conf._res[i]._atom[k],dis,energy,0);
	}
      }
    }
  }

  // energy of res2

  for(i=1;i<=conf._numRes;++i){
    if(i==p2) continue;
    if(res2._center.dis(conf._res[i]._center)<CC_DIS_CUT){
      for(j=0;j<res2._numAtom;++j){
	if(res2._atom[j]._type==UNDEF || res2._atom[j]._type>=22) continue;
	for(k=0;k<conf._res[i]._numAtom;++k){
	  if(conf._res[i]._atom[k]._type==UNDEF || conf._res[i]._atom[k]._type>=22) continue;
	  // for residue at p1 position in conf, conf._res[p1], not res1
	  // only calculate energy of N, CA of conf._res[p1] and N, CA of res2
	  if(i==p1 && j==ATM_C || j==ATM_O || k==ATM_C || k==ATM_O) continue;
	  // for residue p2+1, only calculate energy of C, O of res2 with atoms of p2+1
	  if(i==(p2+1) && j!=ATM_C && j!=ATM_O) continue; 
	  dis = res2._atom[j].dis(conf._res[i]._atom[k]);
	  if(dis>PF_DIS_CUT) continue;
	  aa_int_en(res2, conf._res[i], res2._atom[j],conf._res[i]._atom[k],dis,energy,0);
	}
      }
    }
  }
}

// calculate the energy of a fragment with multiple residues, where length is the size of the fragment
// this function is for general use not for use in FRESS since fress cut residues in the middle
// for fress use one_res_en() and two_res_en()
void frag_en(Structure& conf, Residue* res, int start, int end, int length, double* energy)
{
  int i,j,k,position,len;
  double dis,e_tmp;
  for(len=0;len<length;++len){
    position = res[len]._posn; // get the position of res to check for adjacent residues
    for(i=start;i<=end;++i){
      if(res[len]._center.dis(conf._res[i]._center)<CC_DIS_CUT){
	for(j=0;j<=res[len]._numAtom;++j){
	  if(res[len]._atom[j]._type==UNDEF) continue;
	  for(k=0;k<=conf._res[i]._numAtom;++k){
	    if(conf._res[i]._atom[k]._type==UNDEF) continue;
	    dis = res[len]._atom[j].dis(conf._res[i]._atom[k]);
	    if(dis>PF_DIS_CUT) continue;
	    aa_int_en(res[len], conf._res[i], res[len]._atom[j],conf._res[i]._atom[k],dis,energy,1);
	  }
	}
      }
    }
  }
}

// atom-atom interaction energy, the arguements are distance, atoms, and the energy 
// type == 0: calculating backbone energy during fragment growth, and SIMPL is used
// type == 1: calculating side chain energy after fragment growth, and SIMPL is not used
void aa_int_en(Residue& res1, Residue& res2, Atom& atom1, Atom& atom2, double dis, double* energy, int type) {
  double r, eps; // constants when two atoms are known
  double e_tmp,e_aap,slope,y_intercept;
  // van der waals interaction energy
  double tmp_vdwa1 =0, tmp_vdwa2 =0, tmp_vdwr1 =0, tmp_vdwr2 =0;
  int bin;

  r = Atom::radius[atom1._type]+Atom::radius[atom2._type];

if (PF::cal[EM_VDW]) {
    eps=sqrt((Atom::welldepth[atom1._type])*(Atom::welldepth[atom2._type]));
    // attraction part
    if((r/dis)<1.12){ 
      energy[E_VDWA] += ((pow(r/dis,12)-2*pow(r/dis,6))*eps);
    }
    
    // repulsion part
    else if((r/dis)<1.33){
      energy[E_VDWR] += ((pow(r/dis,12)-2*pow(r/dis,6))*eps);
   }
    else {
       slope=-12*eps*(33.383)*(1/r);  // 33.383=(1.33^13-1.33^7)
      y_intercept=-1*slope*(r/1.33)+eps*19.565; // 19.565=(1.33^12-2(1.33^6))
      energy[E_VDWR] += (y_intercept + dis*slope);
    }
  }
}

// function for calculating vdw and solvation energy, just for surface atoms

// VDW energy within a residue is calculated in calE_surface 
void vdw_sol_list(Structure& Conf, double& e_atr, double& e_rep, double& e_sol,
	     int Start, int End, int StartAtom, double& e_halp, double& e_vaa,
	     double& e_cnt, int* List, int List_size, bool debug) {
  int i,j,k,l,tp1,tp2,disInd,p,q;
  double dis;
  double r_ij,e_ij,slope,y_intercept,e_tmp; // for vdw calculation
  ILIST::iterator ilItr;
  int resType1,resType2,atomIndex1,atomIndex2, isHB;
  string atomName1,atomName2;
  SIMAP::iterator aiItr;
  IIMAP::iterator iiItr;
  int atomPairInd; // index for atom pairs
  e_atr = 0;
  e_rep = 0;
  e_sol = 0;
  Conf._AATermCnt.clear();
  // calculate VDW energy within a residue
  for (int i = Start; i <= End; i++)
    vdw_within_aa(Conf, Conf._res[i], e_vaa);

  // calculate VDW energy between residues including adjacent residues
  for (i = Start; i <= End; i++) {
	for( p = 0; p < List_size; p++) {
          j = List[p];
      if (!((i >= Start && i <= End) || (j >= Start && j <= End)))
	  continue;
      // Move on unless i < j, avoid repeated calculation
      if (j <= i && (j <= End && j >= Start) )   
	  continue;
      // Move on if i == j
      if (i == j)
	  continue;

      if (Conf._res[i]._center.dis(Conf._res[j]._center) < CC_DIS_CUT) {
         
	for (k = 0; k < Conf._res[i]._numAtom; k++) {
	  // All atoms in the first residue
	  // Skip atoms of unknown type and H's
	  if (Conf._res[i]._atom[k]._type == UNDEF)
	    continue;
	  if (Conf._res[i]._atom[k]._type > 21)
	    continue;

	  for (l = 0; l < Conf._res[j]._numAtom; l++) {
       	     if (debug)
	      cout << endl << i << " " << j << " " << k << " " << l << " ";

	    // Skip atoms of unknown type and H's
	    if (Conf._res[j]._atom[l]._type == UNDEF)
	      continue;
	    if (Conf._res[j]._atom[l]._type > 21)
	      continue;

	    if (abs(i - j) == 1) {
	      // adjacent residues, do not calculate VDW for atoms
	      // separated by three bonds or less
	      if (j - i == 1) {
		// most of the time it is j - 1 == 1
		if (k == ATM_N && l == ATM_N)
		  continue;
		if (k == ATM_CA && (l == ATM_N || l == ATM_CA))
		  continue;
		if (k == ATM_C &&
		    (l == ATM_N || l == ATM_CA || l == ATM_C || l == ATM_CB))
		  continue;
		if (k == ATM_O && (l == ATM_N || l == ATM_CA))
		  continue;
		if (k == ATM_CB && l == ATM_N)
		  continue;
	      } else if (j - i == -1) {
		// only when i == Start and j == Start - 1
		if (l == ATM_N && k == ATM_N)
		  continue;
		if (l == ATM_CA && (k == ATM_N || k == ATM_CA))
		  continue;
		if (l == ATM_C && (k == ATM_N || k == ATM_CA || k == ATM_C))
		  continue;
		if (l == ATM_O && (k == ATM_N || k == ATM_CA))
		  continue;
		if (l == ATM_CB && k == ATM_N)
		  continue;
	      }
	    }
	    // Exclude atoms too far from one another
	    dis = Conf._res[i]._atom[k].dis(Conf._res[j]._atom[l]);
	    if (dis > PF_DIS_CUT)
	      continue;
	    // Get HALP energy
	    if (PF::cal[EM_HALP]) {
	      resType1 = Conf._res[i]._type;
	      resType2 = Conf._res[j]._type;
	      atomName1 = Residue::Name1[resType1] +
		  Conf._res[i]._atom[k]._name;
	      atomName2 = Residue::Name1[resType2] +
		  Conf._res[j]._atom[l]._name;

	      // Sanity check
	      if ((aiItr = Residue::AtomIndexMap.find(atomName1)) !=
		  Residue::AtomIndexMap.end())
		atomIndex1 = aiItr->second;
	      else {
		cout << "cal_halp(): atom type error! " << atomName1 << " "
		     << i << " " << k << endl;
		continue;
		//return;
	      }
	      if ((aiItr = Residue::AtomIndexMap.find(atomName2)) !=
		  Residue::AtomIndexMap.end())
		  atomIndex2 = aiItr->second;
	      else {
		  cout << "cal_halp(): atom type error! " << atomName2 << " "
		     << j << " " << l << endl;
		  continue;
	      }

	      // starting distance is 1.5 A
	      disInd = (int)(dis / H_INTV) - (int)(START_DIS / H_INTV);
	      if (disInd >= 10.0 / H_INTV)
		  disInd = (int)(10.0 / H_INTV) - 1;
	      if (disInd < 0)
		  disInd = 0;
	      atomPairInd = PF::APIL[atomIndex1][atomIndex2];
	      e_halp += PF::HADIP[atomPairInd][disInd];
	      if (Conf._checkEn >= 0) {
		// add this into LLIFMAP _EnProfile
		long long int tmpI;
		tmpI = i * 10000000 + k * 100000 + j * 100 + l;
		Conf._EnProfile[E_HALP].insert(LLIFMAP::value_type(tmpI,PF::HADIP[atomPairInd][disInd]*PF::Parameter["E"][E_HALP]));
	      }
	    }
	    if (PF::cal[EM_VDW]) {
	      // Van der waals interaction, atraction and repulsion part
	      e_ij = sqrt(Atom::welldepth[Conf._res[i]._atom[k]._type] *
			  Atom::welldepth[Conf._res[j]._atom[l]._type]);
		  // check for Hbond
	      isHB = 0;
	      if (isHB == 0)
		  r_ij = Atom::radius[Conf._res[i]._atom[k]._type] +
		  Atom::radius[Conf._res[j]._atom[l]._type];
	      else if (isHB == 1)
		  r_ij = 2.95;
	      else if (isHB == 2)
		  r_ij = 2.95;
	      if (r_ij / dis < 1.12) {
		 // attraction part
		  e_tmp  = (pow(r_ij / dis, 12) - 2 * pow(r_ij / dis, 6)) * e_ij;
		  if (Conf._checkEn >= 0) {
		  // add this into SFMAP _EnProfile
		  long long int tmpI;
		  tmpI = i * 10000000 + k * 100000 + j * 100 + l;
		  Conf._EnProfile[E_VDWA].insert(LLIFMAP::value_type(tmpI, e_tmp * PF::Parameter["E"][E_VDWA]));
		}
		e_atr += e_tmp;
	    } else if (r_ij / dis < 1.33) {
		// repulsion part
		e_tmp  = (pow(r_ij / dis, 12) - 2 * pow(r_ij / dis, 6)) * e_ij;
		if (debug)
		  cout << e_tmp;
		e_rep += e_tmp;
	    } else {
		slope = -12 * e_ij * 33.383 / r_ij;
		y_intercept = -slope * r_ij / 1.33 + e_ij * 19.565;
		e_tmp  = y_intercept + dis * slope;
		if (debug)
		  cout << e_tmp;
		if(Conf._checkEn >= 0) { // add this into SFMAP _EnProfile
		  long long int tmpI;
		  tmpI = i * 10000000 + k * 100000 + j * 100 + l;
		  Conf._EnProfile[E_VDWR].insert(LLIFMAP::value_type(tmpI,e_tmp*PF::Parameter["E"][E_VDWR]));
		}
		e_rep += e_tmp;
	    }
	   }
	  // Solvation
	    if (PF::cal[EM_SOL]) {
	      tp1 = Conf._res[i]._atom[k]._type;
	      tp2 = Conf._res[j]._atom[l]._type;
	      r_ij = Atom::radius[tp1]+Atom::radius[tp2];
	      e_tmp = -(Atom::s_dgfree[tp1] * pow(EXPO, -dis * dis) *
			Atom::s_volume[tp2] / (2 * pow(PI,1.5) *
					       Atom::s_lambda[tp1] * r_ij *
					       r_ij) +
			Atom::s_dgfree[tp2] * pow(EXPO, -dis * dis) *
			Atom::s_volume[tp1] / (2 * pow(PI, 1.5) *
					       Atom::s_lambda[tp2] * r_ij *
					       r_ij));
	      if (Conf._checkEn >= 0) { // add this into SFMAP _EnProfile
		long long int tmpI;
		tmpI = i * 10000000 + k * 100000 + j * 100 + l;
		Conf._EnProfile[E_SOL].insert(LLIFMAP::value_type(tmpI,
								  e_tmp *
								  PF::Parameter["E"][E_SOL]));
	      }
	      e_sol += e_tmp;
	    }
	  }
	}
      }
    }
  }
}


// function for calculating vdw and solvation energy,
// the function can also calculate the energy between a fragment with
// the rest of the protein (including intra-fragment energy) by giving
// Start and End different values than 1 and _numRes
// VDW energy are calculated for all atoms separated by more than
// three chemical bonds
// VDW energy within a residue is calculated in calE 
void vdw_sol(Structure& Conf, double& e_atr, double& e_rep, double& e_sol,
	     int Start, int End, int StartAtom, double& e_halp, double& e_vaa,
	     double& e_cnt, bool debug) {
  int i,j,k,l,tp1,tp2,disInd;
  double dis;
  double r_ij,e_ij,slope,y_intercept,e_tmp; // for vdw calculation
  ILIST::iterator ilItr;
  int resType1,resType2,atomIndex1,atomIndex2, isHB;
  string atomName1,atomName2;
  SIMAP::iterator aiItr;
  IIMAP::iterator iiItr;
  int atomPairInd; // index for atom pairs
  e_atr = 0;
  e_rep = 0;
  e_sol = 0;
  Conf._AATermCnt.clear();
  // calculate VDW energy within a residue
  for (i = Start; i <= End; i++)
    vdw_within_aa(Conf, Conf._res[i], e_vaa);
  // calculate VDW energy between residues including adjacent residues
  for (i = Start; i <= End; i++) {
     for (j = 1; j <= Conf._numRes; j++) {
      // All residues in the structure
      // Move on if i == j
      if (i == j)
	  continue;
      if (Conf._res[i]._center.dis(Conf._res[j]._center) < CC_DIS_CUT) {
        for (k = 0; k < Conf._res[i]._numAtom; k++) {
	  // All atoms in the first residue
	  // Skip atoms of unknown type and H's
	  if (Conf._res[i]._atom[k]._type == UNDEF)
	    continue;
	  if (Conf._res[i]._atom[k]._type > 21)
	    continue;

	  for (l = 0; l < Conf._res[j]._numAtom; l++) {
	    // All atoms in the second residue
	    if (debug)
	      cout << endl << i << " " << j << " " << k << " " << l << " ";

	    // Skip atoms of unknown type and H's
	    if (Conf._res[j]._atom[l]._type == UNDEF)
	      continue;
	    if (Conf._res[j]._atom[l]._type > 21)
	      continue;
	    if (abs(i - j) == 1) {
	      // adjacent residues, do not calculate VDW for atoms
	      // separated by three bonds or less
	      if (j - i == 1) {
		// most of the time it is j - 1 == 1
		if (k == ATM_N && l == ATM_N)
		  continue;
		if (k == ATM_CA && (l == ATM_N || l == ATM_CA))
		  continue;
		if (k == ATM_C &&
		    (l == ATM_N || l == ATM_CA || l == ATM_C || l == ATM_CB))
		  continue;
		if (k == ATM_O && (l == ATM_N || l == ATM_CA))
		  continue;
		if (k == ATM_CB && l == ATM_N)
		  continue;
	      } else if (j - i == -1) {
		// only when i == Start and j == Start - 1
		if (l == ATM_N && k == ATM_N)
		  continue;
		if (l == ATM_CA && (k == ATM_N || k == ATM_CA))
		  continue;
		if (l == ATM_C &&  (k == ATM_N || k == ATM_CA || k == ATM_C || k==ATM_CB))
		  continue;
		if (l == ATM_O && (k == ATM_N || k == ATM_CA))
		  continue;
		if (l == ATM_CB && k == ATM_N)
		  continue;
	      }
	    }

	    // Exclude atoms too far from one another
	    dis = Conf._res[i]._atom[k].dis(Conf._res[j]._atom[l]);

	    if (dis > PF_DIS_CUT)
	      continue;

	    if (PF::cal[EM_VDW]) {
	      // Van der waals interaction, atraction and repulsion part
	      e_ij = sqrt(Atom::welldepth[Conf._res[i]._atom[k]._type] *
			  Atom::welldepth[Conf._res[j]._atom[l]._type]);

	      // check for Hbond
	      isHB = 0;
	      if (isHB == 0)
		r_ij = Atom::radius[Conf._res[i]._atom[k]._type] +
		  Atom::radius[Conf._res[j]._atom[l]._type];
	      else if (isHB == 1)
		r_ij = 2.95;
	      else if (isHB == 2)
		r_ij = 2.95;
	      if (r_ij / dis < 1.12) {
		// attraction part
		e_tmp  = (pow(r_ij / dis, 12) - 2 * pow(r_ij / dis, 6)) * e_ij;
		if (Conf._checkEn >= 0) {
		  // add this into SFMAP _EnProfile
		  long long int tmpI;
		  tmpI = i * 10000000 + k * 100000 + j * 100 + l;
		  //Conf._EnProfile[E_VDWA].insert(LLIFMAP::value_type(tmpI, e_tmp * PF::Parameter["E"][E_VDWA]));
		  Conf._EnProfile[E_VDWA].insert(LLIFMAP::value_type(tmpI, e_tmp * 0.6));
		  }
		  e_atr += e_tmp;
	      } else if (r_ij / dis < 1.33) {
		 // repulsion part
		  e_tmp  = (pow(r_ij / dis, 12) - 2 * pow(r_ij / dis, 6)) * e_ij;
		  if (debug)
		  cout << e_tmp;
		  if(Conf._checkEn >= 0) { // add this into SFMAP _EnProfile
		  long long int tmpI;
		  tmpI = i * 10000000 + k * 100000 + j * 100 + l;
		  Conf._EnProfile[E_VDWR].insert(LLIFMAP::value_type(tmpI,e_tmp*0.6));
		}
		e_rep += e_tmp;
	      } else {
		slope = -12 * e_ij * 33.383 / r_ij;
		y_intercept = -slope * r_ij / 1.33 + e_ij * 19.565;
		e_tmp  = y_intercept + dis * slope;
		if (debug)
		  cout << e_tmp;
		if(Conf._checkEn >= 0) { // add this into SFMAP _EnProfile
		  long long int tmpI;
		  tmpI = i * 10000000 + k * 100000 + j * 100 + l;
		  Conf._EnProfile[E_VDWR].insert(LLIFMAP::value_type(tmpI,e_tmp*0.6));
		}
		e_rep += e_tmp;
	      }

	    }
	  }
	}
      }
    }
  }
 // cout << "ATOMNUM   "<<Ela<<"   "<<Elb<<"   "<<Elc<<"   "<<Count1<<"   "<<Eld<<"   "<<Count2<<"   "<<e_atr<<"   "<<e_rep<<"   "<<e_atr+e_rep-Ela<<endl;
}






// this function calculate VDW energy for atoms within an amino acids
void vdw_within_aa(Structure& Conf, Residue& res, double& energy) {
  if (res._type == GLY || res._type == ALA)
    return;
  int type = res._type;
  int atom1, atom2;
  double dis;
  double r_ij, e_ij, slope, y_intercept, e_tmp; // for vdw calculation
  for (int i = 0; i < 7; i++) {
    if (PF::VDW_AA[type][i].size() != 0) {
      atom1 = i;
      for (int j = 0; j < PF::VDW_AA[type][i].size(); j++) {
        
	atom2 = PF::VDW_AA[type][i][j];
	dis = res._atom[atom1].dis(res._atom[atom2]);	
	if(res._atom[atom1]._type >=22 || res._atom[atom2]._type >= 22) continue;
	if (res._atom[atom1]._type == UNDEF || res._atom[atom2]._type == UNDEF)         {
           continue;
	}
	if(dis<PF_DIS_CUT) {
	  e_ij = sqrt((Atom::welldepth[res._atom[atom1]._type])
		      *(Atom::welldepth[res._atom[atom2]._type]));
	  r_ij = Atom::radius[res._atom[atom1]._type]+Atom::radius[res._atom[atom2]._type];
		
	  if((r_ij/dis)<1.12){ // attraction part		  
	    e_tmp = (pow(r_ij/dis,12)-2*pow(r_ij/dis,6))*e_ij;	      
	    //energy += e_tmp;
	    energy += e_tmp*5;
	    //energy += e_tmp*(PF::Parameter["E"][E_VAA]);
	  }
		
	  // repulsion part
	  else if((r_ij/dis)<1.33){	      
	    e_tmp = (pow(r_ij/dis,12)-2*pow(r_ij/dis,6))*e_ij;	      
	    energy += e_tmp*5;	
	    //energy += e_tmp*(PF::Parameter["E"][E_VAA]);	
	  }           	
	  else {          
	    slope = -12*e_ij*(33.383)*(1/r_ij);  // 33.383=(1.33^13-1.33^7)             
	    y_intercept = -1*slope*(r_ij/1.33)+e_ij*19.565; // 19.565=(1.33^12-2(1.33^6))             				
	    e_tmp = (y_intercept + dis*slope);	      				 
	   // energy += e_tmp*0.4;            	  
	    energy += e_tmp*5;            	  
	    //energy += e_tmp*(PF::Parameter["E"][E_VAA]);            	  
	  }
	  if(Conf._checkEn >= 0) { // add this into SFMAP _EnProfile
	    long long int tmpI;
	    tmpI = res._posn * 10000 + atom1 * 100 + atom2;
	  //  res._parent->_EnProfile[E_VAA].insert(LLIFMAP::value_type(tmpI,e_tmp*PF::Parameter["E"][E_VAA]));
	    res._parent->_EnProfile[E_VAA].insert(LLIFMAP::value_type(tmpI,e_tmp*0.6));
	  }
	 
	}	// end of if(dis<PF_DIS_CUT) {
      }	// end of for(int j = 0; j<PF::VDW_AA[type][i].size();++j){  

    } // end of if(PF::VDW_AA[type][i].size() != 0){

  } // end of for(int i = 0; i < 7;++i){
}

// calculate the residue pair term, the SCM of DFIRE
double res_pair_e(Structure& Conf,int Start,int End, int StartAtom)
{
  int i,j,id1,id2,binNum;
  double dis,energy=0;

  for(i=Start;i<=End;i++){
    if(Conf._res[i]._type == GLY) continue;

    for(j=1;j<=Conf._numRes;j++){
      if(j<=i && j>=Start) continue;
      if(abs(j-i)<=1) continue;
      if(Conf._res[j]._type == GLY) continue;
      dis = Conf._res[i]._FG.dis(Conf._res[j]._FG);
      if(dis>15) continue;
      else if(dis<2){
	binNum=0;
      }
      else if(dis>=2 && dis < 8){
	binNum=(int)((dis-2)/0.5)+1;
      }
      else {
	binNum=(int)(dis-8)+13;
      }
      if(Conf._checkEn >= 0){
	long long int tmpI;
	tmpI = i*1000 + j;
	Conf._EnProfile[E_RP].insert(LLIFMAP::value_type(tmpI,
							 PF::ResPair[binNum][Conf._res[i]._type][Conf._res[j]._type]*PF::Parameter["E"][E_RP]));
      }
      energy+=PF::ResPair[binNum][Conf._res[i]._type][Conf._res[j]._type];
    }
  }
  return energy;
}


// Assign secondary structures for a conformation. 1 for helix, 2 for
// sheet and 3 for coils. Following Baker's JMB h-bond paper.
// It calculate the torsion angles first.
// The calculated phi,psi and omega values are verified by torsion
// program. 12-30-2005
// The assignment of secondary structure is verified by eye
// examination of the structure of protein 2ovo and 3ebx.
// Type == 0: calculate phi and psi angles and assign secondary
// structures assuming no information before Start and after End
// Type == 1: do not calculate phi and psi angles since they are
// already calculated previously and assign secondary strutures using
// information before Start and after End if they are available
void assign_ss(Structure& Conf, int Start, int End, int Type) {

  int i,j;
  double phi,psi;

  // calculate phi and psi angles without considering the number of chains
  if (Type == 0) {

    if (Start == 1) {
      Conf._res[1]._phi = 999.99;
      Conf._res[1]._psi =
	torsion(Conf._res[1]._atom[0],
		Conf._res[1]._atom[1],
		Conf._res[1]._atom[2],
		Conf._res[2]._atom[0]);
      Conf._res[1]._omega =
	torsion(Conf._res[1]._atom[1],
		Conf._res[1]._atom[2],
		Conf._res[2]._atom[0],
		Conf._res[2]._atom[1]);
      Conf._res[1]._ss = SS_C;
    }

    for (i = Start; i <= End; i++) {
      if (i == 1 || i == Conf._numRes)
	continue;
      Conf._res[i]._phi =
	torsion(Conf._res[i-1]._atom[2],
		Conf._res[i]._atom[0],
		Conf._res[i]._atom[1],
		Conf._res[i]._atom[2]);
      Conf._res[i]._psi =
	torsion(Conf._res[i]._atom[0],
		Conf._res[i]._atom[1],
		Conf._res[i]._atom[2],
		Conf._res[i+1]._atom[0]);
      Conf._res[i]._omega =
	torsion(Conf._res[i]._atom[1],
		Conf._res[i]._atom[2],
		Conf._res[i+1]._atom[0],
		Conf._res[i+1]._atom[1]);      
    }

    if (End == Conf._numRes) {
      j = Conf._numRes;
      Conf._res[j]._phi =
	torsion(Conf._res[j-1]._atom[2],
		Conf._res[j]._atom[0],
		Conf._res[j]._atom[1],
		Conf._res[j]._atom[2]);
      Conf._res[j]._psi=999.99;
      Conf._res[j]._omega=999.99;
      Conf._res[j]._ss = SS_C;
    }
  }

  // assign secondary structure according to Baker's simple criterion
  // helix 1, strand 2, coil 3.
  // helix: -180<phi<-20, -90<psi<-10; sheet -180<phi<-20, 180>psi>20
  // or -180<psi<-170. Classification of helix or strand required that
  // at least two adjacent residues have the same secondary structure.
  // Two loops are run with the first time just assign the helix or
  // strand types without considering the two adjacent residue
  // criterion, and the second time add this criterion.
  for (i = Start; i <= End; i++) {
    phi = Conf._res[i]._phi;
    psi = Conf._res[i]._psi;
    if (-180 < phi && phi < -20 && -90 < psi && psi < -10)
      Conf._res[i]._ss = SS_H;
    else if (-180 < phi && phi < -20 &&
	     ((180 > psi && psi > 20) || (-180 < psi && psi < -170)))
      Conf._res[i]._ss = SS_E;
    else
      Conf._res[i]._ss = SS_C;
  }

  // assign real secondary structure for all residues
  for (i = Start; i <= End; i++) {
    if (i == 1 || i == Conf._numRes) continue; // they have been assigned as 3
    if (Conf._res[i]._ss == 1 || Conf._res[i]._ss == 2) {
      if (Conf._res[i-1]._ss != Conf._res[i]._ss &&
	  Conf._res[i+1]._ss != Conf._res[i]._ss)
	Conf._res[i]._ss = SS_C;
    }
  }
}

double ss_e(Structure& Conf, int Start, int End) {
  // Compute energy based on secondary structure predicted vs. actual
  // Need to step one outside fragment because SS could've changed there
  int Start2 = (Start == 1) ? 1 : (Start - 1);
  int   End2 = (End == Conf._numRes) ? End : (End + 1);
  // Compare predicted to actual
  double total = 0;
  for (int i = Start2; i <= End2; i++) {
    if(Conf._checkEn >= 0) {
      // add this into SFMAP _EnProfile
      long long int tmpI;
      tmpI = i*100 + Conf._res[i]._ss;
      Conf._EnProfile[E_SS].insert(LLIFMAP::value_type(tmpI,
						       -log((Conf._ssPredProb[i][Conf._res[i]._ss] + 0.5) / 100.)*PF::Parameter["E"][E_SS]));
    }

    double CurProb = Conf._ssPredProb[i][Conf._res[i]._ss];
    if (CurProb <= 0) CurProb = 0.001;
    if (CurProb >= 1) CurProb = 0.999;
    total += -log(CurProb);
  }
  return(total);
}

double contact_num_e(Structure& Conf, int Start, int End) {
  // Compute energy based on contact numbers
  double sa_total = 0;
  double cn_total = 0;

  for (int i = 1; i <= Conf._numRes; i++) {
    cn_total += Conf._res[i]._res_adj.size();
    sa_total += Conf._saPred[i];
  }

  double sa_diff = 0;
  // Sum standardized absolute differences
  for (int i = 1; i <= Conf._numRes; i++)
    sa_diff +=  abs( Conf._saPred[i] / sa_total -
		     Conf._res[i]._res_adj.size() / cn_total );
  return sa_diff;
}

double src_e(Structure& Conf, int Start, int End) {
  // Compute energy based on short-range correlation
  double e_total = 0;

  // Affects residues within a range of 2, so for a regrown fragment
  // recalculate the affected distances
  int calc_start;
  int calc_end;
  if (Start - 2 < 1)
    calc_start = 1;
  else
    calc_start = Start - 2;

  if (End+2 > Conf._numRes)
    calc_end = Conf._numRes-2;
  else
    calc_end = End;
		
  // Sum up the terms
  for (int i = calc_start; i <= calc_end; i++) {
    double distance = Conf._res[i]._center.dis(Conf._res[i+2]._center);
    if (distance < 4)
      e_total += PF::SRC_below4[Conf._res[i]._type][Conf._res[i+2]._type];
    if (distance >= 4)
      e_total += PF::SRC_above4[Conf._res[i]._type][Conf._res[i+2]._type];
  }

  return e_total;
}

double model_e (Structure& Conf,int Start, int End) {
  // energy based on difference of structure from Modeler (or other software)
  // input
  return rmsd(Conf, PF::ModelStart, Start, End, Start, End, 1,0);
}

// rama_e calculate the rama term of the potential function
// it uses the baker potential file name: Rama_smooth_dyn.dat_ss_6.4
double rama_e(Structure& Conf, int Start, int End, double* enArr) {
  double phi,psi,energy=0;

  // Need to step one outside fragment because SS could've changed there
  int Start2 = (Start == 1) ? 1 : (Start - 1);
  int   End2 = (End == Conf._numRes) ? End : (End + 1);

  // the two end residue are omited.
  // assign_ss(Conf, Start, End, 1);
  for (int i = Start2; i <= End2; i++) {
    // no Rama energy for end residues
    if (i == 1 || i == Conf._numRes)
      continue;

    // Get phi in (0, 360)
    double phi = Conf._res[i]._phi;
    if (phi < 0)
      phi += 360;
    double psi = Conf._res[i]._psi;
    if (psi < 0)
      psi += 360;

    // Get Rama energy
    int curType = Conf._res[i]._type;
    int curSS   = Conf._res[i]._ss;
    if (curSS == SS_H)
      enArr[E_RAMA] += PF::Rama[curType][curSS][int(phi/10)][int(psi/10)];
    else if (curSS == SS_E)
      enArr[E_RAMAE] += PF::Rama[curType][curSS][int(phi/10)][int(psi/10)];
    else
      enArr[E_RAMAC] += PF::Rama[curType][curSS][int(phi/10)][int(psi/10)];
    // energy += PF::Rama[curType][curSS][int(phi/10)][int(psi/10)];
  }
  return (energy);
}

void find_hb(Structure& Conf, int StartAtom) {
  // Finds the HBonds in a structure

  // First, clear out all previous H bond information
  for (int i = 1; i <= Conf._numRes; i++)
    for (int j = 0; j <= Conf._res[i]._numAtom; j++) {
      Conf._res[i]._atom[j].HGiven.clear();
      Conf._res[i]._atom[j].HTaken.clear();
    }

  // Next, scan the entire structure for all possible H bonds
  typedef multimap<double, HBond, less<double> > DisHBond;
  DisHBond PossibleHB;
  for (int i = 1; i <= Conf._numRes; i++) {
    // i = Donor residue
    for (int j = StartAtom; j < Conf._res[i]._numAtom; j++) {
      // j = Possible donor atom. Has to be either a backbone or sidechain H.
      // If the type is unknown, forget it
      if (Conf._res[i]._atom[j]._type == UNDEF)
	continue;
      // Proline doesn't have a backbone H
      if (Conf._res[i]._type == PRO && j == ATM_H)                 
	continue;
      // Current atom must be an Hbond donor
      if (Atom::hbondH[Conf._res[i]._atom[j]._type] == 0)
	continue;

      for (ISET::iterator ilItr = Conf._res[i]._res_adj.begin();
	   ilItr != Conf._res[i]._res_adj.end(); ilItr++) {
	// k = current neighbor residue to residue i
	int k = (*ilItr);
	//if(i==48) 
	for (int l = 0; l <= Conf._res[k]._numAtom; l++) {
	  // l = current atom of neighbor residue k
	  // Skip non-acceptors
	  if (Conf._res[k]._atom[l]._type == UNDEF)
	    continue;
	  if (Conf._res[k]._type == PRO && l == ATM_H)
	    continue;
	  if (Atom::acceptor[Conf._res[k]._atom[l]._type] == 0)
	    continue;
	  // h-bond cannot be formed between backbone atoms of adjacent residues
	  if (abs(i-k) < 2 && j <= NUM_BB_ATOM && l <= NUM_BB_ATOM)
	    continue;

	  // This atom pair can be an Hbond!
	  double dis = Conf._res[i]._atom[j].dis(Conf._res[k]._atom[l]);
	  if (dis > MAX_DISTANCE_HBOND)
	    continue;
	  HBond curHB;
	  curHB.DonorAtom    = &Conf._res[i]._atom[j];
	  curHB.AcceptorAtom = &Conf._res[k]._atom[l];
	  if (curHB.DonorAtom->_posn == 32767) {
	    cout << i << " " << j << endl;
	    cout << Conf._res[i]._atom[j]._posn << endl;
	    cout << "DEAD" << endl;
	    exit(0);
	  }
	  PossibleHB.insert(make_pair(dis, curHB));
	}
      }
    }
  }

  // Now that we have assembled all possible Hbonds, we step through the list,
  // which is sorted by distance automatically, and assign the nearest H bonds
  // first, allowing only one bond per donor and two per acceptor
  for (DisHBond::iterator i = PossibleHB.begin(); i != PossibleHB.end(); i++) {
    // Add this HBond if the donor and acceptor both have a free slot
    if (i->second.DonorAtom->HGiven.size() < MAX_HBOND_PER_DONOR &&
	i->second.AcceptorAtom->HTaken.size() < MAX_HBOND_PER_ACCEPTOR) {
      i->second.DonorAtom->HGiven.push_back(i->second.AcceptorAtom);
      i->second.AcceptorAtom->HTaken.push_back(i->second.DonorAtom);
    } else {
    }
  }  
}

// calculate the h-bond energy. All possible hbond have been stored in
// hb_list obtained in function calculating vdw and solvation energy.
void hb_e(Structure& Conf, double& e_hb_bb, double& e_hb_sc, int Start, int End,
	  int StartAtom, double& e_hav, double& e_hb) {
  int type;  // bb (0) or sc (1) 
  int type_ss; // helix (0), sheet (1) or coil (2)
  int type_sp; // sp2 (0) or sp3 (1)
  int type_angle; // chi (0) psi (1) or theta (2)
  int type_dis;  // short (0) or long(1)
  int type_env=0;  // environmental, 0, exposed-exposed and exposed-intermediate
  IIMAP::iterator iiItr;
  // 1, exposed-buried and intermediate-intermediate
  // 2, intermediate-buried and buried-buried

  double hb_bb[3]  = {0, 0, 0};
  double hb_sc[3]  = {0, 0, 0};
  double env_wt[3] = {0.18, 0.28, 0.91};

  ILIST::iterator ilItr;
  double dis, chi, psi, theta, energy_bb = 0, energy_sc = 0;
  cout<<"HB calculation "<<Conf._numRes<<endl;
  find_hb(Conf, StartAtom);
   for (int i = Start; i < End; i++) {
    vector<HBond> AllHB = Conf._res[i].AllHBond(true, true, StartAtom);
    for (int j = 0; j < AllHB.size(); j++) {
      // Loop over HBonds linked to this residue
      // hbondH_res_ind is the residue index of hbondH
      // hbondH_atm_ind is the atom index of the hbondH atom
      // donor_ind is the atom index of the donor atom
      // acc_res_ind is the residue index of the acceptor residue
      // acc_atm_ind is the atom index of the acceptor atom
      // ab_ind is the atom index of the acceptor base atom
      int hbondH_res_ind = AllHB[j].DonorAtom->_parent->_posn;
      int hbondH_atm_ind = AllHB[j].DonorAtom->_posn;
      int donor_ind =
	Residue::prev_atom[Conf._res[hbondH_res_ind]._type][hbondH_atm_ind][2];
      int acc_res_ind    = AllHB[j].AcceptorAtom->_parent->_posn;
      int acc_atm_ind    = AllHB[j].AcceptorAtom->_posn;
      int ab_ind    =
	Residue::prev_atom[Conf._res[acc_res_ind]._type][acc_atm_ind][2];
      int ab_prev_ind =
	Residue::prev_atom[Conf._res[acc_res_ind]._type][acc_atm_ind][1];

      // Determine whether we should skip this HBond to avoid double-counting
      // Rule = if both residues are in the fragment, we only take this HBond
      // if Residue i is the donor
      if (hbondH_res_ind >= Start && hbondH_res_ind <= End &&
	  acc_res_ind >= Start && acc_res_ind <= End &&
	  acc_res_ind == i)
	continue;

      // short, long or ignore
      if ((dis = Conf._res[hbondH_res_ind]._atom[hbondH_atm_ind].dis(Conf._res[acc_res_ind]._atom[acc_atm_ind])) <= 1.4 ||
	  dis >= 3.0)
	// not counted as a hbond
	continue;
      else if (dis < 2.1)
	// short distance, cutoff is 2.1 A
	type_dis = 0;
      else
	// long distance hbond
	type_dis = 1;

      // determine the type of the hbond
      energy_bb = 0;
      energy_sc = 0;
      // bb or sc
      type       = -1;
      type_ss    = -1;
      type_sp    =  0;   // default value is 0
      type_angle = -1;
      type_dis   = -1;
      if (hbondH_atm_ind < 5 && acc_atm_ind < 5)
	type = 0; // bb type
      else
	type = 1;

      // helix, sheet or other
      if (Conf._res[hbondH_res_ind]._ss == SS_H &&
	  Conf._res[acc_res_ind]._ss    == SS_H)
	type_ss = 0; // both residue are helix
      else if (Conf._res[hbondH_res_ind]._ss == SS_E &&
	       Conf._res[acc_res_ind]._ss    == SS_E)
	type_ss = 1; // both residues are sheet
      else
	type_ss = 2; // other type

      // sp2 or sp3 if it is sc (type=1)
      if (type == 1) {
	// one of the atoms should be acceptor, if it is also donor
	// than the type_sp is sp3 since that atom is a hydroxyl
	// oxygen, and hydroxyl oxygen is the only sp3 acceptor
	if(Atom::acceptor[Conf._res[acc_res_ind]._atom[acc_atm_ind]._type] == 1
	   && Atom::donor[Conf._res[acc_res_ind]._atom[acc_atm_ind]._type] == 1)
	  // hbond is sp3 type, since the acceptor is both a donor and
	  // acceptor, hydroxyl.
	  type_sp = 1;
	else
	  // sp2 type
	  type_sp = 0;
      }

      // environment, using the list of adjacent centers, Baker use
      // the criterion defined by the number of Cb atoms within a
      // sphere of 10 A radius of the Cb atom of the residue of
      // interest: exposed 0-14, intermediate 15-20, buried >20.
      // since different criteria are used, the environment weights
      // may need to be recalculated.
      int env_1,env_2;
      if (Conf._res[hbondH_res_ind]._res_adj.size() > 20)
	env_1 = 2;
      else if (Conf._res[hbondH_res_ind]._res_adj.size() > 14)
	env_1 = 1;
      else
	env_1 = 0;

      if (Conf._res[acc_res_ind]._res_adj.size() > 20)
	env_2 = 2;
      else if (Conf._res[acc_res_ind]._res_adj.size() > 14)
	env_2 = 1;
      else
	env_2 = 0;

      if (env_1 + env_2 <= 1)
	type_env = 0;
      else if (env_1 + env_2 == 2)
	type_env = 1;
      else
	type_env = 2;

      // calculate the angles
      chi = torsion(Conf._res[acc_res_ind]._atom[ab_prev_ind],
		    Conf._res[acc_res_ind]._atom[ab_ind],
		    Conf._res[acc_res_ind]._atom[acc_atm_ind],
		    Conf._res[hbondH_res_ind]._atom[hbondH_atm_ind]);
      theta = angle(Conf._res[acc_res_ind]._atom[acc_atm_ind],
		    Conf._res[hbondH_res_ind]._atom[hbondH_atm_ind],
		    Conf._res[hbondH_res_ind]._atom[donor_ind]);
      psi = angle(Conf._res[acc_res_ind]._atom[ab_ind],
		  Conf._res[acc_res_ind]._atom[acc_atm_ind],
		  Conf._res[hbondH_res_ind]._atom[hbondH_atm_ind]);
      if (isnan(chi) || isnan(theta) || isnan(psi) ||
	  Conf._res[acc_res_ind]._atom[acc_atm_ind]._type == UNDEF ||
	  Conf._res[hbondH_res_ind]._atom[donor_ind]._type == UNDEF) {
	cout << "Missing atom involved in H bond!" << endl;
	Conf._res[acc_res_ind]._atom[acc_atm_ind].out(); cout << endl;
	Conf._res[hbondH_res_ind]._atom[hbondH_atm_ind].out(); cout << endl;
	Conf._res[hbondH_res_ind]._atom[donor_ind].out(); cout << endl;
	cout << acc_res_ind    << " " << acc_atm_ind << endl;
	cout << hbondH_res_ind << " " << hbondH_atm_ind << endl;
	cout << hbondH_res_ind << " " << donor_ind << endl;
	continue;
      }

      // update h-bond energy
      // the form of hbond energy is E_h = W_h*(E_d + E_chi + E_theta + E_psi)
      // the hbonds are further divided into three environment categories:
      // exposed, intermediate and buried. 

      // hb_dis
      int dis_label;
      if (type == 0) {
	// backbone-backbone hbond
	dis_label  = type_ss;
	energy_bb += PF::HB_dis[dis_label][int((dis-1.4)/0.05)];
      } else {
	// side-chain hbond
	dis_label  = 3 + type_sp;
	energy_sc += PF::HB_dis[dis_label][int((dis-1.4)/0.05)];
      }

      // angle terms      
      // prevent array subscripts from overflow
      if (chi == 180)
	chi -= 0.0001;
      if (psi == 180)
	psi -= 0.0001;
      if (theta == 180)
	theta -= 0.0001;

      if (type == 0) {
	// backbone-backbone hbond
	// chi
	energy_bb += PF::HB_angle_bb[type_ss*3+0][int((chi+180)/5)];
	// psi
	energy_bb += PF::HB_angle_bb[type_ss*3+1][int((psi+180)/5)];
	// theta
	energy_bb += PF::HB_angle_bb[type_ss*3+2][int((theta+180)/5)];
      } else {
	// side-chain hbond

	if (type_sp == 0) {
	  // sp2

	  if (type_dis == 0) {
	    // short
	    // chi
	    energy_sc += PF::HB_angle_sc[0][int(abs(chi)/10)];
	    // psi
	    energy_sc += PF::HB_angle_sc[1][int(psi/10)];
	    // theta
	    energy_sc += PF::HB_angle_sc[2][int(theta/10)];	    
	  } else {
	    // long
	    // chi
	    energy_sc += PF::HB_angle_sc[3][int(abs(chi)/10)];
	    // psi
	    energy_sc += PF::HB_angle_sc[4][int(psi/10)];
	    // theta
	    energy_sc += PF::HB_angle_sc[5][int(theta/10)];	    
	  }

	} else {
	  // sp3, no chi

	  if (type_dis==0) {
	    // short
	    // psi
	    energy_sc += PF::HB_angle_sc[6][int(psi/10)];
	    // theta
	    energy_sc += PF::HB_angle_sc[7][int(theta/10)];	  
	  } else {
	    // long
	    // psi
	    energy_sc += PF::HB_angle_sc[8][int(psi/10)];
	    // theta
	    energy_sc += PF::HB_angle_sc[9][int(theta/10)];	  
	  }

	}
      }

      hb_bb[type_env] += energy_bb;
      hb_sc[type_env] += energy_sc;
      
      long long int tmpI;
      if (Conf._checkEn >=0) {
	tmpI = hbondH_res_ind * 10000000 + hbondH_atm_ind * 100000 +
	  acc_res_ind * 100 + acc_atm_ind;
	if (type == 0)
	  Conf._EnProfile[E_HBB].insert(LLIFMAP::value_type(tmpI, energy_bb * env_wt[type_env] * PF::Parameter["E"][E_HBB]));
	else
	  Conf._EnProfile[E_HBS].insert(LLIFMAP::value_type(tmpI, energy_sc * env_wt[type_env] * PF::Parameter["E"][E_HBS]));
      }

      // insert the instance into _HBTermCnt
      // if the first digit is
      // 1: distance term
      // 2: chi angle term
      // 3: psi angle term
      // 4: theta angle term
      // insert distance
      tmpI = 1 * 1000000 + type * 100000 + type_ss * 10000 + type_sp * 1000 +
	type_env * 100 + (int)(dis * 20);
      if (tmpI < 1000000)
	cout << "tmpI error " << tmpI << " " << dis << " " << type << " "
	     << type_ss << " " << type_sp << " " << type_env << endl;
      if (PF::COUNT &&
	  Conf._res[hbondH_res_ind]._atom[hbondH_atm_ind]._state != 1 &&
	  Conf._res[acc_res_ind]._atom[acc_atm_ind]._state != 1) {
	if ((iiItr=Conf._HBTermCnt.find(tmpI)) != Conf._HBTermCnt.end())
	  iiItr->second++;
	else
	  Conf._HBTermCnt.insert(IIMAP::value_type(tmpI, 1));
      }
      // Now add the contribution from this term to the entire energy
      IDMAP::iterator ii;
      if ((ii = PF::Parameter["H"].find(tmpI)) != PF::Parameter["H"].end())
	e_hb += ii->second;

      // insert chi angle
      tmpI = 2 * 1000000 + type * 100000 + type_ss * 10000 + type_sp * 1000 +
	type_env * 100 + (int)((int)(chi + 360) % 360 / 5);
      if (PF::COUNT &&
	  Conf._res[hbondH_res_ind]._atom[hbondH_atm_ind]._state != 1 &&
	  Conf._res[acc_res_ind]._atom[acc_atm_ind]._state != 1) {
	if ((iiItr = Conf._HBTermCnt.find(tmpI)) != Conf._HBTermCnt.end())
	  iiItr->second++;
	else
	  Conf._HBTermCnt.insert(IIMAP::value_type(tmpI, 1));
      }
      if ((ii = PF::Parameter["H"].find(tmpI)) != PF::Parameter["H"].end())
	e_hb += ii->second;

      // insert psi angle
      tmpI = 3 * 1000000 + type * 100000 + type_ss * 10000 + type_sp * 1000 +
	type_env * 100 + (int)((int)(psi + 360) % 360 / 5);
      if (PF::COUNT &&
	  Conf._res[hbondH_res_ind]._atom[hbondH_atm_ind]._state != 1 &&
	  Conf._res[acc_res_ind]._atom[acc_atm_ind]._state != 1) {
	if ((iiItr = Conf._HBTermCnt.find(tmpI)) != Conf._HBTermCnt.end())
	  iiItr->second++;
	else
	  Conf._HBTermCnt.insert(IIMAP::value_type(tmpI, 1));
      }
      if ((ii = PF::Parameter["H"].find(tmpI)) != PF::Parameter["H"].end())
	e_hb += ii->second;

      // insert theta angle
      tmpI = 4 * 1000000 + type * 100000 + type_ss * 10000 + type_sp * 1000 +
	type_env * 100 + (int)((int)(theta + 360) % 360 / 5);
      if (PF::COUNT &&
	  Conf._res[hbondH_res_ind]._atom[hbondH_atm_ind]._state != 1 &&
	  Conf._res[acc_res_ind]._atom[acc_atm_ind]._state != 1) {
	if ((iiItr=Conf._HBTermCnt.find(tmpI)) != Conf._HBTermCnt.end())
	  iiItr->second++;
	else
	  Conf._HBTermCnt.insert(IIMAP::value_type(tmpI, 1));
      }
      if ((ii = PF::Parameter["H"].find(tmpI)) != PF::Parameter["H"].end())
	e_hb += ii->second;
      
      // calculate E_HAV, h-bond adjusted VDW energy
      double r_ij,e_ij,slope,y_intercept,e_tmp;
      e_ij = sqrt((Atom::welldepth[Conf._res[acc_res_ind]._atom[acc_atm_ind]._type])
		  * (Atom::welldepth[Conf._res[hbondH_res_ind]._atom[donor_ind]._type]));
      r_ij = Atom::radius[Conf._res[acc_res_ind]._atom[acc_atm_ind]._type] + 
	Atom::radius[Conf._res[hbondH_res_ind]._atom[donor_ind]._type];
      dis = Conf._res[acc_res_ind]._atom[acc_atm_ind].dis(Conf._res[hbondH_res_ind]._atom[donor_ind]);

      // this is the energy calculated in vdw_sol(),
      // which should be subtracted from the total energy
      if ((r_ij / dis) < 1.12)
	// attraction part
	e_tmp = (pow(r_ij/dis,12)-2*pow(r_ij/dis,6))*e_ij;
      else if (r_ij / dis < 1.33)
	// repulsive
	e_tmp = (pow(r_ij/dis,12)-2*pow(r_ij/dis,6))*e_ij;	      
      else {
	slope = -12*e_ij*(33.383)*(1/r_ij);  // 33.383=(1.33^13-1.33^7)
	y_intercept = -1*slope*(r_ij/1.33)+e_ij*19.565; // 19.565=(1.33^12-2(1.33^6))
	e_tmp = (y_intercept + dis*slope);
      }

      e_tmp = -1 * e_tmp;

      // actual vdw energy by setting r_ij to 2.95 to soften it 
      r_ij = 2.95;

      if ((r_ij / dis) < 1.12)
	// attraction part
	e_tmp += (pow(r_ij/dis,12)-2*pow(r_ij/dis,6))*e_ij;
      else if ((r_ij/dis) < 1.33)
	// repulsive
	e_tmp += (pow(r_ij/dis,12)-2*pow(r_ij/dis,6))*e_ij;
      else {
	slope = -12*e_ij*(33.383)*(1/r_ij);
	y_intercept = -1*slope*(r_ij/1.33)+e_ij*19.565;
	e_tmp += (y_intercept + dis*slope);
      }

      e_hav += e_tmp;
      if (Conf._checkEn >=0) {
	tmpI = hbondH_res_ind * 10000000 + donor_ind * 100000 +
	  acc_res_ind * 100 + acc_atm_ind;
	Conf._EnProfile[E_HAV].insert(LLIFMAP::value_type(tmpI, e_tmp *
							  PF::Parameter["E"][E_HAV]));
      }
    }
  }

  for (int i = 0; i < 3; i++) {
    e_hb_bb += hb_bb[i] * env_wt[i];
    e_hb_sc += hb_sc[i] * env_wt[i];
  }
}

IDMAP calEraw(Structure& Conf) {
  // Computes the full energy of a conformation in raw, unweighted terms
  // Restores weights afterwards and energies afterwards
  Structure ConfCopy(Conf);
  IDMAP oldParameter;
  for (int i = 0; i < ENERGY_TYPES; i++) {
    oldParameter[i]       = PF::Parameter["E"][i];
    PF::Parameter["E"][i] = 1;
  }
  calE(ConfCopy, 1, ConfCopy._numRes, true);
  Conf._AATermCnt  = ConfCopy._AATermCnt;
  Conf._HBTermCnt  = ConfCopy._HBTermCnt;
  Conf._BBTTermCnt = ConfCopy._BBTTermCnt;
  Conf._SCTTermCnt = ConfCopy._SCTTermCnt;
  IDMAP rst;
  for (int i = 0; i < ENERGY_TYPES; i++) {
    PF::Parameter["E"][i] = oldParameter[i];
    rst[i] = ConfCopy._enArr[i];
  }
  return(rst);
}

void printRecord(ostream& out, Structure& Conf, string other) {
  printRecord(out, Conf, other, true, true);
}

void printRecord(ostream& out, Structure& Conf, string other, bool weighted) {
  printRecord(out, Conf, other, weighted, true);
}

void printRecord(ostream& out, Structure& Conf, string other, bool weighted,
		 bool count) {
  // Print a data line for the structure

  // Get the energy, directly from the structure if weighted is true.
  // For unweighted, we currently have to recompute the energy.
  double TotalEnergy;
  IDMAP Energy;
  //cout << "weight		"<<weighted <<endl;
  if (weighted) {
    calE(Conf, 1, Conf._numRes, true);
    TotalEnergy = Conf._energy;
    for (int i = 0; i < ENERGY_TYPES; i++)
      Energy[i] = Conf._enArr[i];
  } else {
    TotalEnergy = UNDEF;
    Energy = calEraw(Conf);
  }
  // Print the metadata first, followed by the total energy
  out << other << " " << TotalEnergy << " ";

  // Print the unweighted energy
  for (int i = 0; i < ENERGY_TYPES; i++)
    out << "E" << i << ":" << Energy[i] << " ";

  // Print out all the counts only if we have been asked.
  if (count) {
    // Print the atom-atom interaction terms
    if (PF::cal[EM_CNT]) {
      IIMAP::iterator iiItr;
      for (iiItr  = Conf._AATermCnt.begin();
	   iiItr != Conf._AATermCnt.end();
	   iiItr++)
	out << "A" << iiItr->first << ":" << iiItr->second << " ";
    }

    if (PF::cal[EM_BBT]) {
      IIMAP::iterator iiItr;
      for (iiItr  = Conf._BBTTermCnt.begin();
	   iiItr != Conf._BBTTermCnt.end();
	   iiItr++)
	out << "B" << iiItr->first << ":" << iiItr->second << " ";
    }

    if (PF::cal[EM_HB]) {
      IIMAP::iterator iiItr;
      for (iiItr  = Conf._HBTermCnt.begin();
	   iiItr != Conf._HBTermCnt.end();
	   iiItr++)
	out << "H" << iiItr->first << ":" << iiItr->second << " ";
    }

    if (PF::cal[EM_SCT]) {
      IIMAP::iterator iiItr;
      for (iiItr  = Conf._SCTTermCnt.begin();
	   iiItr != Conf._SCTTermCnt.end();
	   iiItr++)
	out << "S" << iiItr->first << ":" << iiItr->second << " ";
    }
  }

  // Done
  out << endl;
}


// loodis energy calculation, type = false No side-chain, type = true, side chain added
double loodis_e(Structure& Conf, int Start, int End, bool type)
{
   int len,i,j,k, disInd, numatom1, numatom2;
   double dis,r = 0;
   double energy =0;
   int atomIndex1, atomIndex2;

   for(len = Start; len <End+1; ++len)
   {
	for(i = 1; i <=Conf._numRes;++i)
	{
           if(i < len + 2 && i > len -2 )		   // exclude neighbour residue 
	        continue;
	   if (i >= Start && i <= End) 
	     if(i <= len -2)     //avoid calculate the same energy twice
             continue;
	  
	  if(type)
	   {
             numatom1 = Conf._res[len]._numAtom;
	     numatom2 = Conf._res[i]._numAtom;
	   }
	   else
	   {
	      if(Conf._res[len]._type != 5)
	         numatom1 = 6;
	      else
	         numatom1 = 5;
	      
	      if (i >= Start && i <= End) 
	      {
	          if(Conf._res[i]._type != 5)
	            numatom2 = 6;
	          else
	            numatom2 = 5;
	      }
	      else
	      numatom2 = Conf._res[i]._numAtom;
	   }

	   if(Conf._res[len]._center.dis(Conf._res[i]._center)<CC_DIS_CUT){  
	 	  for(j=0;j < numatom1;++j){
          	  if(Conf._res[len]._atom[j]._type==UNDEF || Conf._res[len]._atom[j]._type >= 22) continue;
          	  if(Conf._res[len]._atom[j].x== 0 && Conf._res[len]._atom[j].y== 0 && Conf._res[len]._atom[j].z== 0) continue;
		  for(k=0; k < numatom2;++k){
           	   if(Conf._res[i]._atom[k]._type==UNDEF || Conf._res[i]._atom[k]._type >= 22) continue;
          	   if(Conf._res[i]._atom[k].x== 0 && Conf._res[i]._atom[k].y== 0 && Conf._res[i]._atom[k].z== 0) continue;
	           
            	 	atomIndex1 = Conf._res[len]._atom[j]._type;
            	 	atomIndex2 = Conf._res[i]._atom[k]._type;
            	 	dis = Conf._res[len]._atom[j].dis(Conf._res[i]._atom[k]);
               if(dis < (H_INLO* LOODIS_DIS_BIN)) 
			   {
			   r = Atom::radius[Conf._res[len]._atom[j]._type] + Atom::radius[Conf._res[i]._atom[k]._type];
                           disInd = (int)(dis/H_INLO);
			   if(dis/r < 0.7)
               {
			   energy += PF::LOODIS[atomIndex1-1][atomIndex2-1][disInd];				                                      
			   }
			   else
			   {
            		energy += PF::LOODIS[atomIndex1-1][atomIndex2-1][disInd];
			   }
             }
		  }
		}
	   }
	}
    }
    return energy;
}


void aa_vdw_enum(short type1, short type2, double dis, ofstream& outfile) {
  double r, eps; // constants when two atoms are known
  double e_tmp,e_aap,slope,y_intercept;
  // van der waals interaction energy
  double tmp_vdwa1 =0, tmp_vdwa2 =0, tmp_vdwr1 =0, tmp_vdwr2 =0;
    r = Atom::radius[type1]+Atom::radius[type2];
    eps=sqrt((Atom::welldepth[type1])*(Atom::welldepth[type2]));
    // attraction part
    if((r/dis)<1.12){ 
      tmp_vdwa1 = ((pow(r/dis,12)-2*pow(r/dis,6))*eps);
    }
    // repulsion part
    else if((r/dis)<1.33){
      tmp_vdwr1 = ((pow(r/dis,12)-2*pow(r/dis,6))*eps);
   }
    else {
      slope=-12*eps*(33.383)*(1/r);  // 33.383=(1.33^13-1.33^7)
      y_intercept=-1*slope*(r/1.33)+eps*19.565; // 19.565=(1.33^12-2(1.33^6))
      tmp_vdwr1 = y_intercept + dis*slope;
    }
    outfile <<type1<<" "<<type2<<" "<<dis<<" "<<tmp_vdwa1<<" "<<tmp_vdwr1<<endl;

}
	  
//Steric Clash detection for single loop + REDCELL
void Clash_detection_list(Structure& conf, int Start, int End, vector<int>& ResIdx, vector<int>& ClashNum, int* List, int List_size )
{
     int i, p,j,l,k,clash_count =0;
     double disquare, dis, r, quot =0;

     for(l = Start; l <= End; l++){
      if( conf._res[l]._type == 0 || conf._res[l]._type == 5)
      continue;
        for (p = 0; p < List_size; p++) {
        i = List[p];
	
       if(i == l)
       {
         for (j = NUM_BB_ATOM; j < conf._res[l]._numAtom; j++) {
	  // skip undefined and H atoms
	    if (conf._res[l]._atom[j]._type == UNDEF || conf._res[l]._atom[j]._type >= 22)
	     continue;
            if(conf._res[l]._atom[j].x== 0 && conf._res[l]._atom[j].y== 0 && conf._res[l]._atom[j].z== 0) continue;
	    for(k = 2; k <= 3; ++k) {
	       if(conf._res[i]._atom[k]._type==UNDEF || conf._res[i]._atom[k]._type >= 22) continue;
               if(conf._res[i]._atom[k].x== 0 && conf._res[i]._atom[k].y== 0 && conf._res[i]._atom[k].z== 0) continue;
	       disquare = (conf._res[l]._atom[j].x-conf._res[i]._atom[k].x)*(conf._res[l]._atom[j].x-conf._res[i]._atom[k].x)
	                 +(conf._res[l]._atom[j].y-conf._res[i]._atom[k].y)*(conf._res[l]._atom[j].y-conf._res[i]._atom[k].y)
		         +(conf._res[l]._atom[j].z-conf._res[i]._atom[k].z)*(conf._res[l]._atom[j].z-conf._res[i]._atom[k].z);
	       if(disquare <= PF_DIS_CUT_SQUARE)
	       {
	          dis = sqrt(disquare);	
                  r = Atom::radius[conf._res[l]._atom[j]._type] + Atom::radius[conf._res[i]._atom[k]._type];
		  quot = dis/r;
                  if(quot <= VDW_CLASH_CUTOFF)
  	          {
	            if(quot <= VDW_CLASH_CUTOFF && quot >= 0.5) 
	            clash_count ++;
		    else
		    clash_count += 5;
	          }
	       }
	    }	
	  }  
          
       }
       else if (i > l && i <= End) continue;
       else if (i >= Start && i < l)
       {
         for (j = NUM_BB_ATOM; j < conf._res[l]._numAtom; j++) {
	  // skip undefined and H atoms
	    if (conf._res[l]._atom[j]._type == UNDEF || conf._res[l]._atom[j]._type >= 22)
	     continue;
            if(conf._res[l]._atom[j].x== 0 && conf._res[l]._atom[j].y== 0 && conf._res[l]._atom[j].z== 0) continue;
	    for(k = 0; k < conf._res[i]._numAtom; ++k) {
	       if(conf._res[i]._atom[k]._type==UNDEF || conf._res[i]._atom[k]._type >= 22) continue;
               if(conf._res[i]._atom[k].x== 0 && conf._res[i]._atom[k].y== 0 && conf._res[i]._atom[k].z== 0) continue;
	       disquare = (conf._res[l]._atom[j].x-conf._res[i]._atom[k].x)*(conf._res[l]._atom[j].x-conf._res[i]._atom[k].x)
	                 +(conf._res[l]._atom[j].y-conf._res[i]._atom[k].y)*(conf._res[l]._atom[j].y-conf._res[i]._atom[k].y)
		         +(conf._res[l]._atom[j].z-conf._res[i]._atom[k].z)*(conf._res[l]._atom[j].z-conf._res[i]._atom[k].z);
	       if(disquare <= PF_DIS_CUT_SQUARE)
	       {
	          dis = sqrt(disquare);	
                  r = Atom::radius[conf._res[l]._atom[j]._type] + Atom::radius[conf._res[i]._atom[k]._type];
		  quot = dis/r;
                  if(quot <= VDW_CLASH_CUTOFF)
  	          {
	            if(quot <= VDW_CLASH_CUTOFF && quot >= 0.5) 
	            clash_count ++;
		    else
		    clash_count += 5;
	          }
	       }
	    }	
	  }  
         
       }
       else
       {
        if(( conf._res[l]._center.dis(conf._res[i]._center) < Residue::size[conf._res[l]._type] + Residue::size[conf._res[i]._type] + CUB_SIZE) ||
	 (conf._res[l]._bbc.dis(conf._res[i]._center) < Residue::bb_size + Residue::size[conf._res[i]._type] + CUB_SIZE)) {
	//if (res._center.dis(conf._res[i]._center) < CC_DIS_CUT) {
   	  for (j = NUM_BB_ATOM; j < conf._res[l]._numAtom; j++) {
	  // skip undefined and H atoms
	    if (conf._res[l]._atom[j]._type == UNDEF || conf._res[l]._atom[j]._type >= 22)
	     continue;
            if(conf._res[l]._atom[j].x== 0 && conf._res[l]._atom[j].y== 0 && conf._res[l]._atom[j].z== 0) continue;
	    for(k = 0; k < conf._res[i]._numAtom; ++k) {
	       if(conf._res[i]._atom[k]._type==UNDEF || conf._res[i]._atom[k]._type >= 22) continue;
               if(conf._res[i]._atom[k].x== 0 && conf._res[i]._atom[k].y== 0 && conf._res[i]._atom[k].z== 0) continue;
	       disquare = (conf._res[l]._atom[j].x-conf._res[i]._atom[k].x)*(conf._res[l]._atom[j].x-conf._res[i]._atom[k].x)
	        +(conf._res[l]._atom[j].y-conf._res[i]._atom[k].y)*(conf._res[l]._atom[j].y-conf._res[i]._atom[k].y)
		+(conf._res[l]._atom[j].z-conf._res[i]._atom[k].z)*(conf._res[l]._atom[j].z-conf._res[i]._atom[k].z);
	       if(disquare <= PF_DIS_CUT_SQUARE)
	       {
	          dis = sqrt(disquare);	
                  r = Atom::radius[conf._res[l]._atom[j]._type] + Atom::radius[conf._res[i]._atom[k]._type];
		  quot = dis/r;
                  if(quot <= VDW_CLASH_CUTOFF)
  	          {
	            if(quot <= VDW_CLASH_CUTOFF && quot >= 0.5) 
	            clash_count ++;
		    else
		    clash_count += 5;
	          }
	       }
	    }		
         }
      }
     }
    }
    if(clash_count != 0)
    {
      ResIdx.push_back(l);
      ClashNum.push_back(clash_count);
    }
    clash_count = 0;
   }
}

//Steric Clash Detection for one residue
int Res_clash_detection_list(Structure& conf,  Residue& res, int Start, int End,  int* List, int List_size)
{
    int p,i,j,k, clash_count = 0;
    double disquare, dis, r, quot = 0;    
    for (p = 0; p < List_size; p++) {
       i = List[p];
       if(i == res._posn)
       {
          for (j = NUM_BB_ATOM; j < res._numAtom; j++) {
	  if (res._atom[j]._type == UNDEF || res._atom[j]._type >= 22)
	    continue;
            if(res._atom[j].x== 0 && res._atom[j].y== 0 && res._atom[j].z== 0) continue;
	    for(k = 2; k <= 3; ++k) {        
	    if(conf._res[i]._atom[k]._type==UNDEF || conf._res[i]._atom[k]._type >= 22) continue;
            if(conf._res[i]._atom[k].x== 0 && conf._res[i]._atom[k].y== 0 && conf._res[i]._atom[k].z== 0) continue;
	    disquare = (res._atom[j].x-conf._res[i]._atom[k].x)*(res._atom[j].x-conf._res[i]._atom[k].x)
	              +(res._atom[j].y-conf._res[i]._atom[k].y)*(res._atom[j].y-conf._res[i]._atom[k].y)
		      +(res._atom[j].z-conf._res[i]._atom[k].z)*(res._atom[j].z-conf._res[i]._atom[k].z);
	       if(disquare <= PF_DIS_CUT_SQUARE)
	       {
	          dis = sqrt(disquare);	
                  r = Atom::radius[res._atom[j]._type] + Atom::radius[conf._res[i]._atom[k]._type];
		  quot = dis/r;
                  if(quot <= VDW_CLASH_CUTOFF)
  	          {
	            if(quot <= VDW_CLASH_CUTOFF && quot > 0.5) 
	            clash_count ++;
		    else
		    clash_count += 5;
	          }
	       }
	  }
        }
       }
       else if (i > res._posn && i <= End) continue;
       else if (i >= Start && i < res._posn)
       {
          for (j = NUM_BB_ATOM; j < res._numAtom; j++) {
	  if (res._atom[j]._type == UNDEF || res._atom[j]._type >= 22)
	    continue;
          if(res._atom[j].x== 0 && res._atom[j].y== 0 && res._atom[j].z== 0) continue;
	  //if(type != 0 && j >= NUM_BB_ATOM ) continue;  // when sampling only backbone atoms, skip side chain atoms if there are any
	  for(k = 0; k < conf._res[i]._numAtom; ++k) {
	    if(conf._res[i]._atom[k]._type==UNDEF || conf._res[i]._atom[k]._type >= 22) continue;
            if(conf._res[i]._atom[k].x== 0 && conf._res[i]._atom[k].y== 0 && conf._res[i]._atom[k].z== 0) continue;
	    disquare = (res._atom[j].x-conf._res[i]._atom[k].x)*(res._atom[j].x-conf._res[i]._atom[k].x)
	        +(res._atom[j].y-conf._res[i]._atom[k].y)*(res._atom[j].y-conf._res[i]._atom[k].y)
		+(res._atom[j].z-conf._res[i]._atom[k].z)*(res._atom[j].z-conf._res[i]._atom[k].z);
	       if(disquare <= PF_DIS_CUT_SQUARE)
	       {
	          dis = sqrt(disquare);	
                  r = Atom::radius[res._atom[j]._type] + Atom::radius[conf._res[i]._atom[k]._type];
		  quot = dis/r;
                  if(quot <= VDW_CLASH_CUTOFF)
  	          {
	            if(quot <= VDW_CLASH_CUTOFF && quot >= 0.5) 
	            clash_count ++;
		    else
		    clash_count += 5;
	          }
	       }
	  }
        }
       }
       else
       { 
       if((res._center.dis(conf._res[i]._center) < Residue::size[res._type] + Residue::size[conf._res[i]._type] + CUB_SIZE) ||
	(res._bbc.dis(conf._res[i]._center) < Residue::bb_size + Residue::size[conf._res[i]._type] + CUB_SIZE)) {
	for (j = NUM_BB_ATOM; j < res._numAtom; j++ ) {
	  // skip undefined and H atoms
	  if (res._atom[j]._type == UNDEF || res._atom[j]._type >= 22)
	    continue;
          if(res._atom[j].x== 0 && res._atom[j].y== 0 && res._atom[j].z== 0) continue;
	  //if(type != 0 && j >= NUM_BB_ATOM ) continue;  // when sampling only backbone atoms, skip side chain atoms if there are any
	  for(k = 0; k < conf._res[i]._numAtom; ++k) {
	    if(conf._res[i]._atom[k]._type==UNDEF || conf._res[i]._atom[k]._type >= 22) continue;
            if(conf._res[i]._atom[k].x== 0 && conf._res[i]._atom[k].y== 0 && conf._res[i]._atom[k].z== 0) continue;
	    disquare = (res._atom[j].x-conf._res[i]._atom[k].x)*(res._atom[j].x-conf._res[i]._atom[k].x)
	        +(res._atom[j].y-conf._res[i]._atom[k].y)*(res._atom[j].y-conf._res[i]._atom[k].y)
		+(res._atom[j].z-conf._res[i]._atom[k].z)*(res._atom[j].z-conf._res[i]._atom[k].z);
	       if(disquare <= PF_DIS_CUT_SQUARE)
	       {
	          dis = sqrt(disquare);	
                  r = Atom::radius[res._atom[j]._type] + Atom::radius[conf._res[i]._atom[k]._type];
		  quot = dis/r;
                  if(quot <= VDW_CLASH_CUTOFF)
  	          {
	            if(quot <= VDW_CLASH_CUTOFF && quot >= 0.5) 
	            clash_count ++;
		    else
		    clash_count += 5;
	          }
	       }
	  }
        }
       }
      }
    }
   return clash_count;
}



// Backbone Clash detection
void BBClash_detection(Structure& conf, int Start, int End, vector<int>& ResIdx, vector<int>& ClashNum)
{
     int i, p,j,l,k, t, clash_count =0;
     double disquare, dis, r, quot =0;


     for(l = Start; l <= End; l++){
      if( conf._res[l]._type == 0 || conf._res[l]._type == 5)
      continue;
        for (i = 1; i <= conf._numRes; i++) {
	 if( i >= l && i <= End) continue;
     if(( conf._res[l]._center.dis(conf._res[i]._center) < Residue::size[conf._res[l]._type] + Residue::size[conf._res[i]._type] + CUB_SIZE) ||
	 (conf._res[l]._bbc.dis(conf._res[i]._center) < Residue::bb_size + Residue::size[conf._res[i]._type] + CUB_SIZE)) {
   	  for (j = 0; j < NUM_BB_ATOM; j++) {
	  // skip undefined and H atoms
	    if (conf._res[l]._atom[j]._type == UNDEF || conf._res[l]._atom[j]._type >= 22)
	     continue;
	    if(conf._res[l]._atom[j].x== 0 && conf._res[l]._atom[j].y== 0 && conf._res[l]._atom[j].z== 0) continue;
	    for(k = 0; k < conf._res[i]._numAtom; ++k) {
	       if(conf._res[i]._atom[k]._type==UNDEF || conf._res[i]._atom[k]._type >= 22) continue;
	       if(conf._res[i]._atom[k].x== 0 && conf._res[i]._atom[k].y== 0 && conf._res[i]._atom[k].z== 0) continue;
	       if(i == l+1)
	       {
	         if (j == ATM_CA && k == ATM_CA)
		    continue;
		 if (j == ATM_C &&(k == ATM_N || k == ATM_CA))
		    continue;
		 if (j == ATM_O && k == ATM_N )
		    continue;
               }else if(i == l - 1)
	       {
	         if (k == ATM_CA && j == ATM_CA)
		    continue;
		 if (k == ATM_C &&(j == ATM_N || j == ATM_CA))
		    continue;
		 if (k == ATM_O && j == ATM_N )									           continue;
	       }

	       disquare = (conf._res[l]._atom[j].x-conf._res[i]._atom[k].x)*(conf._res[l]._atom[j].x-conf._res[i]._atom[k].x)
	        +(conf._res[l]._atom[j].y-conf._res[i]._atom[k].y)*(conf._res[l]._atom[j].y-conf._res[i]._atom[k].y)
		+(conf._res[l]._atom[j].z-conf._res[i]._atom[k].z)*(conf._res[l]._atom[j].z-conf._res[i]._atom[k].z);
	       if(disquare <= PF_DIS_CUT_SQUARE)
	       {
	          dis = sqrt(disquare);	
                  r = Atom::radius[conf._res[l]._atom[j]._type] + Atom::radius[conf._res[i]._atom[k]._type];
		  quot = dis/r;
                  if(quot <= VDW_CLASH_CUTOFF)
  	          {
	            if(quot <= VDW_CLASH_CUTOFF && quot >= 0.5) 
	            clash_count ++;
		    else
		    clash_count += 5;
	          }
	       }
	    }		
         }
      }
     }
    if(clash_count != 0)
    {
      ResIdx.push_back(l);
      ClashNum.push_back(clash_count);
    }
    clash_count = 0;
   }
}

// Backbone Loop
void BBClash_detection_list(Structure& conf, int Start, int End, vector<int>& ResIdx, vector<int>& ClashNum, int* List, int List_size )
{
     int i, p,j,l,k, t, clash_count =0;
     double disquare, dis, r, quot =0;


     for(l = Start; l <= End; l++){
      if( conf._res[l]._type == 0 || conf._res[l]._type == 5)
      continue;
        for (p = 0; p < List_size; p++) {
        i = List[p];
        if( i >= l && i <= End) continue;
     if(( conf._res[l]._center.dis(conf._res[i]._center) < Residue::size[conf._res[l]._type] + Residue::size[conf._res[i]._type] + CUB_SIZE) ||
	 (conf._res[l]._bbc.dis(conf._res[i]._center) < Residue::bb_size + Residue::size[conf._res[i]._type] + CUB_SIZE)) {
   	  for (j = 0; j < NUM_BB_ATOM; j++) {
	  // skip undefined and H atoms
	    if (conf._res[l]._atom[j]._type == UNDEF || conf._res[l]._atom[j]._type >= 22)
	     continue;
	    if(conf._res[l]._atom[j].x== 0 && conf._res[l]._atom[j].y== 0 && conf._res[l]._atom[j].z== 0) continue;
	    for(k = 0; k < conf._res[i]._numAtom; ++k) {
	       if(conf._res[i]._atom[k]._type==UNDEF || conf._res[i]._atom[k]._type >= 22) continue;
	       if(conf._res[i]._atom[k].x== 0 && conf._res[i]._atom[k].y== 0 && conf._res[i]._atom[k].z== 0) continue;
	       if(i == l+1)
	       {
	         if (j == ATM_CA && k == ATM_CA)
		    continue;
		 if (j == ATM_C &&(k == ATM_N || k == ATM_CA))
		    continue;
		 if (j == ATM_O && k == ATM_N )
		    continue;
               }else if(i == l - 1)
	       {
	         if (k == ATM_CA && j == ATM_CA)
		    continue;
		 if (k == ATM_C &&(j == ATM_N || j == ATM_CA))
		    continue;
		 if (k == ATM_O && j == ATM_N )									           continue;
	       }
	       disquare = (conf._res[l]._atom[j].x-conf._res[i]._atom[k].x)*(conf._res[l]._atom[j].x-conf._res[i]._atom[k].x)
	        +(conf._res[l]._atom[j].y-conf._res[i]._atom[k].y)*(conf._res[l]._atom[j].y-conf._res[i]._atom[k].y)
		+(conf._res[l]._atom[j].z-conf._res[i]._atom[k].z)*(conf._res[l]._atom[j].z-conf._res[i]._atom[k].z);
	       if(disquare <= PF_DIS_CUT_SQUARE)
	       {
	          dis = sqrt(disquare);	
                  r = Atom::radius[conf._res[l]._atom[j]._type] + Atom::radius[conf._res[i]._atom[k]._type];
		  quot = dis/r;
                  if(quot <= VDW_CLASH_CUTOFF)
  	          {
	            if(quot <= VDW_CLASH_CUTOFF && quot >= 0.5) 
	            clash_count ++;
		    else
		    clash_count += 5;
	          }
	       }
	    }		
         }
      }
     }
    if(clash_count != 0)
    {
      ResIdx.push_back(l);
      ClashNum.push_back(clash_count);
    }
    clash_count = 0;
   }
}
