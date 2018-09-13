// sample_angles.h

#ifndef _SAMPLE_ANGLES_
#define _SAMPLE_ANGLES_

#include "structure.h"
#include "rotamer.h"
#include "reprst.h"
#include "util.h"

void sample_bb_angles(Structure& conf, int position, double(*angles)[3],
		      int numStates, int);

void sample_sc_angles(Residue& res, double(*angles)[6],
		      int numStates, int);

#endif
