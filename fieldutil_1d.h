#ifndef _FIELDUTIL_1D_H_
#define _FIELDUTIL_1D_H_

static int FLD_XCELLS;

typedef struct coordinate coordinate;
struct coordinate {
	int x;
};

typedef struct scalarVar scalarVar;
struct scalarVar {
	float pressure;
	//float density; ignore for now to save memory
	//float temperature; not simulated at the moment
	//float phi
};

typedef struct boundaryLookup boundaryLookup;
struct boundaryLookup {
	int x;
	int lookup_Loc;
};

typedef struct boundaryCell boundaryCell;
struct boundaryCell {
	int x;
	unsigned char p_isdef;
	float p;
	unsigned char uw_isdef;
	float uw;
	unsigned char ue_isdef;
	float ue;
};


void FIELDUTIL_setup(int XscalarCenterNo);
float FLD_getU(coordinate coord);
void FLD_setU(coordinate coord, float value);
scalarVar FLD_getScalar(coordinate coord);
void FLD_setScalar(coordinate coord, scalarVar value);
int FLD_isBoundary(coordinate coord);
boundaryCell FLD_getBoundary(coordinate coord);


#endif