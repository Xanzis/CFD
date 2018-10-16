#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <math.h>
#include "nrutil.h"
#include "matutil.h"

// globals declaration

//end

//struct definitions
typedef struct ijcoord ijcoord;
struct ijcoord {
	int i;
	int j;
	int iCap;
	int jCap;
}
//end

//helper functions (indexing, retrieving values, etc)

int ij_toMat(ijcoord loc) {
	//converts an ij value to a position in the solution matrix. 
	//Since iCap and jCap define what the type of cell this function handles the fact u values index from 1 in the u-direction etc.
}

ijcoord mat_toIJ(int loc) {
	//reverses the previous function, for filling solutions back into arrays
}

float getF(ijcoord loc) {
	// retrieves F from u_field at the given location, averaging when necessary
}

float getA(ijcoord loc, isvelocity) {
	// returns the value of a at the given position
	// if isvelocicty == 1 it handles diffusion differently
	return 0;
}
//end

//differencing functions (phi + u, pressure)

//end

void setup() {
	//setting constants
	//end

	//allocating memory
	//end

	//setting up boundaries
	//end
}