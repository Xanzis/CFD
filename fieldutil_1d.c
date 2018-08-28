#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include "fieldutil_1d.h"


typedef uint_fast32_t  word;
#define WORD_SIZE_BITS (int) (sizeof(word) * 8)
#define WORD_SIZE_BYTES (int) (sizeof(word))

#define F_END 1

float **u_field;
word *is_boundary; // Will be 1 dimenional even for the 2 dimensional case. Uses memory more efficiently and I only have to interact with it in one function anyway
boundaryCell *boundaries;
int IS_BOUNDARY_ROWLEN;

static int FLD_XCELLS;
static int NUMBEROFBOUNDARIES;

scalarVar **scalar_field;

int didSetup = 0;

/*
Things to add - ideally enough to make the actual programming a breeze. (!!!)
This'll mean I don't have to think about this damn stuff when I'm writing code
Plus it's something nice and testable

Setting u field values
Retrieving u field values
Getting binary 'is a boundary here' answers
Have a .h that defines coordinate, boundary structs etc.
Getting boundary info via returning those structs
That General Sort of Thing
*/

void setupError() {
	fprintf(stderr, "Variables referenced before setup\n");
	exit(1);
}

void genericError(char msg[]) {
	fprintf(stderr, "%s\n", msg);
	exit(1);
}

void setupScalarField() {
	int xRange = FLD_XCELLS;
	int yRange = 1;

	scalar_field = malloc(sizeof(scalarVar *) * yRange);
	//if (!scalar_field) genericError("allocation failure 1 in scalar field");
	scalar_field[0] = malloc(sizeof(scalarVar) * yRange * xRange);
	//if (!scalar_field[0]) genericError("allocation failure 2 in scalar field");
	for (int i = 1; i < yRange; i++) scalar_field[i] = scalar_field[i - 1] + xRange;

	//Will require outside initializiation to real values. This just initializes the pointers to rows and columns
}

void setupUField() {
	int xRange = FLD_XCELLS + 1;
	int yRange = 1;

	u_field = (float **) malloc(sizeof(float *) * yRange);
	//if (!u_field) genericError("allocation failure 1 in u field");
	u_field[0] = malloc(sizeof(float) * yRange * xRange);
	//if (!u_field[0]) genericError("allocation failure 2 in u field");
	for (int i = 1; i < yRange; i++) u_field[i] = u_field[i - 1] + xRange;

}

void addIsBoundary(coordinate coord) {
	int x = coord.x + 1; // because -1 is a valid coordinate for an edge past 0.
	int y = 0;
	if ((x >= (FLD_XCELLS + 2) * WORD_SIZE_BITS) || (x < 0)) {
		genericError("addIsBoundary input out of range");
	}

	int place = (y * IS_BOUNDARY_ROWLEN) + x;
	int wordNo = place / WORD_SIZE_BITS;
	int wordPlace = place % WORD_SIZE_BITS;

	is_boundary[wordNo] |= (word) 1 << wordPlace;
}

void setupBoundary() {
	//Add cushion for a boundary past each side of the scalar centers
	int xRange = FLD_XCELLS + 2;
	int yRange = 1;

	IS_BOUNDARY_ROWLEN = xRange;

	int wordsNecessary = ((yRange * xRange) / WORD_SIZE_BITS) + 1; // buffer bc rounding
	is_boundary = malloc(wordsNecessary * sizeof(word *));
	memset(is_boundary, 0, wordsNecessary * WORD_SIZE_BYTES);

	// ENTER BOUNDARIES HERE. Automate (and add bodies function) later.
	NUMBEROFBOUNDARIES = 2;
	boundaries = malloc(sizeof(boundaryCell) * NUMBEROFBOUNDARIES);

	boundaries[0] = (boundaryCell) {-1, 0, 0, 0, 0, 1, 7};
	boundaries[1] = (boundaryCell) {FLD_XCELLS, 1, 50, 0, 0, 0, 0}; // figure out number for the case to test later

	coordinate temp_coord = {-1};
	addIsBoundary(temp_coord);
	temp_coord.x = FLD_XCELLS;
	addIsBoundary(temp_coord);
	//END
}

void FIELDUTIL_setup(int XscalarCenterNo) {
	FLD_XCELLS = XscalarCenterNo;

	setupUField();
	setupScalarField();
	setupBoundary();
	didSetup = 1;
}

float FLD_getU(coordinate coord) {
	if (!didSetup) {
		setupError();
	}
	//Standard coord system (iJ in this case) (takes a coord as input - that's a struct that's now defined in the .h)
	int x = coord.x;
	int y = 0;

	if ((x >= FLD_XCELLS + 1) || (x < 0)) {
		genericError("FLD_getU input out of range");
	}

	return u_field[y][x];
}

void FLD_setU(coordinate coord, float value) {
	if (!didSetup) {
		setupError();
	}
	//Standard coord system (iJ in this case) (takes a coord as input - that's a struct that's now defined in the .h)
	int x = coord.x;
	int y = 0;

	if ((x >= FLD_XCELLS + 1) || (x < 0)) {
		genericError("FLD_setU input out of range");
	}

	u_field[y][x] = value;
}

scalarVar FLD_getScalar(coordinate coord) {
	if (!didSetup) {
		setupError();
	}

	int x = coord.x;
	int y = 0;

	if ((x >= FLD_XCELLS) || (x < 0)) {
		genericError("FLD_getScalar input out of range");
	}

	return scalar_field[y][x];
}

void FLD_setScalar(coordinate coord, scalarVar value) {
	if (!didSetup) {
		setupError();
	}

	int x = coord.x;
	int y = 0;

	if ((x >= FLD_XCELLS) || (x < 0)) {
		genericError("FLD_getScalar input out of range");
	}

	scalar_field[y][x] = value;
}

int FLD_isBoundary(coordinate coord) {
	if (!didSetup) {
		setupError();
	}
	int x = coord.x + 1; // offset because a -1 input is valid. See setup re: boundaries past scalar centers
	int y = 0;

	if ((x >= FLD_XCELLS + 2) || (x < 0)) {
		genericError("FLD_isBoundary input out of range");
	}

	int place = (y * IS_BOUNDARY_ROWLEN) + x;
	int wordNo = place / WORD_SIZE_BITS;
	int wordPlace = place % WORD_SIZE_BITS;
	word theWord = is_boundary[wordNo];
	return (int) (theWord >> wordPlace) & 1;
}

boundaryCell FLD_getBoundary(coordinate coord) {
	if (!didSetup) {
		setupError();
	}

	int x = coord.x;
	if ((x >= FLD_XCELLS + 2) || (x < -1)) {
		genericError("FLD_getBoundary input out of range");
	}
	//search (via wikipedia) (which was found via search lol)
	int m = 0, L = 0, R = NUMBEROFBOUNDARIES - 1;
	boundaryCell currentBoundary;
	int currentX;
	while (L <= R) {
		m = (L + R) / 2;
		currentBoundary = boundaries[m];
		currentX = currentBoundary.x;
		if (currentX < x) {
			L = m + 1;
		}
		else if (currentX > x) {
			R = m - 1;
		}
		else {
			return currentBoundary;
		}
	}
	// Failure 
	genericError("FLD_getBoundary search failed");
}
