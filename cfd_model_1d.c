#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "nrutil.h"
#include "fieldutil_1d.h"

/*
Using the book for most of the math
but the pressure correction equations are from http://web.cecs.pdx.edu/~gerry/class/ME448/notes/pdf/SIMPLEslides.pdf
- actually it looks like the giant matrix inversion might be the only answer

This is an implementation of the SIMPLE method for incompressible steady-state flow. This particular program is built
to handle a 2-dimensional example but will develop tools for use in a later 3d version. 

All units are SI - kg, m, s, K, etc.
Note to future self - do not be afraid to use linear algebra, especially for the row by row a calculations (and eventually pressure correction)
*/



double CELLSTEP; // dx or dy. sticking to square cells for now.
int CELLCOUNT_X;

float *u_field;

scalarVars *scalar_field;

double DENSITY;
double GAMMA; //Diffusion coefficient. Constant for now. Also: looks like no viscosity.

boundaryLookup *boundary_l; //Just a guess for how to do this. Will contain sorted list of boundary points
boundaryCell *boundaries;
int NUMBEROFBOUNDARIES;

void setupBoundaries() {
	NUMBEROFBOUNDARIES = 2;

	boundaries = (boundaryCell *) malloc(NUMBEROFBOUNDARIES);
	boundary_l = (boundaryLookup *) malloc(NUMBEROFBOUNDARIES);

	boundaries[0] = (boundaryCell) {-1, 0, 0, 0, 0, 1, 5}; //change later based on the test case
	boundaries[1] = (boundaryCell) {5, 0, 0, 1, 3, 0, 0}; // sepcifically pressure is important
	// note that both the defined velocites are already part of the velocity field. Something to consider somewhere in this program

	boundary_l[0] = (boundaryLookup) {-1, 0};
	boundary_l[1] = (boundaryLookup) {5, 1};

	/*
	Need to:
	- Create the outer bounding box (maybe depending on what type of outer boundary is desired)
	- Create any interior features
	- Create any inlets or outlets
	- Sort (Yay writing a sort in C!) the resulting boundaries by coordinate
	- initalize and populate the sorted lookup table

	Note to future self: I think the chosen method of doing boundaries is good. Make sure it's handled
	consistently (Particularly regarding the fact that the boundaries will both contain scalar values and 
	face values) I think it will be necessary to have a way of indicating whether a boundary aspect is 
	defined.

	To override any other corrections, boundary conditions will probably have to be reapplied 
	somewhere during solution. Remember this.

	Final thing about boundaries: how to change whatever matrix thing goes on with the pressure correction
	so that the region we don't care about don't factor in? probably just have populating matrix/u fields 
	switch on or off whenever the program scans over a boundary cell. A little hacky but might work. Hmm, 
	but that doesn't totally cover vertical lines etc. Oh well, can deal with that once fluid dynamics in
	a basic cube works.
	*/
}

void setup() {
	CELLSTEP = 0.1;
	CELLCOUNT_X = 5;

	scalar_field = (scalarVars *) malloc(CELLCOUNT_X);
	u_field = (float *) malloc(CELLCOUNT_X + 1);

	DENSITY = 1000; //water
	GAMMA = 0; // fix this later, I don't even know

	setupBoundaries();
}

int findBoundaryCell(coordinate scalarCenterCoord) {
	/*
	binary search for coordinate

	if it has found a boundary:
		return number of that boundary cell for use with array boundaries
	else:
		return -1

	for now it'll just go through all the boundaries. I can write a sort/search later
	*/
	for (int i = 0; i < NUMBEROFBOUNDARIES; i++) {
		if (boundary_l[i].x == scalarCenterCoord.x) {
			return i;
		}
	}
	return -1;
}

int findNextBoundaryCell(coordinate centerCoord, int xDirection) {
	coordinate currentCoord = centerCoord;
}

float findU(coordinate coord) {
	/* Takes a coordinate in the form of the I, J notation used in the book and returns the corresponding u field value.
	To consider:
	- Make sure the coordinate system is rock solid.
	*/

	coordinate westCell = {coord.x - 1};
	int wBound = findBoundaryCell(westCell);
	coordinate eastCell = {coord.x};
	int eBound = findBoundaryCell(eastCell);

	boundaryCell bound = {0, 0, 0, 0, 0, 0, 0};

	if (wBound != -1) {
		boundaryCell = boundaries[wBound];
		if (boundaryCell.ue_isdef) {
			return boundaryCell.ue;
		}
	}

	if (eBound != -1) {
		boundaryCell = boundaries[eBound];
		if (boundaryCell.uw_isdef) {
			return boundaryCell.uw;
		}
	}

	//No need to worry about boundaries past here

	int i = coord.x;

	//check if value is within the u field.
	if(0 <= i < (CELLCOUNT_X + 1)) {
		return u_field[i]; // for now it's nice and one-dimensional :)
	}
	else {
		printf("Definitely shouldn't have reached here. coord is out of bounds of u-field and isn't defined by a boundary")
		exit(0);
	}
}

float findP(J, I) {
	return 0; //for now
	// essentially the same, but the value is stored at a scalar variable. 
}

float applyCentralDifferenceSchemeFace(coordinate pointcoord, int dirX) {
	// dirX gives the direction for the face. 0 is W, 1 is E
	// In higher dimensional models, [dirX, dirY, dirZ] will be a one-hot array
	//format doesn't matter for CD but QUICK will need forward references to cells
	//IMPORTANT - this code currently has no boundary protections. Assuming this will be handled earlier
	
	float massFlux = u_field[pointcoord.x + dirX] * DENSITY;
	float convectionTerm = massFlux * 0.5 * (dirX ? -1 : 1); // inverts if using east face. See CD definition

	float diffusionTerm = 2 * GAMMA/CELLSTEP; // multiply by 2 bc it's distance to boundary not next point

	return diffusionTerm + massFlux
}

float applyQUICKDifferenceScheme(char pointType, coordinate coord) {
	/*
	Using the QUICK difference scheme, returns the value ->a<- used in later calculations.
	Pay close attention to treatment of boundaries.
	for now this only has to deal with u-values. A future version should consider how to also handle v-values (or just be agnostic)
	they might all have to be treated seperately...
	*/
	float a;

	switch(pointType) {
		case 1:
			//It's a scalar
			a = 0;
			break;

		case 2:
			// It's a u value
			a = 0;
			break;

		case 3:
			// it's a v value
			a = 0;
			break;
	}

	return a;
}

void runLine() {
	int reachedEnd = 0;
	int xSectionBegin = 0;
	int xSectionEnd = 0;

	while(!reachedEnd) {
		xSectionEnd = xSectionBegin
		// Find end x location for the uninterrupted group.

		//Then, populate the stuff

		//Then, either set the beginning of the next block, or reachedEnd
		//Beginning of next block must be a boundary with a defined east face
	}
}















