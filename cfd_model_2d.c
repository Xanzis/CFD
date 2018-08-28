#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/*
Using the book for most of the math
but the pressure correction equations are from http://web.cecs.pdx.edu/~gerry/class/ME448/notes/pdf/SIMPLEslides.pdf
- actually it looks like the giant matrix inversion might be the only answer

This is an implementation of the SIMPLE method for incompressible steady-state flow. This particular program is built
to handle a 2-dimensional example but will develop tools for use in a later 3d version. 

All units are SI - kg, m, s, K, etc.
Note to future self - do not be afraid to use linear algebra, especially for the row by row a calculations (and eventually pressure correction)
*/

// define a coordinate struct
typedef struct coordinate coordinate;
struct coordinate {
	int x;
	int y;
}

typedef struct scalarVars scalarVars;
struct scalarVars {
	float pressure;
	//float density; ignore for now to save memory
	//float temperature; not simulated at the moment
	//float phi
}

typedef struct boundaryLookup boundaryLookup;
struct boundaryLookup {
	int x;
	int y;
	int lookup_Loc;
}

typedef struct boundaryCell boundaryCell;
struct boundaryCell {
	int x;
	int y;
	/* insert relevant boundary condition definitions here. For example, pressure, velocity on faces 
	(or also whether to use the velocity already in the velocity field). That last one will probably
	be important for boundaries where pressure is stated but one wants to calculate velocity etc.
	*/
	// just a guess at what might be necessary
	unsigned char p_isdef; //probably the smallest data type
	float p;
	unsigned char uw_isdef;
	float uw;
	unsigned char ue_isdef;
	float ue;
	unsigned char vb_isdef;
	float vb;
	unsigned char vn_isdef;
	float vn;
}

double CELLSTEP; // dx or dy. sticking to square cells for now.
int CELLCOUNT_X;
int CELLCOUNT_Y;

/*initialize pointers to what will be the u and v fields. They should be roughly the same dimensions as the scalar var
field, but there might be an extra v or u value at each end. This depends on how boundary conditions will be handled.
*/
float *v_field;
float *u_field;

scalarVars *scalar_field;

double DENSITY;
double DIFFUSION; //Diffusion coefficient. Constant for now. Also: looks like no viscosity.

boundaryLookup *boundary_l; //Just a guess for how to do this. Will contain sorted list of boundary points
boundaryCell *boundaries;

void setupBoundaries() {
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
	CELLCOUNT_Y = 5;

	scalar_field = (scalarVars *) malloc(CELLCOUNT_X * CELLCOUNT_Y);
	u_field = (float *) malloc((CELLCOUNT_X + 1) * CELLCOUNT_Y);
	v_field = (float *) malloc(CELLCOUNT_X * (CELLCOUNT_Y + 1));

	DENSITY = 1000; //water
	DIFFUSION = 0; // fix this later, I don't even know

	setupBoundaries();
}

int findBoundaryCell(coordinate scalarCenterCoord) {
	/*
	binary search for coordinate

	if it has found a boundary:
		return number of that boundary cell for use with array boundaries
	else:
		return -1
	*/
}

float findU(coordinate coord) {
	/* Takes a coordinate in the form of the I, J notation used in the book and returns the corresponding u field value.
	To consider:
	- Make sure the coordinate system is rock solid.
	- Is this u constrained by a boundary?
		- need totally certain way of relating face coordinates to cel coords for boundary math.
		- Remember to check to both sides :) could be west of a boundary or east of a boundary to the left.
		- if it is, returns value from constraint, I guess. Look at use cases for this function
			and compare with the book's recommended boundary handling.
	*/

	int i = coord.x;
	int J = coord.y;

	//check if value is within the u field.
	if((0 <= i < (CELLCOUNT_X + 1)) && (0 <= J < CELLCOUNT_Y)) {
		return u_field[(J * CELLCOUNT_Y) + i];
		// Nope, gonna have to also check if there is a boundary that influences this u. Means checking scalar centers xy and x-1, y
	}
	else {
		if(!(0 <= J <= CELLCOUNT_Y)) {
			printf("Something's gone terribly wrong (or maybe it's an issue with the four point averages). U value outside y range called.\n");
			// actually check on those four point averages. That might need to be a thing.
		}
		// if execution reaches here it's too far to the left or right and values need to come from boundary conditions.
		//some handling here. Prob need to check boundaries
	}

	//

	return u_field[0];
}

float findV(I, j) {
	// Same as above.
	return v_field[0];
}

float findP(J, I) {
	// essentially the same, but the value is stored at a scalar variable. 
}

float applyDifferenceScheme(char pointType, coordinate coord) {
	/*
	Using the QUICK difference scheme, returns the value ->a<- used in later calculations.
	Pay close attention to treatment of boundaries.
	pointType specifies the type of point at which a is desired. 
	1 is scalar (I, J)
	2 is u (i, J)
	3 is v (I, j)
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

// And so on. It'll get worse from here, but these functions should be very useful.
