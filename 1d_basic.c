#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include "nrutil.h"
#include "matutil.h"

// All units are SI

float L;
int CELLNUM;
float CELLSTEP;

float DENSITY;
float GAMMA; // Diffusive term

float V; // flow velocity: equal in this case
float *U_FIELD; // might as well keep this around
float PHIA;
float PHIB;

char VERBOSE;

void setup() {
	L = 1.0;
	CELLNUM = 25;
	CELLSTEP = 1 / (float) CELLNUM;

	DENSITY = 1;
	GAMMA = 0.1;

	V = 0.5;
	U_FIELD = malloc(sizeof(float) * (CELLNUM + 1));
	for (int i = 0; i <= CELLNUM; i++) {
		U_FIELD[i] = V;
	}

	PHIA = 1;
	PHIB = 0;

	VERBOSE = 0;
}

float max_3f(float a, float b, float c) {
	// maximimum of three inputs :/
	if (a > b) {
		if (c > a) return c;
		return a;
	}
	if (b > c) return b;
	return c;
}

float max_2f(float a, float b) {
	if (a > b) return a;
	return b;
}

float abs_f(float a) {
	if (a < 0) return -1 * a;
	return a;
}

void applydifference_central(float *ap, float *aw, float *ae, float *sp, float *su, float *boundaryVals) {
	// the as are empty vectors to be filled. Su is a vector of p source terms
	/*
	The source terms su and sp are vectors that start out with any non-boundary related
	phi sources. (sp is for sources dependent on the local phi value, su is other influx)
	If boundaries are the only outside influence on the system, sp and su should start out as 0 vectors.
	This function adds the effect of boundary influx/outflux to sp and su and then 
	incorporates those values in calculations of ap.
	This function adds on boundary effects to su because su is needed outside this function, as the
	right hand side of the phi system of equations: A . phi = su
	This method of assembling the coefficients should be easy to adapt to 2d line-by line solving if
	orthogonal boundaries are incorporated.
	*/
	float convection_e, diffusion_e, convection_w, diffusion_w;
	float con_bound, dif_bound;
	float deltaF;
	
	for (int x = 0; x < CELLNUM; x++) {
		diffusion_w = GAMMA / CELLSTEP;
		convection_w = DENSITY * U_FIELD[x];
		diffusion_e = GAMMA/CELLSTEP;
		convection_e = DENSITY * U_FIELD[x + 1];
		// Sort out aw (or if it's the westmost cell sort out western boundary)
		if (x > 0) {
			aw[x] = diffusion_w + (convection_w / 2);
		}
		else {
			// Here the two sources are constructed for convection/diffusion to the first cell (independent)
			// and c-d from the first cell to the boundary (dependent on local phi value and thus put into sp)
			con_bound = DENSITY * U_FIELD[x];
			dif_bound = 2 * (GAMMA / CELLSTEP); // twice standard because the phi value is at the edge of the boudnary, not the next scalar point
			sp[x] += - (dif_bound + con_bound);
			su[x] += (dif_bound + con_bound) * boundaryVals[0];
		}
		//Sort out ae or eastern boundary
		if (x < CELLNUM - 1) {	
			ae[x] = diffusion_e - (convection_e / 2);
		}
		else {
			con_bound = DENSITY * U_FIELD[x + 1]; // Because this is the east boundary, positive U convects phi away from cells
			dif_bound = 2 * (GAMMA / CELLSTEP);
			sp[x] += - (dif_bound - con_bound);
			su[x] += (dif_bound - con_bound) * boundaryVals[1];
		}
		deltaF = diffusion_e - diffusion_w;
		ap[x] = aw[x] + ae[x] + deltaF - sp[x];
	}
}

void applydifference_hybrid(float *ap, float *aw, float *ae, float *sp, float *su, float *boundaryVals) {
	float convection_e, diffusion_e, convection_w, diffusion_w;
	float con_bound, dif_bound;
	float temp;
	float peclet;
	float deltaF;

	for (int x = 0; x < CELLNUM; x++) {
		convection_w = DENSITY * U_FIELD[x];
		convection_e = DENSITY * U_FIELD[x + 1];
		diffusion_w = GAMMA / CELLSTEP;
		diffusion_e = GAMMA/CELLSTEP;
		// Aw or western boundary
		if (x > 0) {
			temp = diffusion_w + (convection_w / 2);
			aw[x] = max_3f(convection_w, temp, 0);
		}
		else {
			con_bound = DENSITY * U_FIELD[x];
			dif_bound = 2 * (GAMMA / CELLSTEP);
			peclet = con_bound / dif_bound;
			if ((peclet > -2) && (peclet < 2)) {
				sp[x] += - (dif_bound + con_bound);
				su[x] += (dif_bound + con_bound) * boundaryVals[0];
			}
			else {
				sp[x] += - (dif_bound + max_2f(con_bound, 0));
				su[x] += (dif_bound + max_2f(con_bound, 0)) * boundaryVals[0];
			}
		}
		// Ae or eastern boundary
		if (x < CELLNUM - 1) {
			temp = diffusion_e - (convection_e / 2);
			ae[x] = max_3f(- convection_e, temp, 0);
		}
		else {
			con_bound = DENSITY * U_FIELD[x + 1];
			dif_bound = 2 * (GAMMA / CELLSTEP);
			peclet = con_bound / dif_bound;
			if ((peclet > -2) && (peclet < 2)) {
				sp[x] += - (dif_bound - con_bound);
				su[x] += (dif_bound - con_bound) * boundaryVals[1];
			}
			else {
				sp[x] += - (dif_bound - max_2f(- con_bound, 0));
				su[x] += (dif_bound - max_2f(- con_bound, 0)) * boundaryVals[1];
			}
		}
		deltaF = diffusion_e - diffusion_w;
		ap[x] = aw[x] + ae[x] + deltaF - sp[x];
	}
}

int main() {
	setup();

	float *ap = vector(0, CELLNUM - 1);
	float *aw = vector(0, CELLNUM - 1);
	float *ae = vector(0, CELLNUM - 1);
	MAT_zerovector(ap, CELLNUM);
	MAT_zerovector(aw, CELLNUM);
	MAT_zerovector(ae, CELLNUM);

	float *sp = vector(0, CELLNUM - 1);
	MAT_zerovector(sp, CELLNUM);
	float *su = vector(0, CELLNUM - 1);
	MAT_zerovector(su, CELLNUM);

	float bounds[2] = {PHIA, PHIB};

	//Put any exterior sources in su here
	//su[10] = 2;

	//applydifference_central(ap, aw, ae, sp, su, bounds);
	applydifference_hybrid(ap, aw, ae, sp, su, bounds);

	//Construct the solution matrix
	float **coeffs = matrix(0, CELLNUM - 1, 0, CELLNUM - 1);

	for (int x = 0; x < CELLNUM; x++) {
		if (x > 0) {
			coeffs[x][x-1] = -1 * aw[x];
		}
		if (x < CELLNUM - 1) {
			coeffs[x][x+1] = -1 * ae[x];
		}
		coeffs[x][x] = ap[x];
	}

	// Equation to solve is: 
	// Coeff . phi = su
	if (VERBOSE) {
		printf("Coefficients:\n");
		MAT_printmatrix(coeffs, CELLNUM, CELLNUM);
		printf("Corresponding b vector (su):\n");
		MAT_printvector(su, CELLNUM);
	}

	float *phi = MAT_solve_gausselim(coeffs, su, CELLNUM);

	if (VERBOSE) {
		printf("Triangularized matrix:\n");
		MAT_printmatrix(coeffs, CELLNUM, CELLNUM);
		printf("Corresponding b vector:\n");
		MAT_printvector(su, CELLNUM);
	}

	printf("Phi solution:\n");
	MAT_printvector(phi, CELLNUM);

	MAT_graphvector(phi, CELLNUM);

	return 1;
}