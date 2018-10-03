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
	CELLNUM = 100;
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

void applydifference_central(float *ap, float *aw, float *ae, float *sp) {
	// the as are empty vectors to be filled. Su is a vector of p source terms
	float convection_e, diffusion_e, convection_w, diffusion_w;
	
	for (int x = 0; x < CELLNUM; x++) {
		if (x > 0) {
			diffusion_w = GAMMA / CELLSTEP;
			convection_w = DENSITY * U_FIELD[x];
			aw[x] = diffusion_w + (convection_w / 2);
		}
		if (x < CELLNUM - 1) {
			diffusion_e = GAMMA/CELLSTEP;
			convection_e = DENSITY * U_FIELD[x + 1];
			ae[x] = diffusion_e - (convection_e / 2);
		}
		ap[x] = aw[x] + ae[x] - sp[x];
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

	float diffusion_far_e = GAMMA / CELLSTEP;
	float convection_far_e = DENSITY * U_FIELD[0];
	float diffusion_far_w = GAMMA / CELLSTEP;
	float convection_far_w = DENSITY * U_FIELD[CELLNUM];

	sp[0] = -(2 * diffusion_far_e + convection_far_e);
	sp[CELLNUM - 1] = -(2 * diffusion_far_w - convection_far_w);

	applydifference_central(ap, aw, ae, sp);

	float *su = vector(0, CELLNUM - 1);
	MAT_zerovector(su, CELLNUM);
	su[0] = (2 * diffusion_far_e + convection_far_e) * PHIA;
	su[CELLNUM - 1] = (2 * diffusion_far_w + convection_far_w) * PHIB;

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