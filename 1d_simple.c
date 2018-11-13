#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <math.h>
#include "nrutil.h"
#include "matutil.h"

float L;
float AREA;
int N; // number of scalar cells. there are one fewer velocity cells.
float CELLSTEP;

float DENSITY;
float GAMMA; // Diffusive term

float *U_FIELD;
float *U_FIELD_STAR; // U field new values
float *U_FIELD_PRIME; // U field correction

float *P_FIELD;
float *P_FIELD_PRIME; // P field correction

float BIGNUM;
float SMALLNUM;
float ALPHA; // under-relaxation factor

int VERBOSE;

typedef struct icoord icoord;
struct icoord {
	int i;
	int iCap;
}

int i_toMat(icoord loc) {
	// If iCap is 1, don't change. If iCap is 0, subtract one
	int shift = 1 - loc.iCap;
	return loc.i - shift;
}

icoord mat_toI(int loc, int iCap) {
	int shift = 1 - iCap;
	return (icoord) {loc + shift, iCap};
}

void simerror(char *error_text) {
	printf("Top level simulation error. Error message follows:\n%s", error_text);
	printf("\nExiting ...\n");
	exit(1);
}

void setupBoundaries() {
	// intializes P and U values at relevant inletsand outlets to what they will stay at for differencing. A different function is in charge of cancelling changes during differencing
	U_FIELD[1] = 10; // inlet velocity = 10 m/s
	P_FIELD[0] = 10; //inlet pressure = 10Pa
	P_FIELD[N - 1] = 0; // outlet pressure = 0Pa
	// That's all the constants in this 1d case. For 2d, adapt to be procedural 
}

int validUIDX(int idx) {
	// returns 1 if idx is a valid index of U_FIELD
	if (1 <= idx && idx < N) {
		return 1;
	}
	return 0;
}

float getF(icoord loc, int w, int e) {
	// Gets the value of F from the coord in the given direction, averaging when necessary
	// If the value is over the edge of the variable array, return 0. It won't matter, as edge variables will be overwritten by boundary conditions
	float res = 0;
	int idx = 0;
	float temp;
	switch (loc.iCap) {
		case 0:
			//It's a u center
			if (w) {
				idx = loc.i - 1;
				if (validUIDX(idx)) {
					temp = U_FIELD[idx] * 0.5; // Need to follow this up with density averaging when I move to compressible flows
					idx += 1;
				}
				else {
					return 0;
				}
				if (validUIDX(idx)) {
					temp += U_FIELD[idx] * 0.5
				}
				else {
					return 0;
				}
			}
			else if (e) {
				idx = loc.i;
				if (validUIDX(idx)) {
					temp = U_FIELD[idx] * 0.5; // Need to follow this up with density averaging when I move to compressible flows
					idx += 1;
				}
				else {
					return 0;
				}
				if (validUIDX(idx)) {
					temp += U_FIELD[idx] * 0.5
				}
				else {
					return 0;
				}
			}
			break;
		case 1:
			// It's a center point
			if (w) {
				idx = loc.i;
				if (validUIDX(idx)) {
					res = U_FIELD[idx]; // Needed to check if it was within array bounds first. If not it'll just return 0
				}
			}
			else if (e) {
				idx = loc.i + 1;
				if (validUIDX(idx)) {
					res = U_FIELD[idx];
				}
			}
			break;
		default:
			simerror("Invalid iCap value in getF");
			return 0;
	}
	res *= DENSITY; //looks like area doesn't actually go in here, interestingly enough
}

float getD(icoord loc, int w, int e) {
	// Ditto
	return GAMMA / CELLSTEP; // This will need rethinking quite soon, but for now it's best to get the system running
	// by the way I need to figure out changing gamma for momentum differencing. Look at k-E theory soon, keep viscosity at 0 for now.
}

float ap_hybrid(icoord loc) {
	//Applies the hybrid difference scheme to get ap at a single point loc.

	Fw = getF(loc, 1, 0);
	Fe = getF(loc, 0, 1);
	Dw = getD(loc, 1, 0);
	De = getD(loc, 0, 1);

	temp = Dw + (Fw / 2);
	float w = max_3f(Fw, temp, 0);

	temp = De - (Fe / 2);
	float e = max_3f(- Fe, temp, 0);

	deltaF = De - Dw;
	return w[i] + e[i] + deltaF; // ap
}

void applydifference_hybrid(p, w, e, iCap) {
	float Fw, Fe, Dw, De;
	float temp, peclet;
	// Starts at 1-iCap because if I is capitalized it's P, which is 0-indexed
	for (int i = 1 - iCap; i <= N; i++) {
		icoord loc = {i, iCap};

		Fw = getF(loc, 1, 0);
		Fe = getF(loc, 0, 1);
		Dw = getD(loc, 1, 0);
		De = getD(loc, 0, 1);

		temp = Dw + (Fw / 2);
		w[i] = max_3f(Fw, temp, 0);

		temp = De - (Fe / 2);
		e[i] = max_3f(- Fe, temp, 0);

		deltaF = De - Dw;
		p[i] = w[i] + e[i] + deltaF; // since sp hasn't been calculated it still needs to be subtracted
	}
	// This incarnation of the function is source-agnostic. All boundary-related source terms should be added after calling this.
	// If I is capitalized, it's a scalar center. If not, it's u

	//Important! Need separate length as for scalar vs. velocity stuff.
	//Or maybe just zero in between runs and feed the solver smaller dimensions
}

void differencePressure(su) {
	// See p186 of the book. During momentum differencing, pressure difference terms are added onto su
	// iterating 1-N through u. But! For differencing, everything needs to shifted back to 0. So u a vectors run from 0 to N-1
	for (int i = 0; i < N - 1; i++) {
		su[i] += (P_FIELD[i] - P_FIELD[i + 1]) * AREA;
	}
}

void initializeUField() {
	//stuff
	//Maybe Gaussian noise?
	float seq[10] = {3, 5, 7, 2, 4, 6, 3, 4, 6, 8};
	int idx = 0;

	for (int i = 1; i < N; i++) {
		if (idx <= 10) {
			idx = 0;
		}
		U_FIELD[i] = seq[idx];
		idx++;
	}
	// Boundary value setting is doen by a setup function
}

void initializePField() {
	// also stuff
	float seq[10] = {2, 6, 3, 8, 6, 4, 5, 3, 2, 9};
	int idx = 0;

	for (int i = 1; i < N; i++) {
		if (idx <= 10) {
			idx = 0;
		}
		U_FIELD[i] = seq[idx];
		idx++;
	}
}

// Following three function apply different types of boundaries: pressure, velocity, and all other generic scalars
// In the future 2d version, these will end up reading from file to determine the boundary structure;
// this file can then be written by a python program, letting me do the higher-level work away from C.

void applyPressureBoundaries(float *su, float *sp) {
	//Pressure corrections at inlet and outlet boundaries are 0
	su[0] = 0;
	sp[0] = SMALLNUM;
	su[N-1] = 0;
}

void applyVelocityBoundaries(float *su, float *sp) {
	// Applies sources to velocity matrix from walls and inlets.
	// Check indexing to make sure indexing works out with inputs to solveSystem.
}

void applyScalarBoundaries(float *su, float *sp, float scalarNum) {
	// Empty for now. scalarNum gives ID of the scalar variable function is applied to.
}

void setup() {
	BIGNUM = (float) pow(10, 20);
	SMALLNUM = - BIGNUM;

	ALPHA = 0.1;

	L = 1.0;
	N = 5;
	CELLSTEP = L / (float) CELLNUM; //Fine but because of boundary overhang this won't actually be the sim length
	AREA = CELLSTEP * CELLSTEP;

	DENSITY = 1;
	GAMMA = 0.1;

	U_FIELD = vector(1, N);
	U_FIELD_PRIME = vector(1, N); //Figure out how to zero these properly given that zerovector only works with zero-indexed.
	U_FIELD_STAR = vector(1, N);
	// Maybe MAT_zerovector(U_FIELD + 1, N)?

	P_FIELD = vector(0, N);
	P_FIELD_PRIME = vector(0, N);

	MAT_zerovector_range(U_FIELD, 1, N);
	MAT_zerovector_range(U_FIELD_PRIME, 1, N);
	MAT_zerovector_range(U_FIELD_STAR, 1, N);

	MAT_zerovector(P_FIELD, N + 1); //0 to N is N+1 elements
	MAT_zerovector(P_FIELD_PRIME, N + 1);

	VERBOSE = 1;

	initializeUField();
	initializePField();
	setupBoundaries();
}

void solveSystem(float *target, float *cp, float *cw, float *ce, float *sp, float *su, int target_isoneindexed, int n_var) {
	// Solves the form of equation that comes up a few times, where (cp[n] - sp[n])phi[n] = sum(c[nb]phi[nb]) + su
	// or rather (cp[n] - sp[n])phi[n] - sum(c[nb]phi[nb]) = su
	// For n_var variables. I know, this is available from the environment conditions, but ideally this function is nice and agnostic.
	// Places resulting solved phi vector into target. 
	// For most transport work - Fw, Fe differential needs to be incorporated into ap 

	float **coeff = matrix(0, n_var - 1, 0, n_var - 1);
	float *rhs = vector(0, n_var - 1); // right hand side of the equation coeff . phi = rhs

	int i_ref = 0;
	int w_idx, e_idx;
	for (int i = 0; i < n_var; i++) {
		i_ref = i + target_isoneindexed;
		// Hopefully this should expand well to 2d. 
		w_idx = i_ref - 1;
		e_idx = i_ref + 1;

		rhs[i] = su[i_ref];
		coeff[i][i] = cp[i_ref] - sp[i_ref];
		if (i > 0) {
			coeff[i][i - 1] = - cw[i_ref];
		}
		if (i < n_var - 1) {
			coeff[i][i + 1] = - ce[i_ref];
		}
	}

	float *temp = MAT_solve_gausselim(coeff, rhs, n_var)
	for (int i = 0; i < n_var; i++) {
		i_ref = i + target_isoneindexed;
		target[i_ref] = temp[i];
	}
	free_vector(temp, 0, n_var - 1);
	free_vector(rhs, 0, n_var - 1);
	free_matrix(coeff, 0, n_var - 1, 0, n_var - 1); // freeing the pointers which aren't temp may not be necessary
	// The above mess was done because I don't know how c memory management works. Is it as simple as target = temp - isoneindexed? Help
}

int main() {
	setup();
	float *ap, *aw, *ae, *su, *sp; //same a vectors for each differencing operation, but U stuff is 1-indexed. Remember that
	sp = vector(0, N);
	su = vector(0, N);
	ap = vector(0, N);
	aw = vector(0, N);
	ae = vector(0, N);
	// setup has already initialized guesses of p and u. Begin convergence loop
	int notconverged = 1
	while (notconverged) {
		// 1. Solve momentum equations
		// 1a momentum differencing
		applydifference_hybrid(ap, aw, ae, 0); // lowercase i for inter scalar-cell values. 
				// Eventually split this off the generic differencing function to properly treat velocity diffusion different from mass diffusion.
				//  -> -> -> Check that indexing is the same as in solvesystem!
		// 1b applying momentum boundary constraints
		applyVelocityBoundaries(su, sp);
		// 1c solve equations, feeding results into u star- this gives the 'initial guess' of velocities from the stated pressure field
		solveSystem(U_FIELD_STAR, ap, aw, ae, sp, su, 0, N); // there should be N velocity elements, and u cells are not i-capitalized

		MAT_zerovector(ap, N+1);
		MAT_zerovector(aw, N+1);
		MAT_zerovector(ae, N+1);
		MAT_zerovector(sp, N+1);
		MAT_zerovector(su, N+1);
		// Assuming that from here forward it should still usethe 

		// 2. Solve pressure correction equation
		// 2a pressure differencing -- the differencing here is quite different from standard
		// 2b apply pressure differencing boundary influences
		// 2c solve equations, feeding results into pressure correction matrix

		// 3. Apply U, P correction
		// 3a add solved correction values into u and p matrices
		// note: see p189 for details on under-relaxation - it's not what you think

		// 4. Solve all other convection-diffusion problems
		// 4a nothing for this so far - but may want a generic convection-diffuision program for scalar variables

		// 5. Decide whether or not the solution has converged
	}

}
//Keep at it!!!
// Store edge boundaries separately from internal stuff. Should make indexing those simple cases (in a consistent way too) easy.
// Also remember - each edge boundary needs to account for all the possible edge influences
// see book's list of boundary types
// In differencing need to check what kind of boundary everything is -- NOPE! boundary-related source should occur outside of differencing

// By the way, the way the differencing works right now it's not necessarily cutting off links to boundaries - and if it does it's accidental
// This means it may be necessary to zero out some aw, ae etc. when post processing the returned as before solving.