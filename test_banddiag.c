#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include "nrutil.h"
#include "matutil.h"

int main() {
	int maxval = 10;
	float num;
	int N = 100;
	float **mat = matrix(0, N - 1, 0, N - 1);
	float *vec = vector(0, N-1);

	// Start off with just 1d-type bands - adjacent
	for (int i = 0; i < N; i++) {
		for (int band = -1; band <= 1; band++) {
			if (i + band >= 0 && i + band < N) {
				num = (float) rand() / (float) RAND_MAX;
				mat[i][i + band] = num * maxval; // negative numbers shouldn't reallyt be necessary
			}
		}
		num = (float) rand() / (float) RAND_MAX;
		vec[i] = num * maxval;
	}

	//MAT_printmatrix(mat, N, N);
	//MAT_printvector(vec, N);

	float *res = MAT_solve_gausselim(mat, vec, N);
	MAT_printvector(vec, N);
}