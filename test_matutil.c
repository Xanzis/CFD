#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include "nrutil.h"
#include "matutil.h"

int main() {
	float testmat_def[9] = {22, 15, 3, -3, -1, 0, -2, 1, 0};
	float **testmatrix = matrix(0, 2, 0, 2);
	float **backupmat = matrix(0, 2, 0, 2);
	for (int j = 0; j < 3; j++) {
		for (int i = 0; i < 3; i++) {
			testmatrix[j][i] = testmat_def[i + (3 * j)];
			backupmat[j][i] = testmat_def[i + (3 * j)];
		}
	}

	float testvec_def[3] = {8, -11, 3};
	float *testvec = vector(0, 2);
	float *backupvec = vector(0, 2);
	for (int i = 0; i < 3; i++) {
		testvec[i] = testvec_def[i];
		backupvec[i] = testvec_def[i];
	}

	printf("Here goes\n");

	float *res = MAT_solve_gausselim(testmatrix, testvec, 3);

	printf("Solution values\n");
	MAT_printvector(res, 3);

	float *rec_b = MAT_multiply_mv(testmatrix, res, 3, 3, 3);
	printf("Recreated b vector (From A . x)\n");
	MAT_printvector(rec_b, 3);

	rec_b = MAT_multiply_mv(backupmat, res, 3, 3, 3);
	printf("Backup, unaltered A:\n");
	MAT_printmatrix(backupmat, 3, 3);
	printf("Reminder: starting b vector\n");
	MAT_printvector(backupvec, 3);
	printf("Recreated b vector from original A:\n");
	MAT_printvector(rec_b, 3);

	return 0;
}
