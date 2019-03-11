#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "nrutil.h"

// Remember, matrices are referenced like so:
// mat[row][col]

void matutilerror(char *error_text) {
	printf("Matrix Utilities error. Error message follows:\n%s", error_text);
	printf("\nExiting ...\n");
	exit(1);
}

void MAT_printmatrix(float **in_matrix, int rowNum, int colNum) {
	for(int j = 0; j < rowNum; j++) {
		for(int i = 0; i < colNum; i++) {
			printf(" %5f ", in_matrix[j][i]);
		}
		printf("\n");
	}
	printf("\n");
}

void MAT_printvector(float *in_vector, int n) {
	for (int i = 0; i < n; i++) {
		printf(" %5f ", in_vector[i]);
		printf("\n");
	}
	printf("\n");
}

void MAT_printvector_range(float *in_vector, int stid, int endid) {
	for (int i = stid; i <= endid; i++) {
		printf(" %5f ", in_vector[i]);
		printf("\n");
	}
	printf("\n");
}


void MAT_graphvector(float *in_vector, int n) {
	int width = 100;
	float lower = in_vector[0];
	float upper = in_vector[0];
	for (int i = 1; i < n; i++) {
		if (in_vector[i] < lower) {
			lower = in_vector[i];
		}
		if (in_vector[i] > upper) {
			upper = in_vector[i];
		}
	}
	lower = 0; // Remove to allow floating origin
	float range = upper - lower;
	float scale = range / (float) width;
	printf("\n");
	for (int i = 0; i < n; i++) {
		for (int count = 0; count * scale + lower < in_vector[i]; count++) {
			printf(" ");
		}
		printf("X\n");
	}
}

void MAT_zerovector(float *in_vector, int n) {
	//Zero out a 0-indexed vector of size n
	memset(in_vector, 0, sizeof (float) * n);
}

void MAT_zerovector_range(float *in_vector, int stid, int endid) {
	for (int i = stid; i <= endid; i++) {
		in_vector[i] = 0;
	}
}

void MAT_zeromatrix(float **in_matrix, int rowNum, int colNum) {
	for (int row = 0; row < rowNum; row++) {
		for (int col = 0; col < colNum; col++) {
			in_matrix[row][col] = 0;
		}
	}
}

void MAT_setmatrix(float **a, float **b, int rowst, int rowend, int colst, int colend) {
	// Sets a = b within the given boundaries (inclusively)
	for (int col = colst; col <= colend; col++) {
		for (int row = rowst; row <= rowend; row++) {
			a[row][col] = b[row][col];
		}
	}
}

void MAT_addmatrix(float **a, float **b, int rowst, int rowend, int colst, int colend) {
	// Sets a += b within the given boundaries (inclusively)
	for (int col = colst; col <= colend; col++) {
		for (int row = rowst; row <= rowend; row++) {
			a[row][col] += b[row][col];
		}
	}
}

float **MAT_multiply_mm(float **a_matrix, float **b_matrix, int a_rowNum, int a_colNum, int b_rowNum, int b_colNum) {
	if (a_colNum != b_rowNum) { matutilerror("Incompatible matrix dimensions in MAT_multiply_mm"); }

	float **result = matrix(0, a_rowNum - 1, 0, b_colNum - 1); //Indexing problems may occur here. Ensure 0-indexing is consistent
	float temp_sum = 0;
	float temp_a = 0;
	float temp_b = 0;

	for (int res_row = 0; res_row < a_rowNum; res_row++) {
		for (int res_col = 0; res_col < b_colNum; res_col++) {
			temp_sum = 0;
			for (int i = 0; i < a_colNum; i++) {
				temp_a = a_matrix[res_row][i];
				temp_b = b_matrix[i][res_col];
				if (temp_b && temp_a) {
					temp_sum += temp_a * temp_b;
				}
			}
			result[res_row][res_col] = temp_sum;
		}
	}
	return result;
}

float *MAT_multiply_mv(float **a_matrix, float *b_vector, int a_rowNum, int a_colNum, int b_rowNum) {
	if (a_colNum != b_rowNum) { matutilerror("Incompatible matrix dimensions in MAT_multiply_mv"); }

	float *result = vector(0, a_rowNum - 1);
	float temp_sum = 0;

	for (int res_row = 0; res_row < a_rowNum; res_row++) {
		temp_sum = 0;
		for (int i = 0; i < a_colNum; i++) {
			temp_sum += a_matrix[res_row][i] * b_vector[i];
		}
		result[res_row] = temp_sum;
	}
	return result;
}

float *MAT_solve_gausselim(float **in_matrix, float *b, int n) {
	// Gaussian elimination and backsubstitution. In a system A . x = b, This algorithm finds x given A *and* b.
	// For standalone inversion of A, use the Gauss-Jordan method or a different algorithm.
	// Warning: to save time and memory, this function does not preserve in_matrix or b
	// Everything assumed to be 0-indexed: A[0...n-1][0...n-1], b[0...n-1]

	float *x = vector(0, n - 1);

	float highval;
	int max_row;
	float max_value;
	float temp; //generic thing for switching values or whatever
	float scaling;

	for (int col = 0; col < n; col++) {
		//find the row with the greatest value in the current column (past finished rows)
		max_row = -1;
		max_value = 0;
		for (int i = col; i < n; i++) {
			temp = fabsf(in_matrix[i][col]);
			if (temp > max_value) {
				max_row = i;
				max_value = temp;
			}
		}
		if (max_row == -1) { matutilerror("Solve error 1 in MAT_solve_gausselim"); } // This might actually happen even if the system is solvable. See if it's a problem.
		// swap the row that was just found to the current row
		for (int i = 0; i < n; i++) {
			temp = in_matrix[col][i];
			in_matrix[col][i] = in_matrix[max_row][i];
			in_matrix[max_row][i] = temp;
		}
		// also swap the b value
		temp = b[col];
		b[col] = b[max_row];
		b[max_row] = temp;
		//scale the current row to make the [col][col] value 1
		scaling = in_matrix[col][col];
		for (int i = col; i < n; i++) {
			in_matrix[col][i] /= scaling;
		}
		b[col] /= scaling;
		// For all lower rows, subtract enough of the current row (and b) to make the col value 0
		for (int row = col + 1; row < n; row++) {
			scaling = in_matrix[row][col] / in_matrix[col][col];
			for (int i = col; i < n; i++) {
				in_matrix[row][i] -= scaling * in_matrix[col][i];
			}
			b[row] -= scaling * b[col];
		}
	}
	// Now in_matrix should be properly triangular, and b should match
	// Time for backsubstitution
	for (int row = n-1; row >= 0; row--) {
		temp = 0;
		for (int j = row + 1; j < n; j++) {
			temp += in_matrix[row][j] * x[j];
		}
		x[row] = b[row] - temp; // assuming that the leading coefficient is 1. Otherwise divide by leading coefficient
	}

	return x;
}

float **MAT_invert_gaussjordan(float **in_matrix, int n) {
	//For later

	float **result = matrix(0, n - 1, 0, n - 1);
	return result;
}