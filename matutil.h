#ifndef _MATRIXUTIL_1D_H_
#define _MATRIXUTIL_1D_H_

void matutilerror(char *error_text);
void MAT_printmatrix(float **in_matrix, int rowNum, int colNum);
void MAT_printvector(float *in_vector, int n);
void MAT_printvector_range(float *in_vector, int stid, int endid);
void MAT_zerovector(float *in_vector, int n);
void MAT_zerovector_range(float *in_vector, int stid, int endid);
void MAT_zeromatrix(float **in_matrix, int rowNum, int colNum);
void MAT_setmatrix(float **a, float **b, int rowst, int rowend, int colst, int colend);
void MAT_addmatrix(float **a, float **b, int rowst, int rowend, int colst, int colend);
void MAT_graphvector(float *in_vector, int n);
float **MAT_multiply_mm(float **a_matrix, float **b_matrix, int a_rowNum, int a_colNum, int b_rowNum, int b_colNum);
float *MAT_multiply_mv(float **a_matrix, float *b_vector, int a_rowNum, int a_colNum, int b_rowNum);
float *MAT_solve_gausselim(float **in_matrix, float *b, int n);
float **MAT_invert_gaussjordan(float **in_matrix, int n);

#endif