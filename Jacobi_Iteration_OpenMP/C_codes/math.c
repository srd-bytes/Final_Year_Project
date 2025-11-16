#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

double** kronecker(double **A, int a_rows, int a_cols, double **B, int b_rows, int b_cols);
// Kronecker product: C = kron(A, B)
double **kronecker(double **A, int a_rows, int a_cols, double **B, int b_rows, int b_cols) {

    int rows = a_rows * b_rows;
    int cols = a_cols * b_cols;

    double **C = calloc(rows, sizeof(double *));
    for (int i = 0; i < rows; i++) {
        C[i] = calloc(cols, sizeof(double));
    }

    for (int i = 0; i < a_rows; i++) {
        for (int j = 0; j < a_cols; j++) {
            for (int k = 0; k < b_rows; k++) {
                for (int l = 0; l < b_cols; l++) {
                    C[i*b_rows + k][j*b_cols + l] =
                        A[i][j] * B[k][l];
                }
            }
        }
    }

    return C;
} // reference : poisson2d_notes

// Identity matrix
double **identity_matrix(int n) {
    double **I = calloc(n, sizeof(double *));
    for (int i = 0; i < n; i++) {
        I[i] = calloc(n, sizeof(double));
        I[i][i] = 1.0;
    }
    return I;
}