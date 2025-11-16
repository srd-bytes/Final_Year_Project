
#include "Generate_poisson_matrix.c"
#include <time.h>
//Helper function to perform the Jacobi iterative method

// Function prototype
double *jacobi(int N, int max_iter, double tol, double **A, double *b, double *x);

double *jacobi(int N, int max_iter, double tol, double **A, double *b, double*x) { 
    //after passing A it is a contiguous 1D array flattened from a 2D array
    double *x_new = (double *) calloc(N, sizeof(double));

    
    x_new[0]=x[0];
    x_new[N-1]=x[N-1];
    for (int k = 0; k < max_iter; k++) {

        
        for (int i = 1; i < N-1; i++) {
            x_new[i] = b[i];
            for (int j = 0; j < N; j++) {
                if (i != j) {
                    x_new[i] -= A[i][j] * x[j];
                }
            }
            x_new[i] /= A[i][i]; // Diagonal element
        }
        
        // Check for convergence
        double max_diff = 0;

        
        for (int i = 0; i < N; i++) {
            max_diff = fmax(max_diff, fabs(x_new[i] - x[i]));
        }
        if (max_diff < tol) {
            break;
        }
        memcpy(x, x_new, N * sizeof(double));
    }
    free(x_new);
    return x;
}

int main(void) {
   
    clock_t start, end;
    double cpu_time_used;
    
    // jacobi parameters
    int max_iter=1000;
    double tol=1e-20;


    //region
    // system parameters
    double x_L =0.0;
    double x_H=10.0;
    int n=50;
    double delta= (x_H-x_L)/n ;
    //end

    //region
    //1D
    // double **A = generate_poisson_matrix1D(n, delta);
    // double *x = calloc(n, sizeof(double));
    // double* rho=generate_charge_1D(n,delta,2);
    // double *b = b_vec_1D(x_L, delta, n,rho);

    //2D
    double **A = generate_poisson_matrix2D(n,n,delta);
    double *x = calloc(n*n, sizeof(double));
    double** rho=generate_charge_2D(n,delta,2);
    double *b = b_vec_2D(x_L,x_L,delta,n,rho);
    //end

    //region Boundary conditions

    //1D
    // x[0] = 10.0;
    // x[n-1] = -10.0;

    //2D
    int N2 = n*n;

    for (int i = 0; i < N2; i++) {
        int row = i / n;
        int col = i % n;

        if (row == 0) {        // top boundary
            for (int j = 0; j < N2; j++) A[i][j] = 0.0;
            A[i][i] = 1.0;
            b[i] = 10.0;
        }
        if (row == n-1) {      // bottom boundary
            for (int j = 0; j < N2; j++) A[i][j] = 0.0;
            A[i][i] = 1.0;
            b[i] = -10.0;
        }
        if (col == 0 || col == n-1) {
            for (int j = 0; j < N2; j++) A[i][j] = 0.0;
            A[i][i] = 1.0;
            b[i] = 0.0; // left/right
        }
    }




    //end

    //region
    start = clock();

    jacobi(n*n, max_iter, tol, A, b, x);

    end = clock();
    //end

    //region
    cpu_time_used = ((double)(end - start)) / (CLOCKS_PER_SEC);
    printf("Total CPU time: %f seconds\n", cpu_time_used);
    //end

    // printf("Solution:\n");
    // for (int i = 0; i < n; i++) {
        printf("x[%d] = %f\n", 3, x[3]);
    // }

    //region
    //1D

    // FILE *fp1 = fopen("potential_1d.csv", "w");
    // FILE *fp2 = fopen("charge_1d.csv", "w");
    // for (int i = 0; i < n; i++) {
    //     fprintf(fp1, "%d,%lf\n", i, x[i]);
    //     fprintf(fp2, "%d,%lf\n", i, rho[i]);
    // }
    // free(rho);
    // fclose(fp1); fclose(fp2);

    //2D

    FILE *fp1 = fopen("potential_2d.csv", "w");
    FILE *fp2 = fopen("charge_2d.csv", "w");
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            fprintf(fp1, "%d,%d,%lf\n", i, j, x[i*n +j]);
            fprintf(fp2, "%d,%d,%lf\n", i, j, rho[i][j]);
        }
        free(rho[i]);
    }
    free(rho);
    fclose(fp1); fclose(fp2);

    //end
    free(x);
    return 0;
}