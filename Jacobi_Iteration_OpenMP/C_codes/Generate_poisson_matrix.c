#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "math.c"

// const int n=2; // number of discrete points in one direction
const double tot_length=3;
// const double delta= 1;



double **generate_poisson_matrix1D(double n,double delta);
double **generate_poisson_matrix1D(double n,double delta){

    
    double **A = calloc(n , sizeof(double *));
    for (int i = 0; i < n; i++)
    {
        double *row= calloc(n,sizeof(double));
        row[i]=-2/(delta*delta);
        if (i<n-1){ row[i+1]=1/(delta*delta);}
        
        if (i>0){ row[i-1]=1/(delta*delta);}
        A[i]=row;
        
    }
   
    return A;
    
}

double **generate_poisson_matrix2D(int nx, int ny,double delta);
double **generate_poisson_matrix2D(int nx, int ny,double delta) {
    // A_x and A_y
    double **Ax = generate_poisson_matrix1D(nx, delta);
    double **Ay = generate_poisson_matrix1D(ny, delta);

    // Identity matrices
    double **Ix = identity_matrix(nx);
    double **Iy = identity_matrix(ny);

    // Kronecker products
    double **IyAx = kronecker(Iy, ny, ny, Ax, nx, nx);
    double **AyIx = kronecker(Ay, ny, ny, Ix, nx, nx);

    //region
    for(int i=0;i<nx;i++){
        free(Ax[i]);
        free(Ay[i]);
        free(Ix[i]);
        free(Iy[i]);
    }
    free(Ax);
    free(Ay);
    free(Ix);
    free(Iy);
    //end

    int N = nx * ny;
    // Sum them into A2D
    double **A2D = calloc(N, sizeof(double *));
    for (int i = 0; i < N; i++) {
        A2D[i] = calloc(N, sizeof(double));
        for (int j = 0; j < N; j++) {
            A2D[i][j] = IyAx[i][j] + AyIx[i][j];
        }
    }
    //region
    for(int i=0;i<nx;i++){
        free(IyAx[i]);
        free(AyIx[i]);
    }
    free(IyAx);
    free(AyIx);
    //end

    return A2D;
}

double **generate_poisson_matrix3D(int nx, int ny, int nz,double delta);
double **generate_poisson_matrix3D(int nx, int ny, int nz,double delta) {
    // 1D operators
    double **Ax = generate_poisson_matrix1D(nx,delta);
    double **Ay = generate_poisson_matrix1D(ny,delta);
    double **Az = generate_poisson_matrix1D(nz, delta);

    // Identities
    double **Ix = identity_matrix(nx);
    double **Iy = identity_matrix(ny);
    double **Iz = identity_matrix(nz);

    // Kronecker products
    double **IzIyAx = kronecker(Iz, nz, nz, kronecker(Iy, ny, ny, Ax, nx, nx), ny*nx, ny*nx);
    double **IzAyIx = kronecker(Iz, nz, nz, kronecker(Ay, ny, ny, Ix, nx, nx), ny*nx, ny*nx);
    double **AzIyIx = kronecker(Az, nz, nz, kronecker(Iy, ny, ny, Ix, nx, nx), ny*nx, ny*nx);
    //region
    for(int i=0;i<nx;i++){
        free(Ax[i]);
        free(Ay[i]);
        free(Az[i]);
        free(Ix[i]);
        free(Iy[i]);
        free(Iz[i]);
    }
    free(Ax);
    free(Ay);
    free(Az);
    free(Ix);
    free(Iy);
    free(Iz);
    //end
    int N = nx * ny * nz;
    double **A3D = calloc(N, sizeof(double *));
    for (int i = 0; i < N; i++) {
        A3D[i] = calloc(N, sizeof(double));
        for (int j = 0; j < N; j++) {
            A3D[i][j] = IzIyAx[i][j] + IzAyIx[i][j] + AzIyIx[i][j];
        }
    }
    //region
    for(int i=0;i<nx;i++){
        free(IzIyAx[i]);
        free(IzAyIx[i]);
        free(AzIyIx[i]);
    }
    free(IzIyAx);
    free(IzAyIx);
    free(AzIyIx);
    //end

    return A3D;
} 
// reference : poisson2d_notes

//region
double* generate_charge_1D( int n, double delta, int type);
double* generate_charge_1D( int n, double delta, int type) {
    double* rho= calloc(n,sizeof(double));
    double xL = 0, xH = n * delta;
    for (int i = 0; i < n; i++) {
        double x = xL + i * delta;
        if (type == 1)
            rho[i] = (i == n/2) ? 1.0 : 0.0;
        else if (type == 2)
            rho[i] = 1.0;
        else if (type == 3)
            rho[i] = exp(-pow((x - xH/2)/0.2, 2));
        else if (type == 4)
            rho[i] = (i == n/4) ? 1.0 : (i == 3*n/4) ? -1.0 : 0.0;
    }
    return rho;
}

double** generate_charge_2D(int n, double delta, int type);
double** generate_charge_2D(int n, double delta, int type) {
    double** rho= calloc(n,sizeof(double*));

    double xL = 0, xH = n * delta;
    double yL = 0, yH = n * delta;

    for (int i = 0; i < n; i++) {
        double x = xL + i * delta;
        rho[i]=calloc(n,sizeof(double));
        for(int j=0;j<n;j++){
            double y = yL + j * delta;
            
            if (type == 1)
                rho[i][j] = (i == n/2)&&(j==n/2) ? 1.0 : 0.0;
            else if (type == 2)
                rho[i][j] = 1.0; // charge density in pC/m^3 as i have normalized the epsilon as 8.85
            else if (type == 3)
                rho[i][j] = exp(-pow((x - xH/2)*(y-yH/2)/0.2, 2));
            else if (type == 4)
                rho[i][j] = (i == n/4)&&(j==n/4) ? 1.0 : (i == 3*n/4)&&(j==3*n/4) ? -1.0 : 0.0;

        }
    }
    return rho;
}

//end
#define EPS0 8.85

double * b_vec_1D(double x_L, double delta, int n,double* rho);
double * b_vec_1D(double x_L, double delta, int n,double* rho){
    double *b= calloc(n, sizeof(double));

    for(int i=0; i<n; i++){
        // b[i]= -cos(M_PI*(x_L+i*delta)); // RULE
        b[i]= -rho[i]/EPS0; // RULE
    }
    
    return b;
}

double * b_vec_2D(double x_L, double y_L, double delta, int n,double** rho);
double * b_vec_2D(double x_L,  double y_L, double delta, int n, double** rho){
    double *b= calloc(n*n, sizeof(double));

    for(int x=0; x<n; x++){

        for(int y=0; y<n; y++){

            // b[x*n + y]= -cos(M_PI*(x_L+x*delta))*cos(M_PI*(y_L+y*delta)); // RULE
            b[x*n + y]= -rho[x][y]/EPS0; // RULE
        }
        
    }
    
    return b;
}

double * b_vec_3D(double x_L, double y_L, double z_L, double delta, int n);
double * b_vec_3D(double x_L,  double y_L, double z_L, double delta, int n){
    double *b= calloc(n*n*n, sizeof(double));

    for(int x=0; x<n; x++){

        for(int y=0; y<n; y++){

            for(int z=0; z<n; z++){

                b[x*n*n + n*y + z]= -cos(M_PI*(x_L+x*delta))*cos(M_PI*(y_L+y*delta))*cos(M_PI*(z_L+z*delta)); // RULE

            }
        }
    }
    return b;
}
// int main(){


    

//     return 0;
// }