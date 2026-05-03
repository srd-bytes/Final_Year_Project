#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "utility.c"

const int iterations = 1;


void jacobi_1D(float* phi_old, float* phi_new, int N){
    

    float delta = 1.0f;
    float epsilon = 1.0f;
    float charge_density= 1.0f; // pC/m

    float b= charge_density*delta*delta/epsilon;
    
    float left_boundary = 1.0f;
    float right_boundary = -1.0f;
    
    for(int i=0; i<iterations; i++){

        for(int j=0; j<N; j++){
            float left = (j==0) ? left_boundary : phi_old[j-1];
            float right = (j==N-1) ? right_boundary : phi_old[j+1];

            phi_new[j] = (left+right+b)/2;
            
        }

        for(int j=0;j<N;j++){phi_old[j] = phi_new[j];}
    }

}

void jacobi_2D(float* phi_old, float* phi_new, int N){
    

    float delta = 0.2f;
    float epsilon = 8.8f;
    float charge_density= 1.0f; // pC/m^2

    float b= charge_density*delta*delta/epsilon;
    
    float left_boundary = 10.0f;
    float right_boundary = -10.0f;
    float top_boundary = 10.0f;
    float bottom_boundary = -10.0f;
    
    for(int i=0; i<iterations; i++){

        for(int y=0; y<N; y++){
            for(int x=0; x<N; x++){
                float left = (x==0) ? left_boundary : phi_old[y*N+x-1];
                float right = (x==N-1) ? right_boundary : phi_old[y*N+x+1];
                float top = (y==0) ? top_boundary : phi_old[y*N+x -N];
                float bottom = (y==N-1) ? bottom_boundary : phi_old[y*N+x +N];

                phi_new[y*N+x] = (left+right+top+bottom+b)/4;
                
            }
        }

        for(int y=0; y<N; y++){
            for(int x=0; x<N; x++){phi_old[y*N+x] = phi_new[y*N+x];}}
    }
}

void jacobi_3D(float* phi_old, float* phi_new, int N){
    
    float delta = 0.2f;
    float epsilon = 8.8f;
    float charge_density= 1.0f; // pC/m^3

    float b= charge_density*delta*delta/epsilon;
    
    float left_boundary = 10.0f;
    float right_boundary = -10.0f;
    float top_boundary = 0.0f;
    float bottom_boundary = 0.0f;
    float front_boundary = 0.0f;
    float back_boundary = 0.0f;
    
    for(int i=0; i<iterations; i++){

        for(int z=0; z<N; z++){
            for(int y=0; y<N; y++){
                for(int x=0; x<N; x++){
                    float left = (x==0) ? left_boundary : phi_old[z*N*N+y*N+x-1];
                    float right = (x==N-1) ? right_boundary : phi_old[z*N*N+y*N+x+1];
                    float top = (y==0) ? top_boundary : phi_old[z*N*N+y*N+x -N];
                    float bottom = (y==N-1) ? bottom_boundary : phi_old[z*N*N+y*N+x +N];
                    float front = (z==0) ? front_boundary : phi_old[z*N*N+y*N+x -N*N];
                    float back = (z==N-1) ? back_boundary : phi_old[z*N*N+y*N+x +N*N];

                    phi_new[z*N*N+y*N+x] = (left+right+top+bottom+front+back+b)/6;
                    
                }
            }
        }

        for(int z=0; z<N; z++){
            for(int y=0; y<N; y++){
                for(int x=0; x<N; x++){phi_old[z*N*N+y*N+x] = phi_new[z*N*N+y*N+x];}}}
    }
}

void run(int N, int times){
    // ----------------------------Writing Performance--------------------------------------------------
    // FILE* fp = fopen("./result/performance/benchmark_cpu_3d.csv", "w");
    // if (fp == NULL) {
    //     printf("Error opening file!\n");
    //     return;
    // }

    // fprintf(fp, "Power,N,CPU\n");
    // int power = 3;

    // -------------------------------------------------------------------------------------------------
    for(int i=0; i<times; i++){
        float* phi_old;
        float* phi_new;
        
        phi_old = (float*)calloc(N*N,sizeof(float));
        phi_new = (float*)calloc(N*N,sizeof(float));
        
        
        time_t start = clock();
        // jacobi_1D(phi_old, phi_new, N);
        jacobi_2D(phi_old, phi_new, N);
        // jacobi_3D(phi_old, phi_new, N);
        time_t end = clock();
        double time_taken = ((double)end - start) / CLOCKS_PER_SEC;
        printf("Time taken for N=%d : %f seconds\n", N, time_taken);

        for(int j=0; j<N*N; j++){printf("phi[%d] = %f\n",j,phi_new[j]);}

        // --------------------write benchmark--------------------------------------
        // fprintf(fp, "%d,%d,%f\n", power+i,N, time_taken);
        // -------------------------------------------------------------------------
        
        
        // write_to_csv_1d(phi_old, N, "./result/potential_1d.csv");
        // write_to_csv_2d(phi_old, N, "phi_old_2d.csv");
        // write_to_csv_3d(phi_old, N, "phi_old_3d.csv");
        
        free(phi_old);
        free(phi_new);
        // N=N*2;
    }
}
int main(){

    int N=2;
    int times=1;
    run(N,times);

    // float* phi_old;
    // float* phi_new;
    
    // phi_old = (float*)calloc(N,sizeof(float));
    // phi_new = (float*)calloc(N,sizeof(float));
    
    
    // time_t start = clock();
    // jacobi_1D(phi_old, phi_new);
    // // jacobi_2D(phi_old, phi_new);
    // // jacobi_3D(phi_old, phi_new);
    // time_t end = clock();
    // double time_taken = ((double)end - start) / CLOCKS_PER_SEC;
    // printf("Time taken for N=%d : %f seconds\n", N, time_taken);

    // fprintf(fp, "%d,%d,%f\n", power+i,N, milliseconds/1000);
    
    
    
    // // write_to_csv_1d(phi_old, N, "./result/potential_1d.csv");
    // // write_to_csv_2d(phi_old, N, "phi_old_2d.csv");
    // // write_to_csv_3d(phi_old, N, "phi_old_3d.csv");
    
    // free(phi_old);
    // free(phi_new);
    
}
