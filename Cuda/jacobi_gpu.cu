#include <stdio.h>
#include <stdlib.h>
#include <cuda_runtime.h>
#include "utility.cu"
#include "physics.cu"

const int iterations = 50;
// const int Tile_width = 1024;
// const int N = 33554432;

// Try to get a way to store 3D values in csv, calculate time, write a cpu jacobi to compare, and a open mp jacobi

void jacobi_cuda_1D(float* phi_old, float* phi_new, int N, int Tile_width){
    // block = N/Tilewidth, thrd_per_blk = Tilewidth
    

    int blocks = N/Tile_width; 
    int threads_per_blk = Tile_width;

    initialize_vector_by_const<<<blocks, threads_per_blk>>>(phi_old, 0.0f);

    for(int i=0; i<iterations; i++){

        physics_solution_1D<<<blocks, threads_per_blk>>>(phi_old, phi_new);
        equate_vectors<<<blocks, threads_per_blk>>>(phi_old, phi_new);
        
    }  
}

void jacobi_cuda_2D(float* phi_old, float* phi_new, int N, int Tile_width){
    // flatten matrix
    // block = N/Tilewidth, thrd_per_blk = Tilewidth
    

    dim3 blocks(N/Tile_width, N/Tile_width); 
    dim3 threads_per_blk(Tile_width,Tile_width);

    initialize_2Dmatrix_by_const<<<blocks, threads_per_blk>>>(phi_old,0.0f);

    for(int i=0; i<iterations; i++){

        physics_solution_2D<<<blocks, threads_per_blk>>>(phi_old, phi_new);
        equate_flat_matrix2D<<<blocks, threads_per_blk>>>(phi_old, phi_new);
        
    }  
}

void jacobi_cuda_3D(float* phi_old, float* phi_new, int N, int Tile_width){
    // block = N/Tilewidth, thrd_per_blk = Tilewidth
    

    dim3 blocks(N/Tile_width, N/Tile_width, N/Tile_width); 
    dim3 threads_per_blk(Tile_width,Tile_width, Tile_width);

    initialize_3Dmatrix_by_const<<<blocks, threads_per_blk>>>(phi_old, 0.0f);

    for(int i=0; i<iterations; i++){

        physics_solution_3D<<<blocks, threads_per_blk>>>(phi_old, phi_new, N);
        equate_flat_matrix3D<<<blocks, threads_per_blk>>>(phi_old, phi_new);
        
    }  
}

void run(int N, int Tile_width, int times){
    // ----------------------------Writing Performance--------------------------------------------------
    FILE* fp = fopen("benchmark_gpu_3d.csv", "w");
    // if (fp == NULL) {
    //     printf("Error opening file!\n");
    //     return -1;
    // }

    fprintf(fp, "Power,N,GPU\n");
    int power = 3;

    // -------------------------------------------------------------------------------------------------

    for(int i=0; i<=times; i++){
        cudaEvent_t start, stop;
        cudaEventCreate(&start);
        cudaEventCreate(&stop);
        // -----------------------------------Multidimensional memory allocation on GPU-----------------------------------
        float* phi_old;
        cudaMalloc((void**)&phi_old, N*N*N*sizeof(float));
        float* phi_new;
        cudaMalloc((void**)&phi_new, N*N*N*sizeof(float));
        // ----------------------------------------------------------------------------------------------------
        
        // ------------------Jacobi iteration---------------------------------------------------------------
        cudaEventRecord(start,0);
        // jacobi_cuda_1D(phi_old, phi_new,N, Tile_width);
        // jacobi_cuda_2D(phi_old, phi_new,N, Tile_width);
        jacobi_cuda_3D(phi_old, phi_new,N, Tile_width);
        cudaEventRecord(stop,0);

        cudaEventSynchronize(stop);
        float milliseconds = 0;
        cudaEventElapsedTime(&milliseconds, start, stop);
        printf("Time taken for N=%d : %f seconds\n", N ,milliseconds/1000);
        // -------------------------------------------------------------------------------------------------
        fprintf(fp, "%d,%d,%f\n", power+i,N, milliseconds/1000);

        // ---------------------------Writing back to cpu-------------------------------------------------
        // float* result= (float*)malloc(N*N*sizeof(float));
        // cudaMemcpy(result, phi_old, N*N*sizeof(float), cudaMemcpyDeviceToHost);
        
        // ------------------------------------------------------------------------------------------------
        
        // ---------------------------Debug------------------------------
        // for(int i=0; i<N*N; i++){
        //     printf("index=%d, result=%f\n", i,result[i]);
        // }
        // --------------------------------------------------------------

        // Writing to a csv
        // write_to_csv_2d(result, N ,"./result/potential_2D.csv");

        



        // -------------------Freeing------------------------------
        cudaFree(phi_old);
        cudaFree(phi_new);
        // free(result);
        N = N*2;
    }
}


int main(void)
{
    // cudaEvent_t start, stop;
    // cudaEventCreate(&start);
    // cudaEventCreate(&stop);
    // // -----------------------------------Multidimensional memory allocation on GPU-----------------------------------
    // float* phi_old;
    // cudaMalloc((void**)&phi_old, N*sizeof(float));
    // float* phi_new;
    // cudaMalloc((void**)&phi_new, N*sizeof(float));
    // // ----------------------------------------------------------------------------------------------------
    
    // // ------------------Jacobi iteration---------------------------------------------------------------
    // cudaEventRecord(start,0);
    // jacobi_cuda_1D(phi_old, phi_new);
    // cudaEventRecord(stop,0);

    // cudaEventSynchronize(stop);
    // float milliseconds = 0;
    // cudaEventElapsedTime(&milliseconds, start, stop);
    // printf("Time taken for N=%d : %f seconds\n", N ,milliseconds/1000);
    // // -------------------------------------------------------------------------------------------------

    // // ---------------------------Writing back to cpu-------------------------------------------------
    // // float* result= (float*)malloc(N*sizeof(float));
    // // cudaMemcpy(result, phi_old, N*sizeof(float), cudaMemcpyDeviceToHost);
    
    // // ------------------------------------------------------------------------------------------------
    
    // // ---------------------------Debug------------------------------
    // // for(int i=0; i<N*N; i++){
    // //     printf("index=%d, result=%f\n", i,result[i]);
    // // }
    // // --------------------------------------------------------------

    // // Writing to a csv
    // // write_to_csv_1d(result, N ,"./result/potential_1D.csv");

    



    // // -------------------Freeing------------------------------
    // cudaFree(phi_old);
    // cudaFree(phi_new);
    // // free(result);
    
    int N = 8;
    int Tile_width = 8;
    int times = 7;
    run(N, Tile_width, times);
    
    
    return 0;
}