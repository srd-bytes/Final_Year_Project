#include <stdio.h>
#include <stdlib.h>
#include <cuda_runtime.h>
#include "utility.cu"
#include "physics.cu"

const int iterations = 100;
const int Tile_width = 8;
const int N = 128;

// Try to get a way to store 3D values in csv, calculate time, write a cpu jacobi to compare, and a open mp jacobi

void jacobi_cuda_1D(float* phi_old, float* phi_new){
    // block = N/Tilewidth, thrd_per_blk = Tilewidth
    

    int blocks = N/Tile_width; 
    int threads_per_blk = Tile_width;

    initialize_vector_by_const<<<blocks, threads_per_blk>>>(phi_old, 0.0f);

    for(int i=0; i<iterations; i++){

        physics_solution_1D<<<blocks, threads_per_blk>>>(phi_old, phi_new, N);
        equate_vectors<<<blocks, threads_per_blk>>>(phi_old, phi_new);
        
    }  
}

void jacobi_cuda_2D(float** phi_old, float** phi_new){
    // block = N/Tilewidth, thrd_per_blk = Tilewidth
    

    dim3 blocks(N/Tile_width, N/Tile_width); 
    dim3 threads_per_blk(Tile_width,Tile_width);

    initialize_2Dmatrix_by_const<<<blocks, threads_per_blk>>>(phi_old, 0.0f);

    for(int i=0; i<iterations; i++){

        physics_solution_2D<<<blocks, threads_per_blk>>>(phi_old, phi_new, N);
        equate_matrix2D<<<blocks, threads_per_blk>>>(phi_old, phi_new);
        
    }  
}

void jacobi_cuda_3D(float*** phi_old, float*** phi_new){
    // block = N/Tilewidth, thrd_per_blk = Tilewidth
    

    dim3 blocks(N/Tile_width, N/Tile_width, N/Tile_width); 
    dim3 threads_per_blk(Tile_width,Tile_width, Tile_width);

    initialize_3Dmatrix_by_const<<<blocks, threads_per_blk>>>(phi_old, 0.0f);

    for(int i=0; i<iterations; i++){

        physics_solution_3D<<<blocks, threads_per_blk>>>(phi_old, phi_new, N);
        equate_matrix3D<<<blocks, threads_per_blk>>>(phi_old, phi_new);
        
    }  
}




int main(void)
{
    // -----------------------------------Multidimensional memory allocation on GPU-----------------------------------
    float** phi_old_cpu= (float**)malloc(N*sizeof(float*));
    float** phi_new_cpu= (float**)malloc(N*sizeof(float*));


    for(int i=0; i<N; i++){
        cudaMalloc((void**)&phi_new_cpu[i], N*sizeof(float));
        cudaMalloc((void**)&phi_old_cpu[i], N*sizeof(float));
    }
    
    float** phi_old_gpu;
    float** phi_new_gpu;
    cudaMalloc((void**)&phi_new_gpu, N*sizeof(float*));
    cudaMalloc((void**)&phi_old_gpu, N*sizeof(float*));
    cudaMemcpy(phi_old_gpu, phi_old_cpu, N*sizeof(float*), cudaMemcpyHostToDevice);
    cudaMemcpy(phi_new_gpu, phi_new_cpu, N*sizeof(float*), cudaMemcpyHostToDevice);
    
    // ----------------------------------------------------------------------------------------------------
    
    // ------------------Jacobi iteration---------------------------------------------------------------
    jacobi_cuda_2D(phi_old_gpu, phi_new_gpu);
    // -------------------------------------------------------------------------------------------------

    // ---------------------------Writing back to cpu-------------------------------------------------
    float** result= (float**)malloc(N*sizeof(float*));
    for(int i=0; i<N;i++){
        result[i]= (float*)malloc(N*sizeof(float));
        cudaMemcpy(result[i], phi_new_cpu[i], N* sizeof(float), cudaMemcpyDeviceToHost);
    }
    // ------------------------------------------------------------------------------------------------
    
    

    // Writing to a csv
    write_to_csv_2d(result,"./result/potential_2D.csv");

    // // 3D solution
    // FILE* fp = fopen("potential_3D.csv", "w");
    // fprintf(fp, "x,y,z,phi\n");

    // for (int i = 0; i < Nx; i++) {
    //     for (int j = 0; j < Ny; j++) {
    //         for (int k = 0; k < Nz; k++) {
                
    //             fprintf(fp, "%f\n", phi[i][j][k]);
    //         }
    //     }
    // }

    // fclose(fp);



    // -------------------Freeing------------------------------
    for(int i=0; i<N;i++){
        cudaFree(phi_old_cpu[i]);
        cudaFree(phi_new_cpu[i]);
        free(result[i]);
        
    }
    cudaFree(phi_new_gpu);
    cudaFree(phi_old_gpu);
    free(phi_old_cpu);
    free(phi_new_cpu);
    free(result);
    
    return 0;
}