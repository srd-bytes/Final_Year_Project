#include <stdio.h>
#include <stdlib.h>
#include <cuda_runtime.h>
#include "physics.cu"

const int iterations = 100;
const int Tile_width = 8;
const int N = 128;



void jacobi_cuda_1D(float* phi_old, float* phi_new){
    // block = N/Tilewidth, thrd_per_blk = Tilewidth
    

    int blocks = N/Tile_width; 
    int threads_per_blk = Tile_width;

    initialize_vector_by_const<<<blocks, threads_per_blk>>>(phi_old, 0.0f);

    for(int i=0; i<iterations; i++){

        solution_1D<<<blocks, threads_per_blk>>>(phi_old, phi_new);
        equate_vectors<<<blocks, threads_per_blk>>>(phi_old, phi_new);
        
    }  
}




int main(void)
{
    float* phi_old;
    float* phi_new;

    cudaMalloc((void**)&phi_new, N*sizeof(float));
    cudaMalloc((void**)&phi_old, N*sizeof(float));

    

    jacobi_cuda_1D(phi_old, phi_new);

    float* phi_cpu = (float*) malloc(N*sizeof(float));
    cudaMemcpy(phi_cpu, phi_old, N* sizeof(float), cudaMemcpyDeviceToHost);

    // Writing to a csv
    FILE* fp;
    fp = fopen("potential.csv","w");
    for(int i = 0; i < N; i++){
        fprintf(fp,"%f\n",phi_cpu[i]);
    }
    fclose(fp);

    cudaFree(phi_new);
    cudaFree(phi_old);
    free(phi_cpu);
    
    return 0;
}