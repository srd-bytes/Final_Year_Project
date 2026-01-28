#include <stdio.h>
#include <stdlib.h>
#include <cuda_runtime.h>
#include "utility.cu"



void jacobi(float** A, float* phi, float* b){


    dim3 blocks(N/Tile_width,N/Tile_width);
    dim3 threads_per_block(Tile_width,Tile_width);

    int iterations = 1;
    for(int k = 0; k < iterations; k++){
        matrix_vector_product<<<blocks,threads_per_block>>>(N,A,phi);
        vector_addition_and_scaler_mult<<<N/Tile_width,Tile_width>>>(phi,b,1.0f);
        
    }
    
}

int main(void)
{
    //Matrix allocation
    float* A_cpu[N];
    for(int i=0;i<N ; i++){
        cudaMalloc((void**)&A_cpu[i], N* sizeof(float));
    }
    
    float** A_gpu;
    cudaMalloc((void**)&A_gpu, N*sizeof(float*) );
    cudaMemcpy(A_gpu, A_cpu, N*sizeof(float*), cudaMemcpyHostToDevice);
    

    dim3 blocks(N/Tile_width,N/Tile_width);
    dim3 threads_per_block(Tile_width,Tile_width);
    initialize_matrix<<<blocks,threads_per_block>>>(N,A_gpu);

    //Vector part
    float* v_gpu;
    cudaMalloc((void**)&v_gpu, N*sizeof(float));
    initialize_vector<<<N/Tile_width, Tile_width>>>(N,v_gpu);

    //b Vector part
    float* b_gpu;
    cudaMalloc((void**)&b_gpu, N*sizeof(float));
    initialize_vector<<<N/Tile_width, Tile_width>>>(N,b_gpu);

    jacobi(A_gpu,v_gpu,b_gpu);

    //printing result
    float* phi_cpu;
    phi_cpu = (float*)malloc(N*sizeof(float));
    cudaMemcpy(phi_cpu, v_gpu, N*sizeof(float), cudaMemcpyDeviceToHost);
    for(int i = 0; i < N; i++){
        printf("phi[%d] = %f\n",i,phi_cpu[i]);
    }

    // Deallocation
    for(int i = 0; i < N; i++) {
        cudaFree(A_gpu[i]);
    }
    cudaFree(A_gpu);
    free(A_cpu);
    cudaFree(v_gpu);
    cudaFree(b_gpu);
    free(phi_cpu);
    
    return 0;
}