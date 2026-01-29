#include <stdio.h>
#include <stdlib.h>
#include <cuda_runtime.h>
#include "physics.cu"

const int iterations = 100;
__global__ void equate(float* phi_old, float* phi_new){
    int index = blockDim.x*blockIdx.x+threadIdx.x;
    phi_old[index] = phi_new[index];
}
void jacobi(float** A, float* phi_old, float* phi_new, float* b){


    dim3 blocks(N/Tile_width,N/Tile_width);
    dim3 threads_per_block(Tile_width,Tile_width);

    float* inter_phi;
    cudaMalloc((void**)&inter_phi, N*sizeof(float));
    for(int k = 0; k < iterations; k++){
        matrix_vector_product<<<blocks,threads_per_block>>>(N,A,phi_old,inter_phi);
        vector_addition_and_scaler_mult<<<N/Tile_width,Tile_width>>>(inter_phi,b, phi_new,0.5f);
        equate<<<N/Tile_width,Tile_width>>>(phi_old,phi_new);
    }
    
    cudaFree(inter_phi);
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
    

    A_matrix_1D(N,A_gpu);

    // dim3 blocks(N/Tile_width,N/Tile_width);
    // dim3 threads_per_block(Tile_width,Tile_width);
    // printf("A Matrix:\n");
    // print_matrix<<<blocks,threads_per_block>>>(A_gpu);
    
    //Vector part
    float* v_gpu;
    cudaMalloc((void**)&v_gpu, N*sizeof(float));
    initialize_vector_const<<<N/Tile_width, Tile_width>>>(N,v_gpu,0.0f); // initialize to 0

    
    
    
    //b Vector part
    float* b_gpu;
    cudaMalloc((void**)&b_gpu, N*sizeof(float));
    // initialize_vector_const<<<N/Tile_width, Tile_width>>>(N,b_gpu,1.0f); // make a function charge density function 1d, 2d, 3d to calculate this
    b_vector_1D<<<N/Tile_width, Tile_width>>>(b_gpu,10.0f,-10.0f);

   
   float* phi_new;
   cudaMalloc((void**)&phi_new, N*sizeof(float));

    jacobi(A_gpu,v_gpu,phi_new,b_gpu);
    print_vector<<<N/Tile_width, Tile_width>>>(phi_new);

    //printing result
    

    // Writing to a csv
    float* phi_cpu;
    phi_cpu = (float*)malloc(N*sizeof(float));
    cudaMemcpy(phi_cpu, phi_new, N*sizeof(float), cudaMemcpyDeviceToHost);
    FILE* fp;
    fp = fopen("potential.csv","w");
    for(int i = 0; i < N; i++){
        fprintf(fp,"%f\n",phi_cpu[i]);
    }
    fclose(fp);

    // Deallocation
    for(int i = 0; i < N; i++) {
        cudaFree(A_gpu[i]);
    }
    cudaFree(A_gpu);
    free(A_cpu);
    cudaFree(v_gpu);
    cudaFree(b_gpu);
    cudaFree(phi_new);
    cudaFree(phi_cpu);
    
    
    return 0;
}