#include <stdio.h>
#include <stdlib.h>
#include <cuda_runtime.h>
#include <math.h>
#include "utility.cu"




__global__ void solution_1D(float* phi_old, float* phi_new){

    int index= blockDim.x*blockIdx.x + threadIdx.x;

    float delta = 0.2f;
    float epsilon = 8.8f;
    float charge_density= 1.0f;

    float b= charge_density*delta*delta/epsilon;
    float phi_0 = 10.0f;
    float phi_N1 = -10.0f;

    if(index==0){
            phi_new[index] = 0.5*(phi_0+ phi_old[index+1] + b);
        } else if (index == N-1) {
            phi_new[index] = 0.5*(phi_old[index-1]+ phi_N1 + b);
        } else {
            phi_new[index] = 0.5*(phi_old[index-1]+ phi_old[index+1] + b);
        }
}


// // -----------------Test-----------------------
// int main(){
//     // float* b_gpu;
//     // cudaMalloc((void**)&b_gpu, N*sizeof(float));
//     // b_vector_1D<<<N/Tile_width,Tile_width>>>(b_gpu,1,1);
//     // cudaDeviceSynchronize();
//     // cudaFree(b_gpu);

//     float* A_cpu[N];
//     for(int i=0;i<N ; i++){
//         cudaMalloc((void**)&A_cpu[i], N* sizeof(float));
//     }
    
//     float** A_gpu;
//     cudaMalloc((void**)&A_gpu, N*sizeof(float*) );
//     cudaMemcpy(A_gpu, A_cpu, N*sizeof(float*), cudaMemcpyHostToDevice);

//     A_matrix_1D(N,A_gpu);


//     // Deallocation
//     for(int i = 0; i < N; i++) {
//         cudaFree(A_gpu[i]);
//     }
//     cudaFree(A_gpu);
//     free(A_cpu);
//     return 0;
// }