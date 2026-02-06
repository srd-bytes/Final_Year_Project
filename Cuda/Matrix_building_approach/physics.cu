#include <stdio.h>
#include <stdlib.h>
#include <cuda_runtime.h>
#include <math.h>
#include "utility.cu"




__global__ void b_vector_1D(float* b, float phi_0, float phi_N1){
    // 1D blocks and threads

    float delta = 0.2f;
    float epsilon = 8.8f;

    int index = blockDim.x*blockIdx.x+threadIdx.x;

    float charge_density = 1.0f; // interpret charge density in pC / m
    // charge_density = sin(index* PI/8);

    b[index] = charge_density* delta*delta/epsilon; // charge density* delta^2/epsilon; 
    b[0]+=phi_0;
    b[N-1]+=phi_N1;

    __syncthreads();
    // printf("b[%d] = %f\n",index,b[index]);
}

__global__ void circular_shift_rows_right(
    float** A,
    int N)
{
    //dim3 threads(64);
    // dim3 blocks2((N+threads.x-1)/threads.x, N);

    int row = blockIdx.y;                              // row index i
    int col = blockIdx.x * blockDim.x + threadIdx.x;  // column index j

    if (row < N && col < N) {
        int src_col = (col - row + N) % N;
        A[row][col] = A[0][src_col];
    }
}
__global__ void adjustment1(float** A){
    // 1 block 2 thread
    if(threadIdx.x == 0){
        A[0][1] = 1.0f;
    }
    if(threadIdx.x == 1){
        A[0][N-1] = 1.0f;
    }
}
__global__ void adjustment2(float** A){
    // 1 block 2 thread
    if(threadIdx.x == 0){
        A[0][N-1] = 0.0f;
    }
    if(threadIdx.x == 1){
        A[N-1][0] = 0.0f;
    }
}

void A_matrix_1D(int N, float** A){
    
    dim3 blocks(N/Tile_width,N/Tile_width);
    dim3 threads_per_block(Tile_width,Tile_width);

    initialize_matrix_by_const<<<blocks,threads_per_block>>>(N,A,0.0f);
    adjustment1<<<1,2>>>(A);

    dim3 threads(64);
    dim3 blocks2((N+threads.x-1)/threads.x, N);

    circular_shift_rows_right<<<blocks2,threads>>>(A,N);
    adjustment2<<<1,2>>>(A);

    // print_matrix<<<blocks,threads_per_block>>>(A); //perfectly works
    
    
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