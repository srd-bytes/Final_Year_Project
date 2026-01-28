#include <stdio.h>
#include <stdlib.h>
#include <cuda_runtime.h>

//-----------------Global Variables ------------
const int Tile_width = 2;
const int N = 4 ;

//-----------------CPU specific code------------

// ---------------- GPU KERNELS ----------------

__global__ void initialize_matrix(int N, float** A){
    // here using shared memory is useless . Also I need to spawn 2D blocks and threads
    int index_col = blockDim.x*blockIdx.x+threadIdx.x;
    int index_row = blockDim.y*blockIdx.y+threadIdx.y;

    A[index_row][index_col] = N*index_row+index_col +1;
    
    printf("A[%d][%d] = %f\n",index_row,index_col,A[index_row][index_col]);
}

__global__ void initialize_vector(int N, float* v){
    int index = blockDim.x*blockIdx.x+threadIdx.x;
    v[index] = index +1;
    
    printf("v[%d] = %f\n",index,v[index]);
}


__global__ void matrix_vector_product(int N, float** A, float* v, float* result){
    
    __shared__ float shared_A[Tile_width][Tile_width];
    __shared__ float shared_v[Tile_width];

    float value=0;
    for(int phase=0; phase < N/Tile_width; phase++){
        shared_v[threadIdx.x] = v[phase*Tile_width+threadIdx.x];
        for(int k=0; k<Tile_width; k++){
            shared_A[threadIdx.x][k] = A[blockIdx.x*blockDim.x+threadIdx.x][phase*Tile_width + k ];
        
        }
        __syncthreads();
        for(int k=0; k<Tile_width; k++){
            value += shared_A[threadIdx.x][k]*shared_v[k];
        }
    }
    __syncthreads();

    result[blockIdx.x*blockDim.x+threadIdx.x] = value;
    __syncthreads();
    printf("result[%d] = %f\n",blockIdx.x*blockDim.x+threadIdx.x,result[blockIdx.x*blockDim.x+threadIdx.x]);
}


//--------------------Test Code----------------
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

    //result
    float* result_gpu;
    cudaMalloc((void**)&result_gpu, N*sizeof(float));
    matrix_vector_product<<<N/Tile_width,Tile_width>>>(N,A_gpu,v_gpu,result_gpu);

    cudaDeviceSynchronize();


    // Deallocation
    for(int i = 0; i < N; i++) {
        cudaFree(A_gpu[i]);
    }
    cudaFree(A_gpu);
    free(A_cpu);
    cudaFree(v_gpu);
    return 0;
}
