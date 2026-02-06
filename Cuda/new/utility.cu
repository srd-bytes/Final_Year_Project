#include <stdio.h>
#include <stdlib.h>
#include <cuda_runtime.h>

//-----------------Global Variables ------------
// const int Tile_width = 8;
// const int N = 128 ;


// --------------Initialization----------------------------------------------------------------

__global__ void initialize_2Dmatrix_by_const(float** matrix_2D, float value){
    int index_x = blockDim.x*blockIdx.x+threadIdx.x;
    int index_y = blockDim.y*blockIdx.y+threadIdx.y;
    // v[index] = index +1;
    matrix_2D[index_y][index_x] = value;
    // printf("v[%d] = %f\n",index,v[index]);
}

__global__ void initialize_3Dmatrix_by_const(float*** matrix_3D, float value){
    int index_x = blockDim.x*blockIdx.x+threadIdx.x;
    int index_y = blockDim.y*blockIdx.y+threadIdx.y;
    int index_z = blockDim.z*blockIdx.z+threadIdx.z;
    // v[index] = index +1;
    matrix_3D[index_z][index_y][index_x] = value;
    // printf("v[%d] = %f\n",index,v[index]);
}
__global__ void initialize_vector_by_const(float* v, float value){
    int index = blockDim.x*blockIdx.x+threadIdx.x;
    // v[index] = index +1;
    v[index] = value;
    // printf("v[%d] = %f\n",index,v[index]);
}




// -------------------------Equate------------------------------------------------------------

__global__ void equate_vectors(float* phi_old, float* phi_new){

    int index= blockDim.x*blockIdx.x + threadIdx.x;

    phi_old[index]= phi_new[index];

}
__global__ void equate_matrix2D(float** phi_old, float** phi_new){

    int index_x = blockDim.x*blockIdx.x+threadIdx.x;
    int index_y = blockDim.y*blockIdx.y+threadIdx.y;

    phi_old[index_y][index_x]= phi_new[index_y][index_x];

}
__global__ void equate_matrix3D(float*** phi_old, float*** phi_new){

    int index_x = blockDim.x*blockIdx.x+threadIdx.x;
    int index_y = blockDim.y*blockIdx.y+threadIdx.y;
    int index_z = blockDim.z*blockIdx.z+threadIdx.z;

    phi_old[index_z][index_y][index_x]= phi_new[index_z][index_y][index_x];

}


//--------------------Test Code----------------
// int main(void)

// {
//     //Matrix allocation
//     float* A_cpu[N];
//     for(int i=0;i<N ; i++){
//         cudaMalloc((void**)&A_cpu[i], N* sizeof(float));
//     }
    
//     float** A_gpu;
//     cudaMalloc((void**)&A_gpu, N*sizeof(float*) );
//     cudaMemcpy(A_gpu, A_cpu, N*sizeof(float*), cudaMemcpyHostToDevice);
    

//     dim3 blocks(N/Tile_width,N/Tile_width);
//     dim3 threads_per_block(Tile_width,Tile_width);
//     initialize_matrix<<<blocks,threads_per_block>>>(N,A_gpu);

//     //Vector part
//     float* v_gpu;
//     cudaMalloc((void**)&v_gpu, N*sizeof(float));
//     initialize_vector<<<N/Tile_width, Tile_width>>>(N,v_gpu);

//     //Matrix Vector Product
//     matrix_vector_product<<<N/Tile_width,Tile_width>>>(N,A_gpu,v_gpu);

//     cudaDeviceSynchronize();


//     // Deallocation
//     for(int i = 0; i < N; i++) {
//         cudaFree(A_gpu[i]);
//     }
//     cudaFree(A_gpu);
//     free(A_cpu);
//     cudaFree(v_gpu);
    
//     return 0;
// }
