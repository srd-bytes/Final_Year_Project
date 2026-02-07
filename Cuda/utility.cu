#include <stdio.h>
#include <stdlib.h>
#include <cuda_runtime.h>

//-----------------Global Variables ------------
// const int Tile_width = 8;
// const int N = 128 ;

// ----------------------------------GPU kernels----------------------------------------------------------------------

// --------------Initialization----------------------------------------------------------------

__global__ void initialize_2Dmatrix_by_const(float* matrix_2D,float value){ // z,y,x for y=0 each x , flat nomenclature stacking the rows
    int index_x = blockDim.x*blockIdx.x+threadIdx.x;
    int index_y = blockDim.y*blockIdx.y+threadIdx.y;
    int n= gridDim.x*blockDim.x;
    // v[index] = index +1;
    int flat_index= n*index_y + index_x;
    matrix_2D[flat_index] = value;
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
__global__ void equate_flat_matrix2D(float* phi_old, float* phi_new){

    int index_x = blockDim.x*blockIdx.x+threadIdx.x;
    int index_y = blockDim.y*blockIdx.y+threadIdx.y;
    int n= gridDim.x*blockDim.x;
    int flat_index= n*index_y + index_x;

    phi_old[flat_index]= phi_new[flat_index];

}
__global__ void equate_matrix3D(float*** phi_old, float*** phi_new){

    int index_x = blockDim.x*blockIdx.x+threadIdx.x;
    int index_y = blockDim.y*blockIdx.y+threadIdx.y;
    int index_z = blockDim.z*blockIdx.z+threadIdx.z;

    phi_old[index_z][index_y][index_x]= phi_new[index_z][index_y][index_x];

}
// ---------------------------------------------------------------------------------------------------------------

// ----------------------------------CPU functions--------------------------------------------------------
void  write_to_csv_1d(float* data, int N,const char *path){
    FILE* fp = fopen(path, "w");

    for (int i = 0; i < N; i++) {
        
        fprintf(fp, "%f", data[i]);
        fprintf(fp, "\n");
    }

    fclose(fp);
}
void  write_to_csv_2d(float* data, int N,const char *path){
    FILE* fp = fopen(path, "w");

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            fprintf(fp, "%f", data[N*i+j]);

            if (j < N - 1)
                fprintf(fp, ",");
        }
        fprintf(fp, "\n");
    }

    fclose(fp);
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
