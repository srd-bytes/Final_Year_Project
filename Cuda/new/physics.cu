#include <stdio.h>
#include <stdlib.h>
#include <cuda_runtime.h>
#include <math.h>




__global__ void physics_solution_1D(float* phi_old, float* phi_new, int N){

    int index= blockDim.x*blockIdx.x + threadIdx.x;

    float delta = 0.2f;
    float epsilon = 8.8f;
    float charge_density= 1.0f;

    float b= charge_density*delta*delta/epsilon;

    float phi_0 = 10.0f;
    float phi_N1 = -10.0f;

    float left= (index==0) ? phi_0: phi_old[index-1];   // condition ? if true : else
    float right= (index==N-1) ? phi_N1 :phi_old[index+1];

    
    phi_new[index] = (left + right + b)/2;
        
}

__global__ void physics_solution_2D(float** phi_old, float** phi_new, int N){

    int index_x = blockDim.x*blockIdx.x+threadIdx.x; //column
    int index_y = blockDim.y*blockIdx.y+threadIdx.y; // rows

    float delta = 0.2f;
    float epsilon = 8.8f;
    float charge_density= 1.0f;

    float b= charge_density*delta*delta/epsilon;

    // --------------------Boundary Values------------------------------------------------------
    float phi_left_boundary = 10.0f;
    float phi_right_boundary = -10.0f;
    float phi_up_boundary = 0.0f;
    float phi_down_boundary = 0.0f;
    // -------------------------------------------------------------------------------------------

    float left= (index_x==0) ? phi_left_boundary : phi_old[index_y][index_x-1];
    float right = (index_x == N-1)? phi_right_boundary :phi_old[index_y][index_x+1];
    float up = (index_y == 0)? phi_up_boundary :phi_old[index_y-1][index_x];
    float down= (index_y == N-1) ? phi_down_boundary :phi_old[index_y+1][index_x];

    phi_new[index_y][index_x] = (left + right + up + down + b)/4 ;
}

__global__ void physics_solution_3D(float*** phi_old, float*** phi_new, int N){

    int index_x = blockDim.x*blockIdx.x+threadIdx.x; //column
    int index_y = blockDim.y*blockIdx.y+threadIdx.y; // rows
    int index_z = blockDim.z*blockIdx.z+threadIdx.z;

    float delta = 0.2f;
    float epsilon = 8.8f;
    float charge_density= 1.0f;

    float b= charge_density*delta*delta/epsilon;

    // ------------------Boundary Values-------------------------------------------------------
    float phi_left_boundary = 10.0f;
    float phi_right_boundary = -10.0f;

    float phi_back_boundary = 0.0f; // y+1
    float phi_front_boundary = 0.0f; // y-1

    float phi_top_boundary = 0.0f;
    float phi_bottom_boundary = 0.0f;
    // --------------------------------------------------------------------------------------

    // ---------------------------------Assignment--------------------------------------------------------

    float left= (index_x==0) ? phi_left_boundary : phi_old[index_z][index_y][index_x-1];
    float right = (index_x == N-1)? phi_right_boundary :phi_old[index_z][index_y][index_x+1];

    float front = (index_y == 0)? phi_front_boundary :phi_old[index_z][index_y-1][index_x];
    float back= (index_y == N-1) ? phi_back_boundary :phi_old[index_z][index_y+1][index_x];


    float bottom = (index_z == 0)? phi_bottom_boundary :phi_old[index_z-1][index_y][index_x];
    float top= (index_z == N-1) ? phi_top_boundary :phi_old[index_z+1][index_y][index_x];
    // ---------------------------------------------------------------------------------------------------

    phi_new[index_z][index_y][index_x] = (left + right + back + front + bottom + top + b)/6 ;
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