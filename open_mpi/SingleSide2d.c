#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mpi.h"
#define N 2 // size of array
#include "structures.c"
#include <time.h>

void jacobi_2d_parallel(float phi[N*N], Physics physics, Simulation sim , Boundary boundary ,int argc, char **argv){
    
    

    MPI_Init(&argc, &argv);

    int np; // Must be power of 4 to avoid any complications
    MPI_Comm_size(MPI_COMM_WORLD, &np);
    

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // Generalizing for any number of processors
    int dims[2];
    MPI_Dims_create(np,2,dims);

     
    Cartesian_int  grid = {dims[1], dims[0]};
    // printf("grid.x = %d, grid.y = %d\n", grid.x, grid.y);
    Cartesian_int block= {N/grid.x, N/grid.y}; // defined in structures.c
    // printf("block.x = %d, block.y = %d\n", block.x, block.y);

    // Blocks and grids are now structs

    // ---------------physics parameter------------------------

    float b= physics.uniform_charge_density*sim.step_size*sim.step_size/physics.epsilon;

    // ---------------------------------------------------------

    // ---------------Main Calculation---------------------------
    
    float* partial_result;
    MPI_Alloc_mem(block.x*block.y*sizeof(float), MPI_INFO_NULL, &partial_result);
    for(int i=0;i<block.x*block.y;i++){partial_result[i]=0;} // initialization

    float* partial_result_new; // this were initialized in every grid
    MPI_Alloc_mem(block.x*block.y*sizeof(float), MPI_INFO_NULL, &partial_result_new); // exposable 
    //free this using MPI_Free_mem(a);

    // Exposing everything
    MPI_Win win;
    MPI_Win_create(partial_result, block.x*block.y*sizeof(float), sizeof(float),MPI_INFO_NULL, MPI_COMM_WORLD, &win);

    //--------------Measuring performance--------------
    clock_t start, end;
    double cpu_time_used;

    start = clock();
    //-------------------------------------------------
    float left;
    float right;
    float up;
    float down;

    float* buffer; // this were initialized in every grid
    MPI_Alloc_mem(block.x*block.y*sizeof(float), MPI_INFO_NULL, &buffer); // exposable 

    for(int iter=0; iter<sim.iteration; iter++){

        
        for(int y=0; y<block.y; y++){
            for(int x=0; x< block.x; x++){

                MPI_Win_fence(0, win);
                
                if(x==0){
                    if(rank % grid.x == 0){
                        left=boundary.left;
                    } else {
                        MPI_Get(buffer, block.x*block.y, MPI_FLOAT, rank-1, 0, block.x*block.y, MPI_FLOAT, win);
                        left = buffer[y*block.x + block.x-1]; // last portion of each row
                    }
                } else {left = partial_result[y*block.x+x-1];}

                if(x==block.x-1){
                    if(rank % grid.x == grid.x - 1){
                        right=boundary.right;
                    } else {
                        MPI_Get(buffer, block.x*block.y, MPI_FLOAT, rank+1, 0, block.x*block.y, MPI_FLOAT, win);
                        right = buffer[y*block.x]; // first portion of each row
                    }
                } else {right = partial_result[y*block.x+x+1];}
                

                

                if(y==0){
                    if(rank < grid.x){
                        up=boundary.up;
                    } else {
                        MPI_Get(buffer, block.x*block.y, MPI_FLOAT, rank-grid.x, 0, block.x*block.y, MPI_FLOAT, win);
                        up = buffer[block.x*(block.y-1) + x]; // last row
                    }
                } else {up = partial_result[y*block.x + x-block.x];}

                if(y==block.y-1){
                    if(rank >= np - grid.x){
                        down=boundary.down;
                    } else {
                        MPI_Get(buffer, block.x*block.y, MPI_FLOAT, rank+grid.x, 0, block.x*block.y, MPI_FLOAT, win);
                        down = buffer[x]; // first row
                    }
                } else {down = partial_result[y*block.x + x + block.x];}
                

                MPI_Win_fence(0, win);
                partial_result_new[y*block.x+x] = (left+right+ up+down+b)/4;
                
            }
        }
        
        for(int u=0;u<block.y;u++){for(int v=0;v<block.x;v++){partial_result[u*block.x+v]=partial_result_new[u*block.x+v];}}
    }
    
    end = clock();
    cpu_time_used= ((double) (end - start))/CLOCKS_PER_SEC;
    
    if(rank==0){printf("Execution time = %f seconds\n", cpu_time_used);}

    MPI_Win_fence(0, win);

    if(rank==0){
        for(int i = 0; i < block.y; i++){
                for(int j = 0; j < block.x; j++){
                    phi[i*N + j] = partial_result_new[i*block.x + j];
                }
            }

        float* buf;
        MPI_Alloc_mem(block.x*block.y*sizeof(float), MPI_INFO_NULL, &buf);

        for(int k=1;k<np;k++){
            MPI_Get(buf, block.x*block.y, MPI_FLOAT, k, 0, block.x*block.y, MPI_FLOAT, win); 

            // compute 2D position
            int row = k / grid.x;
            int col = k % grid.x;

            // place block correctly
            for(int i = 0; i < block.y; i++){
                for(int j = 0; j < block.x; j++){

                    int global_i = row * block.y + i;
                    int global_j = col * block.x + j;

                    phi[global_i * N + global_j] = buf[i * block.x + j];
                }
            }
        }
        MPI_Free_mem(buf);
    
    }
    MPI_Win_fence(0, win);

    //-----------------------Testing-----------------------
    if(rank==0){
        for(int i=0; i<N*N;i++){
            printf("phi[%d] = %f\n",i,phi[i]);
        }
    }
    //-----------------------------------------------------

    MPI_Win_free(&win);
    MPI_Free_mem(partial_result);
    MPI_Free_mem(buffer);
    MPI_Free_mem(partial_result_new);
    // -------------------------------------------------------------------------------------------------------

    MPI_Finalize();


}

int main(int argc, char **argv) {
    
    // printf("hey");
    // fflush(stdout);
    
    static float phi[N*N];

    Physics physics = {1, 1};
    Simulation simulation = {1,1};
    Boundary boundary = {1.0f,-1.0f,0,0 };

    jacobi_2d_parallel(phi, physics, simulation, boundary ,argc, argv); // I am getting wrong answer for num processor < 4
    // Try using Sendrecv instead of Send and Recv
    return 0;
    
}