#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mpi.h"
#define N 1024 // size of array (np must be <= N)
#include "structures.c"
#include <time.h>

// simple for 1D but becomes combination of loops etc and complicated for multidimensions
// transmitter and receiver functions are more like organize and send or receive and decode
// we dont require to use to wait here while sending so request is only needed for validating
// send - receive - wait cycle



void jacobi_1d_parallel(float* phi , Physics physics, Simulation sim , Boundary boundary ,int argc, char **argv){

    clock_t start, end;
    double cpu_time_used;

    start = clock();

    MPI_Init(&argc, &argv);

    int np;
    MPI_Comm_size(MPI_COMM_WORLD, &np);
    

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    

    Cartesian_int grid;
    Cartesian_int block;

    grid.x=np;
    block.x = N/grid.x;

    // ---------------physics parameters------------------------

    float b= physics.uniform_charge_density*sim.step_size*sim.step_size/physics.epsilon;

    // ---------------------------------------------------------

    // ---------------Main Calculation---------------------------
    float* partial_result;
    MPI_Alloc_mem(block.x*sizeof(float), MPI_INFO_NULL, &partial_result);
    for(int i=0;i<block.x;i++){partial_result[i]=0;} // initialization

    float* partial_result_new; // this were initialized in every grid
    MPI_Alloc_mem(block.x*sizeof(float), MPI_INFO_NULL, &partial_result_new); // exposable 
    //free this using MPI_Free_mem(a);

    MPI_Win win; // you can do better here
    // free this using MPI_Win_free(&win);

    // Exposing partial result
    
    MPI_Win_create(partial_result, block.x*sizeof(float), sizeof(float),MPI_INFO_NULL, MPI_COMM_WORLD, &win);

    
    for(int iter=0; iter< sim.iteration; iter++){


        
        for(int j=0; j< block.x; j++){

            float left;
            float right;
            //  code
            MPI_Win_fence(0, win);
            // using get you can get the data of exposed parts
            if(j==block.x-1){
                if(rank!= grid.x-1){
                    
                    MPI_Get(&right, 1, MPI_FLOAT, rank+1, 0, 1, MPI_FLOAT, win);
                   

                } else { right=boundary.right;}
                
            } else { right = partial_result[j+1];}

            if(j==0){
                if(rank!= 0){
                    MPI_Get(&left, 1, MPI_FLOAT, rank-1, block.x-1, 1, MPI_FLOAT, win);
                    
                } else { left=boundary.left; }
                
            } else {left = partial_result[j-1];}

            MPI_Win_fence(0, win);

            partial_result_new[j] = (left+right+b)/2;


            
            
        }

        for(int k=0;k<block.x;k++){partial_result[k]=partial_result_new[k];}

    }
    
    
    end = clock();
    cpu_time_used= ((double) (end - start))/CLOCKS_PER_SEC;
    
    if(rank==0){printf("Execution time = %f seconds\n", cpu_time_used);}

    MPI_Win_fence(0, win);
    if(rank==0){
        for(int u=0;u<block.x;u++){
            phi[u]=partial_result_new[u];
        }

        
        for(int k=1; k<np; k++){
            MPI_Get(&phi[k*block.x], block.x, MPI_FLOAT, k, 0, block.x, MPI_FLOAT, win);   
        }
    
    }
    MPI_Win_fence(0, win);

    //-----------------------Testing-----------------------
    // if(rank==0){
    //     for(int i=0; i<N;i++){
    //         printf("phi[%d] = %f\n",i,phi[i]);
    //     }
    // }
    //-----------------------------------------------------

    MPI_Win_free(&win);
    MPI_Free_mem(partial_result);
    MPI_Free_mem(partial_result_new);
    // -------------------------------------------------------------------------------------------------------

    MPI_Finalize();
}

// Think about 2D case. Later you do the optimization

int main(int argc, char **argv) {
    
    // printf("hey");
    // fflush(stdout);
    
    static float phi[N];

    Physics physics = {1, 1};
    Simulation simulation = {1,10};
    Boundary boundary = {1.0f,-1.0f,0,0 };

    jacobi_1d_parallel(phi, physics, simulation, boundary ,argc, argv); // I am getting wrong answer for num processor < 4
    // Try using Sendrecv instead of Send and Recv
    return 0;
    
}
