#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mpi.h"
#define N 2 // size of array (np must be <= N)
#include "structures.c"

// simple for 1D but becomes combination of loops etc and complicated for multidimensions
// transmitter and receiver functions are more like organize and send or receive and decode
// we dont require to use to wait here while sending so request is only needed for validating
// send - receive - wait cycle



void jacobi_1d_parallel(float phi[N], Physics physics, Simulation sim , Boundary boundary ,int argc, char **argv){

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
    float* partial_result = (float*) calloc(block.x,sizeof(float));
    float* partial_result_new = (float*) calloc(block.x,sizeof(float)); // this were initialized in every grid

    
    

    for(int i=0; i< sim.iteration; i++){


        // Receive
        for(int j=0; j< block.x; j++){

            float left;
            float right;
            
            // MPI_Request request[2];

            if(rank==0){
                
                left= (j==0) ? boundary.left : partial_result[j-1];
                // communicate the right boundary from rank =1
                if(j==block.x-1){
                    if(grid.x>1){
                        
                        MPI_Sendrecv(&partial_result[j], 1, MPI_FLOAT, rank+1, 51, &right, 1, MPI_FLOAT, rank+1, 49, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                    } else {right = boundary.right;}
                
                } else {right = partial_result[j+1];}

            } else if(rank==grid.x-1){
                
                
                right= (j==block.x-1) ? boundary.right : partial_result[j+1];
                // communicate the left boundary from rank=2
                if(j==0){
                    if(grid.x>1){
                        
                        MPI_Sendrecv(&partial_result[j], 1, MPI_FLOAT, rank-1, 49, &left, 1, MPI_FLOAT, rank-1, 51, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                    } else {left = boundary.left;}
                } else {left = partial_result[j-1];}

            } else {
                // communicate the left and right if j==0 and j==process_size-1
                if(grid.x>1){
                    // communicate the right boundary
                    if(j==block.x-1){
                        
                        MPI_Sendrecv(&partial_result[j], 1, MPI_FLOAT, rank+1, 51, &right, 1, MPI_FLOAT, rank+1, 49, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                    } else {right = partial_result[j+1];}

                    // communicate the left boundary
                    if(j==0){
                        
                        MPI_Sendrecv(&partial_result[j], 1, MPI_FLOAT, rank-1, 49, &left, 1, MPI_FLOAT, rank-1, 51, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                    } else {left = partial_result[j-1];}
                }
            }
            partial_result_new[j] = (left + right + b)/2;

            

        }
        // ------------------------------------------------
        for(int y=0; y<block.x; y++){partial_result[y]=partial_result_new[y];}
        // ------------------------------------------------
    }
    
    // ---------------Communicating the partial results (Blocking)-------------------------------------------
    if(rank!=0){
        MPI_Send(partial_result_new, block.x, MPI_FLOAT, 0, 99, MPI_COMM_WORLD);
    } else {
        for(int u=0;u<block.x;u++){
            phi[u]=partial_result_new[u];
        }
        for(int k=1; k<np; k++){
            
            float* buf = (float*) calloc(block.x,sizeof(float));
            MPI_Recv(buf,block.x, MPI_FLOAT,k,99,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            for(int u=0;u<block.x;u++){
                phi[k*block.x+u]=buf[u];
            }
        }
    }
    free(partial_result);
    free(partial_result_new);
    // -------------------------------------------------------------------------------------------------------

    // --------------------------Testing---------------------------------------------
    if(rank ==0){

        
        for(int i=0; i<N;i++){
            printf("phi[%d] = %f\n",i,phi[i]);
        }
        
    }
    // ------------------------------------------------------------------------------

    MPI_Finalize();
}


int main(int argc, char **argv) {
    
    // printf("hey");
    // fflush(stdout);
    float phi[N];

    Physics physics = {1, 1};
    Simulation simulation = {1,1};
    Boundary boundary = {1.0f,-1.0f,0,0 };

    jacobi_1d_parallel(phi, physics, simulation, boundary ,argc, argv); // I am getting wrong answer for num processor < 4
    // Try using Sendrecv instead of Send and Recv
    return 0;
    
}
