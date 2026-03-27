#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mpi.h"
#define N 2 // size of array
#include "structures.c"

// simple for 1D but becomes combination of loops etc and complicated for multidimensions
// transmitter and receiver functions are more like organize and send or receive and decode
// we dont require to use to wait here while sending so request is only needed for validating
// send - receive - wait cycle

void left_end_transmitter_1D(float* partial_result , int process_size, int rank){
    MPI_Request request = MPI_REQUEST_NULL;
    MPI_Isend(&partial_result[process_size-1], 1,MPI_FLOAT,rank+1,51,MPI_COMM_WORLD, &request);
}
void right_end_transmitter_1D(float* partial_result , int process_size, int rank){
    MPI_Request request = MPI_REQUEST_NULL;
    MPI_Isend(&partial_result[0], 1,MPI_FLOAT,rank-1,49,MPI_COMM_WORLD, &request);
}
void mids_transmitter_1D(float* partial_result , int process_size, int rank){
    MPI_Request request = MPI_REQUEST_NULL;
    MPI_Isend(&partial_result[0], 1,MPI_FLOAT,rank-1,49,MPI_COMM_WORLD, &request);
    MPI_Isend(&partial_result[process_size-1], 1,MPI_FLOAT,rank+1,51,MPI_COMM_WORLD, &request);
}



void jacobi_1d_parallel(float phi[N], Physics physics, Simulation sim , Boundary boundary ,int argc, char **argv){

    MPI_Init(&argc, &argv);

    int np;
    MPI_Comm_size(MPI_COMM_WORLD, &np);
    

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int process_size = N/np;

    // ---------------physics parameters------------------------

    float b= physics.uniform_charge_density*sim.step_size*sim.step_size/physics.epsilon;

    // ---------------------------------------------------------

    // ---------------Main Calculation---------------------------
    float* partial_result = (float*) calloc(process_size,sizeof(float));
    float* partial_result_new = (float*) calloc(process_size,sizeof(float));

    
    

    for(int i=0; i< sim.iteration; i++){

        
        // non blocking send all the necessary details
        if(rank==0){
            
            left_end_transmitter_1D(partial_result, process_size, rank);

        } else if(rank==np-1){
        
            right_end_transmitter_1D(partial_result, process_size, rank);

        } else {
            
            mids_transmitter_1D(partial_result, process_size, rank);
            
        }

        // Receive
        for(int j=0; j< process_size; j++){

            float left;
            float right;
            
            MPI_Request request[2];

            if(rank==0){
                
                left= (j==0) ? boundary.left : partial_result[j-1];
                // communicate the right boundary from rank =1
                if(j==process_size-1){
                    MPI_Irecv(&right, 1,MPI_FLOAT, rank+1, 49, MPI_COMM_WORLD, &request[0]);
                    MPI_Wait(&request[0],MPI_STATUS_IGNORE);

                } else {right = partial_result[j+1];}
                
                
            } else if(rank==np-1){
                
                right= (j==process_size-1) ? boundary.right : partial_result[j+1];
                // communicate the left boundary from rank=2
                if(j==0){
                    
                    MPI_Irecv(&left, 1,MPI_FLOAT, rank-1, 51, MPI_COMM_WORLD, &request[0]);
                    MPI_Wait(&request[0],MPI_STATUS_IGNORE);

                } else {left = partial_result[j-1];}

            } else {
                // communicate the left and right if j==0 and j==process_size-1
                
                if(j==0){
                    MPI_Irecv(&left, 1,MPI_FLOAT, rank-1, 51, MPI_COMM_WORLD, &request[0]);
                    MPI_Wait(&request[0],MPI_STATUS_IGNORE);
                } else {left = partial_result[j-1];}

                if(j==process_size-1){
                    MPI_Irecv(&right, 1,MPI_FLOAT, rank+1, 49, MPI_COMM_WORLD, &request[1]);
                    MPI_Wait(&request[1],MPI_STATUS_IGNORE);
                } else {right = partial_result[j+1];}
                
            }
            partial_result_new[j] = (left + right + b)/2;

            // ------------------------------------------------
            for(int y=0; y<process_size; y++){partial_result[y]=partial_result_new[y];}
            // ------------------------------------------------

        }
    }
    
    // ---------------Communicating the partial results (Blocking)-------------------------------------------
    if(rank!=0){
        MPI_Send(partial_result_new, process_size, MPI_FLOAT, 0, 99, MPI_COMM_WORLD);
    } else {
        for(int u=0;u<process_size;u++){
            phi[u]=partial_result_new[u];
        }
        for(int k=1; k<np; k++){
            
            float* buf = (float*) calloc(process_size,sizeof(float));
            MPI_Recv(buf,process_size, MPI_FLOAT,k,99,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            for(int u=0;u<process_size;u++){
                phi[k*process_size+u]=buf[u];
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



//---------------------------2 DIMENSIONS-------------------------------------------------------------------

void transmit_corners(float* partial_result, float* buffer_vertical, float* buffer_horizontal ,int rank, int block_size, int grid_size){
    // for corners you need 2 buffers of size process_size
    //---left corner-----------------
    MPI_Request request = MPI_REQUEST_NULL;

    if(rank==0){
        
        for(int i=0; i< block_size; i++){
            // down
            buffer_vertical[i] = partial_result[block_size*(block_size-1)+i];
            // right
            buffer_horizontal[i] = partial_result[i*block_size + block_size-1];
        }
        
        MPI_Isend(buffer_vertical,block_size, MPI_FLOAT, rank + grid_size ,51, MPI_COMM_WORLD, &request);
        MPI_Isend(buffer_horizontal,block_size, MPI_FLOAT, rank+1 ,11, MPI_COMM_WORLD, &request);

    } else if(rank==grid_size -1){
        
        // right corner
        for(int i=0; i< block_size; i++){
            // down
            buffer_vertical[i] = partial_result[block_size*(block_size-1)+i];
            // left
            buffer_horizontal[i] = partial_result[i*block_size];
        }
        
        MPI_Isend(buffer_vertical,block_size, MPI_FLOAT, rank + grid_size ,51, MPI_COMM_WORLD, &request);
        MPI_Isend(buffer_horizontal,block_size, MPI_FLOAT, rank-1 ,9, MPI_COMM_WORLD, &request);

    } else if(rank==grid_size*(grid_size-1)){
        
        // bottom left corner
        for(int i=0; i< block_size; i++){
            // up
            buffer_vertical[i] = partial_result[i];
            // right
            buffer_horizontal[i] = partial_result[i*block_size + block_size-1];
        }
        
        MPI_Isend(buffer_vertical,block_size, MPI_FLOAT, rank - grid_size ,49, MPI_COMM_WORLD, &request);
        MPI_Isend(buffer_horizontal,block_size, MPI_FLOAT, rank+1 ,11, MPI_COMM_WORLD, &request);
    } else if(rank==grid_size*(grid_size)-1){
        
        // bottom right corner
        for(int i=0; i< block_size; i++){
            // up
            buffer_vertical[i] = partial_result[i];
            // left
            buffer_horizontal[i] = partial_result[i*block_size];
        }
        
        MPI_Isend(buffer_vertical,block_size, MPI_FLOAT, rank - grid_size ,49, MPI_COMM_WORLD, &request);
        MPI_Isend(buffer_horizontal,block_size, MPI_FLOAT, rank-1 ,9, MPI_COMM_WORLD, &request);
    }
}

void transmit_edges(float* partial_result, float* buffer_vertical_1, float* buffer_vertical_2, float* buffer_horizontal_1, float* buffer_horizontal_2,int rank, int block_size, int grid_size){
    
    MPI_Request request = MPI_REQUEST_NULL;

    for(int k=1; k< grid_size - 1; k++ ){

        // upper edge
        if(rank == k){
            for(int i=0; i< block_size; i++){
                // down
                buffer_vertical_2[i] = partial_result[block_size*(block_size-1)+i];
                // left
                buffer_horizontal_1[i] = partial_result[i*block_size];
                // right
                buffer_horizontal_2[i] = partial_result[i*block_size + block_size-1];
            }
        
            MPI_Isend(buffer_vertical_2, block_size, MPI_FLOAT, rank + grid_size ,51, MPI_COMM_WORLD, &request); //down
            MPI_Isend(buffer_horizontal_1,block_size, MPI_FLOAT, rank-1 ,9, MPI_COMM_WORLD, &request); //left
            MPI_Isend(buffer_horizontal_2,block_size, MPI_FLOAT, rank+1 ,11, MPI_COMM_WORLD, &request); //right
        }

        // bottom edge
        if(rank == grid_size*(grid_size-1) +k){
            for(int i=0; i< block_size; i++){
                // up
                buffer_vertical_1[i] = partial_result[i];
                // left
                buffer_horizontal_1[i] = partial_result[i*block_size];
                // right
                buffer_horizontal_2[i] = partial_result[i*block_size + block_size-1];
            }
        
            MPI_Isend(buffer_vertical_1, block_size, MPI_FLOAT, rank - grid_size ,49, MPI_COMM_WORLD, &request); //up
            MPI_Isend(buffer_horizontal_1,block_size, MPI_FLOAT, rank-1 ,9, MPI_COMM_WORLD, &request); //left
            MPI_Isend(buffer_horizontal_2,block_size, MPI_FLOAT, rank+1 ,11, MPI_COMM_WORLD, &request); //right
        }

        // left edge
        if(rank == k*grid_size){
            for(int i=0; i< block_size; i++){
                // up
                buffer_vertical_1[i] = partial_result[i];
                // down
                buffer_vertical_2[i] = partial_result[block_size*(block_size-1)+i];
                // right
                buffer_horizontal_2[i] = partial_result[i*block_size + block_size-1];
            }
        
            MPI_Isend(buffer_vertical_1, block_size, MPI_FLOAT, rank - grid_size ,49, MPI_COMM_WORLD, &request); //up
            MPI_Isend(buffer_vertical_2, block_size, MPI_FLOAT, rank + grid_size ,51, MPI_COMM_WORLD, &request); //down
            MPI_Isend(buffer_horizontal_2,block_size, MPI_FLOAT, rank+1 ,11, MPI_COMM_WORLD, &request); //right
        }

        // right edge
        if(rank == k*grid_size -1){
            for(int i=0; i< block_size; i++){
                // up
                buffer_vertical_1[i] = partial_result[i];
                // down
                buffer_vertical_2[i] = partial_result[block_size*(block_size-1)+i];
                // left
                buffer_horizontal_1[i] = partial_result[i*block_size];
            }
        
            MPI_Isend(buffer_vertical_1, block_size, MPI_FLOAT, rank - grid_size ,49, MPI_COMM_WORLD, &request); //up
            MPI_Isend(buffer_vertical_2, block_size, MPI_FLOAT, rank + grid_size ,51, MPI_COMM_WORLD, &request); //down
            MPI_Isend(buffer_horizontal_1,block_size, MPI_FLOAT, rank-1 ,9, MPI_COMM_WORLD, &request); //left
        }
    }
}

void transmit_rest(float* partial_result, float* buffer_vertical_1, float* buffer_vertical_2, float* buffer_horizontal_1, float* buffer_horizontal_2,int rank, int block_size, int grid_size){

    MPI_Request request = MPI_REQUEST_NULL;

    for(int k=grid_size+1; k< 2*grid_size-1; k++){
        for(int j=k; j< k+ grid_size*(grid_size-2); j=j+grid_size){
            if(rank==j){

                for(int i=0; i< block_size; i++){
                    // up
                    buffer_vertical_1[i] = partial_result[i];
                    // down
                    buffer_vertical_2[i] = partial_result[block_size*(block_size-1)+i];
                    // left
                    buffer_horizontal_1[i] = partial_result[i*block_size];
                    // right
                    buffer_horizontal_2[i] = partial_result[i*block_size + block_size-1];
                }
            
                MPI_Isend(buffer_vertical_1, block_size, MPI_FLOAT, rank - grid_size ,49, MPI_COMM_WORLD, &request); //up
                MPI_Isend(buffer_vertical_2, block_size, MPI_FLOAT, rank + grid_size ,51, MPI_COMM_WORLD, &request); //down
                MPI_Isend(buffer_horizontal_1,block_size, MPI_FLOAT, rank-1 ,9, MPI_COMM_WORLD, &request); //left
                MPI_Isend(buffer_horizontal_2,block_size, MPI_FLOAT, rank+1 ,11, MPI_COMM_WORLD, &request); //right
            }
        }
    }
}


void jacobi_2d_parallel(float phi[N*N], Physics physics, Simulation sim , Boundary boundary ,int argc, char **argv){

    MPI_Init(&argc, &argv);

    int np; // Must be power of 4
    MPI_Comm_size(MPI_COMM_WORLD, &np);
    

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int grid_size = (int) sqrt(np);
    int block_size = N/grid_size;

    // ---------------physics parameter------------------------

    float b= physics.uniform_charge_density*sim.step_size*sim.step_size/physics.epsilon;

    // ---------------------------------------------------------

    // ---------------Main Calculation---------------------------
    float* partial_result = (float*) calloc(block_size*block_size,sizeof(float));
    float* partial_result_new = (float*) calloc(block_size*block_size,sizeof(float));
    
    
    

    for(int iter=0; iter< sim.iteration; iter++){

        
        // non blocking send all the necessary details
        float* buffer_vertical_1 = (float*) calloc(block_size,sizeof(float));
        float* buffer_vertical_2 = (float*) calloc(block_size,sizeof(float));
        float* buffer_horizontal_1 = (float*) calloc(block_size,sizeof(float));
        float* buffer_horizontal_2 = (float*) calloc(block_size,sizeof(float));

        // unoptimized some buffers are unused but eating up space in some ranks. Optimize later

        //------------------sending--------------------------------
        transmit_corners(partial_result, buffer_vertical_1, buffer_horizontal_1, rank, block_size, grid_size);
        transmit_edges(partial_result, buffer_vertical_1, buffer_vertical_2, buffer_horizontal_1, buffer_horizontal_2, rank, block_size, grid_size);
        transmit_rest (partial_result, buffer_vertical_1, buffer_vertical_2, buffer_horizontal_1, buffer_horizontal_2, rank, block_size, grid_size);
        
        
        
    //-----------------------Receiving--------------------------------------
        

        float* rec_buffer_vertical_1 = (float*) calloc(block_size,sizeof(float));
        float* rec_buffer_vertical_2 = (float*) calloc(block_size,sizeof(float));
        float* rec_buffer_horizontal_1 = (float*) calloc(block_size,sizeof(float));
        float* rec_buffer_horizontal_2 = (float*) calloc(block_size,sizeof(float));
        
        MPI_Request request[4];



        // ---------------------------------Corners receive-------------------------------------------------------

        if(rank== 0 || rank== grid_size*(grid_size-1)){ // left corners horizontal
            
            MPI_Irecv(rec_buffer_horizontal_2, block_size, MPI_FLOAT, rank+1, 9, MPI_COMM_WORLD, &request[0]);
            MPI_Wait(&request[0], MPI_STATUS_IGNORE);
            
            free(buffer_horizontal_1);
            free(buffer_horizontal_2);

        } else if(rank== grid_size-1 || rank== grid_size*grid_size-1){ // right corners horizontal

            
            MPI_Irecv(rec_buffer_horizontal_1, block_size, MPI_FLOAT, rank-1, 11, MPI_COMM_WORLD, &request[0]);
            MPI_Wait(&request[0], MPI_STATUS_IGNORE);
            free(buffer_horizontal_1);
            free(buffer_horizontal_2);

        } else if(rank== 0 || rank== grid_size-1){ // up corners vertical

            
            MPI_Irecv(rec_buffer_vertical_2, block_size, MPI_FLOAT, rank+ grid_size, 49, MPI_COMM_WORLD, &request[1]);
            MPI_Wait(&request[1], MPI_STATUS_IGNORE);
            free(buffer_vertical_1);
            free(buffer_vertical_2);

        } else if(rank== grid_size*(grid_size-1) || rank== grid_size*grid_size-1){ // down corners vertical

            
            MPI_Irecv(rec_buffer_vertical_1, block_size, MPI_FLOAT, rank- grid_size, 51, MPI_COMM_WORLD, &request[1]);
            MPI_Wait(&request[1], MPI_STATUS_IGNORE);
            free(buffer_vertical_1);
            free(buffer_vertical_2);
        }
        // -------------------------------------------------------------------------------------------------------------

        // -------------------------------------------------------------------------------------------------------------
        //                                               Edges Receive
        // -------------------------------------------------------------------------------------------------------------

        for(int k=1; k< grid_size - 1; k++ ){

            
            // upper edge
            if(rank == k){
                
                
                MPI_Irecv(rec_buffer_vertical_2, block_size, MPI_FLOAT, rank + grid_size ,49, MPI_COMM_WORLD, &request[0]); //down
                MPI_Irecv(rec_buffer_horizontal_1,block_size, MPI_FLOAT, rank-1 ,11, MPI_COMM_WORLD, &request[1]); //left (data coming from previous rank right edge(tag=11))
                MPI_Irecv(rec_buffer_horizontal_2,block_size, MPI_FLOAT, rank+1 ,9, MPI_COMM_WORLD, &request[2]); //right

                MPI_Wait(&request[0], MPI_STATUS_IGNORE);
                MPI_Wait(&request[1], MPI_STATUS_IGNORE);
                MPI_Wait(&request[2], MPI_STATUS_IGNORE);

                free(buffer_horizontal_1);
                free(buffer_horizontal_2);
                free(buffer_vertical_1);
                free(buffer_vertical_2);
            }

            // bottom edge
            if(rank == grid_size*(grid_size-1) +k){
            
                MPI_Irecv(rec_buffer_vertical_1, block_size, MPI_FLOAT, rank - grid_size ,51, MPI_COMM_WORLD, &request[0]); //up
                MPI_Irecv(rec_buffer_horizontal_1,block_size, MPI_FLOAT, rank-1 ,11, MPI_COMM_WORLD, &request[1]); //left
                MPI_Irecv(rec_buffer_horizontal_2,block_size, MPI_FLOAT, rank+1 ,9, MPI_COMM_WORLD, &request[2]); //right

                MPI_Wait(&request[0], MPI_STATUS_IGNORE);
                MPI_Wait(&request[1], MPI_STATUS_IGNORE);
                MPI_Wait(&request[2], MPI_STATUS_IGNORE);

                free(buffer_horizontal_1);
                free(buffer_horizontal_2);
                free(buffer_vertical_1);
                free(buffer_vertical_2);
            }

            // left edge
            if(rank == k*grid_size){
            
                MPI_Irecv(rec_buffer_vertical_1, block_size, MPI_FLOAT, rank - grid_size ,51, MPI_COMM_WORLD, &request[0]); //up
                MPI_Irecv(rec_buffer_vertical_2, block_size, MPI_FLOAT, rank + grid_size ,49, MPI_COMM_WORLD, &request[1]); //down
                MPI_Irecv(rec_buffer_horizontal_2,block_size, MPI_FLOAT, rank+1 ,9, MPI_COMM_WORLD, &request[2]); //right

                MPI_Wait(&request[0], MPI_STATUS_IGNORE);
                MPI_Wait(&request[1], MPI_STATUS_IGNORE);
                MPI_Wait(&request[2], MPI_STATUS_IGNORE);

                free(buffer_horizontal_1);
                free(buffer_horizontal_2);
                free(buffer_vertical_1);
                free(buffer_vertical_2);
            }

            // right edge
            if(rank == k*grid_size -1){

                MPI_Irecv(rec_buffer_vertical_1, block_size, MPI_FLOAT, rank - grid_size ,51, MPI_COMM_WORLD, &request[0]); //up
                MPI_Irecv(rec_buffer_vertical_2, block_size, MPI_FLOAT, rank + grid_size ,49, MPI_COMM_WORLD, &request[1]); //down
                MPI_Irecv(rec_buffer_horizontal_1,block_size, MPI_FLOAT, rank-1 ,11, MPI_COMM_WORLD, &request[2]); //left
                MPI_Wait(&request[0], MPI_STATUS_IGNORE);
                MPI_Wait(&request[1], MPI_STATUS_IGNORE);
                MPI_Wait(&request[2], MPI_STATUS_IGNORE);

                free(buffer_horizontal_1);
                free(buffer_horizontal_2);
                free(buffer_vertical_1);
                free(buffer_vertical_2);
            }
        }
        //-------------------------------------------------------------------------------------------------------------

        // -------------------------------------------------------------------------------------------------------------
        //                                            Rest Receive
        // -------------------------------------------------------------------------------------------------------------

        for(int k=grid_size+1; k< 2*grid_size-1; k++){
            for(int j=k; j< k+ grid_size*(grid_size-2); j=j+grid_size){
                if(rank==j){
                    
                    MPI_Irecv(rec_buffer_vertical_1, block_size, MPI_FLOAT, rank - grid_size ,51, MPI_COMM_WORLD, &request[0]); //up
                    MPI_Irecv(rec_buffer_vertical_2, block_size, MPI_FLOAT, rank + grid_size ,49, MPI_COMM_WORLD, &request[1]); //down
                    MPI_Irecv(rec_buffer_horizontal_1,block_size, MPI_FLOAT, rank-1 ,11, MPI_COMM_WORLD, &request[2]); //left
                    MPI_Irecv(rec_buffer_horizontal_2,block_size, MPI_FLOAT, rank+1 ,9, MPI_COMM_WORLD, &request[3]); //right

                    MPI_Wait(&request[0], MPI_STATUS_IGNORE);
                    MPI_Wait(&request[1], MPI_STATUS_IGNORE);
                    MPI_Wait(&request[2], MPI_STATUS_IGNORE);
                    MPI_Wait(&request[3], MPI_STATUS_IGNORE);

                    free(buffer_horizontal_1);
                    free(buffer_horizontal_2);
                    free(buffer_vertical_1);
                    free(buffer_vertical_2);
                }
            }
        }
        // ------------------------------------------------------------------------------------------------------------
        
        // -------------------------------------------------------------------------------------------------------------
        //                                            Corners Calculate
        // -------------------------------------------------------------------------------------------------------------

        float left;
        float right;
        float up;
        float down;


        
        for(int i=0; i< block_size; i++){
            for(int j=0; j<block_size; j++){


                // -----------------------Handling Corners----------------------------------------------------------
                if(rank==0){

                    if(j==0){left=boundary.left;} else {left = partial_result[i*block_size +j -1];}
                    if(i==0){up=boundary.up;} else {up = partial_result[i*block_size +j -block_size];}
                    if(j==block_size-1){right=rec_buffer_horizontal_2[i]; } else {right = partial_result[i*block_size +j +1];}
                    if(i==block_size-1){down=rec_buffer_vertical_2[j]; } else {down = partial_result[i*block_size +j + block_size];}

                } else if (rank == grid_size-1){ // top right

                    if(j==0){left=rec_buffer_horizontal_1[i];} else {left = partial_result[i*block_size +j -1];} // left (j==0)
                    if(i==0){up=boundary.up;} else {up = partial_result[i*block_size +j -block_size];} // Up(i==0)
                    if(j==block_size-1){right=boundary.right; } else {right = partial_result[i*block_size +j +1];} // right(j==block_size-1)
                    if(i==block_size-1){down=rec_buffer_vertical_2[j]; } else {down = partial_result[i*block_size +j + block_size];} // down (i==block_size-1)

                } else if (rank == grid_size*(grid_size-1)){ // bottom left

                    if(j==0){left=boundary.left;} else {left = partial_result[i*block_size +j -1];} // left (j==0)
                    if(i==0){up=rec_buffer_vertical_1[j];} else {up = partial_result[i*block_size +j -block_size];} // Up(i==0)
                    if(j==block_size-1){right=rec_buffer_horizontal_2[i]; } else {right = partial_result[i*block_size +j +1];} // right(j==block_size-1)
                    if(i==block_size-1){down=boundary.down; } else {down = partial_result[i*block_size +j + block_size];} // down (i==block_size-1)

                } else if (rank == grid_size*grid_size-1){ // bottom right

                    if(j==0){left=rec_buffer_horizontal_1[i];} else {left = partial_result[i*block_size +j -1];} // left (j==0)
                    if(i==0){up=rec_buffer_vertical_1[j];} else {up = partial_result[i*block_size +j -block_size];} // Up(i==0)
                    if(j==block_size-1){right=boundary.right; } else {right = partial_result[i*block_size +j +1];} // right(j==block_size-1)
                    if(i==block_size-1){down=boundary.down; } else {down = partial_result[i*block_size +j + block_size];} // down (i==block_size-1)
                }
                // -------------------------------------------------------------------------------------------------------------------------------------------

                // -------------------------------EDGES ----------------------------------------------------------------------------------------------
                for(int k=1; k< grid_size - 1; k++ ){
                    
                    // upper edge
                    if(rank == k){
                        if(j==0){left=rec_buffer_horizontal_1[i];} else {left = partial_result[i*block_size +j -1];} // left (j==0)
                        if(i==0){up=boundary.up;} else {up = partial_result[i*block_size +j -block_size];} // Up(i==0)
                        if(j==block_size-1){right=rec_buffer_horizontal_2[i]; } else {right = partial_result[i*block_size +j +1];} // right(j==block_size-1)
                        if(i==block_size-1){down=rec_buffer_vertical_2[j]; } else {down = partial_result[i*block_size +j + block_size];} // down (i==block_size-1)
                    }

                    // down edge
                    if(rank == grid_size*(grid_size-1) +k){

                        if(j==0){left=rec_buffer_horizontal_1[i];} else {left = partial_result[i*block_size +j -1];} // left (j==0)
                        if(i==0){up=rec_buffer_vertical_1[j];} else {up = partial_result[i*block_size +j -block_size];} // Up(i==0)
                        if(j==block_size-1){right=rec_buffer_horizontal_2[i]; } else {right = partial_result[i*block_size +j +1];} // right(j==block_size-1)
                        if(i==block_size-1){down=boundary.down; } else {down = partial_result[i*block_size +j + block_size];} // down (i==block_size-1)
                    
                    }

                    // left edge
                    if(rank == k*grid_size){

                        if(j==0){left=boundary.left;} else {left = partial_result[i*block_size +j -1];} // left (j==0)
                        if(i==0){up=rec_buffer_vertical_1[j];} else {up = partial_result[i*block_size +j -block_size];} // Up(i==0)
                        if(j==block_size-1){right=rec_buffer_horizontal_2[i]; } else {right = partial_result[i*block_size +j +1];} // right(j==block_size-1)
                        if(i==block_size-1){down=rec_buffer_vertical_2[j]; } else {down = partial_result[i*block_size +j + block_size];} // down (i==block_size-1)
                    
                    }

                    // right edge
                    if(rank == k*grid_size -1){

                        if(j==0){left=rec_buffer_horizontal_1[i];} else {left = partial_result[i*block_size +j -1];} // left (j==0)
                        if(i==0){up=rec_buffer_vertical_1[j];} else {up = partial_result[i*block_size +j -block_size];} // Up(i==0)
                        if(j==block_size-1){right=boundary.right; } else {right = partial_result[i*block_size +j +1];} // right(j==block_size-1)
                        if(i==block_size-1){down=rec_buffer_vertical_2[j]; } else {down = partial_result[i*block_size +j + block_size];} // down (i==block_size-1)
                    
                    }
                }
                //-----------------------------------------------------------------------------------------------------------------------------------------------------

                //------------------------------------Rest-------------------------------------------------------------------------------------------------
                for(int k=grid_size+1; k< 2*grid_size-1; k++){
                    for(int t=k; t< k+ grid_size*(grid_size-2); t=t+grid_size){
                        if(rank==t){
                            if(j==0){left=rec_buffer_horizontal_1[i];} else {left = partial_result[i*block_size +j -1];} // left (j==0)
                            if(i==0){up=rec_buffer_vertical_1[j];} else {up = partial_result[i*block_size +j -block_size];} // Up(i==0)
                            if(j==block_size-1){right= rec_buffer_horizontal_2[i]; } else {right = partial_result[i*block_size +j +1];} // right(j==block_size-1)
                            if(i==block_size-1){down=rec_buffer_vertical_2[j]; } else {down = partial_result[i*block_size +j + block_size];} // down (i==block_size-1)
                            
                            
                        }
                    }
                }
                //--------------------------------------------------------------------------------------------------------------
                //---------------------------------Calculation-----------------------------------------------
                partial_result_new[i*block_size +j] = (left + right+ up + down + b)/4;
                partial_result[i*block_size +j] = partial_result_new[i*block_size +j];
            }
            
        }


    }
    
    // ---------------Communicating the partial results (Blocking)-------------------------------------------
    if(rank!=0){
        MPI_Send(partial_result_new, block_size*block_size, MPI_FLOAT, 0, 99, MPI_COMM_WORLD);
    } else {
        for(int u=0;u<block_size*block_size;u++){
            phi[u]=partial_result_new[u];
        }
        for(int k=1; k<np; k++){
            
            float* buf = (float*) calloc(block_size*block_size,sizeof(float));
            MPI_Recv(buf,block_size*block_size, MPI_FLOAT,k,99,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            for(int u=0;u<block_size*block_size;u++){
                phi[k*block_size*block_size+u]=buf[u];
            }
        }
    }
    free(partial_result);
    free(partial_result_new);
    // -------------------------------------------------------------------------------------------------------

    // --------------------------Testing---------------------------------------------
    if(rank ==0){

        
        for(int p=0; p<N*N;p++){
            printf("phi[%d] = %f\n",p,phi[p]);
        }
        
    }
    // ------------------------------------------------------------------------------

    MPI_Finalize();
}

int main(int argc, char **argv) {
    
    // printf("hey");
    // fflush(stdout);
    float phi[N*N];

    Physics physics = {1, 1};
    Simulation simulation = {1,1};
    Boundary boundary = {1.0f,-1.0f,0,0 };

    jacobi_2d_parallel(phi, physics, simulation, boundary ,argc, argv);
    return 0;
    
}
