#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mpi.h"
#define N 2 // size of array
#include "structures.c"

void jacobi_2d_parallel(float phi[N*N], Physics physics, Simulation sim , Boundary boundary ,int argc, char **argv){

    MPI_Init(&argc, &argv);

    int np; // Must be power of 4
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
    float* partial_result = (float*) calloc(block.x*block.y,sizeof(float));
    float* partial_result_new = (float*) calloc(block.x*block.y,sizeof(float));
    
    
    

    for(int iter=0; iter< sim.iteration; iter++){

        
        // non blocking send all the necessary details

        Buffer buffer; // transmitting buffer

        buffer.vertical_1 = (float*) calloc(block.x,sizeof(float));
        buffer.vertical_2 = (float*) calloc(block.x,sizeof(float));
        buffer.horizontal_1 = (float*) calloc(block.y,sizeof(float));
        buffer.horizontal_2 = (float*) calloc(block.y,sizeof(float));
        
        
        
    
        
        Buffer rec_buffer; // Receiving buffer
        rec_buffer.vertical_1 = (float*) calloc(block.x,sizeof(float));
        rec_buffer.vertical_2 = (float*) calloc(block.x,sizeof(float));
        rec_buffer.horizontal_1 = (float*) calloc(block.y,sizeof(float));
        rec_buffer.horizontal_2 = (float*) calloc(block.y,sizeof(float));
        
        // MPI_Request request[4];



        // ---------------------------------Corners-------------------------------------------------------
        if(rank == 0){
            for(int i=0; i< fmax(block.x,block.y); i++){

                if(i<block.x && grid.y>1){
                    // down
                    buffer.vertical_2[i] = partial_result[block.x*(block.y-1)+i];
                }
                
                if(i<block.y && grid.x>1){
                    // right
                    buffer.horizontal_2[i] = partial_result[i*block.x + block.x-1];
                }
            
            }
            //down
            if(grid.y>1){MPI_Sendrecv(buffer.vertical_2 ,block.x, MPI_FLOAT, rank + grid.x ,51, rec_buffer.vertical_1, block.x,MPI_FLOAT, rank+grid.x, 49,MPI_COMM_WORLD, MPI_STATUS_IGNORE);}
            //right
            if(grid.x>1){MPI_Sendrecv(buffer.horizontal_2,block.y, MPI_FLOAT, rank+1 ,11,  rec_buffer.horizontal_1, block.y,MPI_FLOAT, rank+1, 9,MPI_COMM_WORLD, MPI_STATUS_IGNORE);}

        } else if(rank==grid.x -1){ // upper right corner
        
        
            for(int i=0; i< fmax(block.x,block.y); i++){

                if(i<block.x && grid.y>1){
                    // down
                    buffer.vertical_2[i] = partial_result[block.x*(block.y-1)+i];
                }
                
                if(i<block.y && grid.x>1){
                    // left
                    buffer.horizontal_1[i] = partial_result[i*block.x];
                }
                
            }
            
            if(grid.y>1){MPI_Sendrecv(buffer.vertical_2,block.x, MPI_FLOAT, rank + grid.x ,51, rec_buffer.vertical_1, block.x,MPI_FLOAT, rank+grid.x, 49,MPI_COMM_WORLD, MPI_STATUS_IGNORE);}
            if(grid.x>1){MPI_Sendrecv(buffer.horizontal_1,block.y, MPI_FLOAT, rank -1 ,9, rec_buffer.horizontal_2, block.y,MPI_FLOAT, rank-1, 11,MPI_COMM_WORLD, MPI_STATUS_IGNORE);}

        } else if(rank==grid.x*(grid.y-1)){ // bottom left corner
        
            for(int i=0; i< fmax(block.x, block.y); i++){

                if(i<block.x && grid.y>1){
                    // up
                    buffer.vertical_1[i] = partial_result[i];
                }
                
                if(i<block.y && grid.x>1){
                    // right
                    buffer.horizontal_2[i] = partial_result[i*block.x + block.x-1];
                }
                
                
            }
            
            if(grid.y>1){MPI_Sendrecv(buffer.vertical_1,block.x, MPI_FLOAT, rank -grid.x ,49, rec_buffer.vertical_2, block.x,MPI_FLOAT, rank-grid.x, 51,MPI_COMM_WORLD, MPI_STATUS_IGNORE);}
            if(grid.x>1){MPI_Sendrecv(buffer.horizontal_2,block.y, MPI_FLOAT, rank +1 ,11, rec_buffer.horizontal_1, block.y,MPI_FLOAT, rank+1, 9,MPI_COMM_WORLD, MPI_STATUS_IGNORE);} //right

        } else if(rank==grid.x*(grid.y)-1){ // bottom right corner
        
        
            for(int i=0; i< fmax(block.x, block.y); i++){

                if(i<block.x && grid.y>1){
                    // up
                    buffer.vertical_1[i] = partial_result[i];
                }
                
                if(i<block.y && grid.x>1){
                    // left
                    buffer.horizontal_1[i] = partial_result[i*block.x];
                }
                
                
            }
            
            if(grid.y>1){MPI_Sendrecv(buffer.vertical_1,block.x, MPI_FLOAT, rank -grid.x ,49, rec_buffer.vertical_2, block.x,MPI_FLOAT, rank-grid.x, 51,MPI_COMM_WORLD, MPI_STATUS_IGNORE);} //up
            if(grid.x>1){MPI_Sendrecv(buffer.horizontal_1,block.y, MPI_FLOAT, rank -1 ,9, rec_buffer.horizontal_2, block.y,MPI_FLOAT, rank-1, 11,MPI_COMM_WORLD, MPI_STATUS_IGNORE);} //left
        }

        //---------------------------------------Edges------------------------------------------------------------

        for(int k=1; k< fmax(grid.x-1,grid.y-1); k++ ){ // k is meant to be the increment over initial rank

            // upper edge
            if((k<grid.x-1) && (rank == k) ){
                for(int i=0; i< fmax(block.x,block.y); i++){

                    if(i<block.x && grid.y>1){
                        // down
                        buffer.vertical_2[i] = partial_result[block.x*(block.y-1)+i];
                    }
                    
                    if(i<block.y && grid.x>1){
                        // left
                        buffer.horizontal_1[i] = partial_result[i*block.x];
                        // right
                        buffer.horizontal_2[i] = partial_result[i*block.x + block.x-1];
                    }
                    
                    
                }
                if(grid.y>1){MPI_Sendrecv(buffer.vertical_2,block.x, MPI_FLOAT, rank + grid.x ,51, rec_buffer.vertical_1, block.x,MPI_FLOAT, rank+grid.x, 49,MPI_COMM_WORLD, MPI_STATUS_IGNORE);} //down
                if(grid.x>1){ 
                    MPI_Sendrecv(buffer.horizontal_1,block.y, MPI_FLOAT, rank -1 ,9, rec_buffer.horizontal_2, block.y,MPI_FLOAT, rank-1, 11,MPI_COMM_WORLD, MPI_STATUS_IGNORE); //left
                    MPI_Sendrecv(buffer.horizontal_2,block.y, MPI_FLOAT, rank +1 ,11, rec_buffer.horizontal_1, block.y,MPI_FLOAT, rank+1, 9,MPI_COMM_WORLD, MPI_STATUS_IGNORE); //right
                }

            }
            // bottom edge
            if((k<grid.x-1) && (rank == grid.x*(grid.y-1) +k)){

                for(int i=0; i< fmax(block.x, block.y); i++){

                    if(i<block.x && grid.y>1){
                        // up
                        buffer.vertical_1[i] = partial_result[i];
                        
                    }
                    
                    if(i<block.y && grid.x> 1){
                        // left
                        buffer.horizontal_1[i] = partial_result[i*block.x];
                        // right
                        buffer.horizontal_2[i] = partial_result[i*block.x + block.x-1];
                    }

                    
                    
                }
            
                if(grid.y>1){MPI_Sendrecv(buffer.vertical_1,block.x, MPI_FLOAT, rank -grid.x ,49, rec_buffer.vertical_2, block.x,MPI_FLOAT, rank-grid.x, 51,MPI_COMM_WORLD, MPI_STATUS_IGNORE);} //up
                if(grid.x>1){
                    MPI_Sendrecv(buffer.horizontal_1,block.y, MPI_FLOAT, rank -1 ,9, rec_buffer.horizontal_2, block.y,MPI_FLOAT, rank-1, 11,MPI_COMM_WORLD, MPI_STATUS_IGNORE); //left
                    MPI_Sendrecv(buffer.horizontal_2,block.y, MPI_FLOAT, rank +1 ,11, rec_buffer.horizontal_1, block.y,MPI_FLOAT, rank+1, 9,MPI_COMM_WORLD, MPI_STATUS_IGNORE); //right
                }
            }

            // left edge
            if((k<grid.y-1) && (rank == k*grid.x)){
                for(int i=0; i< fmax(block.x, block.y); i++){

                    if(i<block.x && grid.y>1){
                        // up
                        buffer.vertical_1[i] = partial_result[i];
                        // down
                        buffer.vertical_2[i] = partial_result[block.x*(block.y-1)+i];
                    }
                    
                    if(i<block.y && grid.x>1){
                        // right
                        buffer.horizontal_2[i] = partial_result[i*block.x + block.x-1];
                    }

                    
                }
            
                if(grid.y>1){
                    MPI_Sendrecv(buffer.vertical_1,block.x, MPI_FLOAT, rank -grid.x ,49, rec_buffer.vertical_2, block.x,MPI_FLOAT, rank-grid.x, 51,MPI_COMM_WORLD, MPI_STATUS_IGNORE);//up
                    MPI_Sendrecv(buffer.vertical_2,block.x, MPI_FLOAT, rank + grid.x ,51, rec_buffer.vertical_1, block.x,MPI_FLOAT, rank+grid.x, 49,MPI_COMM_WORLD, MPI_STATUS_IGNORE); //down
                }
                if(grid.x>1){
                    MPI_Sendrecv(buffer.horizontal_2,block.y, MPI_FLOAT, rank +1 ,11, rec_buffer.horizontal_1, block.y,MPI_FLOAT, rank+1, 9,MPI_COMM_WORLD, MPI_STATUS_IGNORE); //right
                }
            }

            // right edge
            if((k<grid.y-1) && (rank == k*grid.x -1)){
                for(int i=0; i< fmax(block.x, block.y); i++){

                    if(i<block.x && grid.y>1){
                        // up
                        buffer.vertical_1[i] = partial_result[i];
                        // down
                        buffer.vertical_2[i] = partial_result[block.x*(block.y-1)+i];
                    }
                    
                    if(i<block.y && grid.x>1){
                        // left
                        buffer.horizontal_1[i] = partial_result[i*block.x];
                    }

                }
            
                if(grid.y>1){
                    MPI_Sendrecv(buffer.vertical_1,block.x, MPI_FLOAT, rank -grid.x ,49, rec_buffer.vertical_2, block.x,MPI_FLOAT, rank-grid.x, 51,MPI_COMM_WORLD, MPI_STATUS_IGNORE);//up
                    MPI_Sendrecv(buffer.vertical_2,block.x, MPI_FLOAT, rank + grid.x ,51, rec_buffer.vertical_1, block.x,MPI_FLOAT, rank+grid.x, 49,MPI_COMM_WORLD, MPI_STATUS_IGNORE); //down
                }
                if(grid.x>1){
                    MPI_Sendrecv(buffer.horizontal_1,block.y, MPI_FLOAT, rank -1 ,9, rec_buffer.horizontal_2, block.y,MPI_FLOAT, rank-1, 11,MPI_COMM_WORLD, MPI_STATUS_IGNORE); //left
                }
            }
        }

        //---------------------------------Rest----------------------------------------------------------------------

        for(int k=grid.x+1; k< 2*grid.x-1; k++){
            for(int j=k; j< k + grid.x*(grid.y-2); j=j+grid.x){
                if(rank==j){

                    for(int i=0; i< fmax(block.x,block.y); i++){

                        if(i<block.x && grid.y>1){
                            // up
                            buffer.vertical_1[i] = partial_result[i];
                            // down
                            buffer.vertical_2[i] = partial_result[block.x*(block.y-1)+i];
                        }
                        
                        if(i<block.y && grid.x>1){
                            // left
                            buffer.horizontal_1[i] = partial_result[i*block.x];
                            // right
                            buffer.horizontal_2[i] = partial_result[i*block.x + block.x-1]; 
                        }
                        
                    }
                
                    if(grid.y>1){
                        MPI_Sendrecv(buffer.vertical_1,block.x, MPI_FLOAT, rank -grid.x ,49, rec_buffer.vertical_2, block.x,MPI_FLOAT, rank-grid.x, 51,MPI_COMM_WORLD, MPI_STATUS_IGNORE);//up
                        MPI_Sendrecv(buffer.vertical_2,block.x, MPI_FLOAT, rank + grid.x ,51, rec_buffer.vertical_1, block.x,MPI_FLOAT, rank+grid.x, 49,MPI_COMM_WORLD, MPI_STATUS_IGNORE); //down
                    }
                    if(grid.x>1){
                        MPI_Sendrecv(buffer.horizontal_1,block.y, MPI_FLOAT, rank -1 ,9, rec_buffer.horizontal_2, block.y,MPI_FLOAT, rank-1, 11,MPI_COMM_WORLD, MPI_STATUS_IGNORE); //left
                        MPI_Sendrecv(buffer.horizontal_2,block.y, MPI_FLOAT, rank +1 ,11, rec_buffer.horizontal_1, block.y,MPI_FLOAT, rank+1, 9,MPI_COMM_WORLD, MPI_STATUS_IGNORE); //right
                    }
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

        for(int i=0; i< block.y; i++){ 
            for(int j=0; j<block.x; j++){


                // -----------------------Handling Corners----------------------------------------------------------
                if(rank==0){ //top left

                    if(j==0){left=boundary.left;} else {left = partial_result[i*block.x +j -1];}
                    if(i==0){up=boundary.up;} else {up = partial_result[i*block.x +j -block.x];}

                    if(j==block.x-1){
                        if(grid.x>1){right=rec_buffer.horizontal_1[i]; } else {right=boundary.right;}
                    } else {right = partial_result[i*block.x +j +1];} //right

                    if(i==block.y-1){
                        if(grid.y>1){down=rec_buffer.vertical_1[j]; } else {down=boundary.down;}
                    } else {down = partial_result[i*block.x +j + block.x];} //down

                } else if (rank == grid.x-1){ // top right
                    
                    if(j==0){
                        if(grid.x>1){left=rec_buffer.horizontal_2[i]; } else {left=boundary.left;}
                    } else {left = partial_result[i*block.x +j -1];} //left

                    if(i==0){up=boundary.up;} else {up = partial_result[i*block.x +j -block.x];} // Up(i==0)

                    if(j==block.x-1){right=boundary.right; } else {right = partial_result[i*block.x +j +1];} // right(j==block.x-1)

                    if(i==block.y-1){
                        if(grid.y>1){down=rec_buffer.vertical_1[j]; } else {down=boundary.down;} // down (i==block.y-1)
                    } else {down = partial_result[i*block.x +j + block.x];}

                } else if (rank == grid.x*(grid.y-1)){ // bottom left

                    if(j==0){left=boundary.left;} else {left = partial_result[i*block.x +j -1];} // left (j==0)
                    // if(i==0){
                    //     up=rec_buffer.vertical_2[j]; // Up(i==0)
                    // } else {up = partial_result[i*block.x +j -block.x];}
                    if(i==0){
                        if(grid.y>1){up=rec_buffer.vertical_2[j]; } else {up=boundary.up;}
                    } else {up = partial_result[i*block.x +j - block.x];} //up

                    if(j==block.x-1){
                        if(grid.x>1){right=rec_buffer.horizontal_2[i]; } else {right=boundary.right;} // right(j==block-1)
                    } else {right = partial_result[i*block.x +j +1];}

                    if(i==block.y-1){down=boundary.down; } else {down = partial_result[i*block.x +j + block.x];} // down (i==block-1)

                } else if (rank == grid.x*grid.y-1){ // bottom right
                    
                    
                    if(j==0){
                        if(grid.x>1){left=rec_buffer.horizontal_2[i]; } else {left=boundary.left;}
                    } else {left = partial_result[i*block.x +j -1];} //left
                
                    if(i==0){
                        if(grid.y>1){up=rec_buffer.vertical_2[j]; } else {up=boundary.up;}
                    } else {up = partial_result[i*block.x +j - block.x];} //up

                    
                    if(j==block.x-1){right=boundary.right; } else {right = partial_result[i*block.x +j +1];} // right(j==block-1)
                    if(i==block.y-1){down=boundary.down; } else {down = partial_result[i*block.x +j + block.x];} // down (i==block-1)
                    
                }
                //------------------------Handling Edges----------------------------------------------------

                for(int k=1; k< fmax(grid.x-1, grid.y-1); k++ ){
                    
                    // upper edge
                    if(k<grid.x-1 && rank == k){
                        if(j==0){
                            if(grid.x>1){left=rec_buffer.horizontal_1[i];} else {left=boundary.left;} // left (j==0)
                        } else {left = partial_result[i*block.x +j -1];}
                        if(i==0 && grid.y>1){up=boundary.up;} else {up = partial_result[i*block.x +j -block.x];} // Up(i==0)
                        if(j==block.x-1){
                            if(grid.x>1){right=rec_buffer.horizontal_2[i]; } else {right=boundary.right;} // right(j==block_size-1)
                        } else {right = partial_result[i*block.x +j +1];}
                        if(i==block.y-1){
                            if(grid.y>1){down=rec_buffer.vertical_2[j]; } else {down=boundary.down;} // down (i==block_size-1)
                        } else {down = partial_result[i*block.x +j + block.x];}
                    }

                    // down edge
                    if(k<grid.x-1 && rank == grid.x*(grid.y-1) +k){
                        if(j==0){
                            if(grid.x>1){left=rec_buffer.horizontal_1[i];} else {left=boundary.left;} // left (j==0)
                        } else {left = partial_result[i*block.x +j -1];}
                        if(i==0){
                            if(grid.y>1){up=rec_buffer.vertical_1[j];} else {up=boundary.up;} // Up(i==0)
                        } else {up = partial_result[i*block.x +j -block.x];}
                        if(j==block.x-1){
                            if(grid.x>1){right=rec_buffer.horizontal_2[i]; } else {right=boundary.right;} // right(j==block_size-1)
                        } else {right = partial_result[i*block.x +j +1];}
                        if(i==block.y-1 && grid.y>1){down=boundary.down; } else {down = partial_result[i*block.x +j + block.x];} // down (i==block_size-1)
                    
                    }

                    // left edge
                    if(k<grid.y-1 && rank == k*grid.x){

                        if(j==0 && grid.x>1){left=boundary.left;} else {left = partial_result[i*block.x +j -1];} // left (j==0)
                        if(i==0){
                            if(grid.y>1){up=rec_buffer.vertical_1[j];} else {up=boundary.up;} // Up(i==0)
                        } else {up = partial_result[i*block.x +j -block.x];}
                        if(j==block.x-1){
                            if(grid.x>1){right=rec_buffer.horizontal_2[i]; } else {right=boundary.right;} // right(j==block_size-1)
                        } else {right = partial_result[i*block.x +j +1];}
                        if(i==block.y-1){
                            if(grid.y>1){down=rec_buffer.vertical_2[j]; } else {down=boundary.down;} // down (i==block_size-1)
                        } else {down = partial_result[i*block.x +j + block.x];}
                    
                    }

                    // right edge
                    if(k<grid.y-1 && rank == k*grid.x -1){
                        if(j==0){
                            if(grid.x>1){left=rec_buffer.horizontal_1[i];} else {left=boundary.left;} // left (j==0)
                        } else {left = partial_result[i*block.x +j -1];}
                        if(i==0){
                            if(grid.y>1){up=rec_buffer.vertical_1[j];} else {up=boundary.up;} // Up(i==0)
                        } else {up = partial_result[i*block.x +j -block.x];}
                        
                        if(j==block.x-1 && grid.x>1){right=boundary.right; } else {right = partial_result[i*block.x +j +1];} // right(j==block_size-1)
                        
                        if(i==block.y-1){
                            if(grid.y>1){down=rec_buffer.vertical_2[j]; } else {down=boundary.down;} // down (i==block_size-1)
                        } else {down = partial_result[i*block.x +j + block.x];}
                    
                    }
                }
                //-----------------------------------------------------------------------------------------------------------------------------------------------------
                //------------------------------------Rest-------------------------------------------------------------------------------------------------
                for(int k=grid.x+1; k< 2*grid.x-1; k++){
                    for(int t=k; t< k+ grid.x*(grid.y-2); t=t+grid.x){
                        if(rank==t){
                            if(j==0){
                                if(grid.x>1){left=rec_buffer.horizontal_1[i];} else {left=boundary.left;} // left (j==0)
                            } else {
                                left = partial_result[i*block.x +j -1];
                            }
                            if(i==0){
                                if(grid.y>1){up=rec_buffer.vertical_1[j];} else {up=boundary.up;} // Up(i==0)
                            } else {
                                up = partial_result[i*block.x +j -block.x];
                            }
                            if(j==block.x-1){
                                if(grid.x>1){right= rec_buffer.horizontal_2[i]; } else {right=boundary.right;} // right(j==block_size-1)
                            } else {
                                right = partial_result[i*block.x +j +1];
                            }
                            if(i==block.y-1){
                                if(grid.y>1){down=rec_buffer.vertical_2[j]; } else {down=boundary.down;} // down (i==block_size-1)
                            } else {
                                down = partial_result[i*block.x +j + block.x];
                            }
                            
                            
                        }
                    }
                }
                //--------------------------------------------------------------------------------------------------------------
                partial_result_new[i*block.x +j] = (left + right+ up + down + b)/4;
            
            }
        }
        
        for(int u=0;u<block.y;u++){for(int v=0;v<block.x;v++){
            partial_result[u*block.x +v] = partial_result_new[u*block.x +v];
        }}
        
        
        
    
    // ---------------Communicating the partial results (Blocking)-------------------------------------------
    // if(rank!=0){
    //     MPI_Send(partial_result_new, block.x*block.y, MPI_FLOAT, 0, 99, MPI_COMM_WORLD);
    // } else {
    //     for(int u=0;u<block.x*block.y;u++){
    //         phi[u]=partial_result_new[u];
    //     }
    //     for(int k=1; k<np; k++){
            
    //         float* buf = (float*) calloc(block.x*block.y,sizeof(float));
    //         MPI_Recv(buf,block.x*block.y, MPI_FLOAT,k,99,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    //         for(int u=0;u<block.x*block.y;u++){
    //             phi[k*block.x*block.y+u]=buf[u];
    //         }
    //     }
    // }
    // Real position logic to compare with cpu code
    if(rank != 0){
        MPI_Send(partial_result_new, block.x * block.y, MPI_FLOAT, 0, 99, MPI_COMM_WORLD);
    } 
    else {

        // there are also gaps b/w
        for(int i = 0; i < block.y; i++){
            for(int j = 0; j < block.x; j++){
                phi[i*N + j] = partial_result_new[i*block.x + j];
            }
        }

        for(int k = 1; k < np; k++){

            float* buf = (float*) calloc(block.x * block.y, sizeof(float));
            MPI_Recv(buf, block.x * block.y, MPI_FLOAT, k, 99, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

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

            free(buf);
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
}
int main(int argc, char **argv) {
    
    // printf("hey");
    // fflush(stdout);
    float phi[N*N];

    Physics physics = {1, 1};
    Simulation simulation = {1,1};
    Boundary boundary = {1.0f,-1.0f,0,0 };

    jacobi_2d_parallel(phi, physics, simulation, boundary ,argc, argv);
    // Test for more examples
    return 0;
    
}