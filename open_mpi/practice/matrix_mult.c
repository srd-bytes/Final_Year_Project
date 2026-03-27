#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mpi.h"

#define N 512 // size of array

// Try with countinuous array

int dot_prod(int v_a[N], int v_b[N]){ 

    int result=0;
    
    for (int i=0; i< N; i++){

        result += v_a[i]*v_b[i];
    }
   
    return result;
}
void matrix_mult_parallel(int A[N][N], int B[N][N], int result[N][N],int argc, char **argv){
   
    // 4 processors
    MPI_Init(&argc, &argv);

    int np;
    MPI_Comm_size(MPI_COMM_WORLD, &np);

    int p_id;
    MPI_Comm_rank(MPI_COMM_WORLD, &p_id);

    int effective_size= N/(int)sqrt(np);

    int block_row_id = p_id/(N/effective_size);
    int block_col_id = p_id - block_row_id*(N/effective_size);

    double start, end;

	MPI_Barrier(MPI_COMM_WORLD); // Important
	start = MPI_Wtime();

    // ---------------------dynamic allocation partial result------------------
    int** partial_result = (int**) calloc(effective_size, sizeof(int*));
    for(int i=0; i< effective_size; i++){
        partial_result[i] = (int*) calloc(effective_size, sizeof(int));
    }
    //-------------------------------------------------------------------------
    
    // -------------------Partial result calculation ----------------------------------------------------------------------
    for(int i=0; i< effective_size; i++){
        for(int j=0; j< effective_size; j++){

            int* col_vector = (int*) calloc(N,sizeof(int));
            for(int k=0; k< N; k++){
                col_vector[k] = B[k][block_col_id*effective_size+j];
            }

            partial_result[i][j] = dot_prod(A[block_row_id*effective_size+i],col_vector);

            free(col_vector);
        }

    }
    // --------------------------------------------------------------------------------------------------------------------
    
    // -------------------------------Blocking Communication---------------------------------------------------------------
    if(p_id !=0){
        for(int i=0; i< effective_size; i++){
            MPI_Send(partial_result[i],effective_size, MPI_INT, 0, 99, MPI_COMM_WORLD);
        }
    } else {

        for(int i=0; i< effective_size; i++){
            for(int j=0; j< effective_size; j++){
                result[i][j] = partial_result[i][j];
            }
        }

        for(int process=1; process < np; process++){

            int row_id = process/(N/effective_size);
            int col_id = process - row_id*(N/effective_size);
            // Receive similar to send else there will be deadlock
            for(int i=0; i< effective_size; i++){
                int* buffer = (int*) calloc(effective_size,sizeof(int));
                MPI_Recv(buffer,effective_size, MPI_INT, process, 99, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                // processing the data in the buffer
                for(int col=0; col< effective_size; col++){
                    result[row_id*effective_size + i][col_id*effective_size + col] = buffer[col];
                }
                free(buffer);
            }
        }
        // -----------------------------------------------------------------------------------------------------------

	    end = MPI_Wtime();
	
	    double elapsed = end - start;
	
        if(p_id == 0)
        {
            printf("Execution time = %f seconds\n", elapsed);
        }
        //---------------Testing (printing result) --------------------

        // for(int row=0; row < N; row++){
        //     for (int col=0; col< N; col++){

        //         printf("result[%d][%d] = %d\n", row,col,result[row][col]);
        //     }
        // }
        // -----------------------------------------
    }
    // --------------------------------------------------------------------------------------------------------------------

    for(int i=0; i<effective_size;i++){
        free(partial_result[i]);
    }
    free(partial_result);

    
    MPI_Finalize();
}



int main(int argc, char **argv) {
    

    //Test for parallel matrix multiplication

    int a[N][N];
    int b[N][N];
    for(int i=0; i<N; i++){
        for(int j=0; j<N; j++){
            a[i][j]= 2;
            b[i][j]=2;
        }
    }
    int r[N][N];
    
    matrix_mult_parallel(a,b, r, argc, argv);

    // TODO: Measure the performance

    
    
}
