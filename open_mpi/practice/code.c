#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define N 512 // size of array

#include <time.h>



int dot_prod(int v_a[N], int v_b[N]){

    int result=0;
    
    for (int i=0; i< N; i++){

        result += v_a[i]*v_b[i];
    }
   
    return result;
}

void matrix_mult(int A[N][N], int B[N][N], int result[N][N]){

    for(int i=0; i< N; i++){
        for(int j=0; j<N; j++){
            int* col_vector = (int*) calloc(N, sizeof(int));
            for(int k=0; k<N; k++){
                col_vector[k] = B[k][j];
            }
            result[i][j]= dot_prod(A[i], col_vector);
            free(col_vector);
        }
    }
}

int main() {
    

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
    
    clock_t start, end;
    double cpu_time_used;

    start = clock();

    matrix_mult(a,b, r);

    end = clock();
    cpu_time_used= ((double) (end - start))/CLOCKS_PER_SEC;
    printf("Execution time = %f seconds\n", cpu_time_used);
    //---------------Testing (printing result) --------------------

    // for(int row=0; row < N; row++){
    //     for (int col=0; col< N; col++){

    //         printf("result[%d][%d] = %d\n", row,col,r[row][col]);
    //     }
    // }
    // -----------------------------------------

    // TODO: Measure the performance
    
}