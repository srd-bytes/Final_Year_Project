#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>
#define N 10

int main() {
    clock_t start, end;
    double cpu_time_used;

    start = clock();
    int *arr= calloc(N,sizeof(int));

    #pragma omp parallel for
    for (int i = 0; i < N; i++) {
        // Loop body (empty)
        arr[i]=i*2*3;

    

    
    // printf("Time taken for iterations: seconds\n");
    
    }
    free(arr);
    end = clock();

    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("Time taken for iterations: %f seconds\n", cpu_time_used);
    
    return 0;
}