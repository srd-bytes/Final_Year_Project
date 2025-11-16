#include <stdio.h>
#include <stdlib.h>
#include<omp.h>

int main(void){
    int sum = 0;
    int N=12;
    int *A = (int*) calloc(N,sizeof(int));
    for (int i = 0; i < N; i++){
        A[i]=1;
    }

    #pragma omp parallel for reduction(+:sum) schedule(static,3)
    for (int i = 0; i < N; i++){
        int tid = omp_get_thread_num();
        sum += A[i];
        printf("Thread %d handling i=%d\n", tid, i);
    }

    free(A);
    return 0;
}