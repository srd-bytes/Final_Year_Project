#include <stdio.h>
#include <stdlib.h>
#include <cuda_runtime.h>



int main(void)
{   
    int A[2][2]={{1,2},{3,4}};
    printf("%d\n",A[0][0]);
    return 0;
}