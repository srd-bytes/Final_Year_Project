#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <time.h>

int main(){
    int s=0;
    int s_p=0;


    #pragma omp parallel firstprivate(s_p) reduction(+:s)
    {   
        
        #pragma omp for
        for(int i=0;i<100;i++){
            s_p+=1;
        }
        s=s_p;
        #pragma omp barrier
        printf("sum = %d\n",s);
    }
    printf("sum = %d\n",s);
    
    return 0;
}