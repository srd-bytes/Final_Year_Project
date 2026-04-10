
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

void  write_to_csv_1d(float* data, int N,const char *path){
    FILE* fp = fopen(path, "w");

    for (int i = 0; i < N; i++) {
        
        fprintf(fp, "%f", data[i]);
        fprintf(fp, "\n");
    }

    fclose(fp);
}
void  write_to_csv_2d(float* data, int N,const char *path){
    FILE* fp = fopen(path, "w");

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            fprintf(fp, "%f", data[N*i+j]);

            if (j < N - 1)
                fprintf(fp, ",");
        }
        fprintf(fp, "\n");
    }

    fclose(fp);
}

void write_to_csv_3d(float* data, int N, const char* path) {
    FILE* fp = fopen(path, "w");

    for (int z = 0; z < N; z++) {

        // Optional: mark slice (CSV comment, safe to ignore in parsing)
        fprintf(fp, "# z = %d\n", z);

        for (int y = 0; y < N; y++) {
            for (int x = 0; x < N; x++) {
                int idx = z*N*N + y*N + x;
                fprintf(fp, "%f", data[idx]);

                if (x < N - 1)
                    fprintf(fp, ",");
            }
            fprintf(fp, "\n");
        }

        
    }

    fclose(fp);
}