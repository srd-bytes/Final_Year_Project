double *generate_poisson_matrix(char path[]){
    FILE *fp = fopen(path, "r");
    if (!fp) {
        perror("File opening failed");
        return 1;
    }
}