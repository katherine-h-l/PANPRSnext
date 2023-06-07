#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

// g++ -fopenmp test.cpp

int main() {
    int P = 4;
    int Q = 5;
    double **jointBmatrix = (double **)malloc(P * sizeof(double *));
    jointBmatrix[0] = (double *)calloc(P * Q, sizeof(double));
    for (int j = 0; j < P; j++)
    {
        jointBmatrix[j] = jointBmatrix[0];
    }

    printf("MALLOCED\n");

    #pragma omp parallel
    {
        int thread_id = omp_get_thread_num();
        int num_threads = omp_get_num_threads();
        printf("Thread %d of %d\n", thread_id, num_threads);
    }

    // #pragma omp parallel for
    for (int i = 0; i < P; i++)
    {
        for (int j = 0; j < Q; j++)
        {
            printf("%f ", jointBmatrix[i][j]);
        }
        printf("\n");
    }

    free(jointBmatrix[0]);
    free(jointBmatrix);

    return 0;
}