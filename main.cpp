/*
 * main.cpp
 *
 *  Created on: Jan 19, 2021
 *      Author: d-w-h
 *
 *      Implementation of the QR algorithm for computing eigen
 *      vectors and eigen values of a given real symmetric matrix
 */

#include <stdio.h>
#include "MatrixOps.hpp"
#include "Memory.hpp"

int main(int argc, char* argv[]) {

    double **A, **EigVectors, ** Vec, *EigValues;
    int N, iterations;

    /* Parameters */
    N = 10;
    iterations = 200;

    /* Allocate memory for data */
    A = matrix2D(N, N);
    EigVectors = matrix2D(N, N);
    Vec = matrix2D(N, N);
    EigValues = matrix1D(N);

    /* Initialize A */
    for(int i = 0; i < N; ++i)
        for(int j = 0; j < N; ++j) {
            if(i == j) {
                A[i][j] = i;
            }
            else {
                A[i][j] = -0.25;
            }
        }

    EigenDecomposition(A, N, EigVectors, EigValues, iterations);

    /* Print results */
    for(int i = 0; i < N; ++i) {
        printf("EigValues[%i]: %f\n", i, EigValues[i]);
    }

    for(int i = 0; i < N; ++i) {
        for(int j = 0; j < N; ++j) {
            printf("EigVectors[%i][%i]: %f, ", i, j, EigVectors[i][j]);
        }
        printf("\n");
    }

    /* Verify computation of eigen vectors */
    MatMult(A, N, N, EigVectors, N, N, Vec);

    printf("Verify computation\n");
    for(int i = 0; i < N; ++i) {
        for(int j = 0; j < N; ++j) {
            printf("ratio: %f, ", Vec[i][j]/(EigVectors[i][j] + 1e-20));
        }
        printf("\n");
    }

    /* Free allocated data */
    free2D(A, N);
    free2D(EigVectors, N);
    free2D(Vec, N);
    delete [] EigValues;

    printf("done!\n");
    return 0;
}
