/*
 * MatrixOps.c
 *
 *  Created on: 15 apr. 2016
 *      Author: dharrison
 */

#include <math.h>
#include <stdlib.h>
#include "Memory.hpp"


void MatMult(double** A, int rows1, int columns1, double** B, int rows2, int columns2, double** result) {
    /*
     * This function performs matrix multiplication of two matrices A and B
     */

    for(int i = 0; i < rows1; ++i)
        for(int j = 0; j < columns2; ++j) {
            result[i][j] = 0.0;
        }

    for(int i = 0; i < rows1; ++i)
        for(int j = 0; j < columns2; ++j) {
            for(int k = 0; k < columns1; ++k)
                result[i][j] += A[i][k]*B[k][j];
        }
}


void MatTranspose(double** A, int rows, int columns, double** result) {
    /*
     * This function computes the transpose of matrix A
     */

    for(int i = 0; i < columns; ++i)
        for(int j = 0; j < rows; ++j)
            result[i][j] = A[j][i];

}


void VectorNormalization(double *wp, int sizeVec) {

    double sqrtwpwp = 0.0;

    for(int i = 0; i < sizeVec; ++i)
        sqrtwpwp += wp[i]*wp[i];

    for(int i = 0; i < sizeVec; ++i)
        wp[i] = wp[i]/sqrt(sqrtwpwp);
}


void EigenDecomposition(double** A, int N, double** EigVectors, double* EigValues, int iterations) {
    /*
     * This function computes the eigenvalues and eigenvectors
     * of a real, symmetric, N x N matrix A, using the QR algorithm.
     */

    int i, j, k, p, it;
    double **EigVecs, **Q, **EigVals, **result, *wp, *dumsum;
    double **R, **QT;
    double f;

    EigVecs = matrix2D(N, N);
    result = matrix2D(N, N);
    Q = matrix2D(N, N);
    EigVals = matrix2D(N, N);
    wp = matrix1D(N);
    dumsum = matrix1D(N);
    R = matrix2D(N, N);
    QT = matrix2D(N, N);

    /*Initializing Ait, matrix containing eigenvalues as the diagonal*/
    for(i = 0; i < N; ++i)
        for(j = 0; j < N; ++j)
            EigVals[i][j] = A[i][j];

    /*Initializing Q and E for computation of eigenvectors*/
    for(i = 0; i < N; ++i) {
        Q[i][i] = 1.0;
        EigVecs[i][i] = 1.0;
    }

    /*Eigen decomposition iterations*/
    for(it = 0; it < iterations; ++it) {
        /*Gram-Schmidt decorrelation*/
        for(p = 0; p < N; ++p) {
            for(i = 0; i < N; ++i)
                wp[i] = EigVals[i][p];

            VectorNormalization(wp, N);

            for(i = 0; i < N; ++i)
                dumsum[i] = 0.0;

            for(i = 0; i < N; ++i) {
                for(j = 0; j < p ; ++j) {
                    f = 0.0;
                    for(k = 0; k < N; ++k)
                        f += wp[k]*Q[k][j];

                    dumsum[i] += f*Q[i][j];
                }
            }

            for(i = 0;  i < N; ++i)
                wp[i] -= dumsum[i];

            VectorNormalization(wp, N);

            /*Storing estimated rows of the inverse of the mixing matrix as columns in W*/
            for(i = 0; i < N; ++i)
                Q[i][p] = wp[i];
        }

        MatTranspose(Q, N, N, QT);

        MatMult(QT, N, N, EigVals, N, N, R);

        MatMult(R, N, N, Q, N, N, EigVals);

        MatMult(EigVecs, N, N, Q, N, N, result);

        for(int i = 0; i < N; ++i) {
            for(int j = 0; j < N; ++j) {
                EigVecs[i][j] = result[i][j];
            }
        }
    }

    /* Set results */
    for(i = 0; i < N; ++i)
        EigValues[i] = EigVals[i][i];

    for(i = 0; i < N; ++i)
        for(j = 0; j < N; ++j)
            EigVectors[i][j] = EigVecs[i][j];

    free2D(Q, N);
    free2D(EigVals, N);
    free2D(R, N);
    free2D(QT, N);
    free2D(EigVecs, N);
    free2D(result, N);
    delete [] wp;
    delete [] dumsum;
}
