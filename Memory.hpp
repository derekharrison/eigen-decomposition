/*
 * Memory.h
 *
 *  Created on: 15 apr. 2016
 *      Author: dharrison
 */

#ifndef MEMORY_H_
#define MEMORY_H_


double *matrix1D(int np);
double **matrix2D(int nm, int np);
void free2D(double** f, int nm);


#endif /* MEMORY_H_ */
