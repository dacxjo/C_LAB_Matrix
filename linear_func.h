//
// Created by David Blandon on 7/10/21.
//

#ifndef LAB_TRIANGLE_LINEAR_FUNC_H
#define LAB_TRIANGLE_LINEAR_FUNC_H

/**
 * Resolve an upper triangular matrix
 * @param n - Matrix dimension
 * @param A - Matrix
 * @param b - Vector
 * @param x - Solutions Vector
 * @param tol - Tolerance
 * @return - 0 if success, 1 if not
 */
int resoltrisup(int n, double **A, double *b, double *x, double tol);

/**
 * Resolve a lower triangular matrix
 * @param n - Matrix dimension
 * @param A - Matrix
 * @param b - Vector
 * @param x - Solutions Vector
 * @param tol - Tolerance
 * @return - 0 if success, 1 if not
 */
int resoltriinf(int n, double **A, double *b, double *x, double tol);

/**
 * Check if a given matrix is upper triangular
 * @param A - Matrix to check
 * @param n - Matrix dimension
 * @return - 0 if matrix is not upper triangular and 1 instead
 */
_Bool checkupptri(double **A, int n);


/**
 * Check if a given matrix is lower triangular
 * @param A - Matrix to check
 * @param n - Matrix dimension
 * @return - 0 if matrix is not lower triangular and 1 instead
 */
_Bool checklowertri(double **A, int n);

/**
 * Calc the scalar product between two given vectors
 * @param n - Vectors dimension
 * @param x - Vector X
 * @param y - Vector Y
 * @return - Scalar product
 */
double prod_esc(int n, double *x, double *y);

double *prodMatVect(double **M, double *x,int n);

#endif //LAB_TRIANGLE_LINEAR_FUNC_H
