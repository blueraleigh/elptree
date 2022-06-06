#ifndef MATRIX_H
#define MATRIX_H

#include <R.h>
#include <R_ext/BLAS.h>

// We use dense matrices to represent the embedding coordinates
// and partition centroids. Each row represents a vertex in the
// graph; each column, a coordinate dimension. Because the graph
// structure is dynamic, we implement routines for inserting and
// removing rows that make use of permutation matrices and the
// level 3 BLAS routine dgemm. This is straightforward and hopefully
// also reasonably efficient.
struct matrix {

    // number of rows
    int n;

    // number of columns
    int p;

    // size of leading dimension
    int LDM;

    // values (malloc'd size is LDM * p * sizeof(double))
    double *x;
};

#define mat_elem(m, i, j) (m)->x[(i)+(j)*(m)->LDM]

#define mat_clear(m) memset((m)->x, 0, (m)->LDM*(m)->p*sizeof(double))

#define mat_values(m) (m)->x

/* constructor */
struct matrix *matrix_create(int n, int p, int LDM, double *x);

struct matrix *matrix_copy(struct matrix *A);

/* destructor */
void matrix_free(struct matrix *m);

/* return C = A*B (C must be free'd later) */
struct matrix *matrix_multiply(
    struct matrix *A, char *transA, struct matrix *B, char *transB);

// return a copy of matrix A with row i deleted
struct matrix *matrix_delete_row(struct matrix *A, int i);

// Creates a copy of matrix A with the new row x inserted
// before row i. If i is negative, the new row is appended.
struct matrix *matrix_insert_row(struct matrix *A, int i, double *x);

void matrix_print(struct matrix *A, int brief);

#endif /* MATRIX_H */
