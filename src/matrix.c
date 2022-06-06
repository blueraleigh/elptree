#include "matrix.h"

#ifdef VDEBUG
#define HERE(msg) Rprintf("%s (%s:%u)\n", msg, __FILE__, __LINE__)
#else
#define HERE(msg) ((void)0)
#endif

struct matrix *matrix_create(int n, int p, int LDM, double *x)
{
    if (n < 0 || p <= 0 || LDM <= 0 || LDM < n)
        return NULL;
    
    struct matrix *m = malloc(sizeof(*m));
    
    if (!m)
        return NULL;
    
    m->x = calloc(LDM * p, sizeof(double));
    
    if (!(m->x))
    {
        free(m);
        return NULL;
    }
    
    m->n = n;
    m->p = p;
    m->LDM = LDM;
    
    if (x)
        memcpy(m->x, x, LDM * p * sizeof(double));
    
    return m;
}


struct matrix *matrix_copy(struct matrix *A)
{
    return matrix_create(A->n, A->p, A->LDM, A->x);
}


void matrix_free(struct matrix *m)
{
    if (!m)
        return;
    free(m->x);
    free(m);
}


struct matrix *matrix_multiply(
    struct matrix *A, char *transA, struct matrix *B, char *transB)
{
    HERE("matrix_multiply");

    if (A->p != B->n)
        return NULL;
    
    struct matrix *C = matrix_create(A->n, B->p, A->n, NULL);
    
    if (!C)
        return NULL;
    
    int m = A->n;
    int n = B->p;
    int k = A->p;
    double alpha = 1;
    double *a = mat_values(A);
    int lda = A->LDM;
    double *b = mat_values(B);
    int ldb = B->LDM;
    double beta = 0;
    double *c = mat_values(C);
    int ldc = C->LDM;
    
    F77_CALL(dgemm)(
        transA,
        transB,
        &m,
        &n,
        &k,
        &alpha,
        a,
        &lda,
        b,
        &ldb,
        &beta,
        c,
        &ldc
    );
    
    return C;
}


// creates a copy of matrix A with row i deleted
struct matrix *matrix_delete_row(struct matrix *A, int i)
{
    HERE("matrix_delete_row");

    int ii;
    int last_row = A->n - 1;

    struct matrix *P;
    struct matrix *C;

    if (i > last_row || i < 0)
        i = last_row;

    if (i < last_row)
    {
        // create a permutation matrix
        P = matrix_create(A->n, A->n, A->n, 0);

        // rows before i stay in place
        for (ii = 0; ii < i; ++ii)
            mat_elem(P, ii, ii) = 1;

        // rows after i shift down one position
        for (ii = i+1; ii <= last_row; ++ii)
            mat_elem(P, ii-1, ii) = 1;

        // row i is moved to the last row
        mat_elem(P, last_row, i) = 1;

        C = matrix_multiply(P, "N", A, "N");

        matrix_free(P);
    }
    else
    {
        // row i is already the last row
        C = matrix_copy(A);
    }

    // decrement the row count so subsequent computations
    // ignore the last row
    C->n -=1;

    return C;
}


// Creates a copy of matrix A with the new row x inserted 
// before row i. If row i is negative, the new row is appended.
struct matrix *matrix_insert_row(struct matrix *A, int i, double *x)
{
    HERE("matrix_insert_row");

    int ii;
    int last_row = A->n - 1;

    struct matrix *P;
    struct matrix *B;
    struct matrix *C;

    if (i < 0 || i > last_row)
        i = A->n;
    
    // create an identity permutation matrix but with an extra row
    P = matrix_create(A->n+1, A->n, A->n+1, 0);

    for (ii = 0; ii <= last_row; ++ii)
        mat_elem(P, ii, ii) = 1;

    // B now has row of zeros at the end
    B = matrix_multiply(P, "N", A, "N");

    matrix_free(P);

    if (i <= last_row)
    {
        // a new permutation matrix to perform the insertion
        P = matrix_create(B->n, B->n, B->n, 0);

        // rows before i stay in place
        for (ii = 0; ii < i; ++ii)
            mat_elem(P, ii, ii) = 1;

        // rows i and up shift up one position
        for (ii = i; ii <= last_row; ++ii)
            mat_elem(P, ii+1, ii) = 1;

        // the new row is moved to position i
        mat_elem(P, i, last_row+1) = 1;     

        C = matrix_multiply(P, "N", B, "N");
        matrix_free(P);
        matrix_free(B);
    }
    else // we're all done
    {
        C = B;
    }

    if (x)
    {
        for (ii = 0; ii < C->p; ++ii)
            mat_elem(C, i, ii) = x[ii];
    }

    return C;
}


void matrix_print(struct matrix *A, int brief)
{
    for (int i = 0; i < A->n; ++i)
    {
        for (int j = 0; j < A->p; ++j)
        {
            Rprintf("%-10f", mat_elem(A, i, j));
        }
        Rprintf("\n");
    }
}
