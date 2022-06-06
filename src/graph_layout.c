#include "graph.h"

/*
** Force-directed graph layout as described in:
**
** Graph Drawing by Stress Majorization
** Emden R. Gansner, Yehuda Koren and Stephen North
** In Proceedings 12th Symposium on Graph Drawing (GD), 
** pages 239–250, 2004
**
** https://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.439.26&rep=rep1&type=pdf
** 
*/


static double stress(struct matrix *X, struct matrix *D)
{
    int i;
    int j;
    int n = X->n;
    double d;
    double dx;
    double dy;
    double dij;
    double s = 0;
    for (i = 0; i < n; ++i)
    {
        for (j = i+1; j < n; ++j)
        {
            dij = mat_elem(D, i, j);
            dx = mat_elem(X, i, 0) - mat_elem(X, j, 0);
            dy = mat_elem(X, i, 1) - mat_elem(X, j, 1);
            d = sqrt(dx*dx + dy*dy);
            s += ( (d - dij) * (d - dij) ) / (dij * dij);
        }
    }
    return s;
}


static void chol(struct matrix *A)
{
    char *uplo = "L";
    int n = A->n;
    int lda = A->LDM;
    int info;

    F77_CALL(dpotrf)(
        uplo
        , &n
        , mat_values(A)
        , &lda
        , &info
    );

    if (info != 0)
        error("matrix is not positive definite");
}


static void update_layout(
    struct matrix *Lchol, struct matrix *x, struct matrix *Lx)
{
    struct matrix *Lz = matrix_multiply(Lx, "N", x, "N");

    struct matrix *z = matrix_delete_row(Lz, 0);

    char *uplo = "L";
    int n = Lchol->n;
    int nrhs = 2;
    int lda = Lchol->LDM;
    int ldb = z->LDM;
    int info;

    F77_CALL(dpotrs)(
        uplo
        , &n
        , &nrhs
        , mat_values(Lchol)
        , &lda
        , mat_values(z)
        , &ldb
        , &info
    );

    if (info != 0)
        error("could not update layout coordinates");

    for (int i = 1; i < x->n; ++i)
    {
        mat_elem(x, i, 0) = mat_elem(z, i-1, 0);
        mat_elem(x, i, 1) = mat_elem(z, i-1, 1);
    }

    matrix_free(Lz);
    matrix_free(z);
}


SEXP C_graph_layout(SEXP layout, SEXP G, SEXP V, SEXP X)
{
    HERE("C_graph_layout");

    struct graph *g = graph_init(
        X, G, V, ScalarReal(0), ScalarReal(0), ScalarReal(0),
        ScalarReal(0), ScalarReal(0));

    if (!g)
        error("could not allocate memory for embedding");

    int i;
    int j;
    int n = graph_nnode(g);

    double d;
    double dx;
    double dy;
    double dij;
    double s0;
    double s1;
    double delta;
    double eps = 0.0001;

    // calculate shortest path distances between all vertices
    struct matrix *D = graph_dist(g);

    // layout coordinates
    struct matrix *x = matrix_create(n, 2, n, REAL(layout));

    // laplacian matrices
    struct matrix *Lw = matrix_create(n, n, n, 0);
    struct matrix *Lx = matrix_create(n, n, n, 0);
    
    for (i = 0; i < n; ++i)
    {
        for (j = i+1; j < n; ++j)
        {
            dij = mat_elem(D, i, j);
            mat_elem(Lw, i, j) = mat_elem(Lw, j, i) = -1 / (dij * dij);
        }
    }

    for (j = 0; j < n; ++j)
    {
        for (i = 0; i < n; ++i)
        {
            if (i == j)
                continue;
            dij = mat_elem(D, i, j);
            mat_elem(Lw, j, j) += 1 / (dij * dij);
        }
    }

    // cholesky factorization
    struct matrix *Lchol = matrix_create(n-1, n-1, n-1, 0);
    for (j = 1; j < n; ++j)
    {
        for (i = 1; i < n; ++i)
            mat_elem(Lchol, i-1, j-1) = mat_elem(Lw, i, j);
    }

    chol(Lchol);

    s0 = stress(x, D);

    do {

        for (i = 0; i < n; ++i)
        {
            mat_elem(Lx, i, i) = 0;
            for (j = i+1; j < n; ++j)
            {
                dx = mat_elem(x, i, 0) - mat_elem(x, j, 0);
                dy = mat_elem(x, i, 1) - mat_elem(x, j, 1);
                dij = mat_elem(D, i, j);
                d = sqrt(dx*dx + dy*dy);
                if (d != 0)
                    mat_elem(Lx, i, j) = mat_elem(Lx, j, i) = -(1/d) * (1/dij);
                else
                    mat_elem(Lx, i, j) = mat_elem(Lx, j, i) = 0;
            }
        }
        for (j = 0; j < n; ++j)
        {
            for (i = 0; i < n; ++i)
            {
                if (i == j)
                    continue;
                mat_elem(Lx, j, j) -= mat_elem(Lx, i, j);
            }
        }
        update_layout(Lchol, x, Lx);
        s1 = stress(x, D);
        delta = (s0 - s1) / s0;
        s0 = s1;

    } while (delta > eps);

    SEXP ans = PROTECT(allocMatrix(REALSXP, n, 2));

    memcpy(REAL(ans), mat_values(x), n*2*sizeof(double));

    graph_free(g);
    matrix_free(x);
    matrix_free(Lw);
    matrix_free(Lx);
    matrix_free(Lchol);
    matrix_free(D);

    UNPROTECT(1);
    return ans;
}

