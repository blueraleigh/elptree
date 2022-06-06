#include "graph.h"

/* Modification of cs_cholsol to solve A * X = B, 
** where X and B are dense matrices and 
** A is a square, symmetric sparse matrix. */
static csi cholsol(const cs *A, struct matrix *X, struct matrix *B)
{
    HERE("cholsol");
    
    double *x;
    double *b;
    css *S;
    csn *N;
    csi order = 1; // Cholesky ordering (see: cs_amd)
    csi ok;
    if (!CS_CSC (A) || !B || !X) return (0);
    csi n = A->n;
    S = cs_schol(order, A);
    N = cs_chol(A, S);
    x = cs_malloc(n, sizeof(double));
    ok = (S && N && x);
    if (ok)
    {
        int i;
        int j;
        int k = X->p;
        int ldb = B->LDM;
        for (i = 0, j = 0; i < k; ++i, j += ldb)
        {
            b = mat_values(B) + j; // increment to the appropriate column
            cs_ipvec(S->pinv, b, x, n);
            cs_lsolve(N->L, x);
            cs_ltsolve(N->L, x);
            cs_pvec(S->pinv, x, mat_values(X) + j, n);
        }
    }
    cs_free(x);
    cs_sfree(S);
    cs_nfree(N);
    return ok;
}


/* Modification of cs_lusol to solve A * X = B, 
** where X and B are dense matrices and 
** A is a square, symmetric sparse matrix. */
static csi lusol(const cs *A, struct matrix *X, struct matrix *B)
{
    HERE("lusol");

    /* tol=1 for  partial pivoting; tol<1 gives preference to diagonal */
    double tol = 1;
    double *x;
    double *b;
    css *S;
    csn *N;
    csi order = 2; // LU ordering (see: cs_amd)
    csi ok;
    if (!CS_CSC (A) || !B || !X) return (0);
    csi n = A->n;
    S = cs_sqr(order, A, 0);
    N = cs_lu(A, S, tol);
    x = cs_malloc(n, sizeof(double));
    ok = (S && N && x);
    if (ok)
    {
        int i;
        int j;
        int k = X->p;
        int ldb = B->LDM;
        for (i = 0, j = 0; i < k; ++i, j += ldb)
        {
            b = mat_values(B) + j; // increment to the appropriate column
            cs_ipvec(N->pinv, b, x, n);
            cs_lsolve(N->L, x);
            cs_usolve(N->U, x);
            cs_ipvec(S->q, x, mat_values(X) + j, n);
        }
    }
    cs_free(x);
    cs_sfree(S);
    cs_nfree(N);
    return ok;
}

// assumes that graph_partition and graph_laplacian
// have already been called
static void update_coordinates(struct graph *g)
{
    // solve g->L * g->V = g->K

    HERE("update_coordinates");

    int i;
    int j;

    /* Scale the partition centroids by the partition weights.
    ** This means that K now holds the rhs of equation (5) in
    **   https://arxiv.org/pdf/0801.0176.pdf 
    ** which we will use to update the graph coordinates
    */
    for (i = 0; i < graph_nnode(g); ++i)
    {
        for (j = 0; j < graph_dim(g); ++j)
            mat_elem(g->K, i, j) *= g->W[i];
    }

    // cs_chol requires pos def L
    int ok = cholsol(g->L, g->V, g->K);
    if (!ok)
    {
        // L likely pos semi-definite
        ok = lusol(g->L, g->V, g->K);
    }
    if (!ok)
    {
        cs_print(g->L, 0);
        error("cannot update embedding coordinates");
    }
    
    graph_partition(g);
    graph_laplacian(g);
}


void graph_optimize(struct graph *g)
{
    HERE("graph_optimize");

    double old_score;
    double new_score;
    int i = 0;
    do {
        old_score = g->energy + g->error;
        update_coordinates(g);
        new_score = g->energy + g->error;
        assert(new_score <= old_score || (new_score - old_score) < 1e-8);
    } while ((old_score - new_score) > g->tol);
}


SEXP C_graph_update(SEXP lambda, SEXP mu, SEXP G, SEXP V, SEXP X)
{
    HERE("C_graph_laplacian");

    struct graph *g = graph_init(
        X, G, V, lambda, mu, ScalarReal(0), 
        ScalarReal(R_PosInf), ScalarReal(0.0001));
    
    if (!g)
        error("could not allocate memory for embedding");
    
    graph_partition(g);
    graph_laplacian(g);
    update_coordinates(g);

    int j;
    int i;
    int p;
    int n = graph_nnode(g);

    SEXP ans = PROTECT(allocMatrix(REALSXP, n, graph_dim(g)));

    memcpy(REAL(ans), mat_values(g->V), n * graph_dim(g) * sizeof(double));

    graph_free(g);

    UNPROTECT(1);

    return ans;
}
