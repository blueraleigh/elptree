#include "graph.h"

struct graph *graph_init(
    SEXP X, SEXP G, SEXP V, 
    SEXP lambda, SEXP mu, 
    SEXP alpha, SEXP R02, SEXP tol)
{
    HERE("graph_init");

    struct graph *g = malloc(sizeof(*g));

    if (!g)
        return NULL;

    g->lambda = REAL(lambda)[0];
    g->mu = REAL(mu)[0];
    g->alpha = REAL(alpha)[0];
    g->R02 = REAL(R02)[0];
    g->tol = REAL(tol)[0];

    int k;
    int nedge = INTEGER(getAttrib(G, R_DimSymbol))[0];
    int n = INTEGER(getAttrib(X, R_DimSymbol))[0];
    int p = INTEGER(getAttrib(X, R_DimSymbol))[1];

    g->X = matrix_create(n, p, n, REAL(X));

    n = INTEGER(getAttrib(V, R_DimSymbol))[0];
    p = INTEGER(getAttrib(V, R_DimSymbol))[1];

    g->V = matrix_create(n, p, n, REAL(V));
    g->K = matrix_copy(g->V);
    g->W = malloc(n * sizeof(double));
    g->L = NULL;

    cs *TG = cs_spalloc(n, n, (csi)(nedge), 1, 1);

    for (k = 0; k < nedge; ++k)
    {
        cs_entry(
            TG
            , (csi)(INTEGER(G)[k + nedge] - 1)
            , (csi)(INTEGER(G)[k] - 1)
            , 1
        );
    }

    g->G = cs_compress(TG);

    cs_spfree(TG);

    return g;
}
