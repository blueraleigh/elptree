#include "graph.h"


static void graph_R_repr(struct graph *g, int *A)
{
    cs *G = g->G;
    csi *Gp = G->p;
    csi *Gi = G->i;

    csi i;
    csi j;
    csi p;

    int n = 0;

    for (j = 0; j < graph_nnode(g); ++j)
    {
        for (p = Gp[j]; p < Gp[j+1]; ++p)
        {
            i = Gi[p];

            A[n + 0*graph_nedge(g)] = (int)(j+1);
            A[n + 1*graph_nedge(g)] = (int)(i+1); 
            ++n;
        }
    }
}


SEXP C_graph_fit(
    SEXP X,         // data for embedding
    SEXP G,         // initial graph topology
    SEXP V,         // initial graph embedding
    SEXP size,      // number of graph vertices
    SEXP lambda,
    SEXP mu,
    SEXP alpha,
    SEXP R02,
    SEXP tol
)
{
    HERE("C_graph_fit");
    
    struct graph *g = graph_init(X, G, V, lambda, mu, alpha, R02, tol);

    if (!g)
        error("could not allocate memory for embedding");

    graph_partition(g);
    graph_laplacian(g);
    graph_optimize(g);

    int i;
    int j;
    int maxnodes = INTEGER(size)[0];
    double delta = R_PosInf;
    double new_score;
    double old_score;
    
    Rprintf("%-14s %-14s %-14s %-14s\n", 
        "size", "energy", "error", "energy + error");
    Rprintf("%-14d %-14f %-14f %-14f\n", 
        graph_nnode(g), g->energy, g->error, g->energy + g->error);

    while (graph_nnode(g) < maxnodes)
    {
        old_score = g->energy + g->error;
        graph_grow(g);
        graph_grow(g);
        graph_shrink(g);
        new_score = g->energy + g->error;
        delta = old_score - new_score;
        Rprintf("%-14d %-14f %-14f %-14f\n", 
            graph_nnode(g), g->energy, g->error, g->energy + g->error);
    };

    // prepare SEXP result 
    SEXP ans = PROTECT(allocVector(VECSXP, 5));

    SET_VECTOR_ELT(ans, 0, allocMatrix(INTSXP, graph_nedge(g), 2));
    SET_VECTOR_ELT(ans, 1, allocMatrix(REALSXP, graph_nnode(g), graph_dim(g)));
    SET_VECTOR_ELT(ans, 2, allocVector(REALSXP, graph_nnode(g)));
    SET_VECTOR_ELT(ans, 3, ScalarReal(g->energy));
    SET_VECTOR_ELT(ans, 4, ScalarReal(g->error));

    for (i = 0; i < graph_nnode(g); ++i)
    {
        for (j = 0; j < graph_dim(g); ++j)
        {
            REAL(VECTOR_ELT(ans, 1))[i+j*graph_nnode(g)] = mat_elem(g->V, i, j);
        }
    }

    for (i = 0; i < graph_nnode(g); ++i)
        REAL(VECTOR_ELT(ans, 2))[i] = g->W[i];

    graph_R_repr(g, INTEGER(VECTOR_ELT(ans, 0)));

    graph_free(g);

    UNPROTECT(1);

    return ans;
}
