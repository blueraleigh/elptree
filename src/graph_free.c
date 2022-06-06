#include "graph.h"

void graph_free(struct graph *g)
{
    if (!g)
        return;
    free(g->W);
    matrix_free(g->X);
    matrix_free(g->K);
    matrix_free(g->V);
    cs_spfree(g->G);
    cs_spfree(g->L);
    free(g);
}