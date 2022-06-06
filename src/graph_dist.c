#include "graph.h"

/* Dijkstra's algorith for shortest path lengths */


static double edgelength(struct graph *g, int i, int j)
{
    int k;
    double dv;
    double d = 0;

    for (k = 0; k < graph_dim(g); ++k)
    {
        dv = mat_elem(g->V, i, k) - mat_elem(g->V, j, k);
        d += dv * dv;
    }

    return sqrt(d);
}


static void graph_dist_(
    struct graph *g, struct matrix *D, int i, int *notvisited)
{
    HERE("graph_dist_");
    
    int k;
    int nQ;
    int n = graph_nnode(g);

    double alt;
    
    for (k = 0; k < n; ++k)
    {
        notvisited[k] = 1;
        mat_elem(D, i, k) = R_PosInf;
    }
    mat_elem(D, i, i) = 0;

    csi u;
    csi v;
    csi p;
    csi *Gp = g->G->p;
    csi *Gi = g->G->i;
    
    u = (csi)i;
    nQ = n - 1;

    do 
    {
        notvisited[u] = 0;
        for (p = Gp[u]; p < Gp[u+1]; ++p)
        {
            v = Gi[p];
            if (notvisited[v] && mat_elem(D, i, u) < R_PosInf)
            {
                alt = mat_elem(D, i, u) + edgelength(g, (int)u, (int)v);
                if (alt < mat_elem(D, i, v))
                    mat_elem(D, i, v) = alt;
            }
        }
        alt = R_PosInf;
        for (k = 0; k < n; ++k)
        {
            if (notvisited[k] && mat_elem(D, i, k) <= alt)
            {
                u = k;
                alt = mat_elem(D, i, k);
            }
        }
    } while (nQ-- > 0);
}



struct matrix *graph_dist(struct graph *g)
{
    HERE("graph_dist");

    int i;
    int n = graph_nnode(g);
    struct matrix *D = matrix_create(n, n, n, 0);
    int *notvisited = malloc(n * sizeof(*notvisited));
    for (i = 0; i < n; ++i)
        graph_dist_(g, D, i, notvisited);
    free(notvisited);
    return D;
}

