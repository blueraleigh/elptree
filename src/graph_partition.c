#include "graph.h"


static int nearest_neighbor(
    struct graph *g,
    int i,                   // query point index
    int *outlier,            // flag to record if query point is outlier 
    double *err              // nearest neighbor approximation error
)
{
    int j;
    int p;
    int nn;
    
    double d;
    double d2;
    double min_d2 = R_PosInf;

    struct matrix *X = graph_data(g);
    struct matrix *V = graph_coordinates(g);

    for (j = 0; j < graph_nnode(g); ++j)
    {
        d2 = 0;
        for (p = 0; p < graph_dim(g); ++p)
        {
            d = mat_elem(X, i, p) - mat_elem(V, j, p);
            d2 += d * d;
        }

        if (d2 < min_d2)
        {
            nn = j;
            min_d2 = d2;
        }
    }

    *err += min_d2 < g->R02 ? min_d2 : g->R02;

    *outlier = min_d2 < g->R02 ? 0 : 1;

    return nn;

}


/* partition the data by graph vertices and compute the centroid and
** weight of each partition */
void graph_partition(struct graph *g)
{
    HERE("graph_partition");
    
    int i;
    int j;
    int nn;
    int outlier;
    double err = 0;

    double *W = graph_partition_weights(g);
    struct matrix *K = graph_partition_centers(g);
    struct matrix *X = graph_data(g);

    mat_clear(K);
    memset(W, 0, graph_nnode(g) * sizeof(double));

    for (i = 0; i < graph_ndata(g); ++i)
    {
        nn = nearest_neighbor(g, i, &outlier, &err);
        
        W[nn] += 1;

        if (W[nn] == 1)
        {
            for (j = 0; j < graph_dim(g); ++j)
                mat_elem(K, nn, j) = outlier ? 0 : mat_elem(X, i, j);
        }
        else
        {
            for (j = 0; j < graph_dim(g); ++j)
            {
                mat_elem(K, nn, j) += outlier ? 0 :
                    (mat_elem(X, i, j) - mat_elem(K, nn, j)) / W[nn];
            }
        }
    }

    for (i = 0; i < graph_nnode(g); ++i)
        W[i] /= graph_ndata(g);

    g->error = err / graph_ndata(g);

}

