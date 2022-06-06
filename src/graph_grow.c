#include "graph.h"


// Returns a copy of graph G with a new terminal vertex added to i
static cs *add_vertex_to_vertex(cs *G, csi i)
{
    HERE("add_vertex_to_vertex");

    csi v;
    csi p;
    csi k = G->n;

    csi *Gp = G->p;
    csi *Gi = G->i;

    cs *D;
    cs *DT = cs_spalloc(k+1, k+1, G->nzmax + 2, 1, 1);

    // the new node will receive index k

    for (v = 0; v < k; ++v)
    {
        for (p = Gp[v]; p < Gp[v+1]; ++p)
            cs_entry(DT, Gi[p], v, 1);
    }

    cs_entry(DT, i, k, 1);
    cs_entry(DT, k, i, 1);

    D = cs_compress(DT);

    cs_spfree(DT);

    return D;
}


// Remove the edge (i, j), create a new vertex k,
// and create the edges (i, k) and (k, j)
cs *bisect_edge(cs *G, csi i, csi j)
{
    HERE("bisect_edge");

    csi v;
    csi p;
    csi k = G->n;

    csi *Gp = G->p;
    csi *Gi = G->i;

    cs *D;
    cs *DT = cs_spalloc(k+1, k+1, G->nzmax + 2, 1, 1);

    // the new node will receive index k

    for (v = 0; v < k; ++v)
    {
        for (p = Gp[v]; p < Gp[v+1]; ++p)
        {
            if (v == i && Gi[p] == j)
                continue; // remove the edge (i,j)
            if (v == j && Gi[p] == i)
                continue; // remove the edge (j,i)
            cs_entry(DT, Gi[p], v, 1);
        }
    }

    cs_entry(DT, i, k, 1);
    cs_entry(DT, k, i, 1);
    cs_entry(DT, j, k, 1);
    cs_entry(DT, k, j, 1);

    D = cs_compress(DT);

    cs_spfree(DT);

    return D;
}


void graph_grow(struct graph *g)
{
    HERE("graph_grow");

    csi i;
    csi j;
    csi k;
    csi p;

    // originals
    int n = graph_nnode(g);
    cs *G = g->G;
    csi *Gp = G->p;
    csi *Gi = G->i;
    struct matrix *V = graph_coordinates(g);
    struct matrix *K = graph_partition_centers(g);
    double *W = graph_partition_weights(g);    

    // best new graph
    cs *Gmin = NULL;
    struct matrix *Vmin = NULL;
    struct matrix *Kmin = NULL;
    double *Wmin = NULL;
    double min_error;
    double min_energy;
    double min_score = R_PosInf;

    // centroid coordinates of k-star leaves
    double y[graph_dim(g)];

    int *visited = calloc(n, sizeof(int));

    for (j = 0; j < n; ++j)
    {
        memset(y, 0, graph_dim(g) * sizeof(double));
        for (p = Gp[j]; p < Gp[j+1]; ++p)
        {
            i = Gi[p];
            if (visited[i] == 0)
            {
                // first visit to edge

                g->G = bisect_edge(G, j, i);
                g->V = matrix_insert_row(V, -1, 0);
                g->K = matrix_insert_row(K, -1, 0);
                g->W = malloc((n+1) * sizeof(double));

                for (k = 0; k < graph_dim(g); ++k)
                {
                    mat_elem(g->V, n, k) = ( 
                        mat_elem(V, i, k) + mat_elem(V, j, k) ) / 2;
                }

                graph_partition(g);
                graph_laplacian(g);
                graph_optimize(g);

                if (min_score > (g->energy + g->error))
                {
                    cs_spfree(Gmin);
                    matrix_free(Vmin);
                    matrix_free(Kmin);
                    free(Wmin);
                    Gmin = g->G;
                    Vmin = g->V;
                    Kmin = g->K;
                    Wmin = g->W;
                    min_energy = g->energy;
                    min_error = g->error;
                    min_score = g->energy + g->error;
                }
                else
                {
                    cs_spfree(g->G);
                    matrix_free(g->V);
                    matrix_free(g->K);
                    free(g->W);
                }

            }
            for (k = 0; k < graph_dim(g); ++k)
                y[k] += mat_elem(V, i, k) / (double)(Gp[j+1] - Gp[j]);
        }
        
        g->G = add_vertex_to_vertex(G, j);
        g->V = matrix_insert_row(V, -1, 0);
        g->K = matrix_insert_row(K, -1, 0);
        g->W = malloc((n+1) * sizeof(double));
        
        if ((Gp[j+1] - Gp[j]) > 1)
        {
            if (W[j] > 0)
            {
                // new node position is partition centroid
                for (k = 0; k < graph_dim(g); ++k)
                    mat_elem(g->V, n, k) = mat_elem(K, j, k);
            }
            else
            {
                // no data belong to parent node, so new
                // position is centroid of leaf positions
                for (k = 0; k < graph_dim(g); ++k)
                    mat_elem(g->V, n, k) = y[k];
            }
        }
        else
        {
            // new position is same distance and direction
            // as original leaf was from its neighbor
            for (k = 0; k < graph_dim(g); ++k)
            {
                mat_elem(g->V, n, k) = (
                    2 * mat_elem(V, j, k) - mat_elem(V, Gi[Gp[j]], k) );
            }
        }

        graph_partition(g);
        graph_laplacian(g);
        graph_optimize(g);

        if (min_score > (g->energy + g->error))
        {
            cs_spfree(Gmin);
            matrix_free(Vmin);
            matrix_free(Kmin);
            free(Wmin);
            Gmin = g->G;
            Vmin = g->V;
            Kmin = g->K;
            Wmin = g->W;
            min_energy = g->energy;
            min_error = g->error;
            min_score = g->energy + g->error;
        }
        else
        {
            cs_spfree(g->G);
            matrix_free(g->V);
            matrix_free(g->K);
            free(g->W);
        }

        visited[j] = 1;
    }

    free(visited);

    cs_spfree(G);
    matrix_free(V);
    matrix_free(K);
    free(W);

    g->G = Gmin;
    g->V = Vmin;
    g->K = Kmin;
    g->W = Wmin;
    g->energy = min_energy;
    g->error = min_error;
}
