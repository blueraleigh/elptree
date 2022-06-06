#include "graph.h"

// Returns a copy of graph G with the terminal vertex i removed
static cs *delete_terminal_vertex(cs *G, csi i)
{
    HERE("delete_terminal_vertex");

    csi v;
    csi p;
    csi j;

    csi *Gp = G->p;
    csi *Gi = G->i;

    assert((Gp[i+1] - Gp[i]) == 1);

    cs *D;
    cs *DT = cs_spalloc(G->n-1, G->n-1, G->nzmax - 2, 1, 1);

    // vertex i is being removed. vertex indices 
    // greater than i must be decremented

    for (v = 0; v < G->n; ++v)
    {
        if (v == i)
            continue;
        j = (v < i) ? v : v - 1;
        for (p = Gp[v]; p < Gp[v+1]; ++p)
        {
            if (Gi[p] < i)
                cs_entry(DT, Gi[p], j, 1);
            else if (Gi[p] > i)
                cs_entry(DT, Gi[p] - 1, j, 1);
            else // Gi[p] == i
                continue;
        }
    }

    D = cs_compress(DT);

    cs_spfree(DT);

    return D;
}


// Returns a copy of graph G with vertex i deleted and all of i's neighbors
// (excluding j) connected to j 
static cs *shrink_internal_edge(cs *G, csi i, csi j)
{
    HERE("shrink_internal_edge");

    csi v;
    csi p;
    csi s;
    csi d;

    csi *Gp = G->p;
    csi *Gi = G->i;

    assert((Gp[i+1] - Gp[i]) > 1);
    assert((Gp[j+1] - Gp[j]) > 1);

    cs *D;
    cs *DT = cs_spalloc(G->n-1, G->n-1, G->nzmax - 2, 1, 1);

    // vertex i is being removed. vertex indices 
    // greater than i must be decremented

    csi jnew = (j < i) ? j : j - 1;

    for (v = 0; v < G->n; ++v)
    {
        if (v == j)
            s = jnew;
        else if (v == i) // remap i's edges to j
            s = jnew;
        else
            s = (v < i) ? v : v - 1;
        for (p = Gp[v]; p < Gp[v+1]; ++p)
        {
            if (Gi[p] == i) // remap i's edges to j
                d = jnew;
            else if (Gi[p] == j)
                d = jnew;
            else
                d = (Gi[p] < i) ? Gi[p] : Gi[p] - 1;

            if (s == jnew && d == jnew)
                continue; // remove the edges (i,j) and (j,i)

            cs_entry(DT, s, d, 1);
        }
    }

    D = cs_compress(DT);

    cs_spfree(DT);

    return D;
}


void graph_shrink(struct graph *g)
{
    HERE("graph_shrink");
    
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

    int *visited = calloc(n, sizeof(int));

    for (j = 0; j < n; ++j)
    {
        for (p = Gp[j]; p < Gp[j+1]; ++p)
        {
            i = Gi[p];
            if (visited[i] == 0 && (Gp[j+1] - Gp[j]) > 1 && (Gp[i+1] - Gp[i]) > 1)
            {
                // first visit to internal edge (j,i)

                g->G = shrink_internal_edge(G, j, i);
                g->V = matrix_delete_row(V, j);
                g->K = matrix_delete_row(K, j);
                g->W = malloc((n-1) * sizeof(double));

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
        }

        if ((Gp[j+1] - Gp[j]) == 1)
        {
            g->G = delete_terminal_vertex(G, j);
            g->V = matrix_delete_row(V, j);
            g->K = matrix_delete_row(K, j);
            g->W = malloc((n-1) * sizeof(double));

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
