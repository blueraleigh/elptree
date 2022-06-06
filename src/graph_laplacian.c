#include "graph.h"


static double stretching_penalty(struct graph *g, int i, int j)
{
    //HERE("enter stretching_penalty");

    int k;

    double d;
    double d2 = 0;

    struct matrix *V = graph_coordinates(g);

    for (k = 0; k < graph_dim(g); ++k)
    {
        d = mat_elem(V, i, k) - mat_elem(V, j, k);
        d2 += d * d;
    }

    int deg1 = (int)(g->G->p[i+1] - g->G->p[i]);
    int deg2 = (int)(g->G->p[j+1] - g->G->p[j]);

    return d2 * (g->lambda + g->alpha * (fmax2(2, fmax2(deg1, deg2)) - 2));
}


static double bending_penalty(struct graph *g, int j)
{
    //HERE("bending_penalty");

    int i;
    int k;
    int p;
    int n;

    double d;
    double d2 = 0;
    
    /* centroid coordinates of k-star leaves */
    double y[graph_dim(g)];

    struct matrix *V = graph_coordinates(g);

    csi *Gp = g->G->p;
    csi *Gi = g->G->i;

    assert((Gp[j+1] - Gp[j]) > 1);

    for (p = Gp[j], n = 1; p < Gp[j+1]; ++p, ++n)
    {
        i = Gi[p];
        if (n == 1)
        {
            for (k = 0; k < graph_dim(g); ++k)
                y[k] = mat_elem(V, i, k);
        }
        else
        {
            for (k = 0; k < graph_dim(g); ++k)
                y[k] += (mat_elem(V, i, k) - y[k]) / (double)n;
        }
    }

    for (k = 0; k < graph_dim(g); ++k)
    {
        d = mat_elem(V, j, k) - y[k];
        d2 += d * d;
    }

    return d2 * g->mu;
}


/* Compute the graph Laplacian using the method described in
**
** Gorban AN and Zinovyev AY. Principal Graphs and Manifolds. 
** In Handbook of Research on Machine Learning Applications 
** and Trends: Algorithms, Methods and Techniques 
** (eds. Olivas E.S., Guererro J.D.M., Sober M.M., Benedito J.R.M., Lopes A.J.S.). 
** Information Science Reference, September 4, 2009.
** 
** https://arxiv.org/ftp/arxiv/papers/0809/0809.0490.pdf
**
** see section "Optimisation of the elastic graph algorithm"
** in the linked PDF
*/
void graph_laplacian(struct graph *g)
{
    HERE("graph_laplacian");

    csi i;
    csi j;
    csi jj;
    csi p;
    csi q;

    cs *G = g->G;

    cs *Te = cs_spalloc(G->n, G->n, G->n * G->n, 1, 1);
    cs *Ts = cs_spalloc(G->n, G->n, G->n * G->n, 1, 1);

    double k;
    double energy = 0;
    double lambda = g->lambda;
    double mu = g->mu;
    double *W = graph_partition_weights(g);

    if (!Te || !Ts)
        return;

    csi *Gp = G->p;
    csi *Gi = G->i;

    for (j = 0; j < G->n; ++j)
    {
        k = (double)(Gp[j+1] - Gp[j]);
        for (p = Gp[j]; p < Gp[j+1]; ++p)
        {
            i = Gi[p];
            cs_entry(Te, i, i, lambda);
            cs_entry(Te, j, i, -lambda);

            // we visit each edge twice, so multiply by 1/2
            energy += 0.5 * stretching_penalty(g, j, i);
        }

        // NOTE: graph_partition should be called before
        // this function is invoked so that W is properly
        // initialized
        cs_entry(Te, j, j, W[j]);
        
        if (k > 1) // j is a k-star
        {
            energy += bending_penalty(g, j);
            
            cs_entry(Ts, j, j, mu);
            
            for (p = Gp[j]; p < Gp[j+1]; ++p)
            {
                i = Gi[p];
                cs_entry(Ts, j, i, -(mu / k));
                cs_entry(Ts, i, j, -(mu / k));
                
                for (q = Gp[j]; q < Gp[j+1]; ++q)
                    cs_entry(Ts, Gi[p], Gi[q], mu / (k * k));
            }
        }
    }

    cs *e = cs_compress(Te);
    cs *s = cs_compress(Ts);

    cs_spfree(Te);
    cs_spfree(Ts);

    if (!e || !s) 
        return;

    /* sum up entries with duplicate (i,j) indices */
    cs_dupl(e);
    cs_dupl(s);
    
    cs *L = cs_add(e, s, 1, 1);
    
    cs_spfree(e);
    cs_spfree(s);

    if (g->L)
        cs_spfree(g->L);

    g->L = L;
    g->energy = energy;
}


SEXP C_graph_laplacian(SEXP lambda, SEXP mu, SEXP G, SEXP V, SEXP X)
{
    HERE("C_graph_laplacian");

    struct graph *g = graph_init(
        X, G, V, lambda, mu, ScalarReal(0), 
        ScalarReal(R_PosInf), ScalarReal(0.0001));
    
    if (!g)
        error("could not allocate memory for embedding");
    
    graph_partition(g);
    graph_laplacian(g);

    int j;
    int i;
    int p;
    int n = graph_nnode(g);

    SEXP ans = PROTECT(allocMatrix(REALSXP, n, n));

    memset(REAL(ans), 0, n * n * sizeof(double));

    for (j = 0; j < n; ++j)
    {
        for (p = g->L->p[j]; p < g->L->p[j+1]; ++p)
        {
            i = g->L->i[p];
            REAL(ans)[i + j*n] = g->L->x[p];
        }
    }

    graph_free(g);

    UNPROTECT(1);

    return ans;
}
