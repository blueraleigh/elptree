#ifndef GRAPH_H
#define GRAPH_H

#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <cs.h>
#include "matrix.h"

/* Modified from assert.h included with max os x */
#undef assert
#ifdef NDEBUG
#define assert(e) ((void)0)
#else
#define elptree_assert(e, file, line) \
    error("%s:%u: failed assertion `%s'\n", file, line, e)
#define assert(e) \
    ((void) ((e) ? ((void)0) : elptree_assert(#e, __FILE__, __LINE__)))
#endif

#ifdef VDEBUG
#define HERE(msg) Rprintf("%s (%s:%u)\n", msg, __FILE__, __LINE__)
#else
#define HERE(msg) ((void)0)
#endif

struct graph {

    // stretching elasticity modulus
    double lambda;

    // bending elasticity modulus
    double mu;

    // complex branching penalty
    double alpha;

    // trimming distance
    double R02;

    double tol;

    double energy;

    double error;

    /* partition weights */
    double *W;

    /* data coordinates */
    struct matrix *X;

    /* centroids of partitioned data space */
    struct matrix *K;

    /* embedding coordinates of graph vertices */
    struct matrix *V;

    /* adjacency matrix in compressed sparse column form */
    cs *G;

    /* Laplacian of G */
    cs *L;
};

#define graph_ndata(g) (g)->X->n
#define graph_nnode(g) (g)->G->n
#define graph_nedge(g) (g)->G->p[(g)->G->n]
#define graph_dim(g) (g)->X->p
#define graph_data(g) (g)->X
#define graph_coordinates(g) (g)->V
#define graph_partition_centers(g) (g)->K
#define graph_partition_weights(g) (g)->W

// spring graph laplacian
void graph_laplacian(struct graph *g);

void graph_partition(struct graph *g);

void graph_optimize(struct graph *g);

void graph_grow(struct graph *g);

void graph_shrink(struct graph *g);

struct graph *graph_init(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

void graph_free(struct graph *g);

struct matrix *graph_dist(struct graph *g);

#endif
