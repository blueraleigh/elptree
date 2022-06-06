#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>

#define CALLDEF(name, n)  {#name, (DL_FUNC) &name, n}

SEXP C_graph_fit(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
SEXP C_graph_laplacian(SEXP,SEXP,SEXP,SEXP,SEXP);
SEXP C_graph_update(SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP C_graph_layout(SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    CALLDEF(C_graph_fit, 9),
    CALLDEF(C_graph_laplacian, 5),
    CALLDEF(C_graph_update, 5),
    CALLDEF(C_graph_layout, 4),
    {NULL, NULL, 0}
};


void attribute_visible R_init_elptree(DllInfo *info)
{
    R_registerRoutines(info, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(info, FALSE);
    R_forceSymbols(info, TRUE);
}
