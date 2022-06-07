#' Layout the elastic principle tree in the plane
#'
#' Uses stress majorization to compute an (x,y) coordinate
#' for each vertex.
#'
#' @references 
#' Graph Drawing by Stress Majorization.
#' Emden R. Gansner, Yehuda Koren and Stephen North.
#' In Proceedings 12th Symposium on Graph Drawing (GD), 
#' pages 239–250, 2004.
graph_layout = function(object)
{
    stopifnot(inherits(object, "elptree"))

    G = object$G
    V = object$V
    X = object$X
    PC = object$PC

    # Initial layout
    x = predict(PC, V)[, 1:2]

    # Fix the first vertex at (0,0)
    x = sweep(x, 2, x[1,])

    .Call(C_graph_layout, x, G, V, X)
}
