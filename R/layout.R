graph_node_degree = function(i, G)
{
    length(which(G[,1L] == i))
}


graph_node_children = function(i, G) 
{
    G[which(G[,1L] == i), 2L]
}


graph_iter = function(i, G) 
{   
    N = graph_size(G)

    visited = logical(N)
    notvisited = integer(N)
    
    n = 1L
    notvisited[n] = i

    function() 
    {
        if (n <= 0)
            return (0L)
        v = notvisited[n]
        n <<- n - 1L
        for (d in graph_node_children(v, G)) 
        {
            if (!visited[d]) 
            {
                n <<- n + 1L
                notvisited[n] <<- d
            }
        }
        visited[v] <<- TRUE
        return (v)
    }
}


graph_size = function(G) 
{
    length(unique(G[, 2]))
}


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
