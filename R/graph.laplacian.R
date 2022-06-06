graph.laplacian = function(lambda, mu, G, V, X) 
{
    .Call(C_graph_laplacian, lambda, mu, G, V, X)
}