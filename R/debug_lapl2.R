if (FALSE)
{


    
graph_energy = function(lambda, mu, G, V, X)
{
    N = nrow(V)
    visited = logical(N)
    notvisited = integer(N)
    
    n = 1L
    notvisited[n] = G[1L,1L]

    e = 0

    while (n > 0L) 
    {
        v = notvisited[n]
        n = n - 1L
        k = graph_node_degree(v, G)
        for (d in graph_node_children(v, G))
        {
            if (!visited[d]) # first vist
            {
                z = V[v, ] - V[d, ]
                e = e + lambda * (t(z) %*% z)[1]
                n = n + 1L
                notvisited[n] = d
            }
        }
        if (k > 1L)
        {
            z = V[v,] - colMeans(V[graph_node_children(v, G), ])
            e = e + mu * (t(z) %*% z)[1]
        }
        visited[v] = TRUE
    }
    e
}


graph_error = function(lambda, mu, G, V, X)
{
    e = 0
    for (i in 1:nrow(X)) 
    {
        d = min(as.matrix(dist(rbind(X[i,], V)))[1,-1])
        e = e + d^2
    }
    e / nrow(X)
}


graph_partition = function(lambda, mu, G, V, X)
{
    N = nrow(V)
    w = numeric(N)
    part = integer(nrow(X))
    K = matrix(0, N, ncol(X))
    for (i in 1:nrow(X)) 
    {
        nn = which.min(as.matrix(dist(rbind(X[i,], V)))[1,-1])
        part[i] = nn
        w[nn] = w[nn] + 1
    }
    w = w / nrow(X)

    for (i in 1:N) 
    {
        if (w[i] > 0)
        {
            idx = part == i
            K[i,] = colMeans(X[idx, , drop=FALSE])
        }
    }
    sweep(K, 1, w, "*")
}


graph_laplR = function(lambda, mu, G, V, X) 
{
    N = nrow(V)
    L = matrix(0, N, N)
    L_star_edge = matrix(0, N, N)
    L_star_leaf = matrix(0, N, N)

    w = numeric(N)
    for (i in 1:nrow(X)) 
    {
        nn = which.min(as.matrix(dist(rbind(X[i,], V)))[1,-1])
        w[nn] = w[nn] + 1
    }
    w = w / nrow(X)


    visited = logical(N)
    notvisited = integer(N)
    
    n = 1L
    notvisited[n] = G[1L,1L]

    while (n > 0L) 
    {
        v = notvisited[n]
        n = n - 1L
        k = graph_node_degree(v, G)
        for (d in graph_node_children(v, G))
        {
            if (!visited[d]) # first vist
            {
                L[v, d] = lambda
                n = n + 1L
                notvisited[n] = d
            }
            else # second visit from other direction
            {
                L[v, d] = lambda
            }
        }
        if (k > 1L)
        {
            for (i in graph_node_children(v, G)) {
                L_star_edge[v,i] = L_star_edge[v,i] + mu / k
                L_star_edge[i,v] = L_star_edge[i,v] + mu / k
            }
            for (i in graph_node_children(v, G)) {
                j = i + 1L
                if (j > k)
                    j = 1L
                L_star_leaf[i,j] = -mu / (k*k)
                L_star_leaf[j,i] = -mu / (k*k)
            }
        }
        visited[v] = TRUE
    }
    
    Lapl = diag(w) + 
           (diag(colSums(L)) - L) +
           (diag(colSums(L_star_edge)) - L_star_edge) +
           (diag(colSums(L_star_leaf)) - L_star_leaf)
    
    return (Lapl)
}

library(elptree)
data(iris)

X = data.matrix(iris[, -5])

V = matrix(rnorm(20), 4, 4)
G = rbind(
    c(4,1),
    c(4,2),
    c(4,3),
    c(1,4),
    c(2,4),
    c(3,4)
)
storage.mode(G) = "integer"

V = matrix(rnorm(12), 3, 4)
G = rbind(
    c(1,2),
    c(2,3),
    c(2,1),
    c(3,2)
)
storage.mode(G) = "integer"

f = make.graph(X, 50L)
G = environment(f)$G
V = environment(f)$V

score = matrix(0, 20, 3)
L = graph_laplR(.1, .1, G, V, X)
K = graph_partition(.1, .1, G, V, X)
for (i in 1:20) 
{
    score[i, ] = c(e1<-graph_energy(.1,.1,G,V,X), e2<-graph_error(.1,.1,G,V,X), e1+e2)
    V = solve(L, K)
    L = graph_laplR(.1, .1, G, V, X)
    K = graph_partition(.1, .1, G, V, X)
    #score[i, 2] = graph_energy(.1,.1,G,V,X) + graph_error(.1,.1,G,V,X)
}




}