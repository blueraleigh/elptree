if (FALSE) {



library(elptree3)

data(iris)

X = data.matrix(iris[, -5])

f = make.graph(X, 50L)
ans = f(100,100)
f = make.graph(X, 50L, ans$G, ans$V)
ans = f(10,10)
f = make.graph(X, 50L, ans$G, ans$V)
ans = f(1,1)
f = make.graph(X, 50L, ans$G, ans$V)
ans = f(.1,.1)
f = make.graph(X, 50L, ans$G, ans$V)
ans = f(.01,.01)
f = make.graph(X, 50L, ans$G, ans$V)
system.time(ans <- f(.01,.1))

plot(ans, cex=80, by=iris[,5], piebg=c(2,3,4))

plot(ans$PC$x, cex=0.8, pch="+", col=iris[,5])
edge = which(ans$G[, 2] < ans$G[, 1])
x = predict(ans$PC, ans$V)
x0 = x[ans$G[edge, 1], 1]
y0 = x[ans$G[edge, 1], 2]
x1 = x[ans$G[edge, 2], 1]
y1 = x[ans$G[edge, 2], 2]
segments(x0, y0, x1, y1, col=8)
points(x, pch=19, col=5, cex=0.4)


graph_plot(ans, cex=50)


V = environment(f)$V
G = environment(f)$G
a = foo(1, 1, G, V, X)
w = table(kmeans(X, V, 1L)$cluster) / 150
K = kmeans(X, V, 1L)$centers
solve(a, sweep(K, 1, w, "*"))


V = matrix(rnorm(20), 4, 5)
G = rbind(
    c(4,1),
    c(4,2),
    c(4,3),
    c(1,4),
    c(2,4),
    c(3,4)
)
storage.mode(G) = "integer"
X = matrix(rnorm(500), 100, 5)

foo = function(lambda, mu, G, V, X) 
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
    #print(L_star_leaf)
    #print(L_star_edge)
    #print((diag(colSums(L_star_edge),N,N) - L_star_edge) +
    #       (diag(colSums(L_star_leaf),N,N) - L_star_leaf))
    Lapl = diag(w,N,N) + 
           (diag(colSums(L),N,N) - L) +
           (diag(colSums(L_star_edge),N,N) - L_star_edge) +
           (diag(colSums(L_star_leaf),N,N) - L_star_leaf)
    return (Lapl)
}

lambda = 1
mu = .1
foo(lambda, mu, G, V, X)
graph.laplacian(lambda, mu, G, V, X)





}