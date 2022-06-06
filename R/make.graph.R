make.graph = function(x, size, G, V, ...)
{    
    PC = prcomp(x, ...)

    if (missing(G) || missing(V) || 
        is.null(G) || is.null(V) ||
        is.na(G) || is.na(V))
    {
        mu = colMeans(x)
        # initial vertex coordinates: move the global
        # mean along the first principal axis 1 standard
        # deviation in each direction
        V = rbind(
                mu - PC$sdev[1] * PC$rotation[, 1],
                mu + PC$sdev[1] * PC$rotation[, 1])

        # intial graph topology (edgelist format)
        G = matrix(c(1L,2L,2L,1L), 2L, 2L)
    }
    else
    {
        # should check that inputs make sense   
    }

    function(lambda=0.01, mu=0.1, alpha=0, R02=Inf, tol=1e-4)
    {
        ans = .Call(C_graph_fit,
            x,
            G,
            V,
            as.integer(size),
            as.numeric(lambda),
            as.numeric(mu),
            max(0, min(1, as.numeric(alpha))),
            as.numeric(R02),
            as.numeric(tol)
        )
        structure(
            list(
                G=ans[[1]],
                V=structure(ans[[2]], dimnames=list(NULL, colnames(x))),
                W=ans[[3]],
                energy=ans[[4]],
                error=ans[[5]],
                X=x,
                PC=PC
            ),
            class="elptree"
        )
    }
}
