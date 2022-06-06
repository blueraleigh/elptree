#' Return the index of the vertex that best approximates each
#' row of data in `newdata`
predict.elptree = function(object, newdata, ...)
{
    V = object$V
    x = data.matrix(newdata)
    stopifnot(ncol(x) == ncol(V))
    apply(x, 1, function(y) {
        M = sweep(V, 2, y)
        which.min(diag(M %*% t(M)))
    })
}