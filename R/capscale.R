`capscale` <-
    function (formula, data, distance = "euclidean", sqrt.dist = FALSE,
              comm = NULL, add = FALSE, dfun = vegdist,
              metaMDSdist = FALSE, na.action = na.fail, subset = NULL, ...)
{
    print('#0#')
    EPS <- sqrt(.Machine$double.eps)
    if (!inherits(formula, "formula"))
        stop("needs a model formula")
    if (missing(data)) {
        data <- parent.frame()
    }
    else {
        data <- eval(match.call()$data, environment(formula),
                     enclos = .GlobalEnv)
    }
    print('#1#')
    print('#1.1#')
    print(class(formula))
    print('#1.2#')
    print(formula)
    print('#1.3#')
    print(class(data))
    print('#1.4#')
    print(colnames(data))
    print(length(data))
    print('#1.5#')
    formula <- formula(terms(formula, data = data))
    print('#2#')
    ## The following line was eval'ed in environment(formula), but
    ## that made update() fail. Rethink the line if capscale() fails
    ## mysteriously at this point.
    X <- eval(formula[[2]], envir=environment(formula),
              enclos = globalenv())
    print('#3#')
    ## see if user supplied dissimilarities as a matrix
    if ((is.matrix(X) || is.data.frame(X)) &&
        isSymmetric(unname(as.matrix(X))))
        X <- as.dist(X)
    if (!inherits(X, "dist")) {
        print('#4#')
        comm <- X
        vdata <- as.character(formula[[2]])
        dfun <- match.fun(dfun)
        if (metaMDSdist) {
            commname <- as.character(formula[[2]])
            X <- metaMDSdist(comm, distance = distance, zerodist = "ignore",
                             commname = commname, distfun = dfun, ...)
            commname <- attr(X, "commname")
            comm <- eval.parent(parse(text=commname))
        } else {
            X <- dfun(X, distance)
        }
        print('#5#')
    } else { # vdata name
        print('#6#')
        if (missing(comm))
            vdata <- NULL
        else
            vdata <- deparse(substitute(comm))
        print('#6#')
    }
    inertia <- attr(X, "method")
    print('#7#')
    if (is.null(inertia))
        inertia <- "unknown"
    inertia <- paste(toupper(substr(inertia, 1, 1)),
                     substring(inertia,  2), sep = "")
    inertia <- paste(inertia, "distance")
    print('#8#')
    if (!sqrt.dist)
        inertia <- paste("squared", inertia)
    print('#9#')
    ## postpone info on euclidification till we have done so

    ## evaluate formula: ordiParseFormula will return dissimilarities
    ## as a symmetric square matrix (except that some rows may be
    ## deleted due to missing values)
    d <- ordiParseFormula(formula,
                          data,
                          na.action = na.action,
                          subset = substitute(subset),
                          X = X)
    print('#10#')
    ## ordiParseFormula subsets rows of dissimilarities: do the same
    ## for columns ('comm' is handled later). ordiParseFormula
    ## returned the original data, but we use instead the potentially
    ## changed X and discard d$X.
    if (!is.null(d$subset)) {
        X <- as.matrix(X)[d$subset, d$subset, drop = FALSE]
    }
    print('#11#')
    ## Delete columns if rows were deleted due to missing values
    if (!is.null(d$na.action)) {
        X <- as.matrix(X)[-d$na.action, -d$na.action, drop = FALSE]
    }
    print('#12#')
    X <- as.dist(X)
    k <- attr(X, "Size") - 1
    if (sqrt.dist)
        X <- sqrt(X)
    print('#13#')
    if (max(X) >= 4 + .Machine$double.eps) {
        inertia <- paste("mean", inertia)
        adjust <- sqrt(k)
        X <- X/adjust
    }
    else {
        adjust <- 1
    }
    print('#14#')
    nm <- attr(X, "Labels")
    ## wcmdscale, optionally with additive adjustment
    X <- wcmdscale(X, x.ret = TRUE, add = add)
    if(any(dim(X$points) == 0)) # there may be no positive dims
        X$points <- matrix(0, NROW(X$points), 1)
    print('#15#')
    ## this may have been euclidified: update inertia
    if (!is.na(X$ac) && X$ac > sqrt(.Machine$double.eps))
        inertia <- paste(paste0(toupper(substring(X$add, 1, 1)),
                                substring(X$add, 2)),
                         "adjusted", inertia)
    if (is.null(rownames(X$points)))
        rownames(X$points) <- nm

    print('#16#')
    sol <- ordConstrained(X$points, d$Y, d$Z, method = "capscale")
    print('#17#')

    ## update for negative eigenvalues
    poseig <- length(sol$CA$eig)
    if (any(X$eig < 0)) {
        print('#18#')
        negax <- X$eig[X$eig < 0]
        sol$CA$imaginary.chi <- sum(negax)
        sol$tot.chi <- sol$tot.chi + sol$CA$imaginary.chi
        sol$CA$imaginary.rank <- length(negax)
        sol$CA$imaginary.u.eig <- X$negaxes
        print('#19#')
    }
    if (!is.null(comm)) {
        print('#20#')
        sol$vdata <- vdata
        comm <- scale(comm, center = TRUE, scale = FALSE)
        sol$colsum <- apply(comm, 2, sd)
        ## take a 'subset' of the community after scale()
        if (!is.null(d$subset))
            comm <- comm[d$subset, , drop = FALSE]
        ## NA action after 'subset'
        if (!is.null(d$na.action))
            comm <- comm[-d$na.action, , drop = FALSE]
        if (!is.null(sol$pCCA) && sol$pCCA$rank > 0)
            comm <- qr.resid(sol$pCCA$QR, comm)
        if (!is.null(sol$CCA) && sol$CCA$rank > 0) {
            v.eig <- t(comm) %*% sol$CCA$u/sqrt(k)
            sol$CCA$v <- decostand(v.eig, "normalize", MARGIN = 2)
            comm <- qr.resid(sol$CCA$QR, comm)
        }
        if (!is.null(sol$CA) && sol$CA$rank > 0) {
            v.eig <- t(comm) %*% sol$CA$u/sqrt(k)
            sol$CA$v <- decostand(v.eig, "normalize", MARGIN = 2)
        }
        print('#21#')
    } else {
        ## input data were dissimilarities, and no 'comm' defined:
        ## species scores make no sense and are made NA
        print('#22#')
        sol$CA$v[] <- NA
        if (!is.null(sol$CCA))
            sol$CCA$v[] <- NA
        sol$colsum <- NA
        print('#23#')
    }
    print('#24#')
    if (!is.null(sol$CCA) && sol$CCA$rank > 0)
        sol$CCA$centroids <- centroids.cca(sol$CCA$wa, d$modelframe)
    if (!is.null(sol$CCA$alias))
        sol$CCA$centroids <- unique(sol$CCA$centroids)
    if (!is.null(sol$CCA$centroids)) {
        rs <- rowSums(sol$CCA$centroids^2)
        sol$CCA$centroids <- sol$CCA$centroids[rs > 1e-04, ,
                                               drop = FALSE]
        if (nrow(sol$CCA$centroids) == 0)
            sol$CCA$centroids <- NULL
    }
    print('#25#')
    sol$call <- match.call()
    sol$terms <- terms(formula, "Condition", data = data)
    sol$terminfo <- ordiTerminfo(d, data)
    sol$call$formula <- formula(d$terms, width.cutoff = 500)
    sol$call$formula[[2]] <- formula[[2]]
    sol$sqrt.dist <- sqrt.dist
    print('#26#')
    if (!is.na(X$ac) && X$ac > 0) {
        sol$ac <- X$ac
        sol$add <- X$add
    }
    print('#27#')
    sol$adjust <- adjust
    sol$inertia <- inertia
    if (metaMDSdist)
        sol$metaMDSdist <- commname
    print('#28#')
    sol$subset <- d$subset
    sol$na.action <- d$na.action
    class(sol) <- c("capscale", "rda", "cca")
    print('#29#')
    if (!is.null(sol$na.action))
        sol <- ordiNAexclude(sol, d$excluded)
    print('#30#')
    sol
}
