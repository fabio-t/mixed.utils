
'%!in%' <- function(x,y)!('%in%'(x,y))

colSd     <- function (x, na.rm=FALSE) apply(X=x, MARGIN=2, FUN=sd,     na.rm=na.rm)
colMax    <- function (x, na.rm=FALSE) apply(X=x, MARGIN=2, FUN=max,    na.rm=na.rm)
colMin    <- function (x, na.rm=FALSE) apply(X=x, MARGIN=2, FUN=min,    na.rm=na.rm)
colMedian <- function (x, na.rm=FALSE) apply(X=x, MARGIN=2, FUN=median, na.rm=na.rm)
rowSd     <- function (x, na.rm=FALSE) apply(X=x, MARGIN=1, FUN=sd,     na.rm=na.rm)
rowMax    <- function (x, na.rm=FALSE) apply(X=x, MARGIN=1, FUN=max,    na.rm=na.rm)
rowMin    <- function (x, na.rm=FALSE) apply(X=x, MARGIN=1, FUN=min,    na.rm=na.rm)
rowMedian <- function (x, na.rm=FALSE) apply(X=x, MARGIN=1, FUN=median, na.rm=na.rm)

# Get lower triangle of the correlation matrix
  get_lower_tri<-function(cormat, diag = F){
    cormat[upper.tri(cormat, !diag)] <- NA
    return(cormat)
  }
  # Get upper triangle of the correlation matrix
  get_upper_tri <- function(cormat, diag = F){
    cormat[lower.tri(cormat, !diag)]<- NA
    return(cormat)
  }

cor.test.p <- function(x, method="pearson") {
    # https://stackoverflow.com/a/13112337
    FUN <- function(x, y) cor.test(x, y, method = method)[["p.value"]]
    z <- outer(
        colnames(x),
        colnames(x),
        Vectorize(function(i, j) FUN(x[, i], x[, j]))
    )
    dimnames(z) <- list(colnames(x), colnames(x))
    z
}

count.pairwise <- function(x, y) {
    # same as count.pairwise(x, y) from psych package
    n <- t(!is.na(x)) %*% (!is.na(y))
    n
}

cor2pvalue <- function(r, n) {
  t <- (r*sqrt(n-2))/sqrt(1-r^2)
  # p <- 2*(1 - pt(abs(t),(n-2)))
  p <- -2 *  expm1(pt(abs(t),(n-2),log.p=TRUE))
  se <- sqrt((1-r*r)/(n-2))
  out <- list(r, n, t, p, se)
  names(out) <- c("r", "n", "t", "p", "se")
  return(out)
}

this.file.name <- function() {
    # https://stackoverflow.com/a/1816487
    frame_files <- lapply(sys.frames(), function(x) x$ofile)
    frame_files <- Filter(Negate(is.null), frame_files)

    if (length(frame_files) > 0) {
        frame_files[[length(frame_files)]]
    } else {
        ""
    }
}

this.dir.name <- function() {
    cmdArgs <- commandArgs(trailingOnly = FALSE)
    needle <- "--file="
    match <- grep(needle, cmdArgs)
    if (length(match) > 0) {
        # Rscript
        return(dirname(normalizePath(sub(needle, "", cmdArgs[match]))))
    } else {
        # 'source'd via R console
        filename <- this.file.name()

        if (filename != "") {
            return(dirname(normalizePath(filename)))
        } else {
            # if we couldn't find any file, we assume this was
            # loaded via Hydrogen (and we also assume that the kernel is
            # set to startup where the initial file is, or this won't work)
            print("Warning: probably using Hydrogen to load file?")
            return(getwd())
        }
    }
}

hclustavg <- function(x) {
    hclust(x, method = "average")
}

hclustward <- function(x) {
    hclust(x, method = "ward.D2")
}

hclustfun <- function(x) {
    hclust(x, method = "complete")
}

cordist <- function(x) {
    dd <- as.dist((1 - cor(t(x), use = "pairwise.complete.obs")) / 2, diag = T)
    dd[which(is.na(dd))] <- 0.5
    dd
}

abscordist <- function(x) {
    dd <- as.dist(1 - abs(cor(t(x), use = "pairwise.complete.obs")), diag = T)
    dd[which(is.na(dd))] <- 0
    dd
}

reorderfun <- function(d, w) reorder(d, w)

r2 <- function(preds, actual) {
    rss <- sum((preds - actual)^2, na.rm = T)
    tss <- sum((actual - mean(actual, na.rm = T))^2, na.rm = T)
    rsq <- 1 - rss / tss

    return(rsq)
}

aic <- function(fit) {
    tLL <- fit$nulldev - deviance(fit)
    k <- fit$df
    n <- fit$nobs
    AICc <- -tLL + 2 * k + 2 * k * (k + 1) / (n - k - 1)

    return(AICc)
}

count_norm <- function(countData, types = NULL) {
    if (is.null(types)) {
        design <- ~1
        colData <- data.frame(type = rep(1, ncol(countData)))
    }
    else {
        design <- ~type
        colData <- data.frame(type = factor(types))
    }

    rownames(colData) <- colnames(countData)
    dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = design)
    featureData <- data.frame(gene = rownames(countData))
    (mcols(dds) <- DataFrame(mcols(dds), featureData))

    # you might want to save the structure away, before you play around with things
    # save(dds, file="dds.rda")

    ## now for the actual normalisation. You'll have to experiment
    ## and evaluate the output to decide which one is best for you.
    ## if you add the parameter blind=FALSE to these three functions,
    ## it should be much faster, but I don't know the consequences :)
    ## always evaluate with boxplots or whatever

    # very slow
    # out <- rlog(dds)
    # a bit faster
    out <- varianceStabilizingTransformation(dds, fitType = "local")
    # much faster - it uses only 1000 rows chosen randomly to fit the internal model, though
    # out <- vst(dds)

    # get the data, same for either of the three above
    data_norm <- assay(out)

    return(data_norm)
}

transf_train_test <- function(data, variables, train) {
    tmp <- data

    data <- data[variables]

    t_fun <- c(identity, function(x) {
        log(x + 1)
    }, sqrt, function(x) {
        x**2
    })
    norm_fun <- function(d, tf) {
        sapply(d, function(x) {
            shapiro.test(tf(x))$p.value
        })
    }
    df <- data.frame(
        Identity = norm_fun(data[train, , drop = F], t_fun[[1]]),
        Log = norm_fun(data[train, , drop = F], t_fun[[2]]),
        Sqrt = norm_fun(data[train, , drop = F], t_fun[[3]]),
        Squared = norm_fun(data[train, , drop = F], t_fun[[4]])
    )
    # print(df)
    norm_best_i <- max.col(replace(df, is.na(df), -Inf), ties.method = "first")
    # print(norm_best_i)
    norm_best_f <- t_fun[norm_best_i]
    # print(norm_best_f)
    data_t <- data.frame(sapply(
        1:ncol(data),
        function(i) norm_best_f[[i]](data[, i, drop = F])
    ))
    rownames(data_t) <- rownames(data)
    colnames(data_t) <- colnames(data)

    tmp[variables] <- data_t

    return(list(pvalues=df, transf_data=tmp))
}

make_path <- function(dirname) {
    if (!is.null(dirname)) {
        dirname <- paste0(dirname, "/")

        if (!dir.exists(dirname)) {
            dir.create(dirname, mode = "0755", recursive = T)
        }
    }

    return(dirname)
}

intersect_dist <- function(df) {
    n <- ncol(df)
    mat <- matrix(0, ncol = n, nrow = n)
    colnames(mat) <- rownames(mat) <- colnames(df)

    for (i in 1:nrow(mat))
    {
        for (j in 1:ncol(mat))
        {
            mat[i, j] <- length(which(df[, i] != 0 & df[, j] != 0))
        }
    }

    return(mat)
}

longitudinal_dist <- function(data, variables, abs = F, dist = T) {
    devtools::install_github("lmarusich/rmcorr")

    # just a hack to have a nice variable-by-variable matrix
    df <- as.data.frame(matrix(, length(variables), length(variables), dimnames = list(variables, variables)))
    diag(df) <- 1

    df_pv <- df

    # pairs of variables
    cmbs <- as.data.frame(combn(variables, 2))

    for (var_pair in cmbs) {
        var_pair <- as.character(var_pair)

        res <- rmcorr(PatientId, var_pair[1], var_pair[2], data)

        df[var_pair[1], var_pair[2]] <- res$r
        df[var_pair[2], var_pair[1]] <- res$r

        df_pv[var_pair[1], var_pair[2]] <- res$p
        df_pv[var_pair[2], var_pair[1]] <- res$p
    }

    l <- list()

    if (dist) {
        if (abs) {
            df <- 1 - abs(df)
        }
        else {
            df <- 1 - ((df + 1) / 2)
        }

        l$d <- as.dist(df)
    }
    else {
        df <- as.matrix(df)

        if (abs) {
            l$d <- abs(df)
        }
        else {
            l$d <- df
        }
    }

    l$pv <- df_pv

    return(l)
}

get_top_elements <- function(data, topN = NULL, topPerc = NULL, min_clones = 3, exclude = NULL, keep_excluded) {
    columns <- setdiff(colnames(data), exclude)

    if (!is.null(topN)) {
        rows <- c()

        topN <- max(topN, min_clones)

        # for each column, take the topN clones (non-null rows)
        # and merge them together.
        for (column in columns)
        {
            df <- data[column]
            df <- subset(df, df > 0)
            df <- df[order(df[, c(1)], decreasing = T), c(1), drop = F]
            df <- df[1:min(nrow(df), topN), , drop = F]

            df_rows <- rownames(df)

	    if (keep_excluded) {

            # fix for issue identified on 29/05/2017:
            # if we take the union, some of the clones left out will
            # essentially "come back" for each organ.
            data[setdiff(rownames(data), df_rows), column] <- 0

            rows <- union(rows, df_rows)
        }
    }
    else if (!is.null(topPerc)) {
        rows <- c()

        # for each column, take the topPerc clones (non-null rows)
        # and merge them together.
        for (column in columns)
        {
            df <- data[column]
            if (sum(df, na.rm=T) == 0) {
                print("WARNING: column is made of zeros")
                next
            }
            df <- subset(df, df > 0)
            df <- df[order(df[, c(1)], decreasing = T), c(1), drop = F]
            abundance <- sum(df) * topPerc
            df_rows <- rownames(df)[cumsum(df) <= abundance]

            if (length(df_rows) < min_clones) {
                # we force a minimum number of clones
                df_rows <- rownames(df)[1:min_clones]
            }

            # fix for issue identified on 29/05/2017:
            # if we take the union, some of the clones left out will
            # essentially "come back" for each organ.
            data[setdiff(rownames(data), df_rows), column] <- 0

            rows <- union(rows, df_rows)
        }
    }
    else {
        rows <- rownames(data)
    }

    return(data[rows, ])
}

make_heatmap <- function(orig_data, dirname = "./", prefix = "heatmap", vars = colnames(orig_data),
                         rows = rownames(orig_data), threshold = NULL, k = 4, dist_method = "pearson",
                         abs = F, abs2 = abs, method = "complete", square = T, k2 = k, output = cairo_pdf,
                         cexRow = 1, cexCol = 1, vlim = c(0, 1), na.fix = T, scale = "none",
                         colRow = NULL, colCol = NULL) {
    if (identical(output, svg)) {
        ext <- ".svg"
    }
    else if (identical(output, pdf) || identical(output, cairo_pdf)) {
        ext <- ".pdf"
    }
    else if (identical(output, postscript)) {
        ext <- ".eps"
    }
    else {
        print("WARNING: output is neither svg nor pdf/cairo_pdf, so we default to png")
        output <- png
        ext <- ".png"
    }

    l <- list()

    path <- make_path(dirname)

    print(paste0(path, prefix))

    data <- orig_data[rows, vars]

    if (na.fix) {
        # We remove all-NA rows and columns

        data <- data[rowSums(is.na(data)) < ncol(data), ]
        data <- data[, colSums(is.na(data)) < nrow(data)]
    }

    # scaling
    if (scale == "row") {
        data_scale <- scale(t(data), center = T, scale = T)
    }
    else if (scale == "column") {
        data_scale <- scale(data, center = T, scale = T)
    }
    else if (scale == "rowperc") {
        data_scale <- sweep(data, 1, rowSums(data), "/")
    }
    else if (scale == "columnperc") {
        data_scale <- sweep(data, 2, colSums(data), "/")
    }
    else {
        data_scale <- data
    }

    tdata <- t(data)
    tdata_scale <- t(data_scale)

    set1_cols <- brewer.pal(9, "Set1")

    output(paste0(path, prefix, ext))

    ## we get a similarity matrix first

    if (is.null(dist_method)) {
        # this assumes that the data already IS a similarity
        # matrix between 0 and 1
        dd <- data

        # we don't assume anything in case of missing value
        dd[which(is.na(dd))] <- 0.5

        dd2 <- as.dist(1 - dd)
    }
    else if (dist_method == "intersect") {
        dd <- as.matrix(intersect_dist(data))

        dd[which(is.na(dd))] <- 0

        dd2 <- as.dist(t(apply(dd, 1, function(x) max(x, na.rm = T) - x)))
    }
    else if (dist_method == "longitudinal") {
        dd <- longitudinal_dist(orig_data[rows, ], vars, abs = F, dist = F)$d

        dd[which(is.na(dd))] <- 0
        diag(dd) <- 1

        dd2 <- as.dist(1 - dd)
    }
    else if (grepl("(pearson|kendall|spearman)_pv", dist_method)) {
        r <- cor(data, use = "pairwise.complete.obs", method = dist_method)
        res <- corr2pvalue(r, count.pairwise(data, data))
        pv <- res$p

        dd <- sign(res$r) * -log10(pv)

        dd[which(is.na(dd))] <- 0
        diag(dd) <- 1

        dd2 <- as.dist(1 - dd)
    }
    else if (grepl("(pearson|kendall|spearman)", dist_method)) {
        dd <- cor(data, use = "pairwise.complete.obs", method = dist_method)

        dd[which(is.na(dd))] <- 0
        diag(dd) <- 1

        dd2 <- as.dist(1 - dd)
    }
    else {
        library(vegan)

        dd <- 1 - as.matrix(vegdist(tdata,
            method = dist_method,
            diag = T, upper = T, binary = F
        ))
        dd[which(is.na(dd))] <- 0
        dd2 <- as.dist(1 - dd)
    }

    print(head(dd))

    l$data <- dd

    ## apply absolute only to distance object

    if (abs) {
        dd2 <- abs(dd)
    }

    l$data2 <- dd2

    ## and finally we cluster

    cl <- hclust(dd2, method = method)
    l$cl <- cl

    if (k > 0) {
        clusters <- cutree(cl, k = k)
        sidecols <- as.character(set1_cols[clusters])
        l$clusters <- clusters
    }

    rowv <- reorderfun(as.dendrogram(cl), colMeans(data, na.rm = T))

    if (square) {
        # we want "co-correlations", eg row-vs-row

        if (is.null(vlim)) {
            vlim <- range(as.matrix(dd), finite = T)
        }

        if (is.na(vlim[1])) {
            vlim[1] <- min(dd, na.rm = T)
        }

        if (is.na(vlim[2])) {
            vlim[2] <- max(dd, na.rm = T)
        }

        if (vlim[1] < 0) {
            colours <- colorRampPalette(c("blue", "ghostwhite", "red"), space = "rgb")(99)
        }
        else {
            colours <- colorRampPalette(c("ghostwhite", "red"), space = "rgb")(99)
        }

        bins <- seq(vlim[1], vlim[2], length = 100)

        if (k > 0) {
            heatmap.2(dd,
                dendrogram = "row", RowSideColors = sidecols, density.info = "none", trace = "none",
                Rowv = rowv, Colv = rowv, cexRow = cexRow, cexCol = cexCol, margins = c(9, 9),
                srtCol = 90, col = colours, breaks = bins, symkey = F, lhei = c(1, 5), scale = scale,
                colRow = colRow, colCol = colCol
            )
        }
        else {
            heatmap.2(dd,
                dendrogram = "row", density.info = "none", trace = "none",
                Rowv = rowv, Colv = rowv, cexRow = cexRow, cexCol = cexCol, margins = c(9, 9),
                srtCol = 90, col = colours, breaks = bins, symkey = F, lhei = c(1, 5), scale = scale,
                colRow = colRow, colCol = colCol
            )
        }
    }
    else {
        if (is.null(vlim)) {
            vlim <- range(as.matrix(data_scale), finite = T)
        }

        if (is.na(vlim[1])) {
            vlim[1] <- min(data_scale, na.rm = T)
        }

        if (is.na(vlim[2])) {
            vlim[2] <- max(data_scale, na.rm = T)
        }

        if (vlim[1] < 0) {
            colours <- colorRampPalette(c("blue", "ghostwhite", "red"), space = "rgb")(99)

            # colours <- rev(brewer.pal(11,"Spectral"))
            # colours <- colorRampPalette(colours, space="rgb")(99)
        }
        else {
            colours <- colorRampPalette(c("ghostwhite", "red"), space = "rgb")(99)
        }

        bins <- seq(vlim[1], vlim[2], length = 100)

        # rows-vs-columns in the heatmap (both with dendrograms)

        if (is.null(dist_method)) {
            # this assumes that the data already IS a similarity matrix between 0 and 1
            dd <- tdata

            dd[which(is.na(dd))] <- 0.5 # we don't assume anything in case of missing value

            dd2 <- as.dist(1 - dd)
        }
        else if (dist_method == "intersect") {
            dd <- as.matrix(intersect_dist(tdata))

            dd[which(is.na(dd))] <- 0

            dd2 <- as.dist(t(apply(dd, 1, function(x) max(x, na.rm = T) - x)))
        }
        else if (dist_method == "longitudinal") {
            stop("ERROR: only square heatmap possible with longitudinal distance")
        }
        else if (grepl("(pearson|kendall|spearman)_pv", dist_method)) {
            r <- cor(data, use = "pairwise.complete.obs", method = dist_method)
            res <- corr2pvalue(r, count.pairwise(data, data))
            pv <- res$p

            dd <- sign(res$r) * -log10(pv)

            dd[which(is.na(dd))] <- 0
            diag(dd) <- 1

            dd2 <- as.dist(1 - dd)
        }
        else if (grepl("(pearson|kendall|spearman)", dist_method)) {
            dd <- cor(tdata, use = "pairwise.complete.obs", method = dist_method)

            dd[which(is.na(dd))] <- 0
            diag(dd) <- 1

            dd2 <- as.dist(1 - dd)
        }
        else {
            library(vegan)

            dd <- 1 - as.matrix(vegdist(data,
                method = dist_method,
                diag = T, upper = T, binary = F
            ))
            dd[which(is.na(dd))] <- 0
            dd2 <- as.dist(1 - dd)
        }

        if (abs2) {
            dd2 <- abs(dd)
        }

        print(head(dd))

        cl <- hclust(dd2, method = method)
        l$cl2 <- cl

        if (k2 > 0) {
            clusters <- cutree(cl, k = k2)
            topcols <- as.character(set1_cols[clusters])
            l$clusters2 <- clusters
        }

        colv <- reorderfun(as.dendrogram(cl), colMeans(tdata, na.rm = T))

        if (!is.null(threshold)) {
            data[which(abs(data) < threshold)] <- 0
        }

        if (!is.null(colCol)) {
            colCol <- colCol[colnames(tdata), ] # in case data was subset
            print(colCol)
        }

        if (k2 > 0) {
            if (k > 0) {
                heatmap.2(tdata_scale,
                    dendrogram = "both", RowSideColors = sidecols, ColSideColors = topcols,
                    density.info = "none", trace = "none", Rowv = rowv, Colv = colv,
                    cexRow = cexRow, cexCol = cexCol, margins = c(9, 9), srtCol = 90,
                    col = colours, breaks = bins, symkey = F, lhei = c(1, 5), scale = "none",
                    colRow = colRow, colCol = colCol
                )
            } else {
                heatmap.2(tdata_scale,
                    dendrogram = "both", ColSideColors = topcols,
                    density.info = "none", trace = "none", Rowv = rowv, Colv = colv,
                    cexRow = cexRow, cexCol = cexCol, margins = c(9, 9), srtCol = 90,
                    col = colours, breaks = bins, symkey = F, lhei = c(1, 5), scale = "none",
                    colRow = colRow, colCol = colCol
                )
            }
        }
        else {
            if (k > 0) {
                heatmap.2(tdata_scale,
                    dendrogram = "both", RowSideColors = sidecols,
                    density.info = "none", trace = "none", Rowv = rowv, Colv = colv,
                    cexRow = cexRow, cexCol = cexCol, margins = c(9, 9), srtCol = 90,
                    col = colours, breaks = bins, symkey = F, lhei = c(1, 5), scale = "none",
                    colRow = colRow, colCol = colCol
                )
            } else {
                heatmap.2(tdata_scale,
                    dendrogram = "both",
                    density.info = "none", trace = "none", Rowv = rowv, Colv = colv,
                    cexRow = cexRow, cexCol = cexCol, margins = c(9, 9), srtCol = 90,
                    col = colours, breaks = bins, symkey = F, lhei = c(1, 5), scale = "none",
                    colRow = colRow, colCol = colCol
                )
            }
        }
    }

    dev.off()

    return(l)
}

normalize.medFC <- function(mat) {
    # Perform median fold change normalisation
    #          X - data set [Variables & Samples]
    medSam <- apply(mat, 1, median)
    medSam[which(medSam == 0)] <- 0.0001
    mat <- apply(mat, 2, function(mat, medSam) {
        medFDiSmpl <- mat / medSam
        vec <- mat / median(medFDiSmpl)
        return(vec)
    }, medSam)
    return(mat)
}

elbow_plot <- function(pca, main = NULL) {
    plot(pca@R2, ylim = c(0, 1), xaxt = "n", xlab = NA, ylab = "R2", main = main)
    lines(pca@R2)
    axis(1, at = 1:length(pca@R2), labels = paste0("PC", 1:length(pca@R2)))
}

make_transform_plots <- function(data, variables, dirname) {
    path <- make_path(dirname)

    for (varname in variables) {
        print(varname)

        # for every column (excluding a few),
        # we create a new data frame made of id and variable

        df <- na.omit(data[, varname, drop = T])

        svg(paste0(path, varname, ".svg"))

        par(mfrow = c(2, 4))

        d <- density(df)
        plot(d, main = "Normal")
        polygon(d, col = adjustcolor("red", alpha.f = 0.7))

        d <- density(sqrt(df))
        plot(d, main = "Sqrt")
        polygon(d, col = adjustcolor("red", alpha.f = 0.7))

        d <- density(log(df + 1))
        plot(d, main = "Log")
        polygon(d, col = adjustcolor("red", alpha.f = 0.7))

        d <- density(df**2)
        plot(d, main = "Square")
        polygon(d, col = adjustcolor("red", alpha.f = 0.7))

        qqnorm(df)
        qqline(df)

        qqnorm(sqrt(df))
        qqline(sqrt(df))

        qqnorm(log(df + 1))
        qqline(log(df + 1))

        qqnorm(df**2)
        qqline(df**2)

        dev.off()
    }
}

make_beeswarm <- function(prefix, dirname, output = svg, ...) {
    library(beeswarm)

    path <- make_path(dirname)

    if (identical(output, svg)) {
        ext <- ".svg"
    }
    else if (identical(output, pdf) || identical(output, cairo_pdf)) {
        ext <- ".pdf"
    }
    else if (identical(output, postscript)) {
        ext <- ".eps"
    }
    else {
        print("WARNING: output is neither svg nor pdf/cairo_pdf, so we default to png")
        output <- png
        ext <- ".png"
    }

    output(paste0(path, prefix, ext))
    beeswarm(...)
    bxplot(..., xlab = "", add = TRUE)
    dev.off()
}

make_scatter <- function(prefix, dirname, output = svg, ...) {
    path <- make_path(dirname)

    if (identical(output, svg)) {
        ext <- ".svg"
    }
    else if (identical(output, pdf) || identical(output, cairo_pdf)) {
        ext <- ".pdf"
    }
    else if (identical(output, postscript)) {
        ext <- ".eps"
    }
    else {
        print("WARNING: output is neither svg nor pdf/cairo_pdf, so we default to png")
        output <- png
        ext <- ".png"
    }

    output(paste0(path, prefix, ext))
    plot(...)
    dev.off()
}
