
library(gplots)
library(randomcoloR)

hclustfun <- function(x){hclust(x, method="complete")}
cordist <- function(x) {dd = as.dist((1-cor(t(x), use="pairwise.complete.obs"))/2, diag=T); dd[which(is.na(dd))]=0.5; dd}
abscordist <- function(x) {dd = as.dist(1-abs(cor(t(x), use="pairwise.complete.obs")), diag=T); dd[which(is.na(dd))]=0; dd}
reorderfun <- function(d, w) reorder(d, w)

r2 <- function(preds, actual) {rss <- sum((preds - actual) ^ 2, na.rm=T); tss <- sum((actual - mean(actual, na.rm=T)) ^ 2, na.rm=T); rsq <- 1 - rss/tss; return(rsq)}

transf_train_test <- function(data, variables, train)
{
	tmp = data

    data = data[variables]

    t_fun = c(identity, function(x){log(x+1)}, sqrt, function(x){x**2})
    norm_fun=function(d,tf){sapply(d, function(x) {shapiro.test(tf(x))$p.value})}
    df=data.frame(Identity=norm_fun(data[train,,drop=F], t_fun[[1]]),
                  Log=norm_fun(data[train,,drop=F], t_fun[[2]]),
                  Sqrt=norm_fun(data[train,,drop=F], t_fun[[3]]),
                  Squared=norm_fun(data[train,,drop=F], t_fun[[4]]))
    norm_best_i = max.col(replace(df, is.na(df), -Inf), ties.method="first")
    norm_best_f = t_fun[norm_best_i]
    data_t = data.frame(sapply(1:ncol(data), function(i) norm_best_f[[i]](data[,i,drop=F])))
    rownames(data_t) = rownames(data)
    colnames(data_t) = colnames(data)
    rm(data)

    tmp[variables] = data_t

    return(tmp)
}

make_path <- function(dirname)
{
    if (!is.null(dirname))
    {
        dirname = paste0(dirname, "/")

        if (!dir.exists(dirname))
        {
            dir.create(dirname, mode="0755", recursive=T)
        }
    }

    return(dirname)
}

make_heatmap_from_corrs <- function(corrs, file_prefix, abs=F, method="complete")
{
    if (is.null(hclustfun))
    {
        hclustfun <- function(x){hclust(x, method=method)}
    }

    if (abs)
    {
        colours = colorRampPalette(c("ghostwhite", "red"), space = "rgb")(100)
        bins = seq(0,1,0.01)
        corrs = abs(corrs)
        file_prefix = paste0(file_prefix, "_abs")
    }
    else
    {
        colours = colorRampPalette(c("blue", "ghostwhite", "red"), space = "rgb")(100)
        bins = seq(-1,1,0.02)
    }

    cairo_pdf(paste0(file_prefix, ".pdf"))
    heatmap.2(corrs, dendrogram="both", density.info="none", trace="none", cexRow=0.3, cexCol=1,
              margins=c(9,9), srtCol=90, col=colours, breaks=bins, symkey=F, lhei = c(1, 5), scale="none", hclustfun=hclustfun)
    dev.off()
}

make_heatmap <- function(clinical_data, dirname="./", prefix="heatmap", vars=colnames(clinical_data), rows=rownames(clinical_data), threshold=0.3, k=4, dist_method="pearson", abs=F, abs2=F, method="complete", square=T, k2=4)
{
    l = list()

    path = make_path(dirname)

    print(paste(path, prefix))

    data = clinical_data[rows, vars]

    # print(data)

    colours = colorRampPalette(c("blue", "ghostwhite", "red"), space = "rgb")(100)
    bins = seq(-1,1,0.02)

    set1_cols = brewer.pal(9, "Set1")

    svg(paste0(path, prefix, ".svg"))

    # first get the correlation value, for colouring only (so not absolute)
    if (dist_method == "longitudinal")
    {
        dd = longitudinal_dist(clinical_data[rows, ], vars, abs=F, dist=F)
    }
    else
    {
        dd = cor(data, use="pairwise.complete.obs", method=dist_method)
    }

    diag(dd) = 1 # sometimes the above utterly fails, putting NA in the diagonal which is absurdly wrong
    dd[which(is.na(dd))] = 0

    # now we must make a distance out of it
    if (abs)
    {
        dd2 = as.dist(1 - abs(dd))
    }
    else
    {
        dd2 = as.dist(1 - dd)
    }

    cl = hclust(dd2, method=method)
    clusters = cutree(cl, k=k)
    sidecols = as.character(set1_cols[clusters])
    rowv = reorderfun(as.dendrogram(cl), colMeans(data, na.rm=T))

    l$clusters = clusters

    if (square)
    {
        # we want "co-correlations", eg row-vs-row

        if (!is.null(threshold))
        {
            dd[which(abs(dd) < threshold)] = 0
        }

        colv = rowv

        heatmap.2(dd, dendrogram="both", RowSideColors=sidecols, density.info="none", trace="none", Rowv=rowv,
                  cexRow=1, cexCol=1, margins=c(9,9), srtCol=90, col=colours, breaks=bins, symkey=F, lhei = c(1, 5), scale="none")
    }
    else
    {
        # rows-vs-columns in the heatmap (both with dendrograms)

        tdata = t(data)

        dd = cor(tdata, use="pairwise.complete.obs", method="pearson")
        diag(dd) = 1 # sometimes the above utterly fails, putting NA in the diagonal which is absurdly wrong
        dd[which(is.na(dd))] = 0

        if (abs2)
        {
            dd2 = as.dist(1 - abs(dd))
        }
        else
        {
            dd2 = as.dist(1 - dd)
        }

        cl = hclust(dd2, method=method)
        clusters = cutree(cl, k=k2)
        topcols = as.character(set1_cols[clusters])
        colv = reorderfun(as.dendrogram(cl), colMeans(tdata, na.rm=T))

        l$clusters2 = clusters

        heatmap.2(tdata, dendrogram="both", RowSideColors=sidecols, ColSideColors=topcols, density.info="none", trace="none", Rowv=rowv, Colv=colv,
                  cexRow=1, cexCol=1, margins=c(9,9), srtCol=90, col=colours, breaks=bins, symkey=F, lhei = c(1, 5), scale="none")
    }

    dev.off()

    return(l)
}

normalize.medFC <- function(mat)
{
	# Perform median fold change normalisation
	#          X - data set [Variables & Samples]
	medSam <- apply(mat, 1, median)
	medSam[which(medSam==0)] <- 0.0001
	mat <- apply(mat, 2, function(mat, medSam){
		medFDiSmpl <- mat/medSam
		vec<-mat/median(medFDiSmpl)
		return(vec)
	}, medSam)
	return (mat)
}

elbow_plot <- function(pca, main=NULL)
{
    plot(pca@R2, ylim=c(0,1), xaxt="n", xlab=NA, ylab="R2", main=main)
    lines(pca@R2)
    axis(1, at=1:length(pca@R2), labels=paste0("PC",1:length(pca@R2)))
}

make_transform_plots <- function(data, variables, dirname)
{
    path = make_path(dirname)

    for (varname in variables)
    {
        print(varname)

        # for every column (excluding a few),
        # we create a new data frame made of id and variable

        df = na.omit(data[,varname,drop=T])

        svg(paste(path, paste0(varname, ".svg"), sep = "/"))

        par(mfrow=c(2, 4))

        d = density(df)
        plot(d, main="Normal")
        polygon(d, col=adjustcolor("red", alpha.f=0.7))

        d = density(sqrt(df))
        plot(d, main="Sqrt")
        polygon(d, col=adjustcolor("red", alpha.f=0.7))

        d = density(log(df))
        plot(d, main="Log")
        polygon(d, col=adjustcolor("red", alpha.f=0.7))

        d = density(df**2)
        plot(d, main="Square")
        polygon(d, col=adjustcolor("red", alpha.f=0.7))

        qqnorm(df)
        qqline(df)

        qqnorm(sqrt(df))
        qqline(sqrt(df))

        qqnorm(log(df))
        qqline(log(df))

        qqnorm(df**2)
        qqline(df**2)

        dev.off()
    }
}
