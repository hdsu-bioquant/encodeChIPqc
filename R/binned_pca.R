#' ChIP-seq binned tag counting tool
#'
#' \code{tagCount} returns a \code{data.frame} of tag counts per bin per sample.
#'
#' This function makes use of the \code{GenomicAlignments} package to count
#' tags on bins. The reference genome is split into bins and a
#' \code{data.frame} with the counts per bin is returned. The \code{data.frame}
#' has one column for each sample.
#'
#' @param samples List/vector of paths or list of GAlignments/BamFile/BamViews objects or GAlignmentsList or BamFileList that will be processed.
#' @param org The organism as described in the BSGenome package: Hsapiens, Mmusculus, ...
#' @param assembly The assembly name: UCSC, NCBI, TAIR, ...
#' @param version The assembly version: hg19, mm9, ...
#' @param binSize The size of the bins.
#' @param cluster A \code{snow} cluster to process samples in parallel: \code{cluster <- snow::makeCluster(cores)}.
#' @return A \code{data.frame} with counts on bins. Each column contains the counts of one sample.
#' @examples
#' counts <- tagCount(
#'     samples=c("ctl1.bam", "ctl2.bam", "chip1.bam", "chip2.bam"),
#'     org="Mmusculus", assembly="UCSC", version="mm9"
#' )
#' nCounts <- normalizeTagCount(counts)
#' topHeatmap(nCounts)
#' samplesHeatmap(nCounts)
#' pcaPlot(nCounts, as.factor(c("ctl", "ctl", "chip", "chip")))
tagCount <- function(samples, org, assembly, version,
                     binSize=1000, cluster=NULL) {

    require(GenomicAlignments)

    # check input parms
    if(!class(samples) %in% c("character", "list", "GAlignments", "GAlignmentsList", "BamFile", "BamFileList", "BamViews"))
        stop("samples must be of type character, list, GAlignments, GAlignmentsList, BamFile, BamFileList or BamViews")
    if(!is.character(org)) stop("org has to be a character vector with a valid BSgenome organism name")
    if(!is.character(assembly)) stop("assembly has to be a character vector with a valid BSgenome assembly name")
    if(!is.character(version)) stop("version has to be a character vector with a valid BSgenome version name")
    if(!is.numeric(binSize) || binSize <= 0) stop("binSize has to be a number > 0")

    # get size of genome and tile it
    require(paste("BSgenome", org, assembly, version, sep="."), character.only=T)
    bins <- tileGenome(seqinfo(get(org)), tilewidth=binSize, cut.last.tile.in.chrom=T)

    counter <- function(alignments, bins) {
        require(GenomicAlignments)
        if (class(alignments) != "GAlignments")
            alignments <- readGAlignments(alignments)
        hits <- countOverlaps(alignments, bins)
        # discard reads hitting more than one bin
        counts <- countOverlaps(bins, alignments[hits==1])
        names(counts) <- names(bins)
        return(counts)
    }
    counts <- .clusterApplyLB(cluster, samples, counter, bins)
    counts <- as.data.frame(counts)

    # restore sample names
    if (is.null(names(samples)) && class(samples) != "GAlignments")
        # use file paths, if no names given
        names(samples) <- sapply(samples, function(x) { ifelse(is.character(x), x, "") })
    colnames(counts) <- names(samples)

    return(counts)
}

#' Normalize tag counts
#'
#' Normalize tag counts produced by \code{tagCount} by total tag count.
#'
#' This function normalizes the tag counts produced by \code{tagCount}.
#' The tag counts per bin are scaled, such that each sample has as many
#' normalized tags as the sample with the most tags.
#'
#' @param counts The output of \code{tagCount}.
#' @param n Limit the output to the top \code{n} bins with the greatest variance. The default value of 500 is usually sufficient for clustering. \code{NULL} means all bins are returned.
#' @return A \code{data.frame} with normalized counts on bins. Each column contains the counts of one sample.
#' @examples
#' counts <- tagCount(
#'     samples=c("ctl1.bam", "ctl2.bam", "chip1.bam", "chip2.bam"),
#'     org="Mmusculus", assembly="UCSC", version="mm9"
#' )
#' nCounts <- normalizeTagCount(counts)
#' topHeatmap(nCounts)
#' samplesHeatmap(nCounts)
#' pcaPlot(nCounts, as.factor(c("ctl", "ctl", "chip", "chip")))
normalizeTagCount <- function(counts, n=500) {

    if(!is.numeric(n) || n <= 0) stop("n has to be a number > 0")

    # normalize bins
    s <- apply(counts, 2, sum, na.rm=T) # total sum of reads
    s <- max(s) / s # scaling factor: sum(counts[,i]) * s[i] is same for all i
    normalizedCounts <- sapply(1:length(s), function(i) { return(counts[,i] * s[i]) })

    # restore sample names
    normalizedCounts <- as.data.frame(normalizedCounts)
    colnames(normalizedCounts) <- colnames(counts)

    # calculate sd to filter for bins with highest sd
    if (!is.null(n)) {
        stddevs <- apply(normalizedCounts, 1, sd, na.rm=T)
        normalizedCounts <- tail(normalizedCounts[order(stddevs),], n)
    }

    return(normalizedCounts)
}

#' Heatmap and cluster of top variant regions
#'
#' Plot a \code{gplots} heatmap with the top variant regions.
#'
#' This function uses the \code{heatmap.2} function of \code{gplots} to assess
#' the relationship between samples.
#'
#' @param normalizedCounts The normalized \code{data.frame} with counts per bin as generated by \code{normalizeTagCount}.
#' @param ... Extra parameters sent to \code{heatmap.2}.
#' @return  The plotted cluster matrix.
#' @examples
#' counts <- tagCount(
#'     samples=c("ctl1.bam", "ctl2.bam", "chip1.bam", "chip2.bam"),
#'     org="Mmusculus", assembly="UCSC", version="mm9"
#' )
#' nCounts <- normalizeTagCount(counts)
#' topHeatmap(nCounts)
topHeatmap <- function(normalizedCounts, ...) {

    require(RColorBrewer)
    require(gplots)

    # Heatmap of the VST count table
    hmcol <- colorRampPalette(brewer.pal(9, "Oranges"))(100)
    mat <- as.matrix(normalizedCounts)
    heatmap.2(mat, col=hmcol, trace="none", margin=c(10,6), scale="row",
              dendrogram="column", labRow="", ...)

    return(mat)
}

#' Heatmap and cluster sample-sample distances
#'
#' Plot a \code{gplots} heatmap with the sample-sample distances.
#'
#' This function uses the \code{heatmap.2} function of \code{gplots} to assess
#' the relationship between samples.
#'
#' @param normalizedCounts The normalized \code{data.frame} with counts per bin as generated by \code{normalizeTagCount}.
#' @param ... Extra params sent to \code{heatmap.2}.
#' @return  The plotted distance matrix.
#' @examples
#' counts <- tagCount(
#'     samples=c("ctl1.bam", "ctl2.bam", "chip1.bam", "chip2.bam"),
#'     org="Mmusculus", assembly="UCSC", version="mm9"
#' )
#' nCounts <- normalizeTagCount(counts)
#' samplesHeatmap(nCounts)
samplesHeatmap <- function(normalizedCounts, ...) {

    require(RColorBrewer)
    require(gplots)

    # Heatmap of sample to sample distances
    hmcol <- colorRampPalette(brewer.pal(9, "Oranges"))(100)
    dists <- dist(t(normalizedCounts))
    mat <- as.matrix(dists)
    heatmap.2(mat, trace="none", col=rev(hmcol),
              margin=c(13, 13), ...)

    return(mat)
}

#' PCA plot of ChIP and input samples
#'
#' Plot the first two components of the PCA analysis on all samples.
#'
#' \code{pcaPlot} runs a PCA analysis on all samples (ChIP and input) and plots
#' the first two components of the analysis in a scatter plot. Each dot in the 
#' plot represents a sample. Samples which the analysis deems closely related
#' cluster together.
#'
#' @param normalizedCounts The normalized \code{data.frame} with counts per bin as generated by \code{normalizeTagCount}.
#' @param condition A factor containing the groups which the samples belong to. Samples are colored by condition.
#' @return The result of the PCA analysis.
#' @examples
#' counts <- tagCount(
#'     samples=c("ctl1.bam", "ctl2.bam", "chip1.bam", "chip2.bam"),
#'     org="Mmusculus", assembly="UCSC", version="mm9"
#' )
#' nCounts <- normalizeTagCount(counts)
#' pcaPlot(nCounts, as.factor(c("ctl", "ctl", "chip", "chip")))
pcaPlot <- function(normalizedCounts, condition=as.factor(colnames(counts))) {

    if(!is.factor(condition)) stop("condition has to be a factor vector assigning samples to condition")
    if(length(condition) != ncol(normalizedCounts)) stop("condition doesn't describe the same number of samples")

    require(RColorBrewer)

    # PCA plot of the samples for the top bins
    x <- prcomp(t(normalizedCounts))

    sdev <- round(100 * x$sdev / sum(x$sdev))
    plot(x=x$x[,1],y=x$x[,2],pch=16,cex=.5,
         xlab=paste0("PC1 ",as.character(sdev[1]),"%"),ylab=paste0("PC2 ",as.character(sdev[2]),"%"),
         col=brewer.pal(max(3,length(levels(condition))),"Set1")[condition])
    text(x=x$x[,1],y=x$x[,2],labels=rownames(x$x),cex=.75,
         col=brewer.pal(max(3,length(levels(condition))),"Set1")[condition])

    return(x)
}

