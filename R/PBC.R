#' PCR bottleneck coefficient
#'
#' Calculate the PCR bottleneck coefficient as described in the ENCODE
#' guidelines.
#'
#' The PCR bottleneck coefficient (PBC) is a measure of library complexity, i.e. how skewed the
#' distribution of read counts per location is towards 1 read per location.
#' 
#' Defined in the ENCODE guidelines (https://genome.ucsc.edu/ENCODE/qualityMetrics.html) as:
#'
#' PBC = N1/Nd
#' 
#' with
#' \itemize{
#' \item{\code{N1}: Number of genomic locations to which EXACTLY one unique mapping read maps.}
#' \item{\code{Nd}: Number of genomic locations to which AT LEAST one unique mapping read maps, i.e.
#' the number of non-redundant, unique mapping reads.}
#' }
#' 
#' PBC is further described on the ENCODE Software Tools page. Provisionally, 0-0.5 is severe
#' bottlenecking, 0.5-0.8 is moderate bottlenecking, 0.8-0.9 is mild bottlenecking, while 0.9-1.0
#' is no bottlenecking. Very low values can indicate a technical problem, such as PCR bias, or a
#' biological finding, such as a very rare genomic feature. Nuclease-based assays (DNase, MNase)
#' detecting features with base-pair resolution (transcription factor footprints, positioned 
#' nucleosomes) are expected to recover the same read multiple times, resulting in a lower PBC
#' score for these assays. Note that the most complex library, random DNA, would approach 1.0, 
#' thus the very highest values can indicate technical problems with libraries. It is the practice
#' for some labs outside of ENCODE to remove redundant reads; after this has been done, the value
#' for this metric is 1.0, and this metric is not meaningful. 82\% of TF ChIP, 89\% of His ChIP, 77\%
#' of DNase, 98\% of FAIRE, and 97\% of control ENCODE datasets have no or mild bottlenecking.
#'
#' @param The path to the \code{.bam} file of a ChIP sample or a \code{GAlignments} object of the ChIP sample.
#' @return The PBC coefficient.
#' @examples
#' pbc <- PBC("IP.bam")
PBC <- function(IP) {

    require(GenomicAlignments)
    require(data.table)

    # load ChIP sample if necessary
    if (is.character(IP)) {
        if (!file.exists(IP))
            stop(paste("File", IP, "does NOT exist."))
        else
            aln <- readGAlignments(IP)
    } else if (class(IP) == "GAlignments") {
        aln <- IP
    } else {
        stop("IP must be a file path or a GAlignments object.")
    }

    # convert GAlignments object to data.table for fast aggregation
    aln <- data.table(
        strand=as.factor(BiocGenerics::as.vector(strand(aln))),
        seqnames=as.factor(BiocGenerics::as.vector(seqnames(aln))),
        pos=ifelse(strand(aln) == "+", start(aln), end(aln))
    )

    # aggregate reads by position and count them
    readsPerPosition <- aln[,list(count=.N), by=list(strand, seqnames, pos)]$count

    # PBC = positions with exactly 1 read / positions with at least 1 read
    PBC <- sum(readsPerPosition == 1) / length(readsPerPosition)

    return(PBC)
}
