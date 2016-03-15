#' Plot IP strength
#'
#' Plot IP strength in the way described in the CHANCE paper.
#'
#' Plot IP strength in the way described in the CHANCE paper:
#'
#' Ref: Diaz, A et al. Normalization, bias correction, and peak calling for ChIP-seq.
#'      Stat Appl Genet Mol Biol. 2012 March 31
#'
#' The process follows these steps:
#' \enumerate{
#' \item Tile the genome into bins and count the tags on each bin.
#' \item Sort bins of IP sample and sort input in the same order.
#' \item Calculate for IP and input the cumulative sum of reads.
#' \item Calculate for IP and input the cumulative percentages.
#' \item Calculate the alpha scaling factor (ratio of cumulative sums between
#'       IP and input at the bin with the highest difference).
#' \item Plot IP strength.
#' }
#'
#' @param IP tag counts of ChIP sample as produced by \code{tagCount}.
#' @param input tag counts of input sample as produced by \code{tagCount}.
#' @return The two metrics, alpha and enrichment.
#' @examples
#' counts <- tagCount(
#'     samples=c("IP.bam", "input.bam"),
#'     org="Hsapiens", assembly="UCSC", version="hg19"
#' )
#' IPstrength(counts["IP.bam"], counts["input.bam"])
IPstrength <- function(IP, input) {

    ##   1-Is done by tagCount()

    ##   2-Sort bins in IP, and sort input in the same order as IP
    o <- order(IP)
    IP <- IP[o,,drop=F]    # sort IP
    input <- input[o,,drop=F]    # sort input

    ##   3-Calculate for IP and input the cumulative sum of reads
    cs <- lapply(c(IP, input), cumsum)

    ##   4-Calculate for IP and input the cumulative percentages
    csp <- lapply(cs, function(x) x / x[length(x)])

    ##   5-Calculate the alpha scaling factor (max diff between he cumulative percentages)
    k <- which.max(abs(csp[[2]] - csp[[1]]))
    alpha <- cs[[1]][k] / cs[[2]][k]    # scaling factor: ratio between the cumulative sums at k
    enrichment <- 1 - (k / nrow(IP))

    ##   6-Plot IPstrength
    x <- 1:nrow(IP)
    x <- x / x[length(x)]

    plot(x=x, y=csp[[1]], col="blue", xlab="% of bins", ylab="% of tags", type='l', lwd=6,
         main=paste0("alpha=",round(alpha,4), ", enrichment=", round(enrichment*100,2), "%")) # IP
    lines(x=x, y=csp[[2]], col="red", lwd=6) # input
    abline(v=x[k], col="green", lty=2, lwd=2) # location of the multiplicative scaling factor alpha

    legend = c(
        ifelse(is.null(colnames(input)), "input", colnames(input)),
        ifelse(is.null(colnames(IP)), "IP", colnames(IP))
    )
    legend("topleft", legend=legend, fill=c("red","blue"))

    return(c(alpha=alpha, enrichment=enrichment))
}
