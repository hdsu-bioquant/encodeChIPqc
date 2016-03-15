#' Convert GAlignments object to SPP tags
#'
#' \code{gAlignmentsToSppTags} converts a \code{GAlignments} object to
#' \code{spp} tags, for the purpose of passing \code{GAlignments} objects
#' to functions of the \code{spp} package.
#'
#' \code{spp} uses a custom format to represent tags (see \code{read.bam.tags}).
#' This helper function converts \code{GAlignments} objects to this custom
#' format. The \code{GAlignments} class is more versatile than the custom tag
#' format of \code{spp}. If several metrics shall be calculated on the same BAM
#' file, then it makes sense to load it only once into memory as a
#' \code{GAlignments} object and convert it to \code{spp} tags for
#' \code{spp}-related metrics, instead of loading the BAM file again just to
#' have it in \code{spp} tag format. Note that the \code{GAlignments} object
#' should contain the \code{NM} flag (number of mismatches). Otherwise, the
#' \code{quality} attribute of \code{spp} tags cannot be populated with values.
#' The \code{spp} function \code{select.informative.tags} simply accepts all
#' tags, when the \code{quality} attribute is missing. Create the
#' \code{GAlignments} object using
#'
#' \code{readGAlignments(..., param=ScanBamParam(tag="NM"))}
#'
#' to include the number of mismatches in the \code{GAlignments} object.
#'
#' @param reads \code{GAlignments} object to converted to \code{spp} tags.
#' @return \code{spp} tags which are compatible with the output of \code{read.bam.tags}
#' @examples
#' # the following three lines essentially perform the same as read.bam.tags()
#' # but a GAlignments object is more versatile and serves as input to many
#' # functions
#' library(GenomicAlignments)
#' chipSample <- readGAlignments("IP.bam", param=ScanBamParam(tag="NM"))
#' sppTags <- gAlignmentsToSppTags(chipSample)
gAlignmentsToSppTags <- function(reads) {
    NMTagPresent <- !is.null(values(reads)$NM)
    if (!NMTagPresent)
        warning("Missing NM tag. The GAlignments object should be created using readGAlignments(..., param=ScanBamParam(tag=\"NM\")).")

    # convert GAlignments object to data.table for faster conversion to SPP tag format
    require(data.table)
    readsDT <- data.table(
        seqnames=as.factor(BiocGenerics::as.vector(seqnames(reads))),
        pos=as.integer(ifelse(strand(reads) == "-", -end(reads), start(reads)-1))
    )

    # SPP stores reads by chromosome => aggregate by chromosome
    if (NMTagPresent) {
        readsDT[,NM:=values(reads)$NM] # add NM column if given
        readsByChromosome <- readsDT[,list(pos=list(pos), NM=list(NM)), by=list(seqnames)]
        sppTags <- list(tags=readsByChromosome$pos, quality=readsByChromosome$NM)
    } else {
        readsByChromosome <- readsDT[,list(pos=list(pos)), by=list(seqnames)]
        sppTags <- list(tags=readsByChromosome$pos)
    }

    # restore chromosome names
    for (i in 1:length(sppTags))
        names(sppTags[[i]]) <- readsByChromosome$seqnames

    return(sppTags)
}

#' Get binding characteristics
#'
#' Get binding characteristics from cross-correlation profile using the
#' \code{spp} package.
#'
#' Types of analysis:
#' \itemize{
#' \item Calculate binding characteristics.
#' \item Report QC metrics, RSC and NSC.
#' \item Calculate genome-wide profiles of smoothed tag density.
#' \item Calculate genome-wide profiles providing conservative statistical estimates of fold enrichment.
#' \item Determine statistically significant point binding positions.
#' \item Assessing saturation properties.
#' }
#'
#' @param IP File name/GAlignments object/\code{spp} tags of IP sample.
#' @param input File name/GAlignments object/\code{spp} tags of input sample.
#' @param read.len Read length. If \code{NULL}, the read length is detected automatically by using the median of the first 10,000 reads in the first ChIP sample.
#' @param min.shift Minimum shift range.
#' @param max.shift Maximum shift range.
#' @param bin.size Bin tags within the specified number of basepairs to speed up calculation.
#' @param cluster A \code{snow} cluster: \code{cluster <- snow::makeCluster(cores)}.
#' @return A list with the ChIP data, the input data and the binding characteristics.
#' @examples
#' # call SPP to get the binding characteristics and plot it using phantomPeak()
#' bc <- calculateBindingCharacteristics("IP.bam", "input.bam", read.len=51)
#' pp <- phantomPeak(bc$chip.data, bc$input.data, bc$bc)
calculateBindingCharacteristics <- function(IP, input, read.len=NULL, min.shift=-500, max.shift=1500, bin.size=5, cluster=NULL) {
    require(spp)

    ##
    ## Loading tag data, selecting choosing alignment quality, removing anomalies
    ##
    toSppTags <- function(sample) {
        if (is.character(sample)) # it is a file name
            return(read.bam.tags(sample))
        else if (class(sample) == "GAlignments")
            return(gAlignmentsToSppTags(sample))
        else if (class(sample) == "list") # check if sample is of type "SPP tags"
            if (!is.null(sample$tags))
                return(sample)
        # we only get here, if the format is unrecognized
        stop("Unrecognized format. IP and input must be a file name or a GAlignments object or SPP tags generated with read.bam.tags().")
    }
    chip.data  <- toSppTags(IP)
    input.data  <- toSppTags(input)

    if (is.null(read.len)) {
                # if no read length is given, automatically determine it
                # by taking the median of the first 10000 alignments
        if (is.character(IP)) {
                    require(Rsamtools)
                    reads <- scanBam(BamFile(IP, yieldSize=10000), param=ScanBamParam(what="seq"))
                    read.len <- round(median(width(reads[[1]]$seq)))
        } else if (class(IP) == "GAlignments") {
            read.len <- round(median(qwidth(head(IP, 10000))))
        } else {
            stop("read.len must be specified, if IP is not a file path or GAlignments object.")
        }
        }

    # get binding info from cross-correlation profile
    # srange gives the possible range for the size of the protected region;
    # srange should be higher than tag length; making the upper boundary too high will increase calculation time
    #
    # bin - bin tags within the specified number of basepairs to speed up calculation;
    # increasing bin size decreases the accuracy of the determined parameters
    #
    # IMPORTANT: add remove.tag.anomalies=F if duplicates were already removed
    #
    binding.characteristics <- get.binding.characteristics(chip.data,srange=c(min.shift,max.shift),bin=bin.size,remove.tag.anomalies=F,cluster=cluster)

    # select informative tags based on the binding characteristics
    chip.data  <- select.informative.tags(chip.data,binding.characteristics)
    input.data <- select.informative.tags(input.data,binding.characteristics)

    # restrict or remove singular positions with very high tag counts
    chip.data  <- remove.local.tag.anomalies(chip.data)
    input.data <- remove.local.tag.anomalies(input.data)

    ## Detect the peak at the fragment length
    cc    <- binding.characteristics$cross.correlation

    # detect all the candidate peaks by comparing cc$y[i] to cc$y[i+/-bw]
    bw    <- ceiling(2/bin.size)    # adjust if bin size was smaller than 2
    slope <- as.numeric(cc$y [(1+bw):length(cc$y )] - cc$y [1:(length(cc$y) -bw)] < 0) # 0:pos, 1:neg
    # 1: a peak (from left to right, shifting from pos to neg slope), -1: a valley shifting from neg to pos slope
    shift <- as.numeric(slope[(1+bw):length(slope)] - slope[1:(length(slope)-bw)]) 
    peaks <- which(shift == 1) + bw    # phantompeak tools look for peak with value -1!!

    # Remove fake peaks from 10:read.len+10 bp
    # Correct the peak detected by get.binding.characteristics putting to one found discarding the area around read.len
    # If the highest peak is within the discarded area, it means a problematic IP, and the max peak doesnt correspond to the fragment size
    peaks <- peaks[(cc$x[peaks] < 10) | cc$x[peaks] > read.len+10]
    binding.characteristics$peak$x <- cc$x[peaks[which.max(cc$y[peaks])]]
    binding.characteristics$peak$y <- cc$y[peaks[which.max(cc$y[peaks])]]

    ## Detect the peak at the read length (phantom peak)
    peaks <- which((cc$x >= (read.len - round(2*bin.size))) &
                   (cc$x <= (read.len + round(2*bin.size))))
    binding.characteristics$phantompeak$x <- cc$x[peaks[which.max(cc$y[peaks])]]
    binding.characteristics$phantompeak$y <- cc$y[peaks[which.max(cc$y[peaks])]]

    ## Detect the minimum correlation within the window, which will be the background cross-correlation
    binding.characteristics$back.cc$x <- cc$x[which.min(cc$y)]
    binding.characteristics$back.cc$y <- cc$y[which.min(cc$y)]

    return(list(chip.data=chip.data,input.data=input.data,bc=binding.characteristics))
}

#' QC metrics RSC and NSC
#'
#' Report the QC metrics, RSC and NSC, and plot the cross-correlation profile.
#'
#' Calculate RSC and NSC metrics and plot the cross-correlation profile.
#' \itemize{
#' \item The relative strand correlation (RSC) is the ratio between the
#' peak of the fragment length and the peak of the read length. ENCODE
#' guidelines recommend that the RSC be > 0.8.
#' \item The normalized strand coefficient, NSC, is the normalized ratio
#' between the peak of the fragment-length cross-correlation and the
#' background cross-correlation. ENCODE guidelines recommend that the NSC be
#' > 1.05.
#' }
#'
#' @param chip.data \code{spp} tags of the ChIP sample calculated with \code{calculateBindingCharacteristics}.
#' @param input.data \code{spp} tags of the input sample calculated \code{calculateBindingCharacteristics}.
#' @param binding.characteristics The binding characteristics calculated with \code{calculateBindingCharacteristics}.
#' @param smoothing.bandwidth Smooth the curve of the cross-correlation plot using a moving average (higher value means more smoothing).
#' @return A vector containing the NSC, RSC, estimated read length and estimated fragment length.
#' @examples
#' # call SPP to get the binding characteristics
#' bc <- calculateBindingCharacteristics("IP.bam", "input.bam", read.len=51)
#' # plot the cross-correlation profile
#' pp <- phantomPeak(bc$chip.data, bc$input.data, bc$bc)
phantomPeak <- function(chip.data,input.data,binding.characteristics,smoothing.bandwidth=5) {

    ## plot cross-correlation profile
    NSC <- (binding.characteristics$peak$y / binding.characteristics$back.cc$y)
    RSC <- (binding.characteristics$peak$y        - binding.characteristics$back.cc$y) /
           (binding.characteristics$phantompeak$y - binding.characteristics$back.cc$y)

    cc <- binding.characteristics$cross.correlation
    plot(binding.characteristics$cross.correlation$x,caTools::runmean(binding.characteristics$cross.correlation$y,smoothing.bandwidth,alg="C"),
         type='l',xlab="strand shift",ylab="cross-correlation",
         main=paste("NSC",round(NSC,2),"RSC",round(RSC,2),
                    "Read",binding.characteristics$phantompeak$x,"bp",
                    "Fragment",binding.characteristics$peak$x,"bp"))
    lines(x=c(binding.characteristics$phantompeak$x,binding.characteristics$phantompeak$x,min(cc$x)),
          y=c(min(cc$y),binding.characteristics$phantompeak$y,binding.characteristics$phantompeak$y),
          lty=2,col="green")
    lines(x=c(binding.characteristics$peak$x,binding.characteristics$peak$x,min(cc$x)),
          y=c(min(cc$y),binding.characteristics$peak$y,binding.characteristics$peak$y),
          lty=2,col="blue")
    lines(x=c(binding.characteristics$back.cc$x,binding.characteristics$back.cc$x,min(cc$x)),
          y=c(min(cc$y),binding.characteristics$back.cc$y,binding.characteristics$back.cc$y),
          lty=2,col="red")
    legend("topright",legend=c("phantom peak at read length","peak at fragment length",
                               "background cross-correlation"),fill=c("green","blue","red"))

    ## return
    return(c(NSC=NSC,RSC=RSC,
             Read_len=binding.characteristics$phantompeak$x,
             Fragment_len=binding.characteristics$peak$x))
}

#' Assess saturation properties
#'
#' Determine the minimal saturated enrichment ratio (MSER) and the predicted
#' sequencing depth (PSD).
#'
#' Determine the minimal saturated enrichment ratio (MSER). The saturation 
#' criteria here is 99\% consistency in the set of binding positions when adding
#' 10,000 tags.
#'
#' Determine the predicted sequencing depth.
#'
#' @param chip.data \code{spp} tags of the ChIP sample calculated with \code{calculateBindingCharacteristics}.
#' @param input.data \code{spp} tags of the input sample calculated \code{calculateBindingCharacteristics}.
#' @param binding.characteristics The binding characteristics calculated with \code{calculateBindingCharacteristics}.
#' @param step.size The amount of tags to add in every step to reach the saturation criteria.
#' @param test.agreement Saturation criteria consistency in the set of binding positions when adding \code{step.size} tags.
#' @param n.chains The number of subsampled chains.
#' @param n.steps The number of steps to convergence.
#' @param fold.enrichment The target fold enrichment.
#' @param fdr FDR cutoff.
#' @param cluster A snow cluster: \code{cluster <- snow::makeCluster(cores)}.
#' @return A vector containing the minimal saturated enrichment ratio (MSER) and the predicted sequencing depth (PSD) in million tags.
#' @examples
#' # call SPP to get the binding characteristics
#' bc <- calculateBindingCharacteristics("IP.bam", "input.bam", read.len=51)
#' # calculate MSER and PSD
#' minser <- mser(bc$chip.data, bc$input.data, bc$bc)
mser <- function(chip.data,input.data,binding.characteristics,step.size=1e5,
                 test.agreement=.99,n.chains=8,n.steps=6,fold.enrichment=2,
                 fdr=.01,cluster=NULL) {

    # determine Minimal saturated enrichment ratio (MSER)
    # note: this will take approximately 10-15x the amount of time the initial binding detection did
    # The saturation criteria here is 99% consistency in the set of binding positions when adding 1e5 tags.
    # To ensure convergence the number of subsampled chains (n.chains) should be higher (80)
    mser <- get.mser(
        chip.data,input.data,
        step.size=step.size,
        test.agreement=test.agreement,
        n.chains=n.chains,
        fdr=fdr,
        method=tag.wtd,
        whs=binding.characteristics$whs,
        cluster=cluster
    )

    # Determine the predicted sequencing depth
    # Interpolate MSER dependency on tag count:
    # Note: this requires considerably more calculations than the previous steps
    # (~ 3x more than the first MSER calculation)
    # Here we interpolate MSER dependency to determine a point at which MSER of 2 is reached
    # The interpolation will be based on the difference in MSER at the current depth, and a
    # depth at 5e5 fewer tags (n.steps=6); evaluation of the intermediate points is omitted
    # here to speed up the calculation (excluded.steps parameter)
    # A total of 7 chains is used here to speed up calculation, whereas a higher number of
    # chains (50) would give good convergence
    msers <- get.mser.interpolation(
        chip.data,input.data,
        step.size=step.size,
        test.agreement=test.agreement,
        target.fold.enrichment=fold.enrichment,
        n.chains=(n.chains-1),
        n.steps=n.steps,
        fdr=fdr,method=tag.wtd,
        whs=binding.characteristics$whs,
        cluster=cluster
    )
    psd <- round(unlist(lapply(msers,function(x) x$prediction)) / 1e6,5)

    ## return
    return(setNames(c(mser, psd), c("mser", "psd")))
}

#' Call peaks using SPP
#'
#' Call peaks using the \code{spp} package.
#'
#' \code{callPeaks()} uses the \code{spp} package to identify regions enriched
#' with reads (peaks) when comparing a ChIP sample to a control sample.
#'
#' @param chip.data \code{spp} tags of the ChIP sample calculated with \code{calculateBindingCharacteristics}.
#' @param input.data \code{spp} tags of the input sample calculated with \code{calculateBindingCharacteristics}.
#' @param binding.characteristics The binding characteristics calculated with \code{calculateBindingCharacteristics}.
#' @param fdr Peaks with a corrected p-value above this threshold are not reported.
#' @param method The peak calling method, either \code{tag.wtd} or \code{tag.lwcc}.
#' @param window.size The window size to be used in calculating enrichment.
#' @param z.thr A z-score corresponding to the Poisson ratio threshold used to flag significantly enriched windows.
#' @param cluster A snow cluster: \code{cluster <- snow::makeCluster(cores)}.
#' @return \code{spp} peak data as generated by \code{find.binding.positions} and \code{add.broad.peak.regions}.
#' @examples
#' bc <- calculateBindingCharacteristics("IP.bam", "input.bam", read.len=51)
#' peaks <- callPeaks(bc$chip.data, bc$input.data, bc$bc)
callPeaks <- function(chip.data, input.data, binding.characteristics, fdr=0.01, method=tag.wtd, window.size=1000, z.thr=3, cluster=NULL) {
    require(spp)
    peaks <- find.binding.positions(
        signal.data=chip.data,
        control.data=input.data,
        fdr=fdr,
        method=method,
        whs=binding.characteristics$whs,
        cluster=cluster
    )
    peaks <- add.broad.peak.regions(
        chip.tags=chip.data,
        input.tags=input.data,
        bp=peaks,
        window.size=window.size,
        z.thr=z.thr
    )
    return(peaks)
}

#' Convert output from SPP to GRanges object
#'
#' Convert the peaks called by \code{spp} to a \code{GRanges} object.
#'
#' \code{sppPeaksToGRanges()} takes the output of
#' \code{find.binding.positions()} from the \code{spp} package and converts it
#' into a \code{GRanges} object.
#'
#' @param peaks Output of \code{find.binding.positions()} from the \code{spp} package.
#' @return \code{spp} data converted to \code{GRanges} object.
#' @examples
#' bc <- calculateBindingCharacteristics("IP.bam", "input.bam")
#' peaks <- callPeaks(bc$chip.data, bc$input.data, bc$bc)
#' sppPeaksToGRanges(peaks)
sppPeaksToGRanges <- function(peaks) {
    require(GenomicRanges)

    # extract chromosome names, start positions and peak widths from SPP peak data
    seqname <- unlist(sapply(
        names(peaks$npl),
        function(x, npl) { rep(x, nrow(npl[[x]])) },
    peaks$npl
   ))
    start <- unlist(sapply(
        names(peaks$npl),
        function(chr, npl, whs) {
            ifelse(
                is.na(npl[[chr]]$rs),
                npl[[chr]]$x-whs/2,
                npl[[chr]]$rs
            )
        },
    peaks$npl, peaks$whs
    ))
    width <- unlist(sapply(
        names(peaks$npl),
        function(chr, npl, whs, start) {
            ifelse(
                is.na(npl[[chr]]$re),
                whs,
                npl[[chr]]$re-npl[[chr]]$rs
            )
        },
    peaks$npl, peaks$whs, start
    ))

    # convert to GRanges object
    features <- GRanges(
        seqnames <- Rle(seqname),
        strand=strand(rep("*", length(start))),
        ranges=IRanges(start=start, width=width)
    )
    return(features)
}

#' Fraction of reads in peaks
#'
#' Calculate the FRiP metric (fraction of reads in peaks).
#'
#' \code{frip} counts the number of mapped reads that fall into regions
#' which were identified as binding positions (peaks) in a ChIP-Seq experiment.
#' It returns the ratio [mapped reads in peaks] / [all mapped reads]. High
#' values are an indicator of a good signal-to-noise ratio.
#'
#' Note that the method used for peak identification influences the FRiP value.
#' A comparison of FRiP values between experiments is only possible, if the
#' same method was used for identifying peaks (i.e., same peak-caller with
#' identical parameters).
#'
#' Ref: Stephen G. Landt, Georgi K. Marinov, Anshul Kundaje, Pouya Kheradpour,
#' Florencia Pauli, Serafim Batzoglou, Bradley E. Bernstein, Peter Bickel,
#' James B. Brown, Philip Cayting, Yiwen Chen, Gilberto DeSalvo,
#' Charles Epstein, Katherine I. Fisher-Aylor, Ghia Euskirchen, Mark Gerstein,
#' Jason Gertz, Alexander J. Hartemink, Michael M. Hoffman,
#' Vishwanath R. Iyer, Youngsook L. Jung, Subhradip Karmakar, Manolis Kellis,
#' Peter V. Kharchenko, Qunhua Li, Tao Liu, X. Shirley Liu, Lijia Ma,
#' Aleksandar Milosavljevic, Richard M. Myers, Peter J. Park,
#' Michael J. Pazin, Marc D. Perry, Debasish Raha, Timothy E. Reddy,
#' Joel Rozowsky, Noam Shoresh, Arend Sidow, Matthew Slattery,
#' John A. Stamatoyannopoulos, Michael Y. Tolstorukov, Kevin P. White,
#' Simon Xi, Peggy J. Farnham, Jason D. Lieb, Barbara J. Wold, Michael Snyder:
#' ChIP-seq guidelines and practices of the ENCODE and modENCODE consortia
#' Genome Res. Sep 2012; 22(9): 1813–1831.
#'
#' @param reads File path, \code{GRanges}, \code{GAlignments}, \code{GAlignmentPairs}, \code{BamViews} or \code{BamFile} object containing reads from ChIP sample.
#' @param peaks \code{GRanges} object indicating positions and widths of peaks or \code{spp} peaks as output by \code{find.binding.positions} or \code{add.broad.peak.regions}.
#' @param singleEnd Only applicable when \code{reads} is a \code{BamViews} or \code{BamFile} object; set to \code{TRUE} for single-end data, \code{FALSE} for paired-end data.
#' @return Fraction of reads in peaks.
#' @examples
#' library(rtracklayer)
#' frip("chip.bam", import.bed("peaks.bed"))
frip <- function(reads, peaks, singleEnd=T) {
    require(GenomicAlignments)

    if (is.character(reads)) {
        require(Rsamtools)
        reads <- BamFile(reads)
    }
    
    # check if peaks are from SPP; if so, convert to GRanges
    if (class(peaks) != "GRanges") # if not GRanges, assume SPP peak format
        peaks <- sppPeaksToGRanges(peaks)

    # find reads in peaks
    overlaps <- summarizeOverlaps(
            peaks,
            reads,
            mode="IntersectionNotEmpty",
            ignore.strand=T,
            singleEnd=singleEnd,
            count.mapped.reads=T
    )

    # sum up all reads in peaks and divide by all mapped reads to get FRiP
    readsInPeaks <- sum(assay(overlaps))
    if (class(reads) %in% c("BamViews", "BamFile")) {
        allReads <- colData(overlaps)$mapped
    } else {
        allReads <- colData(overlaps)$records
    }
    result <- readsInPeaks/allReads

    return(result)
}

