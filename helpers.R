###########################
# Scale Shiny app on .io
## Instance Memory
## small  	256 MB (default)
## medium 	512 MB
## large 	1024 MB
## xlarge 	2048 MB
## xxlarge 	4096 MB
# shinyapps::configureApp("aCGH_viewer", size="medium")

message("Sourcing helpers...")
require(ggplot2)
require(grid)
load("data/hg18.rda")
load("data/hg19.rda")
load("data/hg38.rda")

options(warn=-1)
###########################
.welcome <- function(){
    msg <- "Waiting for a file to load ..."
    plot.new()
    legend("center", legend = msg, bty = "n", cex = 1.75, text.col = "blue")
    return(NULL)
}
.getData <- function(filepath, minSeg, gBuild){
	minSeg <- as.numeric(minSeg)

    if (!is.null(filepath)){
        dat <- try(read.delim(filepath, header=TRUE, sep="\t", comment.char = c("#")))
        if(ncol(dat)<2)
            dat <- try(read.delim(filepath, header=TRUE, sep=";", comment.char = c("#")))
        if(ncol(dat)<2)
            dat <- try(read.delim(filepath, header=TRUE, sep=",", comment.char = c("#")))
        if(ncol(dat)<2)
            dat <- try(read.delim(filepath, header=TRUE, sep=" ", comment.char = c("#")))
        if(ncol(dat)<2){
            cat("format not supported.\n")
            return(1)
        }
        
        # Mandatory columns
        expectedCols <- c("chrom", "loc.start", "loc.end", "seg.mean")
        if(!all(expectedCols %in% colnames(dat)))
            stop("Missing columns: ", setdiff(expectedCols, colnames(dat)))
        
        # Not a big deal if one of those is missing.
        if(!"seg.med" %in% colnames(dat))
            dat <- cbind.data.frame(dat, "seg.med" = dat$seg.mean)
        if(!"probes.Sd" %in% colnames(dat))
            dat <- cbind.data.frame(dat, "probes.Sd" = rep(.8, nrow(dat)))
        if(!"ID" %in% colnames(dat))
            dat <- cbind.data.frame(dat, "ID" = rep(NA, nrow(dat)))
        
        # Estimate num.mark from segment length, when missing
        if(!"num.mark" %in% colnames(dat)){
            d <- abs(dat$loc.end - dat$loc.start)
            nm <- ceiling(11.2 + 3.05e-04*d + rnorm(length(d), 0, 10))
            if(any(nm<8))
                nm[nm<8] <- 8
            dat <- cbind.data.frame(dat, "num.mark" = nm)
        }
        
        dat <- dat[,c("ID", "chrom", "loc.start", "loc.end", "num.mark", "seg.mean", "seg.med", "probes.Sd")]

        if(!is.na(minSeg) && minSeg > 10)
            dat <- .smoothSeg(dat, minSeg)

        HG <- switch(gBuild,
        	hg18 = hg18,
        	hg19 = hg19,
        	hg38 = hg38)

        dat <- .chrToGenomeLoc(dat, HG)
        rownames(dat) <- 1:nrow(dat)
        cat("Returning dat\n")
#        system(sprintf("rm %s", filepath))
        return(dat)
    }
    cat("Returning NULL\n")
    return(NULL)
}
.getName <- function(segTable){
    return( unique(as.character(segTable$ID)) )
}
.selectChr <- function(Choice){
    if(Choice == 'All') return(1:23)
    else return(as.numeric(Choice))
}
.mainPlot <- function(segTable, s=3){

    if(is.null(segTable))
        return(NULL)
    
    if(segTable == 1){
        plot.new()
        legend("center", legend = "File format not supported", cex = 2)
        return(NULL)
    }

    N <- sum(segTable$num.mark, na.rm=TRUE)
    w <- N/20e3

    X <- lapply(1:nrow(segTable), function(i){
        n <- ceiling(segTable$num.mark[i]/w)
        n <- max(50, n)
        x <- seq(segTable$loc.start[i], segTable$loc.end[i], len=n)
        y <- rnorm(n, segTable$seg.med[i], segTable$probes.Sd[i]/s)
        return(cbind(loc=x, l2r=y))
    })
    X <- as.data.frame(do.call(rbind, X))

    gPlot <- ggplot(data=X, aes(x=loc, y=l2r)) +
        geom_point(pch = 19, cex = 0.15, col = 'grey50') +
        geom_hline(yintercept = 0) +
        xlab('Genomic position (bp)') +
        ylab('Log2(Ratio)') +
        theme_bw() +
        theme(
            panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            plot.margin=unit(c(4,4,4,0),"mm"),
            axis.text=element_text(size=12),
            axis.title=element_text(size=20),
            axis.title.x=element_text(vjust=-.5),
            axis.title.y=element_text(hjust=.5, vjust = .75),
            plot.title = element_text(lineheight=1.1, size = 25, vjust = 2, face="bold")
            )
        return(gPlot)
}
.updateScale <- function(gPlot, chr, gBuild, Ymin, Ymax){
    if(is.null(gPlot))
        return(NULL)
    
    HG <- switch(gBuild,
        	hg18 = hg18,
        	hg19 = hg19,
        	hg38 = hg38)

    ymin <- max(-3.5, min(gPlot$data$l2r, na.rm=TRUE) - 0.5)*Ymin
    ymin <- min(-1, ymin)
    ymax <- (max(gPlot$data$l2r, na.rm=TRUE) + .75)*Ymax
    if(chr != "All"){
        chr <- as.numeric(chr)
        cumLen <- cumsum(as.numeric(HG$length))
        xmin <- ifelse(chr==1, 0, cumLen[chr-1])
        xmax <- cumLen[chr]
        gPlot <- gPlot + 
            coord_cartesian(xlim = range(xmin, xmax), ylim=range(ymin, ymax)) +
            scale_y_continuous(breaks = seq(round(ymin), round(ymax), by = 0.5)) +
            geom_point(pch = 19, cex = 1, col = 'grey50')
    }
    else
        gPlot <- gPlot + 
            coord_cartesian(xlim = range(-0.5e8, HG$cumlen[24]+0.5e8),
                ylim=range(ymin, ymax)) +
            scale_y_continuous(breaks = seq(round(ymin), round(ymax), by = 0.5))

    return(gPlot)
}
.addSegments <- function(gPlot, segTable, chr, gain, loss, segLen, GLcols){
    if(is.null(gPlot))
        return(NULL)
    
    if(chr=="All"){
        chr <- 1:23
    } else{
        chr <- as.numeric(chr)
    }

    if(segLen %in% c("All", "")){
        segLen <- Inf
    } else{
        segLen <- as.numeric(segLen)
    }

    gainCol <- rgb(0, 0.45, 1, 1)
    lossCol <- "red3"
    if(grepl("red/blue", GLcols)){
        gainCol <- "red3"
        lossCol <- rgb(0, 0.45, 1, 1)
    }

    segTable <- segTable[which(segTable$chrom %in% chr),]
    L <- abs(segTable$loc.end - segTable$loc.start)/1e6
    idx <- which((segTable$seg.med<= loss | segTable$seg.med>= gain) & 
        L <= segLen)

    if(length(idx)>0){
        subTable <- segTable[idx,]
        GLcolors <- ifelse(subTable$seg.med<= loss, lossCol,
            ifelse(subTable$seg.med>= gain, gainCol, "black")
            )
        gPlot <- gPlot+
            geom_segment(
                data=subTable,
                aes(x=loc.start, xend=loc.end, y=seg.med, yend=seg.med),
                colour=GLcolors, size=2
                )
        }

    return(gPlot)   
}
.addChr <- function(gPlot, chr, gBuild){
    if(is.null(gPlot))
        return(NULL)

    HG <- switch(gBuild,
        	hg18 = hg18,
        	hg19 = hg19,
        	hg38 = hg38)

    if(chr=="All") chr <- as.numeric(1:23)
    else chr <- as.numeric(chr)

    ylim <- gPlot$coordinates$limits$y
    if(is.null(ylim)){
        ylim <- range(gPlot$data$l2r)
    }
    cumCentr <- 1/2*HG$length + HG$cumlen
    gPlot <- gPlot+
        geom_vline(xintercept = HG$cumlen[chr], color = 'grey30',
            linetype = 2, size = 0.25) +
        annotate('text', x=cumCentr[chr], y=rep(max(ylim, na.rm=TRUE)*.95,
            length(chr)), label=chr, size = 4, colour = 'grey40')

    return(gPlot)     
}
.addTitle <- function(gPlot, sampleName, gain, loss){
    if(is.null(gPlot))
        return(NULL)

    Title <- paste(unique(sampleName), '\nGain threshold: ', round(gain, 3),
        ', Loss threshold:', round(loss, 3))
    gPlot <- gPlot + ggtitle(Title)

    return(gPlot)  
}
.addTag <- function(gPlot, geneAnnot){
    if(is.null(gPlot))
        return(NULL)
    
#    myBlue <- "darkblue"
    ylim <- gPlot$coordinates$limits$y
    ymin <- min(ylim); ymax <- max(ylim)
    symbol <- as.character(geneAnnot$symbol)
    x <- geneAnnot$genomeStart
    lr <- geneAnnot$Log2Ratio
    yLabel <- ifelse((lr+.8)<ymax, lr+.8, lr-.8)

    if(is.na(lr)) return(gPlot)

    Col <- "grey25" #ifelse(lr<= loss, 'red3', ifelse(lr>=gain, myBlue, 'grey25'))
    gPlot <- gPlot +
        annotate("text",
            x=max(x, 2e8), y=yLabel,
            label=paste0(symbol, '\n(Log2R = ', round(lr, 3), ')'),
            fontface="bold", colour=Col) +
        geom_point(x=x, y=lr, size=6, pch=19, colour="black") +
        geom_point(x=x, y=lr, size=5, pch=19, colour="antiquewhite") +
        geom_point(x=x, y=lr, size=4, pch=19, colour="darkorchid4") +
        geom_point(x=x, y=lr, size=1, pch=19, colour="cyan")
    
    return(gPlot)
}

###########################
# MERGING SEGMENTS
.getSegLen <- function(seg){
    abs(seg$loc.end - seg$loc.start)/1e3
}
.smoothSeg <- function(segTable, minSeg){
    minSeg <- as.numeric(minSeg)
    splitSegTables <- split(segTable, segTable$chrom)
    adjustedLocs <- lapply(splitSegTables, function(sst){
        if(nrow(sst)<2)
            return(sst)
        L <- .getSegLen(sst)
        while(any(L < minSeg) & nrow(sst)>1){
            i <- which(L < minSeg)[1]
            j <- .getCloser(sst, i)
            sst <- .mergeSegments(sst, i, j)
            sst <- sst[-i,]
            L <- .getSegLen(sst)
            }
        return(sst)
        })

    adjustedLocs <- as.data.frame(do.call(rbind, adjustedLocs))
    rownames(adjustedLocs) <- seq(1, nrow(adjustedLocs))

    return(adjustedLocs)
}
.getCloser <- function(segTable, idx){
    if(idx==1){
        return(idx+1)
    } else if (idx==nrow(segTable)){
        return(idx-1)
    } else {
        delta <- abs(segTable$seg.med[c(idx-1,idx+1)] - segTable$seg.med[idx])
        i <- ifelse(which.min(delta)==1, idx-1, idx+1)
        return(i)
    }
}
.mergeSegments <- function(segTable, i, j){
    if(j<i){
        segTable$loc.end[j] <- segTable$loc.end[i]
    } else {
        segTable$loc.start[j] <- segTable$loc.start[i]
    }
    segTable$num.mark[j] <- segTable$num.mark[j] + segTable$num.mark[i]
    return(segTable)
}


# End helper functions
###########################
###########################
