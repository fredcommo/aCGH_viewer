###########################
# Scale Shiny app on .io
## Instance Memory
## small  	256 MB (default)
## medium 	512 MB
## large 	1024 MB
## xlarge 	2048 MB
## xxlarge 	4096 MB
# shinyapps::configureApp("aCGH_viewer", size="medium")

library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library( org.Hs.eg.db)
require(ggplot2)
require(grid)
load("data/hg19.rda")

###########################
.getData <- function(filepath){
    
    if (!is.null(filepath)){
        dat <- try(read.delim(filepath, header=TRUE, sep="\t"))
        if(ncol(dat)<2)
            dat <- try(read.delim(filepath, header=TRUE, sep=";"))
        if(ncol(dat)<2)
            dat <- try(read.delim(filepath, header=TRUE, sep=","))
        if(ncol(dat)<2){
            cat("format not supported.\n")
            return(1)
        }
        if(!"probes.Sd" %in% colnames(dat))
            dat <- cbind.data.frame(dat, "probes.Sd" = rep(1, nrow(dat)))
        dat <- .chrToGenomeLoc(dat)
        cat("Returning dat\n")
        system(sprintf("rm %s", filepath))
        return(dat)
    }
    cat("Returning NULL\n")
    return(NULL)
}
.chrToGenomeLoc <- function(segTable){
    splitTable <- split(segTable, segTable$chrom)
    newTable <- lapply(splitTable, function(tmp){
        chr <- unique(tmp$chrom)
        tmp$loc.start <- tmp$loc.start + hg19$cumlen[chr]
        tmp$loc.end <- tmp$loc.end + hg19$cumlen[chr]
        tmp
    })
    do.call(rbind.data.frame, newTable)
}
.genometoChrLoc <- function(segTable){
    splitTable <- split(segTable, segTable$chrom)
    newTable <- lapply(splitTable, function(tmp){
        chr <- unique(tmp$chrom)
        tmp$loc.start <- tmp$loc.start - hg19$cumlen[chr]
        tmp$loc.end <- tmp$loc.end - hg19$cumlen[chr]
        tmp
    })
    do.call(rbind.data.frame, newTable)
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
        y <- rnorm(n, segTable$seg.mean[i], segTable$probes.Sd[i]/s)
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
            axis.text=element_text(size=15),
            axis.title=element_text(size=20),
            axis.title.x=element_text(vjust=-.5),
            axis.title.y=element_text(hjust=.5, vjust = .75),
            plot.title = element_text(lineheight=1.1, size = 25, vjust = 2, face="bold")
            )
        return(gPlot)
}
.updateScale <- function(gPlot, chr, hg19, Ymin, Ymax){
    if(is.null(gPlot))
        return(NULL)
    
    ymin <- (min(gPlot$data$l2r, na.rm=TRUE) - .75)*Ymin
    ymin <- max(-2.5, ymin)
    ymax <- (max(gPlot$data$l2r, na.rm=TRUE) + .75)*Ymax
    if(chr != "All"){
        chr <- as.numeric(chr)
        cumLen <- cumsum(as.numeric(hg19$length))
        xmin <- ifelse(chr==1, 0, cumLen[chr-1])
        xmax <- cumLen[chr]
        gPlot <- gPlot + 
            coord_cartesian(xlim = range(xmin, xmax),
                ylim=range(ymin, ymax)) +
            scale_y_continuous(breaks = seq(round(ymin), round(ymax), by = 0.5)) +
            geom_point(pch = 19, cex = 1, col = 'grey50')
    }
    else
        gPlot <- gPlot + 
            coord_cartesian(xlim = range(-0.5e8, hg19$cumlen[24]+0.5e8),
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

    #myBlue <- rgb(0, 0.45, 1, 1)

    segTable <- segTable[which(segTable$chrom %in% chr),]
    L <- abs(segTable$loc.end - segTable$loc.start)/1e6
    idx <- which((segTable$seg.mean<= loss | segTable$seg.mean>= gain) & 
        L <= segLen)

    if(length(idx)>0){
        subTable <- segTable[idx,]
        GLcolors <- ifelse(subTable$seg.mean<= loss, lossCol,
            ifelse(subTable$seg.mean>= gain, gainCol, "black")
            )
        gPlot <- gPlot+
            geom_segment(
                data=subTable,
                aes(x=loc.start, xend=loc.end, y=seg.mean, yend=seg.mean),
                colour=GLcolors, size=2
                )
        }

    return(gPlot)   
}
.addChr <- function(gPlot, chr, hg19){
    if(is.null(gPlot))
        return(NULL)

    if(chr=="All") chr <- as.numeric(1:23)
    else chr <- as.numeric(chr)

    ylim <- gPlot$coordinates$limits$y
    if(is.null(ylim)){
        ylim <- range(gPlot$data$l2r)
    }
    cumCentr <- 1/2*hg19$length+hg19$cumlen
    gPlot <- gPlot+
        geom_vline(xintercept = hg19$cumlen[chr], color = 'grey30',
            linetype = 2, size = 0.25) +
        annotate('text', x=cumCentr[chr], y=rep(max(ylim, na.rm=TRUE)*.95,
            length(chr)), label=chr, size = 4, colour = 'grey40')
    return(gPlot)     
}
.addTitle <- function(gPlot, sampleName, gain, loss){
    if(is.null(gPlot))
        return(NULL)

    Title = paste(unique(sampleName), '\nGain threshold: ', round(gain, 3),
        ' Loss threshold:', round(loss, 3))
    gPlot <- gPlot + ggtitle(Title)
    return(gPlot)  
}
.addTag <- function(gPlot, geneAnnot, Yexpand, gain, loss){
    if(is.null(gPlot))
        return(NULL)
    
    myBlue <- "darkblue"
    ylim <- gPlot$coordinates$limits$y
    ymin <- min(ylim); ymax <- max(ylim)
    symbol <- as.character(geneAnnot$symbol)
    x <- geneAnnot$genomeStart
    lr <- geneAnnot$Log2Ratio
    yLabel <- ifelse((lr+1)<ymax, lr+1, lr-1)

    if(is.na(lr)) return(gPlot)

    Col <- ifelse(lr<= loss, 'red3', ifelse(lr>=gain, myBlue, 'grey25'))
    gPlot <- gPlot +
        annotate("text",
            x=max(x, 2e8), y=yLabel,
            label=paste0(symbol, '\n(Log2R = ', round(lr, 3), ')'),
            cex=7, colour=Col) +
        geom_point(x=x, y=lr, cex=4, pch=19, colour="antiquewhite") +
        geom_point(x=x, y=lr, cex=2, pch=19, colour=ifelse(lr>0, "red3", "darkblue"))

    return(gPlot)
}
.geneOfInt <- function(symbol, geneTable){
    if(is.null(geneTable))
        return(NULL)
    
    symbol <- toupper(symbol)
    tmp <- geneTable[which(geneTable$symbol == symbol),]
    
    if(nrow(tmp)==0)
        return(NULL)
        
    tmp
}
.renderLink <- function(uid){
    sprintf("<a href=\"http://www.ncbi.nlm.nih.gov/gene/?term=%s[uid]\" 
    target=\"_blank\" style=\"font-size:18px; \">%s</a>", uid, uid)
}

# End helper functions
###########################
###########################
