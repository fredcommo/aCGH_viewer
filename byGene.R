
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library( org.Hs.eg.db)

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
geneDB  <- genes(txdb, columns=c("gene_id"))

.cmValues <- function(segTable, hg){
    cmLocs <- .locateCM(segTable, hg)
    out <- lapply(cmLocs, function(locs){
        c(segTable$seg.mean[locs[1]], segTable$seg.mean[locs[2]])
        })
    return(out)
}
.locateCM <- function(segTable, hg){

    chrs <- unique(segTable$chrom)
    cmLocs <- lapply(chrs, function(chr){

        cStart <- hg$centromerStart[hg$chrom==chr]
        cEnd <- hg$centromerEnd[hg$chrom==chr]
        tmp <- segTable[segTable$chrom==chr,]

        locStart <- which(tmp$loc.start<=cStart & cStart<=tmp$loc.end)
        if(length(locStart)==0){
            locStart <- which.min(abs(tmp$loc.end - cStart))
        }

        locEnd <- which(tmp$loc.start<=cEnd & cEnd<=tmp$loc.end)
        if(length(locEnd)==0){
            locEnd <- which.min(abs(tmp$loc.start - cEnd))
        }

        as.numeric(rownames(tmp)[c(locStart, locEnd)])
        })

    return(cmLocs)
}
.relativeLog <- function(bygene, cmValues, hg){
    relativeLog <- lapply(1:length(cmValues), function(chr){
        tmp <- bygene[bygene$chr==chr, ]
        ii <- which(tmp$chrStart < hg$centromerStart[hg$chrom==chr])
        jj <- which(tmp$chrStart > hg$centromerEnd[hg$chrom==chr])
        rl <- rep(NA, nrow(tmp))
        rl[ii] <- tmp$Log2Ratio[ii] - cmValues[[chr]][1]
        rl[jj] <- tmp$Log2Ratio[jj] - cmValues[[chr]][2]
        return(rl)
        })
    return(do.call(c, relativeLog))
}
.getGenesFromSeg <- function(chr, Start, End){
    # chr: a integer, from 1 to 24
    # Start, End: numeric. Start/End segment position (from segmentation table)

    geneDB <- geneDB

    if(chr==23) chr <- "X"
    if(chr==24) chr <- "Y"
    
    chr <- sprintf("chr%s", chr)

    ii <- which(as.vector(seqnames(geneDB)) == chr)
    jj <- intersect(ii, which(Start <= start(geneDB) & start(geneDB) <= End))
    kk <- intersect(ii, which(Start <= end(geneDB) & end(geneDB) <= End))
    idx <- unique(union(jj, kk))

    if(length(idx) == 0)
        return(NULL)

    bySymbol <- select(org.Hs.eg.db,
                        keys=geneDB$gene_id[idx],
                        keytype='ENTREZID',
                        columns=c('SYMBOL', 'GENENAME', 'MAP')
        )
    byRange <- as.data.frame(geneDB[idx])
        
    geneList <- merge(bySymbol, byRange,
                        by.x = "ENTREZID", by.y = "gene_id", all = TRUE)
    colnames(geneList) <- tolower(colnames(geneList))
    oldNames <- c("genename", "map", "seqnames", "start", "end")
    newNames <- c("fullName", "cytoband", "chr", "chrStart", "chrEnd")
    colnames(geneList)[colnames(geneList) %in% oldNames] <- newNames
        
    geneList <- geneList[order(geneList$symbol),]
    geneList$chr <- as.character(geneList$chr)
    .chrAsNum(geneList)
}
.chrAsNum <- function(geneList){
    geneList$chr <- gsub("chr", "", geneList$chr)
    geneList$chr <- gsub("X", "23", geneList$chr)
    geneList$chr <- gsub("Y", "24", geneList$chr)
    geneList$chr <- as.numeric(geneList$chr)
    geneList
}
.addGenomeLoc <- function(bygene){
    hg19 <- hg19
    ss <- split(bygene, bygene$chr)
    bygene <- lapply(ss, function(tmp){
        chr <- unique(tmp$chr)
        tmp$genomeStart <- tmp$chrStart + hg19$cumlen[chr]
        return(tmp)
    })
    bygene <- as.data.frame(do.call(rbind, bygene))
    bygene <- bygene[order(bygene$symbol),]
    rownames(bygene) <- seq_len(nrow(bygene))
    bygene
}
ByGene <- function(st){
    if(is.null(st) || st == 1)
        return(NULL)

    st <- .genometoChrLoc(st)
    hg19 <- hg19
    bygene <- lapply(seq_len(nrow(st)), function(ii){
        g <- .getGenesFromSeg(st$chrom[ii], st$loc.start[ii], st$loc.end[ii])
        if(is.null(g))
            return(NULL)

        cbind.data.frame(g,
                        Log2Ratio=st$seg.mean[ii],
                        num.mark=st$num.mark[ii], segNum=ii,
                        "segLength(kb)"=round(abs(g$chrEnd - g$chrStart)/1e3, 2)
                        )
    })
    bygene <- do.call(rbind, bygene)
    cmValues <- .cmValues(st, hg19)
    bygene$relativeLog <- .relativeLog(bygene, cmValues, hg19)
    .addGenomeLoc(bygene)
}

#.ByGene(segTable)
