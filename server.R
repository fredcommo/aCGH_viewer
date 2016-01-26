###########################
###########################

source('helpers.R')
source('byGene.R')
###########################

shinyServer(function(input, output, session) {

    options(shiny.deprecation.messages=FALSE)
    
    geneDB <- reactiveValues(
    	gdb = GRanges()
    	)
    observe(({
    	geneDB$gdb <- .makeGeneDB(input$gBuild)
    	}))

    Input <- reactiveValues(
    	segTable = data.frame(),
    	geneTable = data.frame()
    	)
    observe({
    	Input$segTable <- .getData(input$file1$datapath, input$minSeg, input$gBuild)
    	Input$geneTable <- ByGene(Input$segTable, input$gBuild, geneDB$gdb)
    	})


    reCenterSeg <- reactive({
        st <- Input$segTable
        if(!is.null(st) && st != 1){
            st$seg.med <- st$seg.med + input$center
            return(st)
            }
        return(NULL)
        })

    reCenterGenes <- reactive({
        if(is.null(Input$geneTable))
            return(NULL)

        geneTable <- Input$geneTable
        geneTable$Log2Ratio <- geneTable$Log2Ratio + input$center
        return(geneTable)
        })

    mainPlot <- reactive({
        seg <- reCenterSeg()
        if(is.null(seg)){
#            .welcome()
            return(NULL)
        }
        gPlot <- .mainPlot(seg)
        gPlot
        })

    segmentPlot <- reactive({
        seg <- reCenterSeg()
        gPlot <- mainPlot()
        if(is.null(gPlot))
            return(NULL)

        gPlot <- .addSegments(gPlot, seg, input$chr, input$gain, input$loss, input$segLen, input$GLcols)
        gPlot <- .updateScale(gPlot, input$chr, input$gBuild, input$Ymin, input$Ymax)
        gPlot <- .addChr(gPlot, input$chr, input$gBuild)
        gPlot <- .addTitle(gPlot, Input$segTable$ID, input$gain, input$loss)

        gPlot
    })

    checkFile <- reactive({
            if(is.null(input$file1$datapath))
                return(0)
            return(1)
        })

    welcomeImage <- reactive({
        return(
            list(
                src="www/images/welcome.png",
                filetype="image/png",
                alt="welcome-image"
                )
            )
    })

    createCGHplot <- reactive({

        gPlot <- segmentPlot()
        if(is.null(gPlot))
            return(NULL)

        geneTable <- reCenterGenes()
        gene <- toupper(input$geneSymbol)

        if(!gene %in% c('NONE', '')){
            if(gene %in% geneTable$symbol){
            	geneAnnot <- geneTable[which(geneTable$symbol == gene),]
                gPlot <- .addTag(gPlot, geneAnnot)
            }
        }

        gPlot        
    })

    createSummary <- reactive({

        gene <- toupper(input$geneSymbol)
        if(input$geneSymbol %in% c('NONE', ''))
            return(NULL)

		geneTable <- reCenterGenes()
		if(!gene %in% geneTable$symbol){
            msg <- sprintf("'%s' may not be a valid HUGO symbol", toupper(input$geneSymbol))
            return( data.frame(Error = msg) )			
		}
		else{
			geneAnnot <- geneTable[which(geneTable$symbol == gene),]
	        chr <- as.integer(geneAnnot$chr)
	        chrStart <- as.integer(geneAnnot$chrStart)
	        chrEnd <- as.integer(geneAnnot$chrEnd)
	        geneAnnot$position <- sprintf("chr%s:%s-%s", chr, chrStart, chrEnd)
	        geneAnnot$segNum <- as.integer(geneAnnot$segNum)
	        geneAnnot$"segLength(kb)" <- as.integer(geneAnnot$"segLength(kb)")
	        geneAnnot <- geneAnnot[,c("symbol", "entrezid", "fullName",
	            "position", "segNum", "segLength(kb)", "Log2Ratio")]

	        return(geneAnnot)
			}

        })

    createFullTable <- reactive({

        if(is.null(Input$geneTable))
            return(NULL)

        geneTable <- reCenterGenes()
        geneTable <- geneTable[,c("symbol", "entrezid", "fullName",
            "chr", "cytoband", "chrStart", "chrEnd",
            "segNum", "segLength(kb)", "Log2Ratio")]
        geneTable$Log2Ratio <- round(geneTable$Log2Ratio, 2)
        geneTable$entrezid <- .renderLink(geneTable$entrezid)
        geneTable
        })

    filterGeneTable <- reactive({
        geneTable <- createFullTable()

        if(input$chr=="All"){
            chr <- 1:23
        } else{
            chr <- as.numeric(input$chr)
        }

        if(input$segLen %in% c("All", "")){
            segLen <- Inf
        } else{
            segLen <- as.numeric(input$segLen)
        }

        greater <- as.numeric(input$gain)
        lower <- (-abs(as.numeric(input$loss)))

        geneTable <- .filterBygene(geneTable, chr, greater, lower, segLen)
        return(geneTable)
        })

    createTitle1 <- reactive({
        if(is.null(Input$segTable) || Input$segTable == 1)
            return(NULL)
        return( unique(as.character(Input$segTable$ID)) )
        })

    createTitle2 <- reactive({
        hi <- input$gain
        lo <- input$loss
        return( sprintf("Gain threshold: %s, Loss threshold: %s", hi, lo) )
        })

    progressText <- reactive({
        if(is.null(input$file1$datapath)){
            giturl <- "https://github.com/fredcommo/aCGH_viewer"
            msg <- sprintf("Get example formats at: %s", giturl)
            return(msg)
            }
        return(NULL)
    })

    # outputs to ui
    output$checkFile <- renderText({ checkFile() })
    output$welcomeImage <- renderImage({ welcomeImage() }, deleteFile = FALSE)
    output$Profile <- renderPlot({ createCGHplot() }, res=120, width=1500, height=650)
    output$progress1 <- renderText( progressText() )
    output$progress2 <- renderText( progressText() )
    output$geneSummary <- renderTable({ createSummary() })
    output$tableTitle1 <- renderText({ createTitle1() })
    output$tableTitle2 <- renderText({ createTitle2() })
    output$fullTable <- renderDataTable({ filterGeneTable() },
        options=list(lengthMenu=c(25, 50, 100)), escape = FALSE)

    # Download functions
    output$downloadPlot <- downloadHandler(
        filename <- function(){
            sprintf("%s_genomic_profile.png", createTitle1())
            },
        content <- function(file) {
            png(file, width=1200, height=500)
            print( createCGHplot() )
            dev.off()
            }, contentType = "image/png"
            )

    output$downloadData <- downloadHandler(
        filename <- function(){
            sprintf("%s_geneTable_chr%s_gain%s_loss%s.xls",
            createTitle1(), input$chr, input$gain, input$loss)
            },
        content <- function(file){
            out <- createFullTable()
            out$entrezid <- gsub(".*term=|\\[uid\\].*", "",
                out$entrezid)
            write.table(out, file, sep="\t", row.names=FALSE)
        }
        )

    session$onSessionEnded(function() { stopApp() })
    })
