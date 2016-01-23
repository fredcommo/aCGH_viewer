
shinyUI(
    pageWithSidebar(

        headerPanel("Interactive aCGH Viewer"),

        sidebarPanel(
        #includeCSS("shinySafir.css"),
        tags$head( tags$link(rel="stylesheet", type="text/css", href="style.css") ),

#        tags$hr(),
        h4("Choose a file", br(),
           div(id="ftype", em("(.csv, .tsv, .txt, .bz2, .gz)"))
           ),
        fileInput('file1', '', accept=c(  'text/csv',
                                              'text/comma-separated-values',
                                              'text/tab-separated-values',
                                              'text/plain',
                                              'application/x-bzip2',
                                              'application/x-gzip',
                                              '.csv',
                                              '.tsv',
                                              '.txt',
                                              '.bz2',
                                              '.gz')
        ),
        tags$hr(),

        h4("Genome build"),
        div(class = "row-fuild", style="margin-top: -20px;",
        	radioButtons("gBuild", "",
        		choices = c("hg18", "hg19", "hg38"), selected = "hg19", inline = TRUE)
        	),
        tags$hr(),
                
#        h4('Genome build'),

        h4('Gene symbol'),
        div(id="geneSymbol", style="margin-top: -20px;", textInput("geneSymbol", '', 'NONE')),
#		textInput("geneSymbol", "", "NONE"),
#        tags$hr(),

        h4('Show chromosome'),
        div(id="showChrom", style="margin-top: -20px;",
            selectInput(inputId = "chr", label = "", choices = c('All', 1:23),
                selected = 'All')),
        tags$hr(),
        
        column(width=12,
               div(class = "form-group",
                   radioButtons("GLcols", "Gain/Loss colors",
                                choices = c("blue/red", "red/blue"), inline = TRUE)
               )
        ),

        textInput("minSeg", "Merging segments shorter than (Kb)", 25),
        # sliderInput("center", "Recenter profile", min=-1.5, max=1.5, value=0,
        #     step = .1),
        sliderInput("Ymax", "Rescale max(y)", min=.1, max=1, value=1, step=.1),
        div(id="last-slider",
            sliderInput("Ymin", "Rescale min(y)", min=.1, max=1, value=1, step=.1)
            ),
        # sliderInput("gain", "Gain threshold (Log2ratio)", min=0, max=2,
        #     value=.5, step = .25),
        # div(id="last-slider",
        #     sliderInput("loss", "Loss threshold (Log2ratio)", min=-2, max=0,
        #         value=-.5, step = .25)
        # ),
        withTags(
            div(class="row",
                div(class='col-md-12', style="padding: 0px 0px;",
                    div(class="col-xs-7", style="float: left;",
                        p("Seg length(<Mb)", style="font-size: 14px; 
                            font-weight: bold;")),
                    div(class='col-xs-5', textInput("segLen", "", "All") )
                    )
                )
        ),

        tags$hr(),

        h4("Download"),
        withTags(
            div(class='row',
                div(class="col-md-12 downloads",
                    div(class='col-xs-4 col-xs-offset-1',
                        downloadButton('downloadPlot', 'Profile', class="btn btn-success")
                        ),
                    div(class='col-xs-4 col-xs-offset-2',
                        downloadButton('downloadData', 'Table', class="btn btn-success")
                        )
                    )
                )
            ),

        tags$hr(),
        withTags(
            div(
            class='row-fluid link',
            a("frederic.commo@gustaveroussy.fr",
            href="mailto:frederic.commo@gustaveroussy.fr?Subject=aCGH%20Viewer",
            target="_top")
            )
            )
        ),

        mainPanel(
            tabsetPanel(
                tabPanel("Genomic profile",
                    plotOutput("Profile", width = "100%", height = "100%"),
                    tags$hr(),
                    div(class="row",
                        div(class="col-md-12 sliders-bottom",
                            div(class="col-xs-4", id="recenter",
                                sliderInput("center", "Recenter profile",
                                    min=-1.5, max=1.5, value=0, step = .1)
                                ),
                            div(class="col-xs-4", id="loss-thres",
                                sliderInput("loss", "Loss threshold (Log2ratio)",
                                    min=-2, max=0, value=-.5, step = .25)
                                ),
                            div(class="col-xs-4", id="gain-thres",
                                sliderInput("gain", "Gain threshold (Log2ratio)",
                                    min=0, max=2, value=.5, step = .25)
                                )
                            )
                        ),
                    textOutput("progress1", inline=TRUE),
                    div(class='row-fluid', tableOutput("geneSummary") )
                    ),
                tabPanel("Genes table",
                    h4(textOutput("tableTitle1"), align="center"),
                    h4(textOutput("tableTitle2"), align="center"),
                    tags$hr(),
                    dataTableOutput("fullTable"),
                    textOutput("progress2", inline=TRUE)
                    )
                )
            )
        )
)
