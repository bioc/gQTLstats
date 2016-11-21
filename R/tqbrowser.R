band2feats = function(cbstruct, bandid, gr, featExtractor = function(x)names(x)) {
 # cbstruct is a cytoband GRanges
 # bandid is the name of the cytoband to use
 # gr is the GRanges instance from which feature ids will be drawn
 # featExtractor is a function on subsetByOverlaps(gr, cbstruct[bandid])
 #  yielding a character vector of features retained
 featExtractor( subsetByOverlaps(gr, cbstruct[bandid]) )
 }

tqbrowser = function( mae, felname, gelname, tiling, tsbra, annovec,
  band.init="6q12",  ...) {
# 
# browse ranked trans-associations in tiles
#
# mae has feature data, felname picks feature component,
# gelname picks genotypes component, tiling is a GRanges,
# tsbra is tsByRankAccum instance  
#
# genotypes will come from a VCF stack in mae
# element names must be compatible with tiling
# server starts by checking this and selecting the
# associated component
#
# first, verify that the MAE supplied includes a VcfStack
#
 stopifnot(inherits(experiments(mae)[[gelname]], "VcfStack"))
#
 requireNamespace("shiny")
 ui = fluidPage(
    titlePanel("cytoband chooser"),
    sidebarPanel(
     selectInput("curband", "cytoband", choices=names(tiling),
        selected=band.init),
     uiOutput("snpSelector"),
     numericInput("rank", "rank", 1, min=1, max=5, step=1)
     ),
    tabsetPanel(
     tabPanel("Manh.", plotOutput("manh")),
     tabPanel("y vs GT", plotOutput("eqbox"))
     )
    )
#   
 server = function(input, output) {
   output$selband = renderTable( input$curband )
   output$snpSelector = renderUI({
     tagList(
       selectInput("cursnp", "SNP", choices=
         band2feats(tiling, input$curband, tsbra))
       )
      })
   output$manh = renderPlot({
      x = band2feats(tiling, input$curband, tsbra, function(x) start(x))
      y = band2feats(tiling, input$curband, tsbra, function(x) 
            x$allscores[,input$rank])
      plot(x,y,xlab=input$cursnp)
      })
   output$eqbox = renderPlot({
     curchrn = sub("[pq].*", "", input$curband) # character
#     print(curchrn)
#     print(experiments(mae)[[gelname]]@files)
     fns = experiments(mae)[[gelname]]@files
     fn = fns[curchrn]
     tf = TabixFile(fn)
     eqBox3( tsbra[input$cursnp,]$allfeats[input$rank], 
         experiments(mae)[[felname]],
         tf, tsbra[input$cursnp,], annovec ) 
    })
 } # end server
shinyApp(ui=ui, server=server)
} # end tqbrowser
