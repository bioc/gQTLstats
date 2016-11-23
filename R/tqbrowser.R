band2feats = function(cbstruct, bandid, gr, featExtractor = function(x)names(x)) {
 # cbstruct is a cytoband GRanges
 # bandid is the name of the cytoband to use
 # gr is the GRanges instance from which feature ids will be drawn
 # featExtractor is a function on subsetByOverlaps(gr, cbstruct[bandid])
 #  yielding a character vector of features retained
 featExtractor( subsetByOverlaps(gr, cbstruct[bandid]) )
 }

tqbrowser = function( mae, felname, gelname, tiling, tsbra, annovec,
  band.init="6q12", ermaset=NULL, ...) {
# 
# browse ranked trans-associations in tiles
#
# optional ermaset will supply information on cell-type
# specific chromatin states of genomic intervals
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
     sidebarLayout(
      sidebarPanel(
       selectInput("curband", "cytoband", choices=names(tiling),
          selected=band.init),
       uiOutput("snpSelector"),
       uiOutput("celltypeSelector"),
       numericInput("rank", "rank", 1, min=1, max=5, step=1),
       numericInput("num2lab", "# to label", 5, min=5, max=50, step=5),
       width=3
       ),
      mainPanel(
       tabsetPanel(
        tabPanel("Manh.", plotOutput("manh")),
        tabPanel("y vs GT", plotOutput("eqbox")),
        tabPanel("Manh2", ggvisOutput("p")),
        tabPanel("Manh3", plotlyOutput("manh3"))
        )
       )
    ) # end layout
   ) # end fluidPage
#   
 server = function(input, output, session) {
   output$selband = renderTable( input$curband )
   output$snpSelector = renderUI({
     tagList(
       selectInput("cursnp", "SNP", choices=
         band2feats(tiling, input$curband, tsbra))
       )
      })
   if (!is.null(ermaset)) 
    output$celltypeSelector = renderUI({
     tagList(
       selectInput("celltype", "celltype", choices=
         cellTypes(ermaset), selected=cellTypes(ermaset)[4])
       )
      })
   output$manh = renderPlot({
      x = band2feats(tiling, input$curband, tsbra, function(x) start(x))
      y = band2feats(tiling, input$curband, tsbra, function(x) 
            x$allscores[,input$rank])
      nm = band2feats(tiling, input$curband, tsbra, function(x) names(x))
      topinds = order(y,decreasing=TRUE)[1:input$num2lab]
      plot(x,y,xlab=input$curband,cex.lab=1.2)
      text(x[topinds], jitter(y[topinds]), nm[topinds], cex=1.1)
      })
   output$manh3 = renderPlotly({
    statevec = NULL
    if (!is.null(ermaset)) {
      curr = tiling[input$curband]
      print(curr)     
      seqlevelsStyle(curr) = "UCSC" # match ermaset
      rowRanges(ermaset) = curr
print(input$celltype)
      ind = which(cellTypes(ermaset) == input$celltype)
print(ind)
      curstates = subsetByRanges(ermaset[,ind], curr)[[1]][[1]] #multiple files, multiple ranges permitted, we are using 1,1
 #     print(curstates)     
      seqlevelsStyle(curstates) = seqlevelsStyle(tsbra)
      fo = findOverlaps(subsetByOverlaps(tsbra, tiling[input$curband]), curstates)
      print(fo)
      statevec = curstates$name[ subjectHits(fo) ]
      }
      x = band2feats(tiling, input$curband, tsbra, function(x) start(x))
      y = band2feats(tiling, input$curband, tsbra, function(x) 
            x$allscores[,input$rank])
      nm = band2feats(tiling, input$curband, tsbra, function(x) names(x))
#      topinds = order(y,decreasing=TRUE)[1:input$num2lab]
#      plot(x,y,xlab=input$curband,cex.lab=1.2)
#      text(x[topinds], jitter(y[topinds]), nm[topinds], cex=1.1)
      curdf = data.frame(pos=x, assoc=y, snp=nm, state=paste0(nm, ":", statevec), stringsAsFactors=FALSE)
      pp = ggplot(curdf, aes(x=pos, y=assoc, text=state)) + geom_point()
      ggplotly(pp)
      })
   tt = reactive({
#
# construct manhattanplot and tooltip
#
      pos = band2feats(tiling, input$curband, tsbra, function(x) start(x))
      y = band2feats(tiling, input$curband, tsbra, function(x) 
            x$allscores[,input$rank])
      nm = band2feats(tiling, input$curband, tsbra, function(x) names(x))
      curdf = data.frame(pos=pos, y=y, nm=nm, rowid=1:length(pos), band=input$curband) 
#
# define helper for tooltip
#
     all_values <- function(x) {
             if(is.null(x)) return(NULL)
             curdf[curdf$rowid == x$rowid, "nm" ]
           }
      curdf %>% ggvis(~pos, ~y, key:=~rowid) %>% layer_points() %>% 
                 add_tooltip(function(x) curdf[curdf$rowid==x$rowid, "nm"]) #all_values, "hover")
      })
   tt %>% bind_shiny("p")
   output$eqbox = renderPlot({
     curchrn = sub("[pq].*", "", input$curband) # character
#     print(curchrn)
#     print(experiments(mae)[[gelname]]@files)
     fns = experiments(mae)[[gelname]]@files
     fn = fns[curchrn]
     tf = TabixFile(fn)
     suppressMessages({
     eqBox3( tsbra[input$cursnp,]$allfeats[input$rank], 
         experiments(mae)[[felname]],
         tf, tsbra[input$cursnp,], annovec ) 
     })
    })
 } # end server
shinyApp(ui=ui, server=server)
} # end tqbrowser
