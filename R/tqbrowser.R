band2feats = function(cbstruct, bandid, gr, featExtractor = function(x)names(x)) {
 # cbstruct is a cytoband GRanges
 # bandid is the name of the cytoband to use
 # gr is the GRanges instance from which feature ids will be drawn
 # featExtractor is a function on subsetByOverlaps(gr, cbstruct[bandid])
 #  yielding a character vector of features retained
 featExtractor( subsetByOverlaps(gr, cbstruct[bandid]) )
 }

tqbrowser = function( mae, felname, gelname, tiling, tsbra, annovec,
  band.init="6q12", ermaset, gwascat, ...) {
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
#       actionButton("act", "Submit"),
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
        tabPanel("Manh.", helpText("plotly Manhattan plot, select subplots by mouse, mouseover for point metadata; note points above y=0 are eQTL association scores, points below y=0 are gwascat findings (-log10 p)"), plotlyOutput("manh")),
        tabPanel("y vs GT", plotOutput("eqbox"))
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
    output$celltypeSelector = renderUI({
     tagList(
       selectInput("celltype", "celltype for chrom. state annotation", choices=
         cellTypes(ermaset), selected=cellTypes(ermaset)[4])
       )
      })
   output$eqbox = renderPlot({
     curchrn = sub("[pq].*", "", input$curband) # character
     fns = experiments(mae)[[gelname]]@files
     fn = fns[curchrn]
     tf = TabixFile(fn)
     if (!is.null(input$cursnp)) {  # delay while renderUI sets up
      tb = tsbra[input$cursnp,]
      tbf = tb$allfeats[input$rank]
      suppressMessages({
       eqBox3( tbf, 
         experiments(mae)[[felname]],
         tf, tb, annovec )
       })
      }
   })
  output$manh = renderPlotly({
        curr = tiling[input$curband]
        seqlevelsStyle(curr) = "UCSC" # match ermaset
        rowRanges(ermaset) = curr
        chk = c(is.null(input$celltype), is.null(input$rank), is.null(input$curband))
        if (!any(chk)) {
          ind = which(cellTypes(ermaset) == input$celltype)
          curstates = subsetByRanges(ermaset[,ind], 
                  curr)[[1]][[1]] #multiple files, multiple ranges permitted, we are using 1,1
          seqlevelsStyle(curstates) = seqlevelsStyle(tsbra)[1]
          fo = findOverlaps(subsetByOverlaps(tsbra, tiling[input$curband]), curstates)
          statevec = curstates$name[ subjectHits(fo) ]
          x = band2feats(tiling, input$curband, tsbra, function(x) start(x))
          y = band2feats(tiling, input$curband, tsbra, function(x) 
              x$allscores[,input$rank])
          nm = band2feats(tiling, input$curband, tsbra, function(x) names(x))
          gwascat = as(gwascat, "GRanges")
          genome(gwascat) = genome(curr)[1]
          seqlevelsStyle(gwascat) = seqlevelsStyle(curr)[1]
          gw = subsetByOverlaps(gwascat, curr)
          mcg = mcols(gw)
          seqlevelsStyle(gw) = seqlevelsStyle(curstates)[1]
          gfo = findOverlaps(gw, curstates)
          gstatevec = curstates$name[ subjectHits(gfo) ]
          gwdf = data.frame(pos=start(gw), assoc=-mcg[,"PVALUE_MLOG"],
                  snp=mcg[,"STRONGEST SNP-RISK ALLELE"],
                  stateid=mcg[,"DISEASE/TRAIT"], state=gstatevec,
                  stringsAsFactors=FALSE)
          curdf = data.frame(pos=x, assoc=y, snp=nm, 
                 stateid=paste0(nm, ":", statevec), 
                 state=statevec, stringsAsFactors=FALSE)
          curdf = rbind(curdf, gwdf)
          ggp = ggplot(curdf, aes(x=pos, y=assoc, text=stateid, 
                colour=state))+geom_point() + labs(x=input$curband, y="<0: -log10p gwas, >0:qtl assoc")
          if (!is.null(ggp)) ggplotly( ggp )
#             ggplot(curdf, aes(x=pos, y=assoc, text=state, 
#                colour=stateOnly))+geom_point()) 
          }
     })
 } # end server
shinyApp(ui=ui, server=server)
} # end tqbrowser
