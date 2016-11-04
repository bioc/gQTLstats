
# transBrowse will accept
# 1) tbg: sorted output of tsByGene -- filter by distance (Inf?), MAF, and
# sort by observed score.  answers $feats, names() returns snp id
# 2) anno: annotation resource like gencodeV12 in geuvPack ... must answer
#    seqnames, $gene_id [matches feats from tsByGene], $gene_name
#    [typically  HUGO symbol]
# 3) tivcf: some tabix-indexed vcf relevant to the tbg
# 4) se: SummarizedExperiment with rownames coincident with tbg$feats and anno$gene_id


transBrowse2 = function(tbga, annovec, tivcf, se, title="trans eQTL",
  maxrank=3) { 
# uses transByRankAccum
#
   ui = fluidPage(
     titlePanel(title),
     sidebarPanel(
      helpText("select SNPs"),
      fluidRow(selectInput("snp", "SNP",
        choices = names(tbga), selected=names(tbga)[1],
        selectize=TRUE)
        ), # end row 1
      fluidRow(numericInput("frank", "featureRank", value=1, min=1, max=maxrank,
          step=1)
        ) # end row 2
      ), # end sidebar
     mainPanel(
      tabsetPanel(
       tabPanel("swarm", plotOutput("swarm")),
       tabPanel("stats", dataTableOutput("tab"))
      ) # end tabset
     ) # end main
   ) # end fluidPage
   
   server = function(input, output) {
     output$swarm = renderPlot({
            IND = which(names(tbga) == input$snp)[1]
            curfeat = tbga[IND]$allfeats[input$frank]
            eqBox3( tbga[IND]$allfeats[input$frank], se, tivcf, 
              tbga[IND], annovec )
            })
     output$tab = renderDataTable({
            IND = which(names(tbga) == input$snp)[1]
            curfeat = tbga[IND]$allfeats[input$frank]
            as.data.frame(tbga)
            })
   }
   
 shinyApp(ui=ui, server=server)
}
