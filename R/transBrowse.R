
# transBrowse will accept
# 1) tbg: sorted output of tsByGene -- filter by distance (Inf?), MAF, and
# sort by observed score.  answers $feats, names() returns snp id
# 2) anno: annotation resource like gencodeV12 in geuvPack ... must answer
#    seqnames, $gene_id [matches feats from tsByGene], $gene_name
#    [typically  HUGO symbol]
# 3) tivcf: some tabix-indexed vcf relevant to the tbg
# 4) se: SummarizedExperiment with rownames coincident with tbg$feats and anno$gene_id


transBrowse = function(tbg, anno, tivcf, se, title="trans eQTL") { 
# add annotation to tbg for symbol and feature chr
   gn = anno$gene_name
   gsn = as.character(seqnames(anno))
   names(gn) = anno$gene_id
   names(gsn) = anno$gene_id
   gnu = gn[!duplicated(names(gn))]
   gsnu = gsn[!duplicated(names(gsn))]
   tbg$sym = gnu[tbg$feats]
   tbg$genechr = gsnu[tbg$feats]
#
   ui = fluidPage(
     titlePanel(title),
     sidebarPanel(
      helpText("select SNPs"),
      fluidRow(selectInput("snp", "SNP",
        choices = names(tbg), selected=names(tbg)[1],
        selectize=TRUE)
        ) # end row 1
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
            IND = which(names(tbg) == input$snp)[1]
            eqBox2( tbg$feats[IND], se, tivcf, 
              tbg[IND], main = paste(tbg$sym[IND], tbg$genechr[IND], sep=": "))
            })
     output$tab = renderDataTable({
            as.data.frame(tbg)
            })
   }
   
 shinyApp(ui=ui, server=server)
}
