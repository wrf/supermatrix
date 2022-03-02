# multi-way interactive plot of orthologs from all-v-all MCL clustering
# supermatrix/shinyapp/app.R
# last modified 2022-03-01

library(shiny)
library(ggplot2)
library(dplyr)


cluster_table_headers = c("clusterName", "num_sequences", "num_taxa",
                          "min_per_taxon", "med_per_taxon", "max_per_taxon",
                          "min_length", "mean_length", "max_length", "taxaList")

### change this to output from makehomologs.py ###
input_cluster_file = "~/project/lactobacillus_genomics/fasta_clusters.lactobacillus_v3.tab"

cluster_data = read.table(input_cluster_file, header=FALSE, sep="\t", stringsAsFactors = FALSE,
                          col.names = cluster_table_headers )

d = mutate(cluster_data[,1:9],
           average_seqs_per_sp = num_sequences/num_taxa,
           median_vs_min = (med_per_taxon - min_per_taxon),
           max_vs_median = (max_per_taxon - med_per_taxon),
           is_avg_one = (num_sequences==num_taxa)
           )

ui <- fluidPage(
    titlePanel("Ortholog clusters all-vs-all"),
    h4(input_cluster_file),
    fluidRow(
        column(4,
        selectInput("x_select", label="X-axis control of top left plot",
                    choices = cluster_table_headers[2:9],
                    selected = "num_sequences"),
        selectInput("y_select", label="Y-axis control of top left plot",
                    choices = cluster_table_headers[2:9], 
                    selected = "num_taxa")
        )
    ),
    fluidRow(
        column(6,
           plotOutput(outputId = "seqTaxaCounts",
                      width="100%",
                      brush = "seqTaxaBrush"
           ),
            plotOutput(outputId = "avgProtperSp",
                       width="100%",
                       brush = "avgProtBrush"
            )
        ), # end column
        column(6,
               plotOutput(outputId = "maxperTaxon",
                          width="100%",
                          brush = "maxseqsBrush"
               ),
               plotOutput(outputId = "minMedMax",
                          width="100%",
                          brush = "MMMBrush"
               )
        ) # end column
    ), # end fluidRow
    fluidRow(
        column(4,
               textOutput("selectCountText"),
               ),
        tableOutput("selectedClusters")
    )
) # end fluidPage

# Define server logic ----
server <- function(input, output) {
    
    # setup for combining brushed points from all 4 plots
    vals = reactiveValues(
        cls = d[NULL,]
    )
    observeEvent( c(input$seqTaxaBrush, input$avgProtBrush, input$maxseqsBrush, input$MMMBrush), {
        vals$cls = rbind(
            brushedPoints(d, input$seqTaxaBrush, xvar = input$x_select, yvar = input$y_select ),
            brushedPoints(d, input$avgProtBrush, xvar = "num_sequences", yvar = "average_seqs_per_sp"),
            brushedPoints(d, input$maxseqsBrush, xvar = "num_sequences", yvar = "mean_length"),
            brushedPoints(d, input$MMMBrush, xvar = "median_vs_min", yvar = "max_vs_median")
        )
    })

    # display the 4 plots, and add black dots of the selected points
    output$seqTaxaCounts <- renderPlot({
        ggplot(d, aes(x=get(input$x_select), y=get(input$y_select), color=is_avg_one)) +
            labs(x=input$x_select, y=input$y_select ) +
            theme(legend.position = "none") +
            geom_point(size=3) +
            geom_point(data=vals$cls, colour="#000000")
    })
    output$avgProtperSp <- renderPlot({
         ggplot(d, aes(x=num_sequences, y=average_seqs_per_sp, color=is_avg_one)) +
            theme(legend.position = "none") +
            geom_point(size=3) +
            geom_point(data=vals$cls, colour="#000000")
    })
    output$maxperTaxon <- renderPlot({
        ggplot(d, aes(x=num_sequences, y=mean_length, color=is_avg_one)) +
            theme(legend.position = "none") +
            geom_point(size=3) +
            geom_point(data=vals$cls, colour="#000000")
    })
    output$minMedMax <- renderPlot({
        ggplot(d, aes(x=median_vs_min, y=max_vs_median, color=is_avg_one)) +
            theme(legend.position = "none") +
            geom_point(size=3) +
            geom_point(data=vals$cls, colour="#000000")
    })
    
    # print row count of selected points
    output$selectCountText <- renderText({
        paste(nrow(vals$cls), "points selected")
    })
    
    # make table of the selected points across all 4 plots
    output$selectedClusters <- renderTable({
        vals$cls
    })
}

shinyApp(ui = ui, server = server)

#