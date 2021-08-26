
#' @rdname tmodPanelPlotServer
#' @export
tmodPanelPlotUI <- function(id, dbs, sorting) {

  ttip <- list(

    gene_pval=paste("Determines which genes are considered significant. Significant genes show as colored",
                    "fragments on the plot."),
    filter_auc=paste("Filter the gene sets by setting a minimal AUC threshold on the maximum AUC",
                     "in any of the contrasts. In other words, remove rows in which no test achieves at",
                     "least the AUC threshold."),
    filter_pval=paste("Filter the gene sets by setting a maximum p-value threshold on the minimum p value",
                     "in any of the contrasts. In other words, remove rows in which no test achieves ",
                     "the p-value threshold.")

    )

    sidebarLayout(
      sidebarPanel(
        fluidRow(
        column(
           fluidRow(selectInput(NS(id, "db"),        label="Database", choices=dbs,         width="100%"),
                    bsTooltip(NS(id, "db"), "Gene set database to be shown")),
           fluidRow(selectInput(NS(id, "sort"),      label="Sorting",  choices=sorting,     width="100%"),
                    bsTooltip(NS(id, "sort"), "Sorting order for the enrichment")),
           fluidRow(popify(numericInput(NS(id, "gene_pval"), label="P-value significance threshold for genes", 
                                 value = 0.05, min=0, step=.01, width="100%"),
                    "P-value significance threshold for genes", ttip$gene_pval)),
           fluidRow(popify(numericInput(NS(id, "gene_lfc"), label="L2FC significance threshold for genes", 
                                 value = 0.5, min=0, step=.1, width="100%"),
                    "L2FC significance threshold for genes", ttip$gene_pval)),
           width=5),
        column(
           fluidRow(numericInput(NS(id, "font_size"), label="Relative font size", value = 1, 
                                 min=.1, step=.1, width="100%"),
                    bsTooltip(NS(id, "font_size"), "Change the font size of plot labels")),
           fluidRow(numericInput(NS(id, "filter_auc"),  label="Filter by AUC (per row)",  value=0.5,
                                 min=.1, max=1, step=.05, width="100%"),
                    bsTooltip(NS(id, "filter_auc"), ttip$filter_auc)),
           fluidRow(numericInput(NS(id, "filter_pval"),  label="Filter by p-value (per row)",  
                                 value=0.05, min=0, max=1, step=0.01, width="100%"),
                    bsTooltip(NS(id, "filter_pval"), ttip$filter_auc)),
           fluidRow(downloadButton(NS(id, "save"), "Save plot to PDF", class="bg-success")),
           width=5, offset=1)),
        width=3),
      mainPanel(
          column(width=12,
            plotOutput(NS(id, "panelPlot"), height="100%"),
           ), width=9)

      )

 
}



#' Shiny module displaying tmod panel plots
#'
#' Shiny module displaying tmod panel plots
#' 
#' Tmod results, mapping and databases.
#'
#' The `tmod_res` object is a nested list with the following levels:
#'
#'  * top level are the contrasts. `names(tmod_res)` must be (set)equal to 
#'    `names(cntr)`.
#'  * next level are the names of tmod databases. `names(tmod_res[[1]])`
#'  must be (set)equal to `names(tmod_dbs)`
#*  * third level is the sorting or some other search parameter
#*  * lowest level is a data frame containing the actual result for the
#*    given contrast, database and sorting.
#' @param id Module ID
#' @param dbs named character vector with the IDs and names of the gene set databases 
#' @param sorting named character vector with the sorting types
#' @param tmod_map mapping between the PrimaryIDs from the contrasts and
#' the gene IDs from the gene set databases.
#' @param cntr list of data frames with the results of differential
#' expression analysis. Rownames must correspond to the 'PrimaryID' column
#' of data frame annot.
#' @param annot data frame containing at least the column 'PrimaryID'
#' @param tmod_res list of tmod gene set enrichment analysis
#' results. See Details.
#' @param tmod_dbs list of 
#' @importFrom shinyBS bsTooltip
#' @export
tmodPanelPlotServer <- function(id, cntr, tmod_res, tmod_dbs, tmod_map, annot=NULL) {

	moduleServer(id, function(input, output, session) {
    message("Launching tmod panelplot server")
    disable("save")

    ## Save figure as a PDF
    output$save <- downloadHandler(
      filename = function() {
        req(res())
        ret <- sprintf("tmod_panel_plot_%s_%s.pdf",
                       input$db, input$sort)
        ret <- gsub("[^0-9a-zA-Z_.-]", "", ret)
        return(ret)
      },
      content = function(file) {
        req(res())
        pdf(file=file, width=10, height=12)
        .evidence_plot(res(), pie=pie(), text.cex=input$font_size,
                     filter.rows.auc=input$filter_auc,
                     filter.rows.pval=input$filter_pval)
        dev.off()
      }
    )


    res <- reactive({
      map(tmod_res, ~ .x[[input$db]][[input$sort]])
    })

    pie <- reactive({
      
      .make_pie(res(), cntr, tmod_dbs[[input$db]][["dbobj"]], input$db, tmod_map,
                     gene_pval=input$gene_pval,
                     gene_lfc=input$gene_lfc,
                     text.cex=input$font_size)
    })

                              

    output$panelPlot <- renderPlot({
      enable("save")

      .evidence_plot(res(), pie=pie(), text.cex=input$font_size,
                     filter.rows.auc=input$filter_auc,
                     filter.rows.pval=input$filter_pval)


    }, width=800,height=800)

	})

}




.make_pie <- function(res, cntr=NULL, tmod_db_obj=NULL, dbname=NULL, tmod_map=NULL, 
                           gene_pval=0.05, gene_lfc=0.5, ...) {
  pie <- NULL

  if(!is.null(cntr)) {
    stopifnot(all(names(res) %in% names(cntr)))
    stopifnot(!is.null(tmod_db_obj) && !is.null(dbname) && !is.null(tmod_map))

    tmod_s <- tmodSummary(res)
    dbobj <- tmod_db_obj[ tmod_s[["ID"]] ]
    genes_s <- unique(unlist(dbobj[["MODULES2GENES"]]))

    mp_id <- tmod_map$dbs[[dbname]]
    mp <- tmod_map$maps[[mp_id]]

    genes_sel <- unique(unlist(map(cntr, ~ .x[["PrimaryID"]])))
    genes_sel <- genes_sel[ mp[ genes_sel ] %in% genes_s ]

    lfcs  <- map_dfc(cntr, ~ .x[ match(genes_sel, .x[["PrimaryID"]]), ][["log2FoldChange"]])
    pvals <- map_dfc(cntr, ~ .x[ match(genes_sel, .x[["PrimaryID"]]), ][["padj"]])
    lfcs[ is.na(lfcs) ] <- 0
    pvals[ is.na(pvals) ] <- 1

    pie <- tmodDecideTests(g = mp[ genes_sel ], lfc = lfcs, pval = pvals, mset=dbobj,
      lfc.thr = gene_lfc, pval.thr = gene_pval)

  }

  return(pie)
}



.evidence_plot <- function(res, pie=NULL, ...) {

  tmodPanelPlot(res, pie=pie, grid="b", ...)
}

