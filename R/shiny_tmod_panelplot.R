


### res <- tmod_res$Nonclassical_vs_classical_ID0$tmod$pval
### res <- list(ID1=res)
### res$ID2 <- tmodCERNOtest(cntr$Nonclassical_vs_classical_ID0$symbol[ order(-abs(cntr$Nonclassical_vs_classical_ID0$log2FoldChange)) ])
### pie <- tmodDecideTests(cntr$Nonclassical_vs_classical_ID0$symbol, cntr$Nonclassical_vs_classical_ID0$log2FoldChange, cntr$Nonclassical_vs_classical_ID0$padj)
### pie <- list(ID1=pie)
### pie$ID1 <- pie$ID1$X.1
### pie$ID2 <- pie$ID1

#' Create a tmod panel plot using ggplot
#'
#' Create a tmod panel plot using ggplot
#' @importFrom tidyr pivot_longer pivot_wider everything
#' @importFrom tibble rownames_to_column
#' @export
gg_panelplot <- function(res, pie, auc_thr=.65, q_thr=.05,
                         filter_row_pval=.01,
                         filter_row_auc=.65,
                         pval_cutoff=1e-12,
                         label_angle=45) {

  label_angle=as.numeric(label_angle)
  resS <- tmodSummary(res) 

  if(any(resS$Title != resS$ID)) {
    modnames <- paste(resS$ID, resS$Title)
  } else {
    modnames <- resS$ID
  }

  names(modnames) <- resS$ID

  resS_l <- resS %>%
    pivot_longer(starts_with(c("AUC", "q")), 
                 names_to=c("Param", "Contrast"), 
                 names_sep="\\.", 
                 values_to="Value") %>% 
    pivot_wider(all_of(c("ID", "Title", "Contrast")), 
                names_from="Param", 
                values_from="Value") %>%
    mutate("q" = ifelse(is.na(.data[["q"]]), 1, .data[["q"]])) %>%
    mutate("AUC" = ifelse(is.na(.data[["AUC"]]), .5, .data[["AUC"]])) %>%
    filter(.data[["AUC"]] > auc_thr & .data[["q"]] < q_thr) %>%
    mutate(q = ifelse(.data[["q"]] < pval_cutoff, pval_cutoff, .data[["q"]])) %>%
    mutate(alpha = -log10(.data[["q"]]) / -log10(pval_cutoff))

  selMod <- resS$ID
  if(!is.na(filter_row_auc)) {
    .s <- resS_l %>% filter(.data[["AUC"]] > filter_row_auc) %>% pull("ID")
    selMod <- intersect(selMod, .s)
  }

  if(!is.na(filter_row_pval)) {
    .s <- resS_l %>% filter(.data[["q"]] < filter_row_pval) %>% pull("ID")
    selMod <- intersect(selMod, .s)
  }

  resS_l <- resS_l %>% filter(.data[["ID"]] %in% selMod)

  pie <- imap(pie, ~ { 
               colnames(.x) <- paste0(.y, '.', colnames(.x))
               .x %>% as.data.frame() %>% rownames_to_column("ID")
                })
  pieS <- Reduce(function(x, y) merge(x, y, all=TRUE), pie) %>%
    pivot_longer(-1, names_to=c("Contrast", "Direction"), 
                 names_sep="\\.", 
                 values_to="Number")

  df <- merge(resS_l, pieS, by=c("ID", "Contrast"), all.x=TRUE) %>%
    group_by(paste(.data[["ID"]], .data[["Contrast"]])) %>%
    mutate(Tot=sum(Number)) %>%
    ungroup() %>%
    mutate(Number = .data[["Number"]] * .data[["AUC"]] / .data[["Tot"]]) %>%
    mutate(Direction = factor(.data[["Direction"]], levels=c("up", "N", "down"))) 

  df <- df %>% group_by(paste0(.data[["Contrast"]], .data[["Direction"]])) %>%
    slice(match(resS$ID, .data[["ID"]])) %>%
    ungroup()

  colors <- c("red", "grey", "blue")
  names(colors) <- levels(df$Direction)

  minq <- max(-log10(df[["q"]]))

  ggplot(df, aes(x=.data[["ID"]], y=.data[["Number"]], 
                 fill=.data[["Direction"]],
                 contrast=.data[["Contrast"]],
                 id=.data[["ID"]],
                 alpha=-log10(.data[["q"]]))) + 
    facet_wrap(~ .data[["Contrast"]], nrow=1) + 
    geom_bar(stat="identity") + 
    coord_flip() +
    scale_fill_manual(values=colors) +
    scale_x_discrete(breaks=names(modnames), labels=modnames) +
    theme(strip.text.x = element_text(angle=label_angle), strip.background=element_blank(),
          axis.text.x=element_text(angle=90), axis.title.y=element_blank()) +
    scale_y_continuous(breaks=c(0, .5, 1), labels=c("0", ".5", "1"), limits=c(0, 1)) +
    guides(alpha = guide_legend(override.aes = list(fill = "grey"))) +
    lims(alpha = c(0, minq)) +
    ylab("Effect size")
    
}



#' @importFrom shinycssloaders withSpinner
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
           fluidRow(numericInput(NS(id, "font_size"), label="Font size", value = 12, 
                                 min=3, step=1, width="100%"),
                    bsTooltip(NS(id, "font_size"), "Change the font size of plot labels")),
           fluidRow(selectizeInput(NS(id, "figure_size"), label="Figure size", 
                                choices=c("800x800", 
                                          "600x800",
                                          "1200x800",
                                          "600x1200",
                                          "800x1200",
                                          "1200x1200"),
                                 options = list(create=TRUE, plugins = list('restore_on_backspace')),
                                 width="100%"),
                bsTooltip(NS(id, "figure_size"), 
                  "Change the figure size (in pixels), width x height. Press backspace to enter your own sizes.")),
           fluidRow(selectizeInput(NS(id, "label_angle"), label="Contrast label", 
                                choices=c(Slanted=45, 
                                          Vertical=90,
                                          Horizontal=0),
                                 width="100%"),
                bsTooltip(NS(id, "figure_size"), 
                  "How the contrast label should be displayed on the image.")),
            
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
          fluidRow( textOutput(NS(id, "hover_pos"))),
          column(width=12,
            withSpinner(plotOutput(NS(id, "panelPlot"), 
                                   hover=hoverOpts(NS(id, "plot_hover"), delay=50, delayType="throttle"),
                                   click=NS(id, "plot_click"),
                                   height="100%")),
           ), 
                width=9)
           



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
#' @return Returns a reactive value which is a list with elements
#' `contrast` and `id`.
#' @importFrom shinyBS bsTooltip addTooltip
#' @export
tmodPanelPlotServer <- function(id, cntr, tmod_res, tmod_dbs, tmod_map, annot=NULL) {

	moduleServer(id, function(input, output, session) {
    message("Launching tmod panelplot server")
    disable("save")

    plot_df   <- reactiveVal()
    selection <- reactiveVal()
    selection(list(click=0))

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
        pdf(file=file, width=fig_width() / 75, height=fig_height() / 75)
        gg_panelplot(res(), pie=pie(), 
                     filter_row_auc=input$filter_auc,
                     filter_row_pval=input$filter_pval,
                     label_angle=input$label_angle) + 
                                   theme(text=element_text(size=input$font_size))
        dev.off()
      }
    )

    output$hover_pos <- renderText({

      if(is.null(input$plot_hover)) { return("Hover over the plot to identify the gene sets.") }

      contrast <- input$plot_hover$panelvar1
      id       <- unlist(input$plot_hover$domain$discrete_limits$y)[
                                                        round(input$plot_hover$y) ]

      sprintf("Contrast %s, ID %s. Click to view in tmod browser panel.", 
              contrast, id)
    })


    res <- reactive({
      map(tmod_res, ~ .x[[input$db]][[input$sort]])
    })

    pie <- reactive({
      
      .make_pie(res(), cntr, tmod_dbs[[input$db]][["dbobj"]], input$db, tmod_map,
                     gene_pval=input$gene_pval,
                     gene_lfc=input$gene_lfc,
                     text.cex=input$font_size)
    })

    observeEvent(input$plot_click, {
      df <- plot_df()
      ret <- list(
                  db   = input$db,
                  sort = input$sort,
                  cntr = input$plot_click$panelvar1,
                  id   = unlist(input$plot_click$domain$discrete_limits$y)[
                                round(input$plot_click$y) ],

                  click = selection()$click + 1
                  )
      selection(ret)
    })

    fig_width <- reactiveVal()
    fig_height <- reactiveVal()

    observeEvent(input$figure_size, {
    fig_width(
      as.numeric(gsub(" *([0-9]+) *x *([0-9]+)", "\\1", input$figure_size))
    )

    fig_height(
      as.numeric(gsub(" *([0-9]+) *x *([0-9]+)", "\\2", input$figure_size))
    )
    })


    observe({ output$panelPlot <- renderPlot({
      enable("save")

      g <- gg_panelplot(res(), pie=pie(), 
                     filter_row_auc=input$filter_auc,
                     filter_row_pval=input$filter_pval,
                     label_angle=input$label_angle) + 
                                   theme(text=element_text(size=input$font_size))


      pg <- ggplot_build(g)
      plot_df(pg$data[[1]])

      g
    }, width=fig_width(), height=fig_height()) })

    return(selection)
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

    ## XXX this is a workaround for a bug in tmod; in new versions it will
    ## not be necessary.
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

