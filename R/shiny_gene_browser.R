## prepare contrasts for the gene browser, adding action button, sorting,
## removing unnecessary columns etc.
.gene_browser_prep_res <- function(cntr, but, annot, annot_linkout) {

  cntr   <- cntr %>% 
    map(~ { .x %>% mutate('>'= sprintf(but, PrimaryID)) }) %>%
    map(~ { .x %>% select(all_of(setdiff(colnames(.x), c("stat", "lfcSE", "symbol", "entrez")))) }) %>%
    map(~ { .x %>% { merge(annot, ., by="PrimaryID", all.x=TRUE) } %>% relocate(all_of(">"), .before=1) %>% arrange(pvalue)})

  if(!is.null(annot_linkout)) {
    cntr <- map(cntr, ~ {
                         for(n in names(annot_linkout)) {
                           fmt <- sprintf('<a href="%s" target="_blank">%%s</a>', annot_linkout[[n]])
                           .x[[n]] <- ifelse(is.na(.x[[n]]) | .x[[n]] == "", .x[[n]], sprintf(fmt, .x[[n]], .x[[n]]))
                         }
                         .x
    })
  }

  return(cntr)
}

## Wrapper around plot_gene, mainly to replace "N/A" with NA
.gene_browser_plot <- function(covar, id, covarName, rld, annot, 
                               groupBy = "N/A", colorBy = "N/A", symbolBy = "N/A", trellisBy="N/A") {
  .args <- list(id=id, xCovar=covarName, covar=covar, exprs=rld, groupBy=groupBy, annot=annot,
                colorBy=colorBy, symbolBy=symbolBy, trellisBy=trellisBy)
  ## weirdly, the line below is really, really slow
  #.args <- map(.args, ~ if(!is.na(.x) && .x == "N/A") { NA } else { .x })
  if(.args$groupBy == "N/A") .args$groupBy <- NA
  if(.args$colorBy == "N/A") .args$colorBy <- NA
  if(.args$symbolBy == "N/A") .args$symbolBy <- NA
  if(.args$trellisBy == "N/A") .args$trellisBy <- NA
  #plot_gene(pip, id, xCovar=covarName, covar=covar, rld=rld, groupBy=groupBy, colorBy=colorBy, symbolBy=symbolBy)
  do.call(plot_gene_generic, .args)
}

.default_covar <- function(covar, all_covars, default="group") {
  interesting_covars <- covar %>% 
      summary_colorDF() %>% 
      filter(unique < n()) %>% 
      pull(.data$Col)

  if(default %in% interesting_covars) {
    default_covar <- default
  } else {
    if(length(interesting_covars) > 0) {
      default_covar <- interesting_covars[1]
    } else {
      default_covar <- all_covars[1]
    }
  }

  return(default_covar)
}

## prepare the additional gene info tab panel
.gene_browser_info_tab <- function(id, x, y, covar) {
     ret <- ""

     if(is.numeric(x)) {
       pearson.test  <- cor.test(x, y, use="p")
       spearman.test <- cor.test(x, y, use="p", method="s")
       ret <- paste0(ret,
         sprintf("Correlation: r=%.2f [p = %s], rho=%.2f [p = %s]",
                 pearson.test$estimate,
                 format.pval(pearson.test$p.value, digits=2),
                 spearman.test$estimate,
                 format.pval(spearman.test$p.value, digits=2)))
     }
     return(ret)
}

#' @rdname geneBrowserTableServer
#' @export
geneBrowserTableUI <- function(id, cntr_titles) {
  but <- actionButton("foo", label=" \U25B6 ", class = "btn-primary btn-sm")
  sidebarLayout(
    sidebarPanel(
        fluidRow(selectInput(NS(id, "contrast"), label = "Contrast", choices = cntr_titles, width="100%")),
        fluidRow(
                 column(6,
                 numericInput(NS(id, "f_lfc"),    label="Filter by abs(LFC)", min=0, value=0.5, step=.1, width="100%"),
                 numericInput(NS(id, "f_pval"),    label="Filter by FDR", min=0, max=1.0, value=0.05, step=.1, width="100%")), column(6,
                 selectInput(NS(id, "f_dir"), label="Direction", choices=c(Any="an", Up="up", "Down"="dw"), 
                             width="100%")
                 )),
      tagList(
        HTML(paste("Click on the", but, "buttons to view an expression profile<br/>"))
      ),
      width=3
    ),
    mainPanel(
      dataTableOutput(NS(id, "result_tbl")),
      width=9
    )
  )
}


#' Shiny Module – gene browser table selection
#'
#' Shiny Module – gene browser table selection
#'
#' The linkout feature (parameter `annot_linkout`) allows to define how the
#' different columns from the annotation data frame are represented as
#' linkouts to external data bases. 
#'
#' The parameter `annot_linkout` is a named list. Names must correspond to
#' columns from the annotation data frame. The elements of the list are
#' character strings containing URLs with the `%s` placeholder. For
#' example, if the column `ENSEMBL` contains ENSEMBL identifiers, you can
#' link out by specifying 
#'
#' ```
#' annot_linkout=list(ENSEMBL="https://www.ensembl.org/id/%s")
#' ```
#' @param cntr a list of data frames containing the DE analysis results
#' @param annot annotation data frame containing column 'PrimaryID'
#'        corresponding to the rownames of the contrast data frames
#' @param id identifier for the namespace of the module
#' @param cntr_titles named character vector for contrast choices
#' @param annot_linkout a list; see Details. 
#' @return reactive value containing the gene ID
#' @export
geneBrowserTableServer <- function(id, cntr, annot, annot_linkout=NULL) {

  if(!is.null(annot_linkout)) {
    stopifnot(all(names(annot_linkout) %in% names(annot)))
  }

  moduleServer(id, function(input, output, session) {


    gene_id <- reactiveVal("")

    but <- actionButton("go_%s", label=" \U25B6 ", 
                         onclick=sprintf('Shiny.onInputChange(\"%s-select_button\",  this.id)', id),  
                         class = "btn-primary btn-sm")

    cntr <- .gene_browser_prep_res(cntr, as.character(but), annot, annot_linkout)

    observeEvent(input$select_button, {
      gene_id(gsub("^go_", "", input$select_button))
    })

    output$id_summary <- renderText({
      .cntr <- input$contrast
      sprintf("Contrast is %s", .cntr)
    })

    output$result_tbl <- DT::renderDataTable({
      res <- cntr[[ input$contrast ]]
      if(input$f_dir == "up") {
        res <- res %>% filter(.data[["log2FoldChange"]] > 0)
      } else if(input$f_dir == "dw") {
        res <- res %>% filter(.data[["log2FoldChange"]] < 0)
      }

      res %>% filter(.data[["padj"]] < input$f_pval & abs(.data[["log2FoldChange"]]) > input$f_lfc) %>%
      datatable(escape=FALSE, selection='none', extensions="Buttons",
                options=list(pageLength=5, dom="Bfrtip", scrollX=TRUE, buttons=c("copy", "csv", "excel"))) %>%
        formatSignif(columns=intersect(colnames(cntr[[ input$contrast ]]), 
                                       c("baseMean", "log2FoldChange", "pvalue", "padj")), digits=2)
    })

    gene_id
  })
}

#' @rdname geneBrowserPlotServer
#' @export
geneBrowserPlotUI <- function(id, covar, contrasts=FALSE) {
  all_covars         <- covar %>% summary_colorDF() %>% filter(unique > 1) %>% pull(.data$Col)
  non_unique         <- covar %>% summary_colorDF() %>% filter(unique < nrow(covar)) %>% pull(.data$Col)
  default_covar <- .default_covar(covar, all_covars, default="group")
  col_control <-  sidebarPanel(
      fluidRow(downloadButton(NS(id, "save"), "Save plot to PDF", class="bg-success")),
      fluidRow(
      column(width=5, 
      fluidRow(selectInput(NS(id, "covarName"), "X covariate", non_unique, selected=default_covar, width="100%")),
      fluidRow(selectInput(NS(id, "colorBy"), "Color by", c("N/A", non_unique), selected="N/A", width="100%")),
      fluidRow(selectInput(NS(id, "symbolBy"), "Symbol by", c("N/A", non_unique), selected="N/A", width="100%")),
      ),
      column(width=5,
      fluidRow(selectInput(NS(id, "groupBy"), "Link data points by", c("N/A", non_unique), selected="N/A", width="100%")),
      fluidRow(selectInput(NS(id, "trellisBy"), "Trellis by", c("N/A", non_unique), selected="N/A", width="100%")),
      fluidRow(numericInput(NS(id, "fontSize"),    label="Font size", min=6, value=14, step=1, width="50%")),
      offset=1)
      ),
      fluidRow(textOutput(NS(id, "addInfo"))),
      fluidRow(verbatimTextOutput(NS(id, "geneData")))
    )
  if(contrasts) {
    return(sidebarLayout(col_control,
      mainPanel(
      column(9, style="padding:20px;", tabsetPanel(
      tabPanel("Plot", fluidRow(br(), plotOutput(NS(id, "countsplot")))),
      tabPanel("Contrast overview", fluidRow(br(), dataTableOutput(NS(id, "contr_sum"))))
      )))))
  } else {
    return(col_control,
      mainPanel(plotOutput(NS(id, "countsplot"))))
  }

}

#' Shiny Module – gene browser expression profile plot
#'
#' Shiny Module – gene browser expression profile plot
#'
#' The `gene_id` parameter must be a reactive value, because that is the
#' whole point of the plotting module: observe changes to the gene ID and
#' update the plot accordingly.
#' 
#' In contrast, other parameters must not be reactive values. This may
#' change in future to allow for dynamic exchange of data sets.
#' @param gene_id primary identifier of the gene to show. This must be a
#'        reactive value
#' @param exprs expression matrix; row names must correspond to the primary identifiers
#' @param contrasts (logical) whether or not create an additional panel
#'        next to the plot which can be used to show detailed contrast
#'        information for a gene
#' @param annot (optional) annotation data frame containing column 'PrimaryID'
#'        corresponding to the rownames of the contrast data frames
#' @param id identifier (same as the one passed to geneBrowserTableUI)
#' @param covar data frame with all covariates
#' @param cntr (optional) list of contrasts
#' @param symbol_col name of the column in `annot` which contains the gene
#'        symbols; use NULL if no such column
#' @param description_col name of the column in `annot` which contains the gene
#'        title / description; use NULL if no such column
#' @return does not return anything useful
#' @examples
#' if(interactive()) {
#'    mtx <- matrix(rnorm(40, mean=rep(c(0, 1), each=20)), nrow=1)
#'    rownames(mtx) <- "MUZG"
#'    covar <- data.frame(
#'                        em=rep(LETTERS[1:2], each=20),
#'                        pstrem=rep(letters[1:20], 2),
#'                        bzdrem=rnorm(40))
#'   
#'    ui  <- fluidPage(geneBrowserPlotUI("gplot", covar))
#'    serv <- function(input, output, session) {
#'      gid <- reactiveVal()
#'      gid("MUZG")
#'      geneBrowserPlotServer("gplot", gid, covar, mtx)
#'    }
#'    shinyApp(ui, serv)
#' }
#' @export
geneBrowserPlotServer <- function(id, gene_id, covar, exprs, annot=NULL, cntr=NULL, 
                                  symbol_col="SYMBOL", description_col="GENENAME") {
  stopifnot(is.reactive(gene_id))
  stopifnot(!is.reactive(covar))
  stopifnot(!is.reactive(exprs))
  stopifnot(!is.reactive(annot))
  stopifnot(is.null(annot) || is.data.frame(annot))
  if(!is.null(annot)) {
    stopifnot(all(c("PrimaryID", symbol_col, description_col) %in% 
                  colnames(annot)))
  }

  moduleServer(id, function(input, output, session) {
    disable("save")

    ## Save figure as a PDF
    output$save <- downloadHandler(
      filename = function() {
        req(gene_id())
        ret <- sprintf("expression_profile_%s_covarX_%s_colorBy_%s_groupBy_%s_symbolBy_%s_trellisBy_%s.pdf",
                       gene_id(), input$covarName, input$colorBy, input$groupBy, input$symbolBy, input$trellisBy)
        ret <- gsub("[^0-9a-zA-Z_.-]", "", ret)
        return(ret)
      },
      content = function(file) {
        req(gene_id())
        pdf(file=file, width=8, height=5)
        g <- .gene_browser_plot(covar, gene_id(), input$covarName, exprs, annot, 
                           input$groupBy, input$colorBy, input$symbolBy, input$trellisBy)
        print(g)
        dev.off()
      }
    )

    # Show a turbo card for a gene
    output$geneData <- renderText({
      req(gene_id())
      ret <- sprintf("PrimaryID: %s", gene_id())
      if(!is.null(annot)) {
        m <- match(gene_id(), annot$PrimaryID)
        if(!is.null(symbol_col)) {
          ret <- paste0(ret,
                       sprintf("\nSymbol: %s",
                               annot[m, ][[symbol_col]])) 
        }
        if(!is.null(description_col)) {
          ret <- paste0(ret, 
                        sprintf("\nDescription: %s",
                                annot[m, ][[description_col]]))
        }
      }
      return(ret)
    })
    if(!is.null(cntr)) {
      output$contr_sum <- DT::renderDataTable({
        req(gene_id())
        cn <- names(cntr)
        res <- imap_dfr(cntr, ~ .x %>% 
                        filter(PrimaryID == gene_id()) %>% 
                    mutate(contrast=.y))

        numcol <- map_lgl(res, is.numeric)
        res %>% datatable(escape=FALSE, selection='none', options=list(pageLength=5)) %>%
          formatSignif(columns=colnames(res)[numcol], digits=2)
      })
    }

    ## Additional information - e.g. correlation coefficient if the
    ## covariate is numeric
    output$addInfo <- renderText({
      req(gene_id())
      .gene_browser_info_tab(gene_id(), covar[[input$covarName]], exprs[gene_id(), ])
    })

    ## The actual plot
    output$countsplot <- renderPlot({
      req(gene_id())
      enable("save")
      .gene_browser_plot(covar, gene_id(), input$covarName, exprs, annot, 
                         input$groupBy, input$colorBy, input$symbolBy, input$trellisBy) +
                                    theme(text=element_text(size=input$fontSize))
    })
  })
}


#' Launch a browser of DE analysis results
#'
#' Launch a shiny-based browser of DE analysis results
#'
#' Launches a shiny app, web based, which allows to show gene expression
#' profiles in a pipeline. 
#'
#' To speed up launching, you can pre-load the contrasts with the
#' `get_contrasts` function.
#' @param pip pipeline object returned by `load_de_pipeline`
#' @param cntr (optional) pre-loaded contrasts (returned by `get_contrasts`)
#' @param annot (optional) pre-loaded annotation table (returned by `get_annot`)
#' @return does not return a value
#' @importFrom rlang .data
#' @importFrom stats cor.test
#' @importFrom bslib bs_theme
#' @examples
#' \dontrun{
#' pip <- load_de_pipeline(config_file="DE_config.yaml")
#' gene_browser(pip, tmod_dbs)
#' }
#' @export
gene_browser <- function(pip, cntr=NULL, annot=NULL) {

  message("preparing...")
  if(is.null(annot)) {
    message(" * Loading Annotation (consider using the annot option to speed this up)")
    annot  <- get_annot(pip)
  }

  # prepare the contrast tables
  if(is.null(cntr)) {
    message(" * Loading contrasts... (consider using the cntr option to speed this up)")
    cntr <- get_contrasts(pip)
  }

  config <- get_config(pip)
  covar  <- get_covariates(pip)

  rld    <- get_object(pip, step="DESeq2", extension="rld.blind.rds")
  rld    <- rld@assays@data@listData[[1]]
  
  cntr_titles <- map_chr(config$contrasts$contrast_list, `[[`, "ID")
  names(cntr_titles) <- map_chr(config$contrasts$contrast_list, `[[`, "title")

  thematic_shiny(font="auto")

  ## prepare the UI
  ui <- fluidPage(
    useShinyjs(),
    theme = bs_theme(bootswatch = "united"),
    fluidRow(titlePanel(h1("Gene browser")), class="bg-primary"),
    fluidRow(HTML("<hr>")),
    geneBrowserTableUI("geneTab", cntr_titles),
    geneBrowserPlotUI("genePlot", covar, contrasts=TRUE)
  )

  ## prepare the server
  server <- function(input, output, session) {
    gene_id <- geneBrowserTableServer("geneTab", cntr, annot)
    geneBrowserPlotServer("genePlot", gene_id, covar, rld, annot)

    message("launching!")
  }

  shinyApp(ui, server)
}

