
#' Launch a browser of tmod gene set enrichment analysis results
#'
#' Launch a shiny-based browser of tmod gene set enrichment analysis results
#' @param pip pipeline object returned by `load_de_pipeline`
#' @param tmod_dbs tmod db object returned by get_tmod_dbs
#' @return does not return a value
#' @import dplyr
#' @importFrom shiny shinyApp renderText renderPlot verbatimTextOutput textOutput
#' @importFrom shiny actionButton column fluidPage fluidRow h1 mainPanel plotOutput reactiveValues selectInput sidebarLayout sidebarPanel titlePanel
#' @importFrom DT datatable formatSignif renderDataTable dataTableOutput
#' @importFrom colorDF summary_colorDF
#' @examples
#' \dontrun{
#' pip <- load_de_pipeline(config_file="DE_config.yaml")
#' tmod_dbs <- get_tmod_dbs(pip)
#' tmod_browser(pip, tmod_dbs)
#' }
#' @export
tmod_browser <- function(pip, tmod_dbs=NULL) {

  message("preparing...")
  but <- as.character(actionButton("go_%s-!-%s-!-%s-!-%s", label=" > ", 
                                   onclick='Shiny.onInputChange(\"select_button\",  this.id)'))
  annot    <- get_annot(pip)
  tmod_map <- get_tmod_mapping(pip)

  config <- get_config(pip)
  cntr_titles <- map_chr(config$contrasts$contrast_list, `[[`, "ID")
  names(cntr_titles) <- map_chr(config$contrasts$contrast_list, `[[`, "title")

  cntr   <- get_contrasts(pip) %>% map(~ .x %>% rownames_to_column("PrimaryID")) 

  # prepare the contrast tables
  tmod_res <- get_tmod_res(pip) %>% imap(~ {
    .cntr <- .y
    imap(.x, ~ {
           .dbname <- .y
           imap(.x, ~ {
                  .sort <- .y
                  .x %>% mutate(">"=sprintf(but, ID, .cntr, .dbname, .sort)) %>%
                    select(-cES, -cerno) %>%
                    arrange(P.Value) %>%
                    relocate(all_of(">"), .before=1)
            })
         })
  })

  if(is.null(tmod_dbs)) {
    tmod_dbs <- get_tmod_dbs(pip)
  }
  dbs <- names(tmod_dbs)
  sorting <- config$tmod$sort_by

  ui <- fluidPage(
    titlePanel(h1("tmod browser")),
    sidebarLayout(
      sidebarPanel(
         fluidRow(selectInput("contrast", label="Contrast", choices=cntr_titles)),
         fluidRow(selectInput("db", label="Database", choices=dbs)),
         fluidRow(selectInput("sort", label="Sorting", choices=sorting))
      ),
      mainPanel(
        dataTableOutput("tmodResTab")
      )
    ),
    sidebarLayout(
      sidebarPanel(
        fluidRow(verbatimTextOutput("modinfo"))
      ),
      mainPanel(
        fluidRow(verbatimTextOutput("cmdline")),
        fluidRow(plotOutput("evidencePlot"))
      )
    )
  )

  server <- function(input, output, session) {
    output$tmodResTab <- renderDataTable({
      res <- tmod_res[[input$contrast]][[input$db]][[input$sort]]
      datatable(res, escape=FALSE, selection='none', options=list(pageLength=5)) %>%
        formatSignif(columns=intersect(colnames(res), 
                                       c("AUC", "cerno", "P.Value", "adj.P.Val")), digits=2)
    })
    
    selvar <- reactiveValues()

    output$modinfo <- renderText({
      if(is.null(selvar$id)) { return("") }
      .db <- tmod_dbs[[selvar$db]][["dbobj"]]
      .id <- selvar$id
      ret <- sprintf("Module ID: %s\nDescription: %s\nContrast: %s\nDatabase: %s\nSorting: %s",
              .id,
              .db[["MODULES"]][.id, ][["Title"]],
              selvar$cntr,
              selvar$db,
              selvar$sort)

      genes <- .db[["MODULES2GENES"]][[.id]]
      genes_pid <- tmod_rev_db_map_ids(pip, ids=genes, dbname=selvar$db, tmod_dbs_mapping_obj=tmod_map)
      genes_sign <- genes_pid %in% (cntr[[selvar$cntr]] %>% filter(.data$padj < 0.05) %>% pull(.data$PrimaryID))
      ret <- paste0(ret, "\nSignificant genes:\n",     paste(genes[genes_sign],  collapse=","), "\n")
      ret <- paste0(ret, "\nNon-significant genes:\n", paste(genes[!genes_sign], collapse=","), "\n")
      return(ret)
    })

    output$evidencePlot <- renderPlot({

      if(length(input$select_button) == 0) {
        return(NULL)
      }
      id <- gsub("go_", "", input$select_button)
      tmp <- unlist(strsplit(gsub("go_", "", id), "-!-"))
      names(tmp) <- c("id", "cntr", "db", "sort")
      selvar$id <- tmp[1]
      selvar$cntr <- tmp[2]
      selvar$db <- tmp[3]
      selvar$sort <- tmp[4]
      output$cmdline <- renderText({

        sprintf('plot_evidence(pip, id="%s", dbname="%s", contrast="%s")', 
                selvar$id, selvar$db, selvar$cntr)

      })
      plot_evidence(pip, id=selvar$id, dbname=selvar$db, contrast=selvar$cntr, tmod_dbs_obj=tmod_dbs)
    })
  }
  shinyApp(ui, server)
}

#' Launch a browser of DE analysis results
#'
#' Launch a shiny-based browser of DE analysis results
#' @param pip pipeline object returned by `load_de_pipeline`
#' @return does not return a value
#' @importFrom rlang .data
#' @importFrom stats cor.test
#' @examples
#' \dontrun{
#' pip <- load_de_pipeline(config_file="DE_config.yaml")
#' gene_browser(pip, tmod_dbs)
#' }
#' @export
gene_browser <- function(pip) {

  message("preparing...")
  but <- as.character(actionButton("go_%s", label=" > ", 
                                   onclick='Shiny.onInputChange(\"select_button\",  this.id)'))
  annot  <- get_annot(pip)

  # prepare the contrast tables
  cntr   <- get_contrasts(pip) %>% map(~ .x %>% rownames_to_column("PrimaryID")) %>% 
    map(~ { .x %>% mutate(go= sprintf(but, PrimaryID)) }) %>%
    map(~ { .x %>% select(all_of(setdiff(colnames(.x), c("symbol", "entrez")))) }) %>%
    map(~ { .x %>% { merge(annot, ., by="PrimaryID", all.x=TRUE) } %>% select(-baseMean, -stat, -lfcSE) %>% arrange(pvalue)})
  config <- get_config(pip)

  covar  <- get_covariates(pip)

  interesting_covars <- covar %>% summary_colorDF() %>% filter(unique < n()) %>% pull(.data$Col)
  all_covars         <- covar %>% summary_colorDF() %>% filter(unique > 1) %>% pull(.data$Col)
  if("group" %in% interesting_covars) {
    default_covar <- "group"
  } else {
    if(length(interesting_covars) > 0) {
      default_covar <- interesting_covars[1]
    } else {
      default_covar <- all_covars[1]
    }
  }

  ds     <- get_deseq2(pip)
  rld    <- get_object(pip, step="DESeq2", extension="rld.blind.rds")
  rld    <- rld@assays@data@listData[[1]]
  
  cntr_titles <- map_chr(config$contrasts$contrast_list, `[[`, "ID")
  names(cntr_titles) <- map_chr(config$contrasts$contrast_list, `[[`, "title")

  ui <- fluidPage(
    titlePanel(h1("Gene browser")),
    fluidRow(
      column(2, selectInput("contrast", label = "Dataset", choices = cntr_titles)),
      column(10, dataTableOutput("result_tbl"))
      ),
    sidebarLayout(
    sidebarPanel(
       fluidRow(selectInput("covarName", "X covariate", all_covars, selected=default_covar)),
       fluidRow(selectInput("groupBy", "Group by", c("N/A", all_covars), selected="N/A")),
       fluidRow(selectInput("colorBy", "Color by", c("N/A", all_covars), selected="N/A")),
       fluidRow(selectInput("symbolBy", "Symbol by", c("N/A", all_covars), selected="N/A")),
       fluidRow(textOutput("addInfo")),
       fluidRow(verbatimTextOutput("geneData"))
    ),

    mainPanel(
      plotOutput("countsplot")
      
    ))
  )

  server <- function(input, output, session) {
   output$id_summary <- renderText({
     .cntr <- input$contrast
     sprintf("Contrast is %s", .cntr)
   })

    output$result_tbl <- DT::renderDataTable({
      message(glue("preparing table contrast {input$contrast}"))
      message(nrow(cntr[[ input$contrast ]]))
      datatable(cntr[[ input$contrast ]], escape=FALSE, selection='none', options=list(pageLength=5)) %>%
        formatSignif(columns=intersect(colnames(cntr[[ input$contrast ]]), 
                                       c("log2FoldChange", "pvalue", "padj")), digits=2)
    })

    output$geneData <- renderText({
      id <- gsub("go_", "", input$select_button)
      if(length(id) == 0) {
        id <- rownames(rld)[1]
      }

      sprintf("PrimaryID: %s\nSymbol: %s\nDescription: %s",
              id,
              annot[match(id, annot$PrimaryID), "SYMBOL"],
              annot[match(id, annot$PrimaryID), "GENENAME"])
    })

   #output$addInfo <- renderText({
   #  message("calculating info")
   #  id <- gsub("go_", "", input$select_button)
   #  if(length(id) == 0) {
   #    id <- rownames(rld)[1]
   #  }
   #  ret <- ""

   #  x <- covar[[input$covarName]]

   #  if(is.numeric(x)) {
   #    y <- rld[id, ]
   #    pearson.test  <- cor.test(x, y, use="p")
   #    spearman.test <- cor.test(x, y, use="p", method="s")
   #    ret <- paste0(ret,
   #      sprintf("Correlation: r=%.2f [p = %s], rho=%.2f [p = %s]",
   #              pearson.test$estimate,
   #              format.pval(pearson.test$p.value, digits=2),
   #              spearman.test$estimate,
   #              format.pval(spearman.test$p.value, digits=2)))
   #  }
   #  return(ret)
   #})

    output$countsplot <- renderPlot({
      id <- gsub("go_", "", input$select_button)
      if(length(id) == 0) {
        id <- rownames(rld)[1]
      }
      .gene_browser_plot(pip, covar, id, input$covarName, rld, annot, input$groupBy, input$colorBy, input$symbolBy)
    })
  }

  message("launching!")

  shinyApp(ui, server)
}

.gene_browser_plot <- function(pip, covar, id, covarName, rld, annot, 
                               groupBy = "N/A", colorBy = "N/A", symbolBy = "N/A") {

  .args <- list(x=pip, id=id, xCovar=covarName, covar=covar, rld=rld, groupBy=groupBy, annot=annot,
                colorBy=colorBy, symbolBy=symbolBy)
  #.args <- map(.args, ~ if(!is.na(.x) && .x == "N/A") { NA } else { .x })
  if(.args$groupBy == "N/A") .args$groupBy <- NA
  if(.args$colorBy == "N/A") .args$colorBy <- NA
  if(.args$symbolBy == "N/A") .args$symbolBy <- NA
  #plot_gene(pip, id, xCovar=covarName, covar=covar, rld=rld, groupBy=groupBy, colorBy=colorBy, symbolBy=symbolBy)
  do.call(plot_gene, .args)
}


