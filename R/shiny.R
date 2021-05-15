.tmod_browser_prepare_res <- function(pip, but, tmod_res=NULL) {

  if(is.null(tmod_res)) {
    message(" * Reading tmod results (consider using it as an argument)...")
    tmod_res <- get_tmod_res(pip)
  }

  # prepare the tmod res
  tmod_res <- tmod_res %>% imap(~ {
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

  return(tmod_res)
}


## given a module, contrast, sorting prepare a module info tab contents,
## including which genes are significant in the given contrast
.tmod_browser_mod_info <- function(pip, id, db_name, cntr_name, sort_name, tmod_dbs, cntr, tmod_map) {

  db <- tmod_dbs[[db_name]][["dbobj"]]
  ret <- sprintf("Module ID: %s\nDescription: %s\nContrast: %s\nDatabase: %s\nSorting: %s",
          id,
          db[["MODULES"]][id, ][["Title"]],
          cntr_name,
          db_name, 
          sort_name)

  genes <- db[["MODULES2GENES"]][[id]]
  genes_pid <- tmod_rev_db_map_ids(pip, ids=genes, dbname=db_name, tmod_dbs_mapping_obj=tmod_map)
  genes_sign <- genes_pid %in% (cntr[[cntr_name]] %>% filter(.data$padj < 0.05) %>% pull(.data$PrimaryID))
  ret <- paste0(ret, "\nSignificant genes:\n",     paste(genes[genes_sign],  collapse=","), "\n")
  ret <- paste0(ret, "\nNon-significant genes:\n", paste(genes[!genes_sign], collapse=","), "\n")
  return(ret)

}

## returns the UI of the tmod browser
tmodBrowserUI <- function(cntr_titles, dbs, sorting, but) {

  ui <- fluidPage(
    theme = bs_theme(bootswatch = "sandstone"),
    fluidRow(titlePanel(h1("tmod browser")), class="bg-primary"),
    fluidRow(HTML("<hr>")),
    sidebarLayout(
      sidebarPanel(
         fluidRow(selectInput("contrast", label="Contrast", choices=cntr_titles, width="100%")),
         fluidRow(selectInput("db", label="Database", choices=dbs, width="100%")),
         fluidRow(selectInput("sort", label="Sorting", choices=sorting, width="100%")),
         HTML(paste("Click on the", but, "buttons to view an expression profile"))
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

  return(ui)
}

#' Launch a browser of tmod gene set enrichment analysis results
#'
#' Launch a shiny-based browser of tmod gene set enrichment analysis results
#'
#' To speed up launching the browser, you can load the tmod_dbs and
#' tmod_res objects first (using the functions get_tmod_dbs and
#' get_tmod_res).
#' @param pip pipeline object returned by `load_de_pipeline`
#' @param tmod_dbs tmod db object returned by get_tmod_dbs
#' @param tmod_res tmod results object returned by get_tmod_res
#' @return does not return a value
#' @import dplyr
#' @importFrom shiny shinyApp renderText renderPlot verbatimTextOutput textOutput renderUI uiOutput
#' @importFrom shiny actionButton column fluidPage fluidRow h1 mainPanel plotOutput reactiveValues selectInput sidebarLayout sidebarPanel titlePanel 
#' @importFrom shiny showNotification removeNotification req
#' @importFrom DT datatable formatSignif renderDataTable dataTableOutput
#' @importFrom htmltools HTML
#' @importFrom colorDF summary_colorDF
#' @importFrom thematic thematic_shiny
#' @examples
#' \dontrun{
#' pip <- load_de_pipeline(config_file="DE_config.yaml")
#' tmod_dbs <- get_tmod_dbs(pip)
#' tmod_browser(pip, tmod_dbs)
#' }
#' @export
tmod_browser <- function(pip, tmod_dbs=NULL, tmod_res=NULL) {

  message("preparing...")
  but <- as.character(actionButton("go_%s-!-%s-!-%s-!-%s", label=" > ", 
                                   onclick='Shiny.onInputChange(\"select_button\",  this.id)',  class = "btn-primary btn-sm"))
  annot    <- get_annot(pip)
  tmod_map <- get_tmod_mapping(pip)

  config <- get_config(pip)
  cntr_titles <- map_chr(config$contrasts$contrast_list, `[[`, "ID")
  names(cntr_titles) <- map_chr(config$contrasts$contrast_list, `[[`, "title")

  cntr   <- get_contrasts(pip) %>% map(~ .x %>% rownames_to_column("PrimaryID")) 
  tmod_res <- .tmod_browser_prepare_res(pip, but, tmod_res)

  if(is.null(tmod_dbs)) {
    message(" * Reading tmod_dbs (consider using it as an argument)...")
    tmod_dbs <- get_tmod_dbs(pip)
  }

  dbs <- names(tmod_dbs)
  sorting <- config$tmod$sort_by

  thematic_shiny(font="auto")

  ## prepare the tmod browser UI
  ui <- tmodBrowserUI(cntr_titles, dbs, sorting, but)

  server <- function(input, output, session) {
    output$tmodResTab <- renderDataTable({
      id <- showNotification("Rendering table...", duration = NULL, closeButton = FALSE)
      on.exit(removeNotification(id), add = TRUE)
      res <- tmod_res[[input$contrast]][[input$db]][[input$sort]]
      datatable(res, escape=FALSE, selection='none', options=list(pageLength=5)) %>%
        formatSignif(columns=intersect(colnames(res), 
                                       c("AUC", "cerno", "P.Value", "adj.P.Val")), digits=2)
    })
    
    selvar <- reactiveValues()

    output$modinfo <- renderText({
      req(selvar$id)
      ret <- .tmod_browser_mod_info(pip, 
                                    selvar$id, selvar$db, selvar$cntr, selvar$cntr, 
                                    tmod_dbs, cntr, tmod_map)
      return(ret)
    })

    output$evidencePlot <- renderPlot({

      req(input$select_button)
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


## prepare contrasts for the gene browser, adding action button, sorting,
## removing unnecessary columns etc.
.gene_browser_prep_res <- function(cntr, but, annot) {

  cntr   <- cntr %>% map(~ .x %>% rownames_to_column("PrimaryID")) %>% 
    map(~ { .x %>% mutate('>'= sprintf(but, PrimaryID)) }) %>%
    map(~ { .x %>% select(all_of(setdiff(colnames(.x), c("stat", "lfcSE", "symbol", "entrez")))) }) %>%
    map(~ { .x %>% { merge(annot, ., by="PrimaryID", all.x=TRUE) } %>% relocate(all_of(">"), .before=1) %>% arrange(pvalue)})

  return(cntr)
}

## prepare the additional gene info tab panel
.gene_browser_info_tab <- function(id, x, y, covar) {
     message("calculating info")
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

## prepare the UI of the gene browser
geneBrowserUI <- function(cntr_titles, all_covars, default_covar, but) {

  ui <- fluidPage(
    theme = bs_theme(bootswatch = "united"),
    fluidRow(titlePanel(h1("Gene browser")), class="bg-primary"),
    fluidRow(HTML("<hr>")),
    fluidRow(
      column(2, selectInput("contrast", label = "Dataset", choices = cntr_titles),
        HTML(paste("Click on the", but, "buttons to view an expression profile"))
             ),
      column(10, dataTableOutput("result_tbl"))
      ),
    sidebarLayout(
    sidebarPanel(
       fluidRow(selectInput("covarName", "X covariate", all_covars, selected=default_covar, width="100%")),
       fluidRow(selectInput("groupBy", "Group by", c("N/A", all_covars), selected="N/A", width="100%")),
       fluidRow(selectInput("colorBy", "Color by", c("N/A", all_covars), selected="N/A", width="100%")),
       fluidRow(selectInput("symbolBy", "Symbol by", c("N/A", all_covars), selected="N/A", width="100%")),
       fluidRow(textOutput("addInfo")),
       fluidRow(verbatimTextOutput("geneData"))
    ),

    mainPanel(
      plotOutput("countsplot")
      
    ))
  )

 

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
  but <- as.character(actionButton("go_%s", label=" > ", 
                                   onclick='Shiny.onInputChange(\"select_button\",  this.id)',  class = "btn-primary btn-sm"))
  if(is.null(annot)) {
    message(" * Loading Annotation (consider using the annot option to speed this up)")
    annot  <- get_annot(pip)
  }

  # prepare the contrast tables
  if(is.null(cntr)) {
    message(" * Loading contrasts... (consider using the cntr option to speed this up)")
    cntr <- get_contrasts(pip)
  }

  cntr <- .gene_browser_prep_res(cntr, but, annot)
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

  thematic_shiny(font="auto")

  ## prepare the UI
  ui <- geneBrowserUI(cntr_titles, all_covars, default_covar, but) 

  ## prepare the server
  server <- function(input, output, session) {
    output$id_summary <- renderText({
      .cntr <- input$contrast
      sprintf("Contrast is %s", .cntr)
    })

    output$result_tbl <- DT::renderDataTable({
      id <- showNotification("Rendering table...", duration = NULL, closeButton = FALSE)
      on.exit(removeNotification(id), add = TRUE)
      datatable(cntr[[ input$contrast ]], escape=FALSE, selection='none', options=list(pageLength=5)) %>%
        formatSignif(columns=intersect(colnames(cntr[[ input$contrast ]]), 
                                       c("baseMean", "log2FoldChange", "pvalue", "padj")), digits=2)
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

    output$addInfo <- renderText({
      req(input$select_button)
      id <- gsub("go_", "", input$select_button)
      ret <- .gene_browser_info_tab(id, covar[[input$covarName]], rld[id, ])
      return(ret)
    })

    output$countsplot <- renderPlot({
      req(input$select_button)
      id <- gsub("go_", "", input$select_button)
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
  ## weirdly, the line below is really, really slow
  #.args <- map(.args, ~ if(!is.na(.x) && .x == "N/A") { NA } else { .x })
  if(.args$groupBy == "N/A") .args$groupBy <- NA
  if(.args$colorBy == "N/A") .args$colorBy <- NA
  if(.args$symbolBy == "N/A") .args$symbolBy <- NA
  #plot_gene(pip, id, xCovar=covarName, covar=covar, rld=rld, groupBy=groupBy, colorBy=colorBy, symbolBy=symbolBy)
  do.call(plot_gene, .args)
}


