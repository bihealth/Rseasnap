#' @rdname discoServer
#' @export
discoUI <- function(id, cntr_titles) {

  if(!length(cntr_titles) > 1) {
    h4("You need at least two contrasts for this plot")
  } else {

  fluidRow(
    column(width=1),
    column(width=2,
        fluidRow(selectInput(NS(id, "contrast1"), label = "Contrast 1", 
                             choices = cntr_titles, width="50%"),
        selectInput(NS(id, "contrast2"), label = "Contrast 2", 
                             choices = cntr_titles, selected=cntr_titles[2], width="50%")),
        fluidRow(checkboxInput(NS(id, "autoscale"), "Automatic scale", value=TRUE)),
        fluidRow(sliderInput(NS(id, "min"), "Min", min=-150, max=0, value=-100, width="80%")),
        fluidRow(sliderInput(NS(id, "max"), "Max", min=0, max=150, value=100, width="80%")),
        fluidRow(downloadButton(NS(id, "save"), "Save plot to PDF", class="bg-success")),
        HTML("<br/>Hover to identify genes, click to select, or click & drag to select an area<br/><br/>"),
        fluidRow(verbatimTextOutput(NS(id, "msg"))),
        fluidRow(tableOutput(NS(id, "point_id")))
    ),

    column(width=6, 
      plotOutput(NS(id, "discoplot"), 
                 hover=hoverOpts(NS(id, "plot_hover"), delay=50, delayType="throttle"),
                 click=NS(id, "plot_click"),
                 brush=NS(id, "plot_brush")

      )
    ),
    column(width=2,
      HTML("Click on the button to view an expression profile"),
      tableOutput(NS(id, "sel_genes"))
      )
  )
  }
}

.get_gene_df <- function(pid, selcols, annot=NULL) {

  if(is.null(annot)) {
    ret <- data.frame(PrimaryID=pid)
  } else {
    ret <- annot %>% filter(.data[["PrimaryID"]] %in% pid) %>%
      select(any_of(selcols))
  }
  return(ret)
}

#' Shiny Module – disco plots
#'
#' Shiny Module – disco plots
#' @param id identifier of the shiny module (character vector)
#' @param cntr list of data frames containing the contrast information.
#'        Data frames must have the columns log2FoldChange and pvalue. Rownames of
#'        the data frames should be IDs.
#' @param annot Annotation data frame. The annotation data frame must have
#'        a column named "PrimaryID" which corresponds to the rownames of the data
#'        frames in the `cntr` list.
#' @param selcols which column in the gene table when genes are selected
#'        from the plot
#' @param cntr_titles character vector containing the IDs of the contrasts
#'        (same as `names(cntr)`).
#' @return Returns a reactive expression returning the ID of the activated gene
#' @examples
#' if(interactive()) {
#'    cntr1 <- data.frame(log2FoldChange=rnorm(5000),
#'                        pvalue=runif(5000))
#'    rownames(cntr1) <- paste0("ID", 1:5000)
#'    cntr2 <- data.frame(log2FoldChange=cntr1$log2FoldChange + 
#'                                       rnorm(5000),
#'                        pvalue=runif(5000) * cntr1$pvalue)
#'    rownames(cntr2) <- paste0("ID", 1:5000)
#'    cntr <- list("Contrast 1"=cntr1, "Contrast 2"=cntr2)
#'    shinyApp(ui=fluidPage(discoUI("disco", names(cntr))),
#'             server=function(input, output, session) {
#'                discoServer("disco", cntr)
#'             })
#' }
#' @export
discoServer <- function(id, cntr, annot=NULL,
    selcols=c("PrimaryID", "ENTREZ", "SYMBOL")) {

  moduleServer(id, function(input, output, session) {
    disable("min")
    disable("max")

    disco          <- reactiveVal()
    current_genes  <- reactiveVal()
    selected_genes <- reactiveVal()

    selcols <- c("PrimaryID", "ENTREZ", "SYMBOL")

    ## enable manual color scale
    observeEvent(input$autoscale, { 
      if(input$autoscale) { 
        disable("min") 
        disable("max")
      }
      else { 
        enable("min") 
        enable("max")
      }
    })

    ## save the disco plot to a PDF file
    output$save <- downloadHandler(
      filename = function() {
        ret <- sprintf("disco_plot_%s_vs_%s.pdf", input$contrast1, input$contrast2)
        ret <- gsub("[^0-9a-zA-Z_.-]", "", ret)
        return(ret)
      },
      content = function(file) {
        c1 <- input$contrast1
        c2 <- input$contrast2
        pdf(file=file, width=8, height=8)
        if(input$autoscale) {
          g <- plot_disco(cntr[[c1]], cntr[[c2]], disco=disco())
        } else {
          g <- plot_disco(cntr[[c1]], cntr[[c2]], lower=input$min, upper=input$max, disco=disco())
        }
        print(g)
        dev.off()
      }
    )

    ## creating the actual plot
    output$discoplot <- renderPlot({
      c1 <- input$contrast1
      c2 <- input$contrast2
      disco(disco_score(cntr[[c1]], cntr[[c2]], by="PrimaryID"))

      if(input$autoscale) {
        g <- plot_disco(cntr[[c1]], cntr[[c2]], disco=disco())
      } else {
        g <- plot_disco(cntr[[c1]], cntr[[c2]], lower=input$min, upper=input$max, disco=disco())
      }
      return(g)
    }, width=600, height=600, res=90)
    
    ## React to clicking on the plot: save the current list of genes as a
    ## table on the output, adding buttons for selecting a gene
    output$sel_genes <- renderTable({
      df <- req(selected_genes())
      #df <- isolate(current_genes())
      link <- actionButton(NS(id, "gene_id_%s"), label="%s",
                           onclick=sprintf('Shiny.onInputChange(\"%s-foobar\",  this.id)', id),
                           class = "btn-primary btn-sm")
      df[["PrimaryID"]] <- sprintf(as.character(link), df$PrimaryID, df$PrimaryID)
      df
    }, sanitize.text.function=function(x) x)

    observeEvent(input$plot_click, {
      selected_genes(current_genes())
    })

    ## react to hover over points: enter the close genes into current list
    observeEvent(input$plot_hover, {
      pid <- disco() %>% nearPoints(input$plot_hover) %>% pull(.data[["PrimaryID"]])
      ret <- .get_gene_df(pid, selcols, annot)
      current_genes(ret)
    })

    ## react to points selected by brush: enter the genes into current list
    observeEvent(input$plot_brush, {
      pid <- disco() %>% brushedPoints(input$plot_brush) %>% pull(.data[["PrimaryID"]])
      ret <- .get_gene_df(pid, selcols, annot)
      selected_genes(ret)
    })

    ## enter current genes into the output table
    output$point_id <- renderTable({ 
      current_genes()
    })

    ## returns a reactive with the gene ID of the clicked gene
    gene_id <- eventReactive(input$foobar, {
      gsub(paste0("^", NS(id, "gene_id_")), "", input$foobar)
    })
  })
}

