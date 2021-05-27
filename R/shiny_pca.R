.gethovertext <- function(covar, ids) {

  ids <- setdiff(colnames(covar), ids)
  if(length(ids) < 1) {
    ids <- colnames(covar)
  }

  tmp <- map_dfc(rlang::set_names(ids), ~ paste0(.x, ": ", covar[[.x]]))
  apply(tmp, 1, function(x) paste(x, collapse="\n"))
}

.getcol <- function(x, df) {
  if(x == "N/A" || ! x %in% colnames(df)) {
    return(NULL)
  }
  return(as.formula(paste("~", x)))
}


#' @rdname pcaServer
#' @export
pcaUI <- function(id, covar, pca_names, colorBy=NULL, symbolBy=NULL) {

  if(is.null(colnames(covar))) {
    colnames(covar) <- paste0("V", 1:ncol(covar))
  }
  interesting_covariates <- covar %>% summary_colorDF() %>%
    dplyr::filter(unique > 1 & (unique < nrow(covar) | .data[["Class"]] == '<dbl>')) %>%
    pull("Col")

  categorical_covs <- covar %>% summary_colorDF() %>%
    dplyr::filter(unique > 1 & unique < nrow(covar) & .data[["Class"]] != '<dbl>') %>%
    pull("Col")

  if(is.null(colorBy) || !colorBy %in% interesting_covariates) {
    colorBy <- interesting_covariates[1]
  }

  if(is.null(symbolBy) || !symbolBy %in% interesting_covariates) {
    symbolBy <- "N/A"
  }


    sidebarLayout(
      sidebarPanel(
          selectInput(NS(id, "color"), "Color by:",   c("N/A", interesting_covariates), selected=colorBy),
          selectInput(NS(id, "symbol"), "Symbol by:", c("N/A", categorical_covs), selected=symbolBy),
          checkboxInput(NS(id, "threeD"), "3D "),
        fluidRow(
          column(width=4, selectInput(NS(id, "x"), "X: ", choices=pca_names, selected=pca_names[1], width="100%")),
          column(width=4, selectInput(NS(id, "y"), "Y: ", choices=pca_names, selected=pca_names[2], width="100%")),
          column(width=4, selectInput(NS(id, "z"), "Z: ", choices=pca_names, selected=pca_names[3], width="100%"))
        ),
      width=3),
      mainPanel(
       plotlyOutput(NS(id, "pca_plot"), width="100%", height=600),
       width=9, useShinyjs()
      )
    )
}

#' Shiny Module – PCA plots
#'
#' Shiny Module – PCA plots
#' 
#' @param id identifier of the shiny module (character vector)
#' @param pca pca matrix – columns correspond to principal components, rows
#' to observations. Rows must be named and must correspond to the ID column 
#' of the covariate data frame. 
#' @param covar data frame containing covariates. The identifiers of the
#' samples in the covariate data frame are taken from the ID column (by
#' default, "ID").
#' @param colorBy selected covariate to use for coloring the plot
#' @param symbolBy selected covariate to use for symbols on the plot
#' @param threeD whether the plot should be three-dimensional by default
#' @param idcol name of the ID column in the covariate data frame.
#' @param pca_names names of the PCA components to include in the plot. If
#'        not specified, all components will be accessible from the interface.
#' @importFrom plotly renderPlotly plotlyOutput
#' @importFrom shiny isTruthy
#' @importFrom stats as.formula prcomp
#' @examples
#' if(interactive()) {
#'   data(iris)
#'   covar <- iris
#'  
#'   pca <- prcomp(iris[,1:4], scale.=TRUE)
#'  
#'   ui <- fluidPage(pcaUI("pca", covar, colnames(pca$x), colorBy="Species"))
#'  
#'   server <- function(input, output, session) {
#'     pcaServer("pca", pca$x, covar)
#'   }
#'  
#'   shinyApp(ui, server)
#' }
#' @export
pcaServer <- function(id, pca, covar, idcol="ID", threeD=FALSE) {
  stopifnot(is.character(id))
  stopifnot(is.data.frame(covar))
  stopifnot(is.matrix(pca) || is.data.frame(pca))
  stopifnot(is.logical(threeD))

  if(is.null(colnames(covar))) {
    colnames(covar) <- paste0("V", 1:ncol(covar))
  }

  if(is.null(colnames(pca))) {
    colnames(pca) <- paste0("PC", 1:ncol(pca))
  }

  if(any(ops <- colnames(pca) %in% colnames(covar))) {
    colnames(covar)[ops] <- paste(colnames(covar), "_covariate")
  }

  df <- cbind(covar, pca)
  if(is.null(df[["hoverinfo"]])) {
    df[["hoverinfo"]] <- .gethovertext(covar, names(covar))
  }

  moduleServer(id, function(input, output, session) {

    output$pca_plot <- renderPlotly({ 
      x <- input$x
      y <- input$y
      z <- input$z
      message(sprintf("x=%s", x))

      message("creating plot")
      symbol <- .getcol(input$symbol, df)
      color  <- .getcol(input$color, df)

      if(isTruthy(input$threeD)) {
        enable("z")
        plot_ly(data=df, type="scatter3d", x=df[[x]], y=df[[y]], z=df[[z]], 
                mode="markers", color=color, symbol=symbol,
                hoverinfo="hoverinfo")
      } else {
        disable("z")
        plot_ly(data=df, type="scatter", x=df[[x]], y=df[[y]], 
                mode="markers", color=color, symbol=symbol,
                hovertext=df[["hoverinfo"]])
      }
    })

  })
}


#' Simple PCA plot in a browser
#'
#' Simple PCA plot in a browser using shiny
#'
#' @param covar covariate file
#' @param x matrix or data frame for PCA
#' @param colorBy name of the covariate column for coloring
#' @param symbolBy name of the covariate column for chosing symbols
#' @examples
#' if(interactive()) {
#'    plot_pca_shiny(iris[,1:4], iris[,5,drop=FALSE])
#' }
#' @export
plot_pca_shiny <- function(x, covar, colorBy=NULL, symbolBy=NULL) {

  pca <- prcomp(x, scale.=TRUE)

  ui <- fluidPage(pcaUI("pca", covar, colnames(pca$x), 
                  colorBy=colorBy, symbolBy=symbolBy))

  server <- function(input, output, session) {
    pcaServer("pca", pca$x, covar)
  }

  shinyApp(ui, server)
}

