#' Visualization of SColoc results
#'
#' @param adj_mat Adjacent matrix
#' @param meta_data Optional, must contain id column match with adj_mat col/rownames
#' @param directed Generate a directed graph? (Default is false)
#'
#' @return shiny app
#' @export
#'
#' @import reshape2
#' @import magrittr
#' @import shiny
#' @import dplyr
#' @import ggpubr
#' @import visNetwork
#'
#' @examples scoloc_vis(scoloc_test)
scoloc_vis <- function(adj_mat, meta_data=NULL, directed=F) {
  mst_cons_am <- adj_mat
  mst_cons_node <- data.frame(id=rownames(mst_cons_am), label=rownames(mst_cons_am))

  if (!directed) mst_cons_am[upper.tri(mst_cons_am, diag = T)] <- NA

  mst_cons_am <- data.frame(id=rownames(mst_cons_am), mst_cons_am, check.names=F)
  mst_cons_edge <- reshape2::melt(mst_cons_am) %>% na.omit() %>% magrittr::set_colnames(c('from', 'to', 'value'))

  app <- list(
    ui=fluidPage(
      sidebarLayout(
        sidebarPanel(
          sliderInput('edge_val', 'Edge Value Cutoff',
                      min=round(range(mst_cons_edge$value)[1], 2), max=round(range(mst_cons_edge$value)[2], 2), value=round(range(mst_cons_edge$value)[1], 2), step=diff(round(range(mst_cons_edge$value), 2))/100),
          selectInput(inputId='node_col', label='Color', choices=c('None', colnames(meta_data)), selected='None'),
          selectInput(inputId='node_size', label='Size', choices=c('None', colnames(meta_data)), selected='None'),
          checkboxInput(inputId='smooth', label='Smooth', value=FALSE),
          checkboxInput(inputId='physics', label='Physics', value=FALSE),
          numericInput('mass', 'Mass', value=.5, step=.01),
          sliderInput('fontsize', 'FontSize', value=15, min=5, max=25),
          tags$hr(),
          actionButton('StopID', 'Stop')
        ),
        mainPanel(visNetworkOutput("network", height = '800px'))
        )
      ),
    server=function(input, output) {

      observeEvent(input$StopID, {
        stopApp()
      })

      output$network <- renderVisNetwork({
        if (!is.null(meta_data)) {
          if (input$node_col=='color') {
            col_colmn <- as.character(input$node_col)
            col_df <- data.frame(id=meta_data$id, color=meta_data[, col_colmn])
            mst_cons_node <- dplyr::left_join(mst_cons_node, col_df, by = "id") %>% data.frame
          } else if (input$node_col!='None') {
            col_colmn <- as.character(input$node_col)
            col_df <- data.frame(id=meta_data$id, col_=meta_data[, col_colmn])
            node_cols <- ggpubr::get_palette('Set1', length(unique(col_df$col_)))
            names(node_cols) <- unique(col_df$col_)
            col_df$color <- node_cols[col_df$col_]
            mst_cons_node <- dplyr::left_join(mst_cons_node, col_df, by = "id") %>% data.frame
          }
          if (input$node_size!='None') {
            size_colmn <- as.character(input$node_size)
            size_df <- data.frame(id=meta_data$id, value=as.numeric(meta_data[, size_colmn])/sum(as.numeric(meta_data[, size_colmn])))
            mst_cons_node <- dplyr::left_join(mst_cons_node, size_df, by = "id") %>% data.frame
          }
        }
        mst_cons_edge <- mst_cons_edge[mst_cons_edge$value > input$edge_val, ]
        visNetwork::visNetwork(mst_cons_node, mst_cons_edge) %>%
          visNetwork::visNodes(mass=input$mass, size=15, font = list(size=input$fontsize)) %>%
          visNetwork::visEdges(smooth=input$smooth, physics=input$physics,
                               arrows=list(to=list(enabled=directed, scaleFactor=.2)),
                               color=list(color='rgba(132, 132, 132, .5)'))
      })
    }
  )
  shiny::runApp(app)
}

