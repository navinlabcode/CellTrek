#' Convert RGB image to gray scale image
#'
#' @param img Image data with RGB-channels
#'
#' @return A gray scale image
#' @export
#'
#' @examples img_gs <- img_gs(img_raw)
img_gs <- function (img) {
  img_out <- img[,,1]+img[,,2]+img[,,3]
  img_out <- img_out/max(img_out)
  return(img_out)
}

#' Visualization of CellTrek results
#'
#' @param celltrek_df CellTrek output meta data, must contain coord_x and coord_y
#' @param img Image data
#' @param scale_fac Scale factor
#'
#' @return Shiny GUI
#' @export
#'
#' @import shiny
#' @import plotly
#' @import RColorBrewer
#'
#' @examples celltrek_vis(test_celltrek@meta.data, test_celltrek@images[[1]]@image, scale_fac=test_celltrek@images[[1]]@scale.factors$lowres)
celltrek_vis <- function(celltrek_df, img, scale_fac) {

  img_fact <- scale_fac
  img_temp <- img
  img_data <- celltrek_df
  img_data$coord_x_new=img_data$coord_y*img_fact
  img_data$coord_y_new=dim(img_temp)[1]-img_data$coord_x*img_fact
  if (!('id_new' %in% colnames(img_data))) img_data$id_new <- rownames(img_data)

  app <- list(
    ui=fluidPage(
      titlePanel('CellTrek visualization'),
      sidebarLayout(sidebarPanel(width=3,
                                 selectInput(inputId='color_inp', label='Color', choices=colnames(img_data), selected='None'),
                                 radioButtons(inputId='color_typ', label='Type', choices=c('Categorical', 'Continuous'), selected='Categorical'),
                                 selectInput(inputId='shape_inp', label='Shape', choices=c('None', colnames(img_data)), selected='None'),
                                 actionButton('Plot', 'Plot'),
                                 tags$hr(),
                                 textInput('colID', 'Add Type:'),
                                 actionButton('AddID', 'Add'),
                                 downloadButton('downloadData', 'Download'),
                                 tags$hr(),
                                 actionButton('StopID', 'Stop')),
                    mainPanel(plotlyOutput('CellTrek', height='1000px', width='1200px'),
                              dataTableOutput('Tab_temp')))
      ),

    server=function(input, output, session) {
      options(warn=-1)
      data_react <- reactiveValues()
      data_react$DF <- data.frame(id_new=character(), coord_x=numeric(), coord_y=numeric(), add_col=character())

      observeEvent(input$StopID, {
        stopApp()
      })

      observeEvent(input$Plot, {
        color_var <- isolate(input$color_inp)
        type_var <- isolate(input$color_typ)
        shape_var <- isolate(input$shape_inp)

        if (type_var=='Categorical') {
          output$CellTrek <- renderPlotly({
            img_data$color_var <- factor(img_data[, color_var])

            if (shape_var=='None') {img_data$shape_var <- ''}
            else {img_data$shape_var <- factor(img_data[, shape_var])}

            if (length(levels(img_data$color_var))<=9) {pnt_colors <- brewer.pal(length(levels(img_data$color_var)), "Set1")}
            else {pnt_colors <- colorRampPalette(brewer.pal(9, "Set1"))(length(levels(img_data$color_var)))}

            plotly::plot_ly(d=img_data, x=~coord_x_new, y=~coord_y_new, customdata=~id_new,
                            color=~color_var, type = 'scatter', mode = 'markers', text = ~color_var, symbol=~shape_var,
                            colors=pnt_colors,
                            marker = list(
                              line = list(color = 'rgb(1, 1, 1)', width = .5),
                              size=8,
                              opacity=.8)) %>%
              plotly::layout(
                xaxis = list(range = c(0, dim(img_temp)[2]), showgrid = FALSE, showline = FALSE),
                yaxis = list(range = c(0, dim(img_temp)[1]), showgrid = FALSE, showline = FALSE),
                images = list(source = plotly::raster2uri(as.raster(img_temp)), x=0, y=0,
                sizex=dim(img_temp)[2], sizey=dim(img_temp)[1],
                xref = "x", yref = "y", xanchor = "left", yanchor = "bottom", layer = "below", sizing = "stretch"))
            })}
        if (type_var=='Continuous') {
          output$CellTrek <- renderPlotly({
            img_data$color_var <- img_data[, color_var]
            if (shape_var=='None') {img_data$shape_var <- ''}
            else {img_data$shape_var <- factor(img_data[, shape_var])}
            plotly::plot_ly(d=img_data, x=~coord_x_new, y=~coord_y_new, customdata=~id_new,
                            color=~color_var, type = 'scatter', mode = 'markers', text=~color_var, symbol=~shape_var,
                            colors=c('#377EB8', 'white', '#E41A1C'),
                            marker = list(line = list(color = 'rgb(1, 1, 1)', width = .5), size=8, opacity=.8)) %>%
              plotly::layout(xaxis = list(range = c(0, dim(img_temp)[2]), showgrid = FALSE, showline = FALSE),
                             yaxis = list(range = c(0, dim(img_temp)[1]), showgrid = FALSE, showline = FALSE),
                             images = list(source = plotly::raster2uri(as.raster(img_temp)), x=0, y=0,
                             sizex=dim(img_temp)[2], sizey=dim(img_temp)[1],
                             xref = "x", yref = "y", xanchor = "left", yanchor = "bottom", layer = "below", sizing = "stretch"))
            })}
        })

      observeEvent(input$AddID, {
        tab_temp <- event_data('plotly_selected')
        if (!is.null(tab_temp)) {
          data_temp <- img_data[match(tab_temp$customdata, img_data$id_new), c('id_new', 'coord_x', 'coord_y')]
          data_temp$add_col <- input$colID
          data_react$DF <- bind_rows(data_react$DF, data_temp)
          output$Tab_temp <- renderDataTable(data_react$DF)
        }
      })
      output$downloadData <- downloadHandler(
        filename = function() {
          paste("changemyname.csv", sep = ",")
        },
        content = function(file) {
          write.csv(data_react$DF, file, row.names = F)
        }
      )
    }
    )
    shiny::runApp(app)
}

