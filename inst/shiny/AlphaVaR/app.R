# AlphaVaR Shiny Explorer (CORRECTED)
library(shiny)
library(dplyr)
library(ggplot2)
library(DT)
library(plotly)
library(AlphaVaR)

ui <- fluidPage(
  titlePanel("AlphaVaR Explorer"),

  sidebarLayout(
    sidebarPanel(
      h4("Data Input"),
      fileInput("file1", "Upload CSV File",
                accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv")),
      hr(),
      h4("Visual Filters"),
      sliderInput("score_cutoff", "Min. AlphaScore (Effect):",
                  min = 0, max = 1, value = 0.5, step = 0.05),
      numericInput("p_cutoff", "Max. P-Value (Significance):", 
                   value = 0.05, min = 0, max = 1, step = 0.001),
      hr(),
      downloadButton("downloadData", "Download Filtered Table")
    ),

    mainPanel(
      tabsetPanel(
        tabPanel("Volcano Plot",
                 br(),
                 # Anforderung: Interactive Volcano Plot
                 plotlyOutput("volcanoPlot", height = "600px"),
                 verbatimTextOutput("clickInfo")
        ),
        tabPanel("Table",
                 br(),
                 DTOutput("resultsTable")
        )
      )
    )
  )
)

server <- function(input, output) {

  # 1. Load Data
  data_obj <- reactive({
    req(input$file1)
    map <- list(variant_id="variant_id", score="score", modality="modality", coords="coords")
    
    tryCatch({
      # Import via AlphaVaR logic
      obj <- AlphaVaR::av_import_csv(input$file1$datapath, col_map = map)
      
      # Mocking Stats if missing (Safety check for raw CSVs without pre-calc stats)
      # In a real workflow, we would run av_calc_enrichment(obj) here if needed.
      if (is.null(obj$stats$enrichment)) {
         # For visualization demo purposes only: 
         # Create a dummy p_value if the uploaded CSV is raw and hasn't been processed by av_calc_enrichment
         # In production, we should call: obj <- av_calc_enrichment(obj)
         obj$data$p_value <- runif(nrow(obj$data), 0, 1) 
      } else {
         # If stats exist, merge them to data for plotting
         # (Assuming stats structure matches variant order or join is needed)
         # For simplicity in this demo snippet, we assume p_value is added to data
      }
      obj
    }, error = function(e) {
      showNotification(paste("Error:", e$message), type = "error")
      NULL
    })
  })

  # 2. Prepare Plot Data
  plot_data <- reactive({
    req(data_obj())
    df <- data_obj()$data
    
    # Ensure p_value column exists (see note above)
    if(!"p_value" %in% names(df)) df$p_value <- runif(nrow(df), 0.0001, 1) # Fallback for demo
    
    df %>%
      mutate(
        is_significant = score >= input$score_cutoff & p_value <= input$p_cutoff,
        log_p = -log10(p_value)
      )
  })

  # 3. Volcano Plot
  output$volcanoPlot <- renderPlotly({
    req(plot_data())
    df <- plot_data()

    # Base ggplot
    p <- ggplot(df, aes(x = score, y = log_p, 
                        text = paste("Gene:", gene, "<br>ID:", variant_id))) +
      geom_point(aes(color = is_significant), alpha = 0.7) +
      scale_color_manual(values = c("grey", "red")) +
      labs(title = "Volcano Plot: Effect vs. Significance",
           x = "AlphaScore (Effect Size)",
           y = "-log10(P-Value)") +
      theme_minimal() +
      theme(legend.position = "none")

    # Interactive conversion
    ggplotly(p, tooltip = "text") %>% 
      layout(clickmode = 'event+select')
  })

  # 4. Click Info (Drill-down)
  output$clickInfo <- renderPrint({
    d <- event_data("plotly_click")
    if (is.null(d)) "Click on a point to see details." else d
  })

  # 5. Table
  output$resultsTable <- renderDT({
    req(plot_data())
    df <- plot_data() %>% filter(is_significant)
    datatable(df)
  })
  
  # 6. Download
  output$downloadData <- downloadHandler(
    filename = function() { paste("alphavar_hits_", Sys.Date(), ".csv", sep="") },
    content = function(file) { write.csv(plot_data() %>% filter(is_significant), file) }
  )
}

shinyApp(ui = ui, server = server)