library(shiny)
library(plotly)

ui <- fluidPage(
  titlePanel("Ammonia-Ammonium Equilibrium with Multiple Samples"),
  
  sidebarLayout(
    sidebarPanel(
      numericInput("NH4_total", "Standard NH4Cl concentration (moles/L)", value = 0.075, min = 0.001, step = 0.1),
      numericInput("V_total", "Sample Volume (mL)", value = 100, min = 0.1, step = 0.1),
      selectInput("num_samples", "Number of Samples", choices = 1:6, selected = 1),
      numericInput("Temp", "Temperature (°C)", value=15, min=3, max=45, step=1),
      numericInput("Sal", "Salinity (PSU)", value=25, min=2, max=35, step=1),
      uiOutput("sample_selector"),
      uiOutput("sample_panel")
    ),
    
    mainPanel(
      tableOutput("resultsTable"),
      br(),
      plotlyOutput("NH3_NH4_plot")
    )
  )
)

server <- function(input, output, session) {
  
  # -----------------------------------------------------------
  # Keep the selector UI but make sure selection is kept valid
  # -----------------------------------------------------------
  output$sample_selector <- renderUI({
    req(input$num_samples)
    selectInput("selected_sample", "Select Sample Panel:",
                choices = paste("Sample", seq_len(as.numeric(input$num_samples))))
  })
  
  # When num_samples changes, update the selectInput choices and keep selection valid
  observeEvent(input$num_samples, {
    n <- as.numeric(input$num_samples)
    if (is.na(n) || n < 1) return()
    choices <- paste("Sample", seq_len(n))
    current <- isolate(input$selected_sample)
    new_sel <- if (!is.null(current) && current %in% choices) current else choices[1]
    updateSelectInput(session, "selected_sample", choices = choices, selected = new_sel)
  })
  
  # -----------------------------------------------------------
  # The single sample panel (only displayed for selected sample)
  # -----------------------------------------------------------
  output$sample_panel <- renderUI({
    req(input$selected_sample)
    i <- as.numeric(gsub("Sample ", "", input$selected_sample))
    
    wellPanel(
      h4(paste("Sample", i)),
      numericInput(paste0("V_added_", i), "Added Volume (mL)", value = 10, min = 0.1, step = 0.1),
      numericInput(paste0("pH_", i), "pH", value = 8, min = 1, max = 13, step = 0.1)
    )
  })
  
  # -----------------------------------------------------------
  # pKa function
  # -----------------------------------------------------------
  calc_pKa <- function(T, S) {
    -0.467 + 0.00113*S + 2887.9/(T+273.15)
  }
  
  # -----------------------------------------------------------
  # Compute sample results (mg/L) — robust to missing inputs
  # -----------------------------------------------------------
  calc_results <- reactive({
    Temp <- input$Temp
    Sal  <- input$Sal
    n    <- as.numeric(input$num_samples)
    if (is.na(n) || n < 1) return(NULL)
    
    res <- data.frame(
      Sample      = 1:n,
      NH3_mg_L    = rep(NA_real_, n),
      NH4_mg_L    = rep(NA_real_, n),
      Percent_NH3 = rep(NA_real_, n),
      pH          = rep(NA_real_, n)
    )
    
    pKa <- calc_pKa(Temp, Sal)
    C_total <- input$NH4_total  # mol/L
    
    for(i in 1:n){
      # safely fetch inputs; if missing, set NA and skip calculations
      pH_i <- input[[paste0("pH_", i)]]
      V_added_i <- input[[paste0("V_added_", i)]]
      
      if (is.null(pH_i) || length(pH_i) == 0) {
        # leave NA entries (panel not yet created / no input)
        next
      }
      
      # ensure numeric
      pH_i <- as.numeric(pH_i)
      
      # fraction and mg/L calculations
      frac_NH3 <- 1 / (1 + 10^(pKa - pH_i))
      frac_NH4 <- 1 - frac_NH3
      
      NH3_mg <- C_total * frac_NH3 * 14.01 * 1000 *input[[paste0("V_added_", i)]]/input[[paste0("V_total")]]
      NH4_mg <- C_total * frac_NH4 * 14.01 * 1000 *input[[paste0("V_added_", i)]]/input[[paste0("V_total")]]
      
      res$NH3_mg_L[i] <- NH3_mg
      res$NH4_mg_L[i] <- NH4_mg
      res$Percent_NH3[i] <- frac_NH3 * 100
      res$pH[i] <- pH_i
    }
    
    res
  })
  
  # -----------------------------------------------------------
  # Equilibrium curve (Percent NH3 vs pH)
  # -----------------------------------------------------------
  calc_curve <- reactive({
    Temp <- input$Temp
    Sal  <- input$Sal
    pKa  <- calc_pKa(Temp, Sal)
    
    pH_seq  <- seq(4, 12, length.out = 150)
    NH3_pct <- 100 / (1 + 10^(pKa - pH_seq))
    
    data.frame(pH_seq, NH3_pct)
  })
  
  # -----------------------------------------------------------
  # Table
  # -----------------------------------------------------------
  output$resultsTable <- renderTable({
    calc_results()
  }, digits = 4)
  
  # -----------------------------------------------------------
  # Dual-Axis Plot
  # -----------------------------------------------------------
  output$NH3_NH4_plot <- renderPlotly({
    df    <- calc_results()
    curve <- calc_curve()
    if (is.null(df)) return(NULL)
    
    # safe y-range: sum of NH3+NH4 for each row, then overall max
    totals <- rowSums(df[, c("NH3_mg_L", "NH4_mg_L")], na.rm = TRUE)
    ymax <- suppressWarnings(max(totals, na.rm = TRUE))
    if (!is.finite(ymax) || ymax <= 0) ymax <- 1
    
    plot_ly() %>%
      add_trace(
        data = df,
        x = ~pH,
        y = ~NH3_mg_L,
        type = "scatter",
        mode = "markers+lines",
        name = "Samples (mg/L)",
        text = ~paste(
          "Sample:", Sample,
          "<br>NH3:", round(NH3_mg_L,2), "mg/L",
          "<br>NH4+:", round(NH4_mg_L,2), "mg/L",
          "<br>% NH3:", round(Percent_NH3,1), "%"
        ),
        marker = list(
          size = 12,
          color = df$Percent_NH3,
          colorscale = "Viridis",
          showscale = TRUE
        ),
        hoverinfo = "text",
        yaxis = "y"
      ) %>%
      add_trace(
        data = curve,
        x = ~pH_seq,
        y = ~NH3_pct,
        type = "scatter",
        mode = "lines",
        name = "% NH3 (Equilibrium)",
        line = list(color = "blue", width = 4),
        yaxis = "y2"
      ) %>%
      add_trace(
        data = curve,
        x = ~pH_seq,
        y = ~(100-NH3_pct),
        type = "scatter",
        mode = "lines",
        name = "% NH4 (Equilibrium)",
        line = list(color = "lightblue", width = 4),
        yaxis = "y2"
      ) %>%
      layout(
        title = "NH₃/NH₄⁺ Samples vs Equilibrium Curve",
        xaxis = list(title = "pH"),
        yaxis = list(
          title = "NH₃ & NH₄⁺ (mg/L)",
          side = "left",
          rangemode = "tozero",
          range = c(0, ymax)
        ),
        yaxis2 = list(
          title = "% NH₃ (Equilibrium)",
          overlaying = "y",
          side = "right",
          range = c(0,100)
        ),
        legend = list(x = 0.02, y = 0.98)
      )
  })
}

shinyApp(ui, server)
