library(shiny)
library(dplyr)
library(DT)
library(httr)
library(jsonlite)
library(purrr)

# Increase the maximum upload size to 500 MB
options(shiny.maxRequestSize = 800 * 1024^2)

# Function to fetch all available panels (handles pagination)
fetch_all_panels <- function() {
  url <- "https://panelapp.genomicsengland.co.uk/api/v1/panels/"
  all_panels <- list()
  current_url <- url  # Start with the base URL
  
  repeat {
    response <- GET(current_url)
    if (status_code(response) == 200) {
      content <- fromJSON(content(response, "text"), flatten = TRUE)
      all_panels <- append(all_panels, list(content$results))  # Ensure it's added as a list
      # Check if there is a next page
      if (is.null(content$`next`)) break
      current_url <- content$`next`
    } else {
      stop("Failed to fetch panels: ", status_code(response))
    }
  }
  
  all_panels <- bind_rows(all_panels)
  return(all_panels)
}

# Function to fetch genes for a specific panel
get_genes_for_panel <- function(panel_id) {
  url <- paste0("https://panelapp.genomicsengland.co.uk/api/v1/panels/", panel_id, "/")
  response <- GET(url)
  
  if (status_code(response) == 200) {
    details <- fromJSON(content(response, "text"), flatten = TRUE)
    genes <- details$genes %>%
      mutate(
        gene_symbol = map_chr(gene_data.gene_symbol, ~paste(.x, collapse = "; "), .default = NA_character_)
      ) %>%
      pull(gene_symbol)  # Extract gene symbols
    return(genes)
  } else {
    stop("Failed to fetch genes: ", status_code(response))
  }
}

ui <- fluidPage(
  titlePanel("T2T Variant Filter Tool with PanelApp"),
  
  sidebarLayout(
    sidebarPanel(
      fileInput("csv_file", "Upload CSV File", accept = c(".csv")),
      
      # PanelApp dropdown
      selectizeInput("panel_select", 
                     "Select Panels from PanelApp:", 
                     choices = NULL, 
                     multiple = TRUE),
      
      # Filters
      textInput("af_filter", "Filter Allele Frequency (AF) ≤", value = ""),
      textInput("ac_filter", "Filter Allele Count (AC) ≤", value = ""),
      textInput("sv_filter", "Filter SV Probability ≤", value = ""),
      
      fluidRow(
        column(
          6,
          checkboxGroupInput(
            "consequence", "Select Consequence",
            choices = NULL, selected = NULL, inline = FALSE
          )
        ),
        column(
          6,
          checkboxGroupInput(
            "variant_class", "Select Variant Class",
            choices = NULL, selected = NULL, inline = FALSE
          )
        )
      ),
      
      selectizeInput("impact_filter", "Filter by IMPACT", choices = NULL, multiple = TRUE),
      selectizeInput("gene_filter", "Filter by SYMBOL", choices = NULL, multiple = FALSE),
      checkboxInput("uniquet2t", "Show Only Variants Unique to T2T:", value = FALSE),
      checkboxInput("remove_empty_genotypes", 
                    "Remove variants with empty values following 0/0", 
                    value = FALSE),
      
      # Dynamic genotype filter
      uiOutput("genotype_filter_ui"),
      
      actionButton("apply_filters", "Apply Filters"),
      downloadButton("download_filtered", "Download Filtered Data")
    ),
    
    mainPanel(
      DTOutput("filtered_table")
    )
  )
)

server <- function(input, output, session) {
  # Fetch all PanelApp panels once the app starts
  all_panels <- fetch_all_panels()
  updateSelectInput(session, "panel_select", choices = setNames(all_panels$id, all_panels$name))
  
  # Reactive for reading and processing the CSV
  raw_data <- reactive({
    req(input$csv_file)
    read.csv(input$csv_file$datapath, stringsAsFactors = FALSE)
  })
  
  # Detect genotype sample columns dynamically
  sample_columns <- reactive({
    req(raw_data())
    colnames(raw_data())[grep("^(D\\d+|DGS.*|X88549.*)", colnames(raw_data()))]
  })
  
  # Update dropdown and checkbox choices dynamically
  observeEvent(raw_data(), {
    data <- raw_data()
    updateCheckboxGroupInput(session, "consequence", 
                             choices = unique(unlist(strsplit(data$Consequence, "&"))),
                             selected = NULL)
    updateCheckboxGroupInput(session, "variant_class", 
                             choices = unique(data$VARIANT_CLASS),
                             selected = NULL)
    updateSelectizeInput(session, "gene_filter", 
                         choices = unique(data$SYMBOL), server = TRUE)
    updateSelectizeInput(session, "impact_filter", 
                         choices = unique(data$IMPACT), server = TRUE)
  })
  
  # Generate UI for genotype filter dynamically
  output$genotype_filter_ui <- renderUI({
    req(sample_columns())
    tagList(
      lapply(sample_columns(), function(sample) {
        selectInput(
          inputId = paste0("genotype_", sample),
          label = paste("Genotype for", sample),
          choices = c("All", "0/0", "0/1", "1/1"),
          selected = "All"
        )
      })
    )
  })
  
  # Reactive for filtering data
  filtered_data <- eventReactive(input$apply_filters, {
    req(raw_data())
    data <- raw_data()
    
    # Apply AF filter to specific fields
    af_columns <- c("AF_HPRC", "AF_1000G", "AF_gnomad", "AF_hgsvc3")
    af_columns <- intersect(af_columns, colnames(data)) # Ensure columns exist in the data
    if (input$af_filter != "" && length(af_columns) > 0) {
      af_threshold <- as.numeric(input$af_filter)
      if (!is.na(af_threshold)) {
        data <- data %>% filter(
          rowSums(sapply(af_columns, function(col) {
            suppressWarnings(as.numeric(data[[col]]) > af_threshold) # Convert to numeric and filter
          }), na.rm = TRUE) == 0 # Keep rows where no column exceeds the threshold
        )
      }
    }
    
    # Apply AC filter to all columns containing "AC"
    ac_columns <- grep("AC", colnames(data), value = TRUE)
    if (input$ac_filter != "" && length(ac_columns) > 0) {
      ac_threshold <- as.numeric(input$ac_filter)
      if (!is.na(ac_threshold)) {
        data <- data %>% filter(
          rowSums(sapply(ac_columns, function(col) {
            suppressWarnings(as.numeric(data[[col]]) > ac_threshold)
          }), na.rm = TRUE) == 0
        )
      }
    }
    
    # Apply SV filter 
    if (input$sv_filter != "") {
      sv_threshold <- suppressWarnings(as.numeric(input$sv_filter))
      if (!is.na(sv_threshold)) {
        data <- data %>% filter(is.na(MeanPROB) | MeanPROB > sv_threshold)
      }
    }
    
    # Inclusive filtering for Consequence and VARIANT_CLASS
    if (!is.null(input$consequence) && length(input$consequence) > 0 &&
        !is.null(input$variant_class) && length(input$variant_class) > 0) {
      data <- data %>% filter(
        sapply(Consequence, function(x) {
          any(input$consequence %in% unlist(strsplit(x, "&")))
        }) |
          VARIANT_CLASS %in% input$variant_class
      )
    } else {
      # Apply Consequence filter independently if no VARIANT_CLASS filter
      if (!is.null(input$consequence) && length(input$consequence) > 0) {
        data <- data %>% filter(
          sapply(Consequence, function(x) {
            any(input$consequence %in% unlist(strsplit(x, "&")))
          })
        )
      }
      # Apply VARIANT_CLASS filter independently if no Consequence filter
      if (!is.null(input$variant_class) && length(input$variant_class) > 0) {
        data <- data %>% filter(
          VARIANT_CLASS %in% input$variant_class
        )
      }
    }
    
    # Apply SYMBOL filter if not "All"
    if (!is.null(input$gene_filter) && input$gene_filter != "") {
      data <- data %>% filter(SYMBOL == input$gene_filter)
    }
    
    # Apply IMPACT filter
    if (!is.null(input$impact_filter) && length(input$impact_filter) > 0) {
      data <- data %>% filter(IMPACT %in% input$impact_filter)
    }
    
    # Apply UniqueT2T filter
    if (input$uniquet2t) {
      data <- data %>% filter(UniqueT2T == "Yes")
    }
    
    # Apply genotype filters
    genotype_inputs <- lapply(sample_columns(), function(sample) {
      input[[paste0("genotype_", sample)]]
    })
    names(genotype_inputs) <- sample_columns()
    
    for (sample in names(genotype_inputs)) {
      genotype <- genotype_inputs[[sample]]
      if (genotype != "All") {
        data <- data %>% filter(grepl(genotype, .data[[sample]]))
      }
    }
    
    # Reorder columns
    selected_columns <- c("SYMBOL", "VARIANT_CLASS", "HGVSc", "Consequence", 
                          "CHROM", "POS", "ID", "QUAL", "FORMAT")
    genotype_columns <- sample_columns()
    remaining_columns <- setdiff(colnames(data), c(selected_columns, genotype_columns))
    data <- data %>% select(all_of(selected_columns), all_of(genotype_columns), all_of(remaining_columns))
    
    # Dynamically remove unfilled columns
    data <- data[, colSums(!is.na(data) & data != "") > 0]
  
    # Apply PanelApp filter for multiple panels
    if (!is.null(input$panel_select) && length(input$panel_select) > 0) {
      panel_genes <- unlist(lapply(input$panel_select, get_genes_for_panel))
      data <- data %>% filter(SYMBOL %in% panel_genes)
    }
    
    # Remove rows where any genotype column has empty values (0/0:0:0... or 0/0:.:.:...)
    if (input$remove_empty_genotypes) {
      genotype_cols <- sample_columns()
      if (length(genotype_cols) > 0) {
        data <- data %>% filter(
          rowSums(sapply(genotype_cols, function(col) {
            grepl("^0/0(:0)+$", data[[col]]) | grepl("^0/0(:\\.)+$", data[[col]])  # Match empty patterns
          })) == 0  # Keep rows where no column matches any of the empty patterns
        )
      }
    }

    return(data)
  })
  
  # Render filtered table
  output$filtered_table <- renderDT({
    req(filtered_data())
    datatable(filtered_data())
  })
  
  # Download handler for filtered data
  output$download_filtered <- downloadHandler(
    filename = function() {
      paste("filtered_data", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      write.csv(filtered_data(), file, row.names = FALSE)
    },
    contentType = "text/csv"
  )
}

shinyApp(ui, server)

