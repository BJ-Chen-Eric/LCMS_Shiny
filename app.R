library(shiny)
library(data.table)
library(dplyr)
library(stringr)
library(ggplot2)
library(cowplot)
library(openxlsx) 
library(plotly)

# ==========================================
# 1. Global Settings
# ==========================================
gg_theme <- theme(
  title = element_text(size = 36),
  axis.title.x = element_text(size = 30),
  axis.text.x = element_text(size = 28), 
  axis.title.y = element_text(size = 30),
  axis.text.y = element_text(size = 28),
  plot.subtitle = element_text(size = 28),
  plot.caption = element_text(size = 32), 
  legend.text = element_text(size = 20), 
  legend.key.size = unit(4, 'lines'), 
  legend.key.height = unit(1, "cm"),
  strip.text = element_text(size = 20), 
  strip.background = element_blank()
)

header_cleaning <- function(header_vector, pattern='/') {
  test <- str_split(header_vector, pattern = pattern) %>% do.call(what=rbind) %>% as.data.frame()
  A <- lapply(str_split(header_vector, pattern = pattern), length) %>% unlist() %>% as.data.frame() %>% 
    mutate(index=seq(1, nrow(test)))
  for(i in unique(A$.)) {
    if(identical(i, ncol(test))) {next}
    test[A[A$. %in% i, 'index'], (i+1):ncol(test)] <- ''
  }
  return(test)
}

# ==========================================
# 2. User Interface (UI)
# ==========================================
ui <- fluidPage(
  titlePanel("LC-MS Automated Analysis Tool"),
  
  sidebarLayout(
    sidebarPanel(
      # File upload
      fileInput("file1", "1. Please upload a CSV file", accept = c(".csv")),
      
      # Standard group setting
      textInput("standard_strain", "2. Set Standard Group", value = ""),
      
      # Minimum number of points for statistical analysis
      numericInput("n_points", "3. Minimum Number of Points", value = 3, min = 2, max = 10),
      
      # Statistical test method selection
      radioButtons("stat_method", "4. Select Statistical Test Method", 
                   choices = list("t-test (Student's t-test)" = "t_test",
                                  "Wilcoxon test (Mann-Whitney U)" = "wilcox_test"), 
                   selected = "t_test"),
      
      # COLOR SETTINGS ADDED HERE
      fluidRow(
        column(6, selectInput("col_time1", "Bar Color 1", choices = colors(), selected = "black")),
        column(6, selectInput("col_time2", "Bar Color 2", choices = colors(), selected = "white"))
      ),

      # Execute button
      actionButton("run_btn", "5. Start Analysis", class = "btn-primary", style = "width: 100%; margin-bottom: 20px;"),
      
      # Download buttons (visible only after analysis is complete)
      uiOutput("download_ui"),
      uiOutput("download_excel_ui"), 
      
      hr(),
      
      # Detailed Data Format Instructions
      tags$div(
        tags$h4(tags$b("📝 Data Format Requirements:")),
        tags$ul(style = "padding-left: 20px; line-height: 1.6;",
                tags$li("The uploaded file must be in ", tags$b(".csv"), " format."),
                tags$li("The file ", tags$b("MUST"), " contain these exact column names (case-sensitive): ", 
                        tags$br(), tags$code("Genotype"), ", ", tags$code("Treatment"), ", ", tags$code("Rep"), ", ", tags$code("pmol [Compound]/ g FW"), ". The Compound name can be any string but must start with 'pmol'."),
                tags$li("Example ", tags$code("Genotype"), ": ", tags$i("B73-W0, double-W2")),
                tags$li("Example ", tags$code("Treatment"), ": ", tags$i("KODA, KODA")),
                tags$li("Example ", tags$code("Compound"), " column name: ", tags$i("pmol 10-HOD/ g FW, pmol JA/ g FW"))
        )
      )
    ),
    
    mainPanel(
      # Status message
      textOutput("status_msg"),
      tags$head(tags$style("#status_msg{color: blue; font-size: 16px; font-weight: bold;}")),
      
      br(),
      # Chart preview selector
      uiOutput("preview_selector_ui"),
      
      # Display interactive preview plot
      plotlyOutput("preview_plot", height = "700px")  
    )
  )
)

# ==========================================
# 3. Server Logic
# ==========================================
server <- function(input, output, session) {
  
  # Store a list of all drawn ggplot objects
  stored_plots <- reactiveVal(list())
  # Store a list of dataframes for Excel Export
  stored_excel_data <- reactiveVal(list())
  
  # Triggered when "Start Analysis" is clicked
  observeEvent(input$run_btn, {
    
    # Fail-safe: Ensure file is uploaded
    req(input$file1)
    
    # Clear previous plots and status
    stored_plots(list())
    stored_excel_data(list())
    output$status_msg <- renderText("Reading data and processing...")
    
    tryCatch({
      # Start progress bar
      withProgress(message = 'Analysis in progress...', value = 0, {
        
        # --- Read and clean data (original logic) ---
        hplc_koda <- fread(input$file1$datapath) %>% as.data.frame()
        colnames(hplc_koda) <- str_replace_all(colnames(hplc_koda), ' ', '|')
        hplc_koda$p <- paste(hplc_koda$Genotype, hplc_koda$Treatment, hplc_koda$Rep, sep = '_') %>% 
          str_replace_all(pattern = '-', replacement = '_')
        
        meta_koda <- hplc_koda$p %>% header_cleaning('_') %>% 
          mutate(strain = V1, time = V2, treat = V3, rep= V4) %>% 
          select(-c(V1, V2, V3, V4))
        meta_koda$group <- paste(meta_koda$strain, meta_koda$time, meta_koda$treat, sep = '_')
        meta_koda$group <- str_remove(meta_koda$group, '_NA')
        cols <- c(grep(colnames(hplc_koda), pattern = 'pmol'))
        hplc_koda <- cbind(meta_koda, hplc_koda[, cols])
        
        n_points <- input$n_points
        standard <- input$standard_strain
        if(!any(hplc_koda$strain %in% standard)) {stop("Standard group setting error. Please ensure the input standard name matches the strain name in the data.")}
        
        # Prepare empty lists for plots and Excel
        temp_plot_list <- list()
        temp_excel_list <- list()
        
        # Get the columns to loop through (Compound)
        cols_to_run <- 6:(ncol(hplc_koda)-1)
        total_cols <- length(cols_to_run)
        
        # --- Start looping through each Compound ---
        for(idx in seq_along(cols_to_run)) {
          col <- cols_to_run[idx]
          compond <- colnames(hplc_koda)[col]
          
          # Update progress bar
          incProgress(1/total_cols, detail = paste("Processing:", compond))
          
          sub <- hplc_koda %>% select(strain, time, treat, rep, group, all_of(compond)) 
          colnames(sub)[6] <- 'value'
          
          # 10 times fail-safe removal of extreme values
          for(i in 1:10)  {
            sub <- sub %>% group_by(group) %>%
              mutate(mean = mean(value, na.rm = TRUE), 
                     sd = sd(value, na.rm = TRUE), 
                     z_score = (value - mean) / sd, 
                     n = sum(!is.na(value)), 
                     all_na_z = all(is.na(z_score))) %>%
              slice(
                if (
                  first(all_na_z) ||        
                  first(n) <= n_points      
                ) {
                  row_number()
                } else {
                  -which.max(z_score)       
                }
              ) %>%
              as.data.frame()
          }
          
          plot_df <- sub %>% group_by(group) %>%
            summarise(n = sum(!is.na(value)), 
                      mean_value = mean(value, na.rm = TRUE), 
                      sem_value = sd(value, na.rm = TRUE) / sqrt(n), 
                      .groups = "drop") %>%
            mutate(strain = header_cleaning(group, '_')[, 1], 
                   time = header_cleaning(group, '_')[, 2])
          
          if(any(is.na(plot_df$sem_value))) next
          
          strain_list <- plot_df$strain %>% unique()
          time_list <- plot_df$time %>% unique()
          
          plot_df$p_value <- NA
          plot_df$stars <- NA
          
          for(k in strain_list) {
            for(i in time_list) {
              if(k == standard) next
              
              # Filter data for specific time and strain
              sub_test <- sub %>% filter(time %in% i, strain %in% c(k, standard))
              
              # Fail-safe: Skip test if a group is empty or has only one unique value
              if(length(unique(sub_test$value)) < 2) next
              
              stat_choice <- input$stat_method
              pval <- NA
              
              tryCatch({
                if (stat_choice == "t_test") {
                  tt <- t.test(value ~ strain, data = sub_test, var.equal = FALSE, alternative = "two.sided")
                  pval <- tt$p.value
                } else if (stat_choice == "wilcox_test") {
                  wt <- wilcox.test(value ~ strain, data = sub_test, alternative = "two.sided", exact = FALSE)
                  pval <- wt$p.value
                }
              }, error = function(e) {
                pval <<- NA 
              })
              
              if (!is.na(pval)) {
                stars_str <- case_when(
                  pval <= 0.001 ~ "***", 
                  pval <= 0.01 ~ "**", 
                  pval <= 0.05 ~ "*", 
                  TRUE ~ "ns"
                )
                
                plot_df[plot_df$strain == k & plot_df$time == i, 'p_value'] <- round(pval, digits = 3)
                plot_df[plot_df$strain == k & plot_df$time == i, 'stars'] <- stars_str
              }
            }
          }
          
          plot_df$p_value <- ifelse(plot_df$p_value == 0, "<0.001", plot_df$p_value)
          plot_df$label_text <- paste0(plot_df$p_value, '\n', plot_df$stars)
          
          facet_top <- plot_df %>% filter(!is.na(p_value)) %>% group_by(group) %>%
            summarise(top = max(mean_value + sem_value, na.rm = TRUE), .groups = "drop") %>% 
            mutate(pad = pmax(0.075 * top, 0.5))
          
          label_df <- plot_df %>% filter(!is.na(p_value)) %>% 
            left_join(facet_top, by = "group") %>% 
            mutate(y = top + pad)
          
          final_levels <- plot_df %>%
            mutate(
              time_num = str_extract(time, "[0-9]+") %>% as.numeric(), 
              is_std = if_else(str_detect(group, standard), 1, 2)
            ) %>%
            arrange(time_num, is_std) %>%
            pull(group)
          
          plot_df$group <- factor(plot_df$group, levels = final_levels)
          if(nrow(label_df) > 0) label_df$group <- factor(label_df$group, levels = final_levels)
          
          point_DF <- sub
          point_DF$group <- factor(point_DF$group, levels = final_levels)
          
          compond_clean <- compond %>% str_remove('pmol\\|') %>% str_remove('\\/\\|g\\|FW') %>% str_remove_all(pattern = '\\*|\\:') 
          
          # ==============================================================
          # Prepare and store clean dataframe for Excel
          # ==============================================================
          excel_df <- sub %>% 
            select(strain, time, treat, rep, group, value) %>% 
            group_by(group) %>%
            mutate(
              Mean_in_group = mean(value, na.rm = TRUE),
              SD_in_group = sd(value, na.rm = TRUE)
            ) %>%
            ungroup() %>%
            left_join(plot_df %>% select(group, p_value, stars), by = "group") %>%
            arrange(match(group, final_levels), rep) 
          
          temp_excel_list[[compond_clean]] <- excel_df
          
          # --- Plotting ---
          # Update scale_fill_manual to use user inputs
          p <- ggplot(plot_df, aes(x = group, y = mean_value, fill = time)) +
            geom_col(width = 0.6, color = 'black') +
            geom_errorbar(aes(ymin = mean_value - sem_value, ymax = mean_value + sem_value), width = 0.15, linewidth = 1.2, color = "grey65") +
            geom_point(data = point_DF, aes(x = group, y = value), inherit.aes = FALSE, position = position_jitter(width = 0.08, height = 0), size = 3.5, alpha = 1, color = "grey50") +
            scale_fill_manual(values = c(input$col_time1, input$col_time2)) +
            labs(x = 'Strain', y = "pmol / g FW", title = compond_clean, fill='') +
            theme_minimal() + gg_theme +
            theme(plot.title = element_text(hjust = 0.5), axis.line.x = element_line(color = "black", linewidth = 0.6), axis.line.y = element_line(color = "black", linewidth = 0.6), plot.margin = margin(t = 0.5, r = 1.5, b = 1.5, l = 3.5, unit = "lines"), panel.grid.major.x = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))
          
          if(nrow(label_df) > 0) {
            y_limit_max <- max(max(label_df$y, na.rm=TRUE), max(point_DF$value, na.rm=TRUE)) * 1.2 + 0.5
            p <- p + geom_text(data = label_df, aes(x = group, y = y, label = label_text), inherit.aes = FALSE, size = 6, vjust = -1) +
              scale_y_continuous(limits = c(0, y_limit_max))
          } else {
            p <- p + scale_y_continuous(limits = c(0, max(point_DF$value, na.rm=TRUE) * 1.2 + 0.5))
          }
          
          temp_plot_list[[compond_clean]] <- p
        }
        
        stored_plots(temp_plot_list)
        stored_excel_data(temp_excel_list)
        output$status_msg <- renderText(paste("Analysis complete! Generated", length(temp_plot_list), "compounds in total."))
      })
      
    }, error = function(e) {
      output$status_msg <- renderText(paste("Analysis failed. Error message:", e$message))
    })
  })
  
  # --- UI update: Dropdown ---
  output$preview_selector_ui <- renderUI({
    plots <- stored_plots()
    if(length(plots) == 0) return(NULL)
    selectInput("preview_choice", "Select a metabolite plot to preview:", choices = names(plots), width = "100%")
  })
  
  # --- Render Plotly ---
  output$preview_plot <- renderPlotly({
    req(input$preview_choice)
    p <- stored_plots()[[input$preview_choice]]
    ggplotly(p, tooltip = c("x", "y", "fill")) %>% layout(hovermode = "closest")
  })
  
  # --- Download ZIP ---
  output$download_ui <- renderUI({
    if(length(stored_plots()) > 0) {
      downloadButton("download_plots", "Download all plots (ZIP file)", class = "btn-success", style = "width: 100%; margin-bottom: 10px;")
    }
  })
  
  output$download_plots <- downloadHandler(
    filename = function() { paste("LCMS_Plots_", format(Sys.time(), "%Y%m%d_%H%M"), ".zip", sep = "") },
    content = function(file) {
      temp_dir <- tempdir()
      setwd(temp_dir)
      plots <- stored_plots()
      file_paths <- c()
      for (name in names(plots)) {
        file_name <- paste0(name, "_combined.png")
        num_groups <- length(unique(plots[[name]]$data$group))
        ggsave(filename = file_name, plot = plots[[name]], width = 2 * num_groups, height = 9, dpi = 300, bg = "white")
        file_paths <- c(file_paths, file_name)
      }
      zip(zipfile = file, files = file_paths)
    }
  )
  
  # --- Download Excel ---
  output$download_excel_ui <- renderUI({
    if(length(stored_excel_data()) > 0) {
      downloadButton("download_excel", "Download Data Summary (Excel)", class = "btn-info", style = "width: 100%; margin-bottom: 10px;")
    }
  })
  
  output$download_excel <- downloadHandler(
    filename = function() { paste("LCMS_Data_Summary_", format(Sys.time(), "%Y%m%d_%H%M"), ".xlsx", sep = "") },
    content = function(file) {
      wb <- createWorkbook()
      excel_data <- stored_excel_data()
      sheet_names <- make.unique(substr(names(excel_data), 1, 31))
      for (i in seq_along(excel_data)) {
        addWorksheet(wb, sheet_names[i])
        writeData(wb, sheet_names[i], excel_data[[i]], rowNames = FALSE)
      }
      saveWorkbook(wb, file, overwrite = TRUE)
    }
  )
}

shinyApp(ui = ui, server = server)

# shiny::runApp('/Users/minieric/Desktop/projects/Jean_maize/Jean/LCMS_Shiny') 
# rsconnect::deployApp("/Users/minieric/Desktop/projects/Jean_maize/Jean/LCMS_Shiny")