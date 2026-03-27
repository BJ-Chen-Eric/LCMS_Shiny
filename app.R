library(shiny)
library(data.table)
library(dplyr)
library(stringr)
library(ggplot2)
library(cowplot)

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
      textInput("standard_strain", "2. Set Standard Group", value = "B73"),
      
       # Statistical test method selection
      radioButtons("stat_method", "3. Select Statistical Test Method", 
                   choices = list("t-test (Student's t-test)" = "t_test",
                                  "Wilcoxon test (Mann-Whitney U)" = "wilcox_test"), 
                   selected = "t_test"),

      # Execute button
      actionButton("run_btn", "4. Start Analysis", class = "btn-primary", style = "width: 100%; margin-bottom: 20px;"),
      
      # Download button (visible only after analysis is complete)
      uiOutput("download_ui"),
      
      hr(),
      
      # Added: Detailed Data Format Instructions
      tags$div(
        tags$h4(tags$b("📝 Data Format Requirements:")),
        tags$ul(style = "padding-left: 20px; line-height: 1.6;",
          tags$li("The uploaded file must be in ", tags$b(".csv"), " format."),
          tags$li("The file ", tags$b("MUST"), " contain these exact column names (case-sensitive): ", 
                  tags$br(), tags$code("Genotype"), ", ", tags$code("Treatment"), ", ", tags$code("Rep"), "."),
          # tags$li("Metabolite compounds should start from the 6th column."),
          tags$li("Example ", tags$code("Genotype"), ": ", tags$i("B73-W0, double-W2")),
          tags$li("Example ", tags$code("Treatment"), ": ", tags$i("KODA, KODA"))
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
      
      # Display preview plot
      plotOutput("preview_plot", height = "700px")
    )
  )
)

# ==========================================
# 3. Server Logic
# ==========================================
server <- function(input, output, session) {
  
  # Store a list of all drawn ggplot objects
  stored_plots <- reactiveVal(list())
  
  # Triggered when "Start Analysis" is clicked
  observeEvent(input$run_btn, {
    
    # Fail-safe: Ensure file is uploaded
    req(input$file1)
    
    # Clear previous plots and status
    stored_plots(list())
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


        standard <- input$standard_strain
        if(!any(hplc_koda$strain %in% standard)) {stop("Standard group setting error. Please ensure the input standard name matches the strain name in the data.")}
        # Prepare an empty list for plots
        temp_plot_list <- list()
        
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
          
          # 10 times fail-safe removal of extreme values (fixed NA error)
          for(i in 1:10)  {
            sub <- sub %>% group_by(group) %>%
              mutate(mean = mean(value, na.rm = TRUE), 
                     sd = sd(value, na.rm = TRUE), 
                     z_score = (value - mean) / sd, 
                     n = sum(!is.na(value)), 
                     all_na_z = all(is.na(z_score))) %>%
              slice(
                if (isTRUE(first(all_na_z)) || coalesce(first(n), 0) <= 3) {
                  row_number()
                } else {
                  idx_to_remove <- which.max(z_score)
                  if(length(idx_to_remove) > 0) -idx_to_remove else row_number()
                }
              ) %>%
              ungroup() %>% as.data.frame()
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
              
              # Read user selected statistical method
              stat_choice <- input$stat_method
              pval <- NA # Initialize p-value
              
              # Try executing statistical test (with tryCatch to prevent crashing from extreme data)
              tryCatch({
                if (stat_choice == "t_test") {
                  # Execute t-test
                  tt <- t.test(value ~ strain, data = sub_test, var.equal = FALSE, alternative = "two.sided")
                  pval <- tt$p.value
                } else if (stat_choice == "wilcox_test") {
                  # Execute Wilcoxon test (add exact = FALSE to avoid ties warning)
                  wt <- wilcox.test(value ~ strain, data = sub_test, alternative = "two.sided", exact = FALSE)
                  pval <- wt$p.value
                }
              }, error = function(e) {
                # If test fails (e.g., a group is all NA), keep pval as NA
                pval <<- NA 
              })
              
              # If p-value is successfully calculated, assign stars
              if (!is.na(pval)) {
                stars_str <- case_when(
                  pval <= 0.001 ~ "***", 
                  pval <= 0.01 ~ "**", 
                  pval <= 0.05 ~ "*", 
                  TRUE ~ "ns"
                )
                
                # Write results to plot_df
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
          
          # Dynamic sorting logic
          final_levels <- plot_df %>%
            mutate(
              time = str_extract(time, "[0-9]+") %>% as.numeric(), 
              is_std = if_else(str_detect(group, standard), 1, 2)
            ) %>%
            arrange(time, is_std) %>%
            pull(group)
          
          plot_df$group <- factor(plot_df$group, levels = final_levels)
          if(nrow(label_df) > 0) label_df$group <- factor(label_df$group, levels = final_levels)
          
          point_DF <- sub
          point_DF$group <- factor(point_DF$group, levels = final_levels)
          
          compond_clean <- compond %>% str_remove('pmol\\|') %>% str_remove('\\/\\|g\\|FW') %>% str_remove_all(pattern = '\\*|\\:') 
          
          # --- Plotting ---
          p <- ggplot(plot_df, aes(x = group, y = mean_value, fill = time)) +
            geom_col(width = 0.6, color = 'black') +
            geom_errorbar(aes(ymin = mean_value - sem_value, ymax = mean_value + sem_value), width = 0.15, linewidth = 1.2, color = "grey65") +
            geom_point(data = point_DF, aes(x = group, y = value), inherit.aes = FALSE, position = position_jitter(width = 0.08, height = 0), size = 3.5, alpha = 1, color = "grey50") +
            scale_fill_manual(values = c('black', 'white')) +
            labs(x = 'Strain', y = "pmol / g FW", title = compond_clean, fill='') +
            theme_minimal() + gg_theme +
            theme(plot.title = element_text(hjust = 0.5), axis.line.x = element_line(color = "black", linewidth = 0.6), axis.line.y = element_line(color = "black", linewidth = 0.6), plot.margin = margin(t = 0.5, r = 1.5, b = 1.5, l = 3.5, unit = "lines"), panel.grid.major.x = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))
          
          # Add geom_text only when label_df has data (to avoid error when all NA)
          if(nrow(label_df) > 0) {
            y_limit_max <- max(max(label_df$y, na.rm=TRUE), max(point_DF$value, na.rm=TRUE)) * 1.2 + 0.5
            p <- p + geom_text(data = label_df, aes(x = group, y = y, label = label_text), inherit.aes = FALSE, size = 6, vjust = -1) +
                     scale_y_continuous(limits = c(0, y_limit_max))
          } else {
             p <- p + scale_y_continuous(limits = c(0, max(point_DF$value, na.rm=TRUE) * 1.2 + 0.5))
          }
          
          # Save the generated plot into the List
          temp_plot_list[[compond_clean]] <- p
        }
        
        # Save the List into reactiveVal for UI usage
        stored_plots(temp_plot_list)
        output$status_msg <- renderText(paste("Analysis complete! Generated", length(temp_plot_list), "plots in total."))
      })
      
    }, error = function(e) {
      output$status_msg <- renderText(paste("Analysis failed, please check the data format. Error message:", e$message))
    })
  })
  
  # --- UI update: Generate dropdown menu ---
  output$preview_selector_ui <- renderUI({
    plots <- stored_plots()
    if(length(plots) == 0) return(NULL)
    selectInput("preview_choice", "Select a metabolite plot to preview:", choices = names(plots), width = "100%")
  })
  
  # --- Render preview plot ---
  output$preview_plot <- renderPlot({
    req(input$preview_choice)
    plots <- stored_plots()
    plots[[input$preview_choice]]
  })
  
  # --- Generate download button ---
  output$download_ui <- renderUI({
    if(length(stored_plots()) > 0) {
      downloadButton("download_plots", "Download all plots (ZIP file)", class = "btn-success", style = "width: 100%;")
    }
  })
  
  # --- Download logic (Generate ZIP) ---
  output$download_plots <- downloadHandler(
    filename = function() {
      paste("LCMS_Plots_", format(Sys.time(), "%Y%m%d_%H%M"), ".zip", sep = "")
    },
    content = function(file) {
      # Create a temporary directory
      temp_dir <- tempdir()
      setwd(temp_dir)
      
      plots <- stored_plots()
      file_paths <- c()
      
      length(final_levels)
      # Loop to save all plots as png
      for (name in names(plots)) {
        file_name <- paste0(name, "_combined.png")
        ggsave(filename = file_name, plot = plots[[name]], width = 3*length(final_levels), height = 9, dpi = 300, bg = "white")
        file_paths <- c(file_paths, file_name)
      }
      
      # Pack into ZIP
      zip(zipfile = file, files = file_paths)
    }
  )
}

shinyApp(ui = ui, server = server)

# shiny::runApp('/Users/minieric/Desktop/projects/Jean_maize/Jean/LCMS_Shiny') 
# rsconnect::deployApp("/Users/minieric/Desktop/projects/Jean_maize/Jean/LCMS_Shiny")