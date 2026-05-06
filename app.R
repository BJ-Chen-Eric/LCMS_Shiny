library(shiny)
library(data.table)
library(dplyr)
library(stringr)
library(ggplot2)
library(cowplot)
library(openxlsx) 
library(plotly)

# ==========================================
# 1. Global Settings & Theme
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

# ==========================================
# 2. User Interface (UI)
# ==========================================
ui <- fluidPage(
  titlePanel("LC-MS Automated Analysis Tool v3.3"),
  
  sidebarLayout(
    sidebarPanel(
      # --- INSTRUCTION SECTION ---
      tags$div(
        style = "background-color: #f8f9fa; padding: 15px; border-radius: 8px; border: 1px solid #dee2e6; margin-bottom: 20px;",
        tags$h4(tags$b("📋 Quick Start Guide")),
        tags$p(tags$b("1. Required Data Columns:"), "Your CSV must contain:"),
        tags$ul(
          tags$li(tags$code("Genotype"), ": Plant variety (e.g., B73, Mutant)"),
          tags$li(tags$code("Time"), ": Development stage (e.g., W0, W2)"),
          tags$li(tags$code("Treatment"), ": Condition applied (e.g., KODA)"),
          tags$li(tags$code("Rep"), ": Biological replicate number")
        ),
        hr(),
        tags$p(tags$b("2. Comparison Logic:")),
        tags$ul(
          tags$li(tags$i("Within-strain:"), " Compares time points (W0 vs. Wx) for each genotype separately. Ideal for growth dynamics."),
          tags$li(tags$i("Cross-strain:"), " Compares genotypes (B73 vs. Mutant) at every specific time point. Ideal for trait discovery.")
        )
      ),
      
      fileInput("file1", "1. Upload CSV File", accept = c(".csv")),
      
      radioButtons("comp_mode", "2. Comparison Strategy",
                   choices = list("Within-strain (Time Course)" = "within_strain",
                                  "Cross-strain (Genotype Comparison)" = "cross_strain"),
                   selected = "within_strain"),
      
      uiOutput("mode_params_ui"),
      
      numericInput("n_points", "3. Min Points for Outlier Filter", value = 3, min = 2, max = 10),
      
      radioButtons("stat_method", "4. Statistical Test", 
                   choices = list("t-test" = "t_test", "Wilcoxon" = "wilcox_test"), 
                   selected = "t_test"),
      
      uiOutput("color_assign_title"),
      uiOutput("dynamic_colors"),
      br(),

      actionButton("run_btn", "Start Analysis", class = "btn-primary", style = "width: 100%; font-weight: bold;"),
      br(), br(),
      
      uiOutput("download_ui"),
      uiOutput("download_excel_ui")
    ),
    
    mainPanel(
      textOutput("status_msg"),
      tags$head(tags$style("#status_msg{color: #007bff; font-size: 18px; font-weight: bold; padding: 10px;}")),
      br(),
      uiOutput("preview_selector_ui"),
      plotlyOutput("preview_plot", height = "750px")  
    )
  )
)

# ==========================================
# 3. Server Logic
# ==========================================
server <- function(input, output, session) {
  
  stored_plots <- reactiveVal(list())
  stored_excel_data <- reactiveVal(list())
  
  raw_data <- reactive({
    req(input$file1)
    df <- fread(input$file1$datapath) %>% as.data.frame()
    # Case-insensitive mapping for user convenience
    colnames(df)[tolower(colnames(df)) == "genotype"] <- "strain"
    colnames(df)[tolower(colnames(df)) == "time"] <- "time"
    colnames(df)[tolower(colnames(df)) == "treatment"] <- "treat"
    colnames(df)[tolower(colnames(df)) == "rep"] <- "rep"
    return(df)
  })

  output$mode_params_ui <- renderUI({
    df <- raw_data()
    if(input$comp_mode == "within_strain") {
      time_pts <- unique(df$time)
      selectInput("baseline_val", "Select Baseline Time (Control Group)", choices = time_pts)
    } else {
      strains <- unique(df$strain)
      selectInput("baseline_val", "Select Standard Genotype (Control Group)", choices = strains)
    }
  })
  
  output$color_assign_title <- renderUI({
    if(input$comp_mode == "within_strain") tags$b("5. Assign Colors (Time Points)")
    else tags$b("5. Assign Colors (Genotypes)")
  })

  output$dynamic_colors <- renderUI({
    df <- raw_data()
    req(input$comp_mode)
    
    items <- if(input$comp_mode == "within_strain") {
      pts <- unique(df$time)
      pts[order(as.numeric(str_extract(pts, "\\d+")))]
    } else {
      unique(df$strain)
    }
    
    default_colors <- c("black", "white", "grey60", "steelblue", "firebrick", "forestgreen", "gold", "purple")
    
    lapply(seq_along(items), function(i) {
      column(6, selectInput(paste0("color_item_", items[i]), items[i], 
                            choices = colors(), selected = default_colors[(i-1) %% length(default_colors) + 1]))
    }) %>% do.call(fluidRow, .)
  })

  observeEvent(input$run_btn, {
    req(input$file1, input$baseline_val)
    
    withProgress(message = 'Crunching Data...', value = 0, {
      hplc_koda <- raw_data()
      colnames(hplc_koda) <- str_replace_all(colnames(hplc_koda), ' ', '|')
      hplc_koda$group <- paste(hplc_koda$strain, hplc_koda$time, hplc_koda$treat, sep = '_')
      
      control_ref <- input$baseline_val
      compound_cols <- grep("pmol", colnames(hplc_koda), value = TRUE)
      
      temp_plot_list <- list()
      temp_excel_list <- list()
      
      # Determine mapping for legend and colors
      items_to_color <- if(input$comp_mode == "within_strain") unique(hplc_koda$time) else unique(hplc_koda$strain)
      user_color_mapping <- sapply(items_to_color, function(x) input[[paste0("color_item_", x)]])
      names(user_color_mapping) <- items_to_color
      
      for(idx in seq_along(compound_cols)) {
        compond <- compound_cols[idx]
        incProgress(1/length(compound_cols), detail = paste("Processing:", compond))
        
        sub <- hplc_koda %>% select(strain, time, treat, rep, group, all_of(compond)) 
        colnames(sub)[6] <- 'value'
        
        # Iterative Outlier Filter
        for(i in 1:10) {
          sub <- sub %>% group_by(group) %>%
            mutate(mean = mean(value, na.rm=T), sd = sd(value, na.rm=T), n = sum(!is.na(value))) %>%
            slice(if (first(n) <= input$n_points) row_number() else -which.max((value-mean)/sd)) %>%
            ungroup() %>% as.data.frame()
        }
        
        plot_df <- sub %>% group_by(group, strain, time) %>%
          summarise(mean_val = mean(value, na.rm=T), sem_val = sd(value, na.rm=T)/sqrt(n()), .groups="drop")
        
        plot_df$p_value <- NA
        
        # Comparison logic
        if(input$comp_mode == "within_strain") {
          for(s in unique(plot_df$strain)) {
            ctrl_vals <- sub %>% filter(strain == s, time == control_ref) %>% pull(value)
            if(length(ctrl_vals) < 2) next
            for(t in setdiff(unique(plot_df$time), control_ref)) {
              exp_vals <- sub %>% filter(strain == s, time == t) %>% pull(value)
              if(length(exp_vals) < 2) next
              pval <- if(input$stat_method == "t_test") t.test(exp_vals, ctrl_vals)$p.value else wilcox.test(exp_vals, ctrl_vals, exact=F)$p.value
              plot_df[plot_df$strain == s & plot_df$time == t, "p_value"] <- round(pval, 3)
            }
          }
        } else {
          for(t in unique(plot_df$time)) {
            ctrl_vals <- sub %>% filter(time == t, strain == control_ref) %>% pull(value)
            if(length(ctrl_vals) < 2) next
            for(s in setdiff(unique(plot_df$strain), control_ref)) {
              exp_vals <- sub %>% filter(strain == s, time == t) %>% pull(value)
              if(length(exp_vals) < 2) next
              pval <- if(input$stat_method == "t_test") t.test(exp_vals, ctrl_vals)$p.value else wilcox.test(exp_vals, ctrl_vals, exact=F)$p.value
              plot_df[plot_df$strain == s & plot_df$time == t, "p_value"] <- round(pval, 3)
            }
          }
        }
        
        plot_df$stars <- case_when(plot_df$p_value <= 0.001 ~ "***", plot_df$p_value <= 0.01 ~ "**", 
                                   plot_df$p_value <= 0.05 ~ "*", is.na(plot_df$p_value) ~ "", TRUE ~ "ns")
        
        # Dynamic Fill and Sorting
        if(input$comp_mode == "within_strain") {
          final_levels <- plot_df %>% mutate(t_num = as.numeric(str_extract(time, "\\d+"))) %>%
            arrange(desc(strain == "B73"), strain, t_num) %>% pull(group) # B73 often prioritized visually
          fill_var <- "time"
        } else {
          final_levels <- plot_df %>% mutate(t_num = as.numeric(str_extract(time, "\\d+"))) %>%
            arrange(t_num, desc(strain == control_ref), strain) %>% pull(group)
          fill_var <- "strain"
        }
        
        plot_df$group <- factor(plot_df$group, levels = final_levels)
        plot_df$label_text <- paste0(ifelse(plot_df$p_value == 0, "<.001", as.character(plot_df$p_value)), "\n", plot_df$stars)
        
        # Final Plotting Logic
        compond_name <- compond %>% str_remove_all("pmol\\||/\\|g\\|FW")
        max_y_data <- max(c(plot_df$mean_val + plot_df$sem_val, sub$value), na.rm=T)
        label_df <- plot_df %>% filter(stars != "") %>% mutate(y_pos = mean_val + sem_val + (max_y_data * 0.05))
        
        p <- ggplot(plot_df, aes_string(x="group", y="mean_val", fill=fill_var)) +
          geom_col(width=0.6, color="black") +
          geom_errorbar(aes(ymin=mean_val-sem_val, ymax=mean_val+sem_val), width=0.15, linewidth=1) +
          geom_point(data=sub %>% mutate(group=factor(group, levels=final_levels)), aes_string(x="group", y="value"), 
                     position=position_jitter(0.08), color="grey40", alpha=0.6) +
          scale_y_continuous(limits = c(0, (max(c(label_df$y_pos, plot_df$mean_val + plot_df$sem_val), na.rm=T)*1.35))) +
          scale_fill_manual(values=user_color_mapping) +
          labs(title=compond_name, x="", y="pmol / g FW", fill=str_to_title(fill_var)) +
          theme_minimal() + gg_theme +
          theme(axis.text.x = element_text(angle=45, hjust=1))
        
        if(nrow(label_df) > 0) p <- p + geom_text(data=label_df, aes(x=group, y=y_pos+2, label=label_text), size=7, vjust=0, fontface="bold")

        temp_plot_list[[compond_name]] <- p
        temp_excel_list[[compond_name]] <- plot_df
      }
      stored_plots(temp_plot_list)
      stored_excel_data(temp_excel_list)
      output$status_msg <- renderText("Analysis Complete. Ready for Export.")
    })
  })
  
  # --- UI Renderers ---
  output$preview_selector_ui <- renderUI({
    req(stored_plots())
    selectInput("preview_choice", "🔬 Select Metabolite for Preview:", choices = names(stored_plots()), width = "100%")
  })
  
  output$preview_plot <- renderPlotly({
    req(input$preview_choice)
    ggplotly(stored_plots()[[input$preview_choice]])
  })

  output$download_ui <- renderUI({
    if(length(stored_plots()) > 0) downloadButton("download_plots", "🖼️ Download All Plots (ZIP)", class = "btn-success", style="width:100%; margin-bottom: 5px;")
  })

  output$download_plots <- downloadHandler(
    filename = function() { paste0("LCMS_Batch_Plots_", Sys.Date(), ".zip") },
    content = function(file) {
      tmpdir <- tempdir()
      setwd(tmpdir)
      fs <- c()
      for(n in names(stored_plots())) {
        path <- paste0(n, ".png")
        ggsave(path, stored_plots()[[n]], width=14, height=10, bg="white", dpi=300)
        fs <- c(fs, path)
      }
      zip(file, fs)
    }
  )
}

shinyApp(ui, server)

# shiny::runApp('/Users/minieric/Desktop/projects/Jean_maize/Jean/LCMS_Shiny') 
# rsconnect::deployApp("/Users/minieric/Desktop/projects/Jean_maize/Jean/LCMS_Shiny")