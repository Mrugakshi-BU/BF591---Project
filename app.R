## Author: Mrugakshi Chidrawar
## mruga77@bu.edu
## BU BF591
## Project

# Welcome to R Shiny. All that glitters is not gold.
library(shiny)
library(ggplot2)
library(colourpicker)
library(dplyr)
library(tidyverse)
library('RColorBrewer')
library(scales)
library("gplots")
library(pheatmap)
library(ggbeeswarm)

options(shiny.maxRequestSize=30*1024^2)

# Define UI for application that draws a histogram
ui <- fluidPage(
    titlePanel("Mrugakshi Chidrawar : BF591 - Project"),
    tabsetPanel(
        tabPanel("Sample Information",
                 p('\n'),
                 sidebarLayout(
                     sidebarPanel(
                         fileInput("sample_file", label="Load Sample Data:", accept = c(".csv")), 
                         width = 3
                     ),
                     mainPanel(
                         tabsetPanel(
                             tabPanel("Summary",
                                      p('\n'),
                                      verbatimTextOutput("rows_cols"),
                                      p('\n'),
                                      textOutput('cont_summary'),
                                      p('\n'),
                                      tableOutput('sample_continuous_summary'),
                                      p('\n'),
                                      textOutput('cat_summary'),
                                      p('\n'),
                                      tableOutput('sample_categorical_summary')
                             ),
                             tabPanel("Table",
                                      p('\n'),
                                      dataTableOutput('sample_table')
                             ),
                             tabPanel("Plots",
                                      p('\n'),
                                      plotOutput('sample_plot1'),
                                      plotOutput('sample_plot2'),
                                      plotOutput('sample_plot3'),
                                      plotOutput('sample_plot4')
                             )
                         ),
                         width = 9
                     )
                 ),
        ),
        tabPanel("Counts Matrix",
                 p('\n'),
                 sidebarLayout(
                     sidebarPanel(
                         fileInput("counts_file", label="Load Counts Matrix:", accept = c(".csv")), 
                         sliderInput(inputId = "variance_percentile", min = 0, max = 1,
                                     label = "Select percentile(X) to display genes with at least X percentile of variance:", value = 0.25, step = 0.01),
                         sliderInput(inputId = "non_zero_samples", min = 0, max = 69,
                                     label = "Select X to display genes with at least X samples that are non-zero:", value = 15, step = 1),
                         submitButton("Submit", icon("check")),
                         width = 3
                     ),
                     mainPanel(
                         tabsetPanel(
                             tabPanel("Summary",
                                      p('\n'),
                                      verbatimTextOutput("counts_data_summary"),
                             ),
                             tabPanel("Diagnostic Scatter plots",
                                      p('\n'),
                                      plotOutput('scatter_plot_1'),
                                      plotOutput('scatter_plot_2'),
                             ),
                             tabPanel("Clustered Heatmap",
                                      p('\n'),
                                      plotOutput('clustered_heatmap'),
                             ),
                             tabPanel("PCA",
                                      p('\n'),
                                      sidebarPanel(
                                          selectInput("x_pc", "Select 1st Principle Component to plot on x axis:", 'PC1', choices = c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8")),
                                          selectInput("y_pc", "Select 2nd Principle Component to plot on y axis:", 'PC2', choices = c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8")),
                                          submitButton("Submit", icon("check")),
                                      ),
                                      mainPanel(
                                          p('\n'),
                                          plotOutput('pca_plot'),
                                      )
                             )
                         ),
                         width = 9
                     )
                 )
        ),
        tabPanel("Differential Expression",
                 p('\n'),
                 sidebarLayout(
                     sidebarPanel(
                         fileInput("diff_exp_file", label="Load Differential Expression results:", accept = c(".csv")), 
                         width = 3
                     ),
                     mainPanel(
                         tabsetPanel(
                             tabPanel("Differential Expression Results",
                                      p('\n'),
                                      dataTableOutput('diff_exp_table'),
                             ),
                             tabPanel("Visualization",
                                      p('\n'),
                                 sidebarPanel(
                                     p('A volcano plot can be generated with "log2 fold-change" on the x-axis and "p-adjusted" on the y-axis.'),
                                     
                                     radioButtons("x_axis", "Choose the column for the X-axis:", 'log2FoldChange',
                                                  choices = c('baseMean',
                                                              'log2FoldChange',
                                                              'lfcSE',
                                                              'stat',
                                                              'pvalue',
                                                              'padj')),
                                     
                                     radioButtons("y_axis", "Choose the column for the Y-axis:", 'padj',
                                                  choices = c('baseMean',
                                                              'log2FoldChange',
                                                              'lfcSE',
                                                              'stat',
                                                              'pvalue',
                                                              'padj')),
                                     
                                     colourInput("base_point_col", 'Base point color', "#FFCF56"),
                                     
                                     colourInput("highlight_col", 'Highlight point color', "#8DC0E3"),
                                     
                                     sliderInput(inputId = "p_magnitude", min = -200, max = 0,
                                                 label = "Select the magnitude of the p adjusted coloring:", value = -15, step = 1),
                                     
                                     submitButton("Submit", icon("check")),
                                     width = 3,
                                 ),
                                 mainPanel(
                                     tabsetPanel(
                                         tabPanel("Filtered Table",
                                                  p('\n'),
                                                  tableOutput('filtered_diff_exp_table')
                                         ),
                                         tabPanel("Plot",
                                                  p('\n'),
                                                  plotOutput('volcano_plot')
                                         ),
                                     ),
                                     width = 9,
                                 )
                             )
                         ),
                         width = 9
                     )
                 )
        ),
        tabPanel("Visualization of Individual Gene Expression(s)",
                 p('\n'),
                 sidebarLayout(
                     sidebarPanel(
                         fileInput("sample_file_2", label="Load Sample data:", accept = c()), 
                         fileInput("counts_file_2", label="Load counts matrix:", accept = c()),
                         textInput("search_gene", label="Enter a gene:", placeholder = "ENSG00000000003.10"),
                         selectInput("sample_column", "Select a categorical column from the sample data: ", 'Diagnosis', choices = c("PMI", "Age_of_death", "RIN", "Seq_reads", "Diagnosis")),
                         selectInput("plot_type", "Select a typy of plot to visualize: ", 'Bar plot', choices = c("Bar plot", "Violin plot", "Boxplot", "Beeswarm plot")),
                         submitButton("Submit", icon("check")),
                         width = 3
                     ),
                     mainPanel(
                         p('\n'),
                         plotOutput('gene_plot')
                     )
                 )
        ),
    )
)

# Define server logic required to draw a histogram
server <- function(input, output, session) {
    
    load_sample_data <- reactive({
        file <- input$sample_file
        if(is.null(file)) {     
            return(NULL) 
        }
        dataf <- read.csv(file$datapath, header = TRUE)
        dataf$Diagnosis<-as.factor(dataf$Diagnosis)
        return(dataf)
    })
    
    load_counts_data <- reactive({
        file <- input$counts_file
        if(is.null(file)) {     
            return(NULL) 
        }
        return(read.csv(file$datapath, header = TRUE))
    })
    
    load_sample_data_2 <- reactive({
        file <- input$sample_file_2
        if(is.null(file)) {     
            return(NULL) 
        }
        dataf <- read.csv(file$datapath, header = TRUE)
        dataf$Diagnosis<-as.factor(dataf$Diagnosis)
        return(dataf)
    })
    
    load_counts_data_2 <- reactive({
        file <- input$counts_file_2
        if(is.null(file)) {     
            return(NULL) 
        }
        return(read.csv(file$datapath, header = TRUE))
    })
    
    load_diff_exp_data <- reactive({
        file <- input$diff_exp_file
        if(is.null(file)) {     
            return(NULL) 
        }
        return(read.csv(file$datapath, header = TRUE))
    })

    draw_table <- function(dataf, file, flag) {
        if(is.null(file)) {     
            return(NULL) 
        }
        if(flag == 1) {
            return(dataf)
        }
        if(flag == 'continuous_columns_summary') {
            num_cols <- unlist(lapply(dataf, is.numeric))
            x <- colMeans(dataf[, num_cols], na.rm = TRUE)
            y <- sapply(dataf[, num_cols], sd, na.rm = TRUE)
            summary <- cbind(colnames(dataf[, num_cols]), x, y)
            colnames(summary) <- c("Column name", "Mean", "Std Dev")
            return(summary)
        }
        if(flag == 'categorical_columns_summary') {
            str <- paste(levels(dataf$Diagnosis), collapse = ", ")
            summary <- cbind("Diagnosis", str)
            colnames(summary) <- c("Column name", "Levels")
            return(summary)
        }
        if(flag == 3) {
            slider <- input$p_magnitude
            dataf <- dplyr::filter(dataf, dataf$padj < (1*(10**slider)))
            dataf$pvalue <- lapply(dataf$pvalue, formatC, digits = 9)
            dataf$padj <- lapply(dataf$padj, formatC, digits = 9)
            return(dataf)
        }
    }
    
    draw_plot <- function(dataf, file, flag) {
        if(is.null(file)) {     
            return(NULL) 
        }
        if(flag == 'variance' || flag == 'count') {
            var <- apply(dataf[-1],1,var)
            zeros <- rowSums(dataf==0)
            med <- apply(dataf[-1],1,median)
            
            if(flag=='variance') {
                col_label = sprintf("variance >= %s", quantile(var, input$variance_percentile))
                plot <-  ggplot(dataf, aes(y = log10(var), x = log10(med), color = (var >= quantile(var, input$variance_percentile)) & (zeros >= input$non_zero_samples))) +
                    ylab('log10(row_variance)')
            }
            if(flag=='count') {
                col_label = sprintf("count of zeros >= %s", input$non_zero_samples)
                plot <-  ggplot(dataf, aes(y = zeros, x = log10(med), color = (var >= quantile(var, input$variance_percentile)) & (zeros >= input$non_zero_samples))) +
                    ylab('Count of zeros')
            }
            plot <-  plot +
                geom_point() +
                scale_colour_manual(name = col_label, values = setNames(c('aquamarine4', 'coral2', 'grey'),c(T, F, NA))) +
                xlab('log10(row_median)') +
                theme_linedraw()
            return(plot)
        }
        if (flag == 'heatmap') {
            filtered_data <- get_counts_data_summary(load_counts_data(), input$variance_percentile, input$non_zero_samples, 'plot')
            heatmap <- pheatmap(as.matrix(filtered_data[,-1]), scale = "row", show_rownames = F,  col = brewer.pal(9, 'YlOrRd'))
            return(heatmap)
        }
        if(flag == 'filter') {
            x_name <- input$x_axis
            y_name <- input$y_axis
            slider <- input$p_magnitude
            color1 <- input$base_point_col
            color2 <- input$highlight_col
            
            y_label = sprintf("-log10(%s)", y_name)
            col_label = sprintf("padj < (1*(10^%s))", slider)
            
            vol_plot <-  ggplot(dataf, aes(x = get(x_name), y = -log10(get(y_name)), color = padj < (1*(10**slider)))) +
                geom_point() +
                scale_colour_manual(name = col_label, values = setNames(c(color1, color2, 'grey'),c(T, F, NA))) +
                xlab(x_name) +
                ylab(y_label) +
                theme_linedraw()
            
            return(vol_plot)
        }
        if(flag == 'PCA') {
            x_pc <- input$x_pc
            y_pc <- input$y_pc
            dataf[ ,c(1,2)] <- list(NULL)
            pca_results <- prcomp(scale(t(dataf)), center=FALSE, scale=FALSE)
            
            var = (pca_results$sdev)**2
            prop_var = var/sum(var)
            
            x_label = sprintf("%s (Variance: %s)", x_pc, prop_var[strtoi(substr(x_pc, 3, 3))])
            y_label = sprintf("%s (Variance: %s)", y_pc, prop_var[strtoi(substr(y_pc, 3, 3))])
            pca_plot <- as_tibble(pca_results$x) %>%
                ggplot(aes(x=get(x_pc),y=get(y_pc))) +
                geom_point() +
                xlab(x_label) +
                ylab(y_label)
            return(pca_plot)
        }
        else {
            mv_plot <- ggplot(dataf, aes(x=get(flag))) +
                geom_histogram(alpha=0.5, bins = 40, fill='maroon') +
                xlab(flag) +
                ylab("Count") +
                ggtitle(paste0("Histogram representing the distribution of ", flag))
            return(mv_plot)
        }   
    }
    
    get_rows_cols <- function(dataf, flag) {
        file <- input$sample_file
        if(is.null(file)) {     
            return(invisible(NULL))
        }
        if(flag == 'cont_summary') {
            return("Summary of the continuous variables in the dataset: ")
        }
        if(flag == 'cat_summary') {
            return("Summary of the categorical variables in the dataset: ")
        }
        
        one <- paste('Number of rows: ', nrow(dataf))
        two <- paste('Number of columns: ', ncol(dataf))
        return(paste(one, two, sep = '\n'))
    }
    
    get_counts_data_summary <- function(dataf, variance_percentile, non_zero_samples, flag) {
        file <- input$counts_file
        if(is.null(file)) {     
            return(invisible(NULL))
        }
        rows <- nrow(dataf)
        cols <- ncol(dataf)
        dataf <- dataf[rowSums(dataf == 0) >= non_zero_samples, ]
        var <- apply(dataf[-1],1,var)
        filtered_dataf <- dataf[var >= quantile(var, variance_percentile), ]
        
        if (flag == 'plot') {
            return(filtered_dataf)
        }
        
        filtered_genes <- nrow(filtered_dataf)
        filtered_percent <- (filtered_genes/rows)*100
        
        one <- paste('Total number of Genes: ', rows, '\n')
        two <- paste('Total number of Samples: ', cols, '\n')
        three <- paste('Number and % of genes passing current filter: ', filtered_genes, ' ; ', filtered_percent, '%', '\n')
        four <- paste('Number and % of genes not passing current filter: ', (rows-filtered_genes), ' ; ',  (100-filtered_percent), '%', '\n')
        return(paste(one, two, three, four, sep = '\n'))
    }
    
    draw_gene_plot <- function(sample_dataf, counts_dataf, gene, column, gene_plot_type) {
        sample_file <- input$sample_file_2
        counts_file <- input$counts_file_2
        if(is.null(sample_file) || is.null(counts_file) || gene=='') {     
            return(invisible(NULL))
        }
        gene_row <- counts_dataf[counts_dataf['X'] == gene, ]
        print(sample_dataf)
        sample_dataf$Age_of_death   <- cut(sample_dataf$Age_of_death,seq(40,100,10),right=FALSE,labels=c("Age group:40-49","Age group:50-59","Age group:60-69","Age group:70-79","Age group:80-89","Age group:90-99"))
        sample_dataf$PMI <- cut(sample_dataf$PMI,seq(0,35,7),right=FALSE,labels=c("PMI:0-6","PMI:7-14","PMI:15-21","PMI:22-28","PMI:29-35"))
        sample_dataf$RIN <- cut(sample_dataf$RIN, seq(6,10,1), right=FALSE,labels=c("RIN:6-6.9","RIN:7-7.9","RIN:8-8.9","RIN:9-10"))
        sample_dataf$Seq_reads <- cut(sample_dataf$Seq_reads, seq(38420004,167044880,20000000), right=FALSE, labels=c("Group1", "Group2", "Group3", "Group4", "Group5","Group6"))
        
        xcol <- unlist(select(sample_dataf, column))
        ycol <- unlist(gene_row[-1]) 
        gene_plot <- ggplot(mapping=aes(x=xcol,y=ycol))+xlab('Categories')+ylab('Count')
        
        if(gene_plot_type == "Bar plot"){
            gene_plot <- gene_plot + geom_bar(stat = 'identity', fill = 'blue')
        }
        if(gene_plot_type == "Boxplot"){
            gene_plot <- gene_plot + geom_boxplot(fill = 'pink', color = 'blue')
        }
        if(gene_plot_type == "Violin plot"){
            gene_plot <- gene_plot + geom_violin(fill="#E69F00", col = "blue")
        }
        if(gene_plot_type == "Beeswarm plot"){
            gene_plot <- gene_plot + geom_beeswarm(col = "blue")
        }
        return(gene_plot)
    }
    
    #tab 1 : Sample Information
    output$sample_table <- renderDataTable(draw_table(load_sample_data(), input$sample_file, 1))
    output$cont_summary <- renderText(get_rows_cols(load_sample_data(), 'cont_summary'))
    output$sample_continuous_summary <- renderTable(draw_table(load_sample_data(), input$sample_file, 'continuous_columns_summary'))
    output$cat_summary <- renderText(get_rows_cols(load_sample_data(), 'cat_summary'))
    output$sample_categorical_summary <- renderTable(draw_table(load_sample_data(), input$sample_file, 'categorical_columns_summary'))
    output$rows_cols <- renderText(get_rows_cols(load_sample_data(), 'summary'))
    output$sample_plot1 <- renderPlot(draw_plot(load_sample_data(), input$sample_file, 'PMI'))
    output$sample_plot2 <- renderPlot(draw_plot(load_sample_data(), input$sample_file, 'Age_of_death'))
    output$sample_plot3 <- renderPlot(draw_plot(load_sample_data(), input$sample_file, 'RIN'))
    output$sample_plot4 <- renderPlot(draw_plot(load_sample_data(), input$sample_file, 'Seq_reads'))
    
    #tab 2 : Counts Matrix
    output$counts_data_summary <- renderText(get_counts_data_summary(load_counts_data(), input$variance_percentile, input$non_zero_samples, 'summary'))
    output$scatter_plot_1 <- renderPlot(draw_plot(load_counts_data(), input$counts_file, 'variance'))
    output$scatter_plot_2 <- renderPlot(draw_plot(load_counts_data(), input$counts_file, 'count'))
    output$clustered_heatmap <- renderPlot(draw_plot(load_counts_data(), input$counts_file, 'heatmap'), height = 700, width = 700)
    output$pca_plot <- renderPlot(draw_plot(load_counts_data(), input$counts_file, 'PCA'))
    
    #tab3 : Differential Expression
    output$diff_exp_table <- renderDataTable(draw_table(load_diff_exp_data(), input$diff_exp_file, 1))
    output$filtered_diff_exp_table <- renderTable(draw_table(load_diff_exp_data(), input$diff_exp_file, 3))
    output$volcano_plot <- renderPlot(draw_plot(load_diff_exp_data(), input$diff_exp_file, 'filter'))
    
    #tab4 : Visualization of Individual Gene Expression(s)
    output$gene_plot <- renderPlot(draw_gene_plot(load_sample_data_2(), load_counts_data_2(), input$search_gene, input$sample_column, input$plot_type))
}

# Run the application
shinyApp(ui = ui, server = server)
