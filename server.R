library("Rsamtools")
library("jsonlite")
source("backend.R")
library("reshape2")
library("ggplot2")
library("genomeIntervals")
library("shinyURL")
# library('BSgenome.Hsapiens.UCSC.hg19')
setwd("/data/home/apattison/ShinyApps/andrew/")

shinyServer(function(input, output, session) {
  shinyURL.server(session)
  output$select_file_path <- renderUI({
    selectizeInput("file_path", label = ("Select a dataset"), 
                   choices = list.dirs(full.names=F, recursive =F), 
                   selected =  list.dirs(full.names=F, recursive =F)[1])   
  })
  
  
  found_gff_files <- reactive({
    find_gff_files(paste0(input$file_path))
    
  })
  
  # Being lazy and using a gff parser to to get possible inputs here
  possible_inputs <- reactive({
    if(length(found_gff_files())==0){
      return(NULL)
    }
    withProgress(message = 'Grabbing read names from the gff file',
                 detail = 'This only happens once.', value = 0,{  
                   
                   if (file.exists(paste0("/data/home/apattison/ShinyApps/dev/PAT-seq-explorer/gff_name_dump/",input$file_path, ".txt"))){
                     namess <- read.csv(paste0("/data/home/apattison/ShinyApps/dev/PAT-seq-explorer/gff_name_dump/",input$file_path, ".txt"),  stringsAsFactors=FALSE)
                     return(as.character(namess[,1]))
                   }
                   
                   if (substring(found_gff_files()[1],1,1)=="/"){
                     parsed_gff <- readGff3(paste(found_gff_files()[1]))
                     
                   }
                   else{
                     parsed_gff <- readGff3(paste("./", input$file_path,"/", found_gff_files()[1],sep=""))
                     
                   }
                   names_list <- as.character(getGffAttribute(parsed_gff, "Name"))
                   if (sum (is.na(names_list)) == length(names_list)){
                     names_list <- as.character(getGffAttribute(parsed_gff, "id"))
                   }     
                   if (sum (is.na(names_list)) == length(names_list)){
                     names_list <- as.character(getGffAttribute(parsed_gff, "ID"))                     
                   }
                   if (substring(names_list[1],1,4)=="peak"){
                     gff <- read.delim(paste(found_gff_files()[1]), header=FALSE,
                                       comment.char="",stringsAsFactors=F)
                     gff <- gff[-1,]
                     new_gff<-gff
                     rm <- regmatches(gff[,9], regexpr("Name=[^;\\s]+",gff[,9] , perl=T))
                     names_list <- gsub(x = rm, pattern = "(Name=)", replacement = "", perl = T)
                   }
                   write.csv(names_list, quote=F, row.names=F, paste0("/data/home/apattison/ShinyApps/dev/PAT-seq-explorer/gff_name_dump/",input$file_path, ".txt"))
                   return(names_list)
                 }
    )
    
  })
  
  output$gene_list <- renderUI({
    if(length(possible_inputs())==0){
      return(NULL)
    }
    
    selectizeInput("select_genes", label = ("Select a gene"),
                   choices = possible_inputs(), selected = possible_inputs()[16] , multiple = T,
                   options = list(maxOptions = 50))
    
  })
  
  select_gene_peak <- reactive({
    if(input$gene_or_peak == 1){
      input$select_genes
    } else {
      input$select_peak
    }
  })
  
  
  found_bam_files <- reactive({
    find_bam_files(paste0(input$file_path, "/"))
  })
  output$bam_files <- renderUI({
    
    if (class(found_bam_files())=='data.frame'){
      selectizeInput(inputId = "select_bam_files", label = ("Select samples"),
                     choices = found_bam_files()[[1]], 
                     selected = found_bam_files()[[1]][1], multiple =T) 
      #Goes into the data frame and gets the file paths corresponding 
      #to the selected BAM files      
    }
    else{
      selectizeInput("select_bam_files", label = ("Select samples"),
                     choices = found_bam_files(), 
                     selected = found_bam_files()[1], multiple =T)  
    }
    
  })
  processed_gff <- reactive({
    #     if(length(found_gff_files()) == 0){
    #       return (NULL)
    #     }
    withProgress(message = 'Processing the gff file, this may take a few seconds.',
                 detail = 'This only happens once.', value = 0,{   
                   if (substring(found_gff_files()[1],1,1)=="/"){
                     modify_gff_inplace(paste(found_gff_files()[1]), input$file_path)
                   }
                   else{
                     modify_gff_inplace(paste("./", input$file_path,"/", found_gff_files()[1],sep=""), input$file_path)
                   }
                   
                 }
    )
  })
  output$gff_check <- renderPrint({
    cat("The gff file may take a moment to be processed, the plot will appear when finished")
  })
  
  gffInput <- reactive({
    #     if (length(processed_gff())== 0){
    #       return(NULL)
    #     }
    filter_gff_for_rows(processed_gff(), select_gene_peak())
  })
  output$gff_rows<- renderDataTable({
    gffInput () 
  })
  
  select_group_fun <- reactive({
    count <- 1
    lapply(input$select_bam_files, 
           function(i) {          
             selectInput(paste0('snumber', i),              
                         h5(paste0('Select a group for ', i)),
                         choices = 1:length(input$select_bam_files))
           }
    )
  })
  output$select_group <- renderUI({
    select_group_fun()
  })    
  group_list <- reactive({
    res <- lapply(input$select_bam_files, 
                  function(i) { 
                    input[[paste0('snumber', i)]]
                  })
  })
  
  
  poly_a_counts<- reactive({
    if (class(found_bam_files())=='data.frame'){
      bam_files <- found_bam_files()$bam [found_bam_files()$name %in% input$select_bam_files]
    }
    else{
      bam_files <- input$select_bam_files
    }
    withProgress(message = 'Processing the bam files.',
                 detail = 'This will take longer if there are a lot of reads.', 
                 value = 0,{  
                   initial_table <- get_a_counts (input$file_path, gffInput(),bam_files,
                                                  group_list(),found_bam_files())
                   initial_table[is.na(initial_table)]<-0
                   
                   if (input$all_reads ==F){  
                     initial_table <- initial_table[initial_table$number_of_as > 0,]
                   }
                   subsetted_by_sliders <- initial_table[initial_table$number_of_ad_bases >=            
                                                           input$ad_slider&
                                                           initial_table$width >=
                                                           input$al_length[1]&
                                                           initial_table$width <=
                                                           input$al_length[2],]  
                 })
  })
  output$means_frame <- renderDataTable({
    make_means_and_meds_frame(poly_a_counts())
  })
  
  output$print_poly_a_counts <- renderDataTable({
    poly_a_counts()    
  })
  output$n_reps<- renderText ({
    input$n_replicates   
  })
  
  output$gene_info <- renderText({
    full_df <- poly_a_counts()  
    split_frame <- split(full_df, full_df$sample)
    names_string (split_frame, input$merge, input$all_reads) 
  })
  gene_expression_plot_calcs <-
    reactive({
      gene_expression_plot(poly_a_counts())
    })
  output$gene_expression_plot <-renderPlot({
    gene_expression_plot_calcs()
  }) 
  plot_calcs <- reactive({
    poly_a_plot(poly_a_counts(), input$xslider,select_gene_peak(), input$legend, 
                input$merge) 
  })
  pilup_plot_calcs <- reactive({
    pileup_plot(poly_a_counts(), input$xslider,select_gene_peak(), input$legend,
                group = F, input$order_alt, alt_cumu_dis = F,show_poly_a =F, input$poly_a_pileup )
  })
  coverage_plot_calcs <- reactive({
    
    igv_plot (poly_a_counts(), input$xslider,select_gene_peak(),input$legend,group = input$merge, 
              input$order_alt, alt_cumu_dis =F,show_poly_a =input$spa, poly_a_pileup=T,gffInput ())
  })
  
  output$igv_plot<- renderPlot({ 
    if(input$recalc == 0){
      return()
      #coverage_plot_calcs()
    }
    isolate(coverage_plot_calcs())
  })
  
  
  output$pilup_plot<- renderPlot({  
    pilup_plot_calcs()
  })
  
  
  output$scp_plot<- renderPlot({  
    plot_calcs()
  })
  #Workaround for a shiny bug thatdoesn't handle reactive plots well. 
  plot_calcs2 <- function(){
    poly_a_plot(poly_a_counts(), input$xslider,select_gene_peak(), input$legend, 
                input$merge)    
  }
  selected_plot_points <- reactive({
    input$plot_brush
  })
  output$seq_plot <- renderDataTable({
    min_points <-selected_plot_points()$xmin
    max_points <-selected_plot_points()$xmax
    if (input$alt_plot ==F){
      sequence <- poly_a_counts()$sequence[poly_a_counts()$number_of_as > min_points & poly_a_counts()$number_of_as < max_points]  
      name <- poly_a_counts()$gene_or_peak_name[poly_a_counts()$number_of_as > min_points & poly_a_counts()$number_of_as < max_points] 
      sample <- poly_a_counts()$sample[poly_a_counts()$number_of_as > min_points & poly_a_counts()$number_of_as < max_points]
      
    }
    else{
      sequence <- poly_a_counts()$sequence[poly_a_counts()$width > min_points & poly_a_counts()$width < max_points]  
      name <- poly_a_counts()$gene_or_peak_name[poly_a_counts()$width > min_points & poly_a_counts()$width < max_points] 
      sample <- poly_a_counts()$sample[poly_a_counts()$width > min_points & poly_a_counts()$width < max_points]
    }
    data.frame(sample, name, sequence)
  })
  
  output$frame_type <-renderText({
    if (input$alt_plot ==F){
      return("Reads That Fall Within The Selected Poly (A)-Tail Lengths")
    }
    else{
      return("Reads That Fall Within The Selected Read Lengths")      
    }
  }) 
  
  output$downloadPlot <- downloadHandler(
    filename = function(){
      paste(trim(select_gene_peak()), '.eps', sep='')
    },
    content = function(file){
      setEPS(width = 10)
      postscript(file)
      plot_calcs2()     
      dev.off()       
    })  
  output$dlplot1pdf <- downloadHandler(
    filename = function() { paste0('.pdf', sep='') },
    content = function(file) {
      ggsave(file,coverage_plot_calcs())
    })
  
  output$dlplot1eps <- downloadHandler(
    filename = function(){
      paste(trim(select_gene_peak()), '.eps', sep='')
    },
    content = function(file){
      setEPS(width = 10)
      postscript(file)
      coverage_plot_calcs()     
      dev.off() 
    })
  output$dlplot2pdf <- downloadHandler(
    filename = function() { paste0('.pdf', sep='') },
    content = function(file) {
      ggsave(file,gene_expression_plot(poly_a_counts()))
    })
  
  output$dlplot2eps <- downloadHandler(
    filename = function(){
      paste(trim(select_gene_peak()), '.eps', sep='')
    },
    content = function(file){
      setEPS(width = 10)
      postscript(file)
      gene_expression_plot(poly_a_counts())    
      dev.off() 
    })
})