#library(markdown)

shinyUI(navbarPage("PAT-seq explorer",
                   tabPanel("Genome Coverage",
                            sidebarLayout(
                              sidebarPanel(
                                
                                uiOutput("select_file_path"),
                                
                                uiOutput("gene_list"),
                                
                                checkboxInput("all_reads", label = "Include reads that do not have a 
                  poly (A)-tail", value = F),
                                
                                checkboxInput("merge", label = "Combine samples", value = F),
                                
                                conditionalPanel(
                                  condition = "input.merge == false",
                                  uiOutput("bam_files")                                   
                                )
                                
                              ),
                              mainPanel(
                                plotOutput("igv_plot")
                               
 
                              )
                            )
                   ),
                   tabPanel("Gene Expression",
                            sidebarLayout(
                              sidebarPanel(
                              ),
                              mainPanel(
                                plotOutput("gene_expression_plot")
                                
                              )
                            )
                   ),
                   tabPanel("Poly (A)-tail cumulative distribution",
                            sidebarLayout(
                              sidebarPanel(
                                radioButtons("plotType", "Plot type",
                                             c("Scatter"="p", "Line"="l")
                                )
                              ),
                              mainPanel(
                                plotOutput('scp_plot'),
                                checkboxInput("legend", label = "Display a legend on the plot", value = T),
                                sliderInput("xslider", label= 'x axis slider', min=0, max=400,
                                            value =c(0, 300), step = 25,ticks = TRUE, 
                                            sep = ","),
                                sliderInput("ad_slider", label= 'number of adapter bases (filters out reads that
                were not sequenced completely to the end of the poly (A)-tail', min=0, max=23,
                                            value =0, step = 1,ticks = TRUE, 
                                            sep = ",")
                              )
                            )
                   ),
                   tabPanel("Additional summary statistics",
                            
                            mainPanel(
                              dataTableOutput("means_frame"),
                              dataTableOutput("gff_rows"),
                              dataTableOutput('print_poly_a_counts')
                            )
                   ),
                   tabPanel("Help",
                            sidebarLayout(
                              sidebarPanel(
                                checkboxInput("order_alt", label = 
                                                "Order reads by alignment length", value = T),
                                checkboxInput("poly_a_pileup", label = 
                                                "Show raw read and poly (A)-tail lengths", value = F)
                              ),
                              mainPanel(
                                plotOutput("pilup_plot"),
                                sliderInput("al_length", label= 'alignment length range (gets reads that had a 
                genomic aligment within the selected length)', min=0, max=400,
                                            value =c(0,400))
                              )
                            )
                   )
))