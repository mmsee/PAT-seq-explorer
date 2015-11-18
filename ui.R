#library(markdown)

shinyUI(navbarPage("PAT-seq explorer",
                   tabPanel("Genome Coverage",
                            sidebarLayout(
                              sidebarPanel(
                                
                                uiOutput("select_file_path"), 

                                conditionalPanel(
                                  condition = "input.merge == false",
                                  uiOutput("bam_files")                                   
                                ),
                                conditionalPanel(
                                  condition = "input.merge == true",
                                  uiOutput("select_group")                                
                                ),
                                checkboxInput("merge", label = "Combine samples", value = F),
                                uiOutput("gene_list")
                                
                                
                              ),
                              mainPanel(
                                plotOutput("igv_plot"),
                                checkboxInput("spa", label = "Show the Poly (A) tail", value = T),
                                checkboxInput("all_reads", label = "Include reads that do not have a 
                  poly (A)-tail", value = F),
                                sliderInput("al_length", label= 'Aligned reads length range', min=0, max=400,
                                            value =c(0,400)),
                                sliderInput("ad_slider", label= "Number of sequenced adpater bases", min=0, max=23,
                                            value =0, step = 1,ticks = TRUE, 
                                            sep = ",")
 
                              )
                            )
                   ),
                   tabPanel("Gene Expression",                            
                              mainPanel(
                                plotOutput("gene_expression_plot")                                
                            )
                   ),
                   tabPanel("Poly (A)-tail cumulative distribution",
                            sidebarLayout(
                              sidebarPanel(
                                downloadButton("downloadPlot", label = "Download EPS")
                              ),
                              mainPanel(
                                plotOutput('scp_plot'),
                                checkboxInput("legend", label = "Display a legend on the plot", value = T),
                                sliderInput("xslider", label= 'x axis slider', min=0, max=400,
                                            value =c(0, 300), step = 25,ticks = TRUE, 
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
                                plotOutput("pilup_plot")

                              )
                            )
                   )
))