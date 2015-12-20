#library(markdown)

shinyUI(navbarPage("PAT-seq explorer",
                   tabPanel("Genome Coverage",
                            sidebarLayout(
                              sidebarPanel(
                                
                                uiOutput("select_file_path"), 

                                conditionalPanel(
                                  condition = "input.merge == false",
                                  uiOutput("bam_files"),
                                  actionButton("recalc", "Re-calculate")
                                ),
                                conditionalPanel(
                                  condition = "input.merge == true",
                                  uiOutput("select_group")                                
                                ),
                                radioButtons("gene_or_peak", "Find gene or peak", choices=list("Gene"=1, "Peak"=2), selected=1, inline=T),
                                conditionalPanel(
                                    condition= "input.gene_or_peak == 1",
                                    checkboxInput("merge", label = "Combine samples", value = F),
                                    uiOutput("gene_list")
                                ),
                                conditionalPanel(
                                    condition= "input.gene_or_peak == 2",
                                    checkboxInput("merge", label = "Combine samples", value = F),
                                    textInput("select_peak", "Please enter a peak number (as PeakNUM)")
                                    )
                                
                              ),
                              mainPanel(
#                                     tags$style(type="text/css",
#                                                ".shiny-output-error { visibility: hidden; }",
#                                                ".shiny-output-error:before { visibility: hidden; }"
#                                     ),
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