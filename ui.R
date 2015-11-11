#library(markdown)

shinyUI(navbarPage("PAT-seq explorer",
                   tabPanel("Gene expression",
                            sidebarLayout(
                              sidebarPanel(
                                radioButtons("plotType", "Plot type",
                                             c("Scatter"="p", "Line"="l")
                                )
                              ),
                              mainPanel(
                                plotOutput("plot")
                              )
                            )
                   ),
                   tabPanel("Coverage plot",
                            sidebarLayout(
                              sidebarPanel(
                                radioButtons("plotType", "Plot type",
                                             c("Scatter"="p", "Line"="l")
                                )
                              ),
                              mainPanel(
                                plotOutput("plot")
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
                                plotOutput("plot")
                              )
                            )
                   )
))