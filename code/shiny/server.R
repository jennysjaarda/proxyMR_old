#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#
# Define server logic required to draw a histogram ----

library(shiny)

# source( 'shiny_input.R' )
# loadd(gwas_results)
# loadd(gwas_results)
# loadd(gwas_results)


server <- function(input, output) {

    output$distPlot <- renderPlot({

        x    <- faithful$waiting
        bins <- seq(min(x), max(x), length.out = input$bins + 1)

        hist(x, breaks = bins, col = "#75AADB", border = "white",
             xlab = "Waiting time to next eruption (in mins)",
             main = "Histogram of waiting times")

    })

    output$simple_xy <- renderPlot({
        plot(1:10,1:10)
    })

    output$view <- renderTable({
        head(faithful)
    })

    output$view2 <- renderTable({
        head(faithful)
    })

}
