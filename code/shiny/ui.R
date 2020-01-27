#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#
library(shiny)
source( 'shiny_input.R' )

trait_list <- c("a", "b", "c")
phenotype_list <- c("raw", "phesant")

# Define UI for app that draws a histogram ----
ui <- fluidPage(

    # App title ----
    titlePanel("Assortative Mating MR in the UKBiobank"),

    # Sidebar layout with input and output definitions ----
    sidebarLayout(

        # Sidebar panel for inputs ----
        sidebarPanel(

            # Input: Select trait ---
            selectInput(inputId = "trait",
                        label = "Select a trait",
                        choices = trait_list),

            # Input: Select phenotype ---
            radioButtons(inputId = "phenotype",
                        label = "Select a phenotype",
                        choices = phenotype_list),



            # Input: Slider for the number of bins ----
            sliderInput(inputId = "bins",
                        label = "Number of bins:",
                        min = 1,
                        max = 50,
                        value = 30)

        ),

        # Main panel for displaying outputs ----
        mainPanel(

            # Output: HTML table with requested number of observations ----
            #tableOutput("view"),
            fluidRow(
                splitLayout(cellWidths = c("50%", "50%"), tableOutput("view"), tableOutput("view2"))
            ),

            fluidRow(
                splitLayout(cellWidths = c("50%", "50%"), plotOutput("distPlot"), plotOutput("simple_xy"))
            ),


            # Output: Histogram ----
            #plotOutput(outputId = "distPlot"),

            # Output: Histogram ----
            #plotOutput(outputId = "simple_xy"),






        ),


    )

)
