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

            # tags$head(tags$style(type="text/css", ".dataTables_filter {display: none;    }"
            #     ))


            # Output: HTML table with requested number of observations ----
            #tableOutput("view"),
            tabsetPanel(
                id = 'summary',
                tabPanel("Trait description", tableOutput("trait_description")),
                tabPanel("MR summary", DT::dataTableOutput("mr_summary"))
            ),

            fluidRow(
                splitLayout(cellWidths = c("50%", "50%"), plotOutput("male_female_mr"), plotOutput("female_male_mr"))
            ),

            tabsetPanel(
                id = 'binned_results',
                tabPanel("Male > female, age", DT::dataTableOutput("male_female_age")),
                tabPanel("Female > male, age", DT::dataTableOutput("female_male_age")),
                tabPanel("Male > female, time together", DT::dataTableOutput("male_female_tt")),
                tabPanel("Female > male, time together", DT::dataTableOutput("female_male_tt"))
            ),

            tabsetPanel(
                id = 'binned_figures',
                tabPanel("Binned by age", splitLayout(cellWidths = c("50%", "50%"), plotOutput("male_female_age_fig"), plotOutput("female_male_age_fig"))),
                tabPanel("Binned by time together", splitLayout(cellWidths = c("50%", "50%"), plotOutput("male_female_tt_fig"), plotOutput("female_male_tt_fig")))
            ),

        ),


    )

)
