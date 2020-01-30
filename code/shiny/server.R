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

source('shiny_input.R')



server <- function(input, output) {

    trait_descript <- reactive(input$trait)
    i <- reactive(which(traits[["shiny_description"]]==trait_descript())) # i represents the trait
    test_vars <- reactive(expand_grid(tibble(household_GWAS_result = household_GWAS_result, i = i(), grouping_var = grouping_var), exposure_sex = sex))
    mr_results <- reactive(purrr::pmap(test_vars(), household_MR))
    bin_input <- reactive(tibble(household_mr_output = mr_results(), i = i(), grouping_var = test_vars()$grouping_var))
    bin_plots <- reactive(purrr::pmap(bin_input(), bin_plot))




    output$trait_description <-renderTable({
        output <- shiny_data[[i()]]$trait_info
        colnames(output) <- "Description"
        output
    }, rownames = TRUE)

    output$view <- renderTable({
        #trait_descript <- input$trait
        #j <- which(traits$shiny_description==trait_descript)
        shiny_data[[i()]]$trait_info
    })

    output$mr_summary <- DT::renderDataTable({
        DT::datatable(rbind(mr_results()[[1]]$full_MR_summary, mr_results()[[2]]$full_MR_summary),
                      options = list(scrollX = TRUE,
                                     bFilter=0,                                    # global search box on/off
                                     bInfo=0                                      # information on/off (how many records filtered, etc
                                     ))

    })

    output$male_female_age <- DT::renderDataTable({

        output <- mr_results()[[1]]$bin_summary
        DT::datatable(output,
                      options = list(scrollX = TRUE,
                                     bFilter=0,                                   # global search box on/off
                                     bInfo=0                                      # information on/off (how many records filtered, etc
                      ))


    })
    output$female_male_age <- DT::renderDataTable({
        output <- mr_results()[[2]]$bin_summary
        DT::datatable(output,
                      options = list(scrollX = TRUE,
                                     bFilter=0,                                   # global search box on/off
                                     bInfo=0                                      # information on/off (how many records filtered, etc
                      ))

    })
    output$male_female_tt <- DT::renderDataTable({
        output <- mr_results()[[3]]$bin_summary
        DT::datatable(output,
                      options = list(scrollX = TRUE,
                                     bFilter=0,                                   # global search box on/off
                                     bInfo=0                                      # information on/off (how many records filtered, etc
                      ))

    })
    output$female_male_tt <- DT::renderDataTable({
        output <- mr_results()[[4]]$bin_summary
        DT::datatable(output,
                      options = list(scrollX = TRUE,
                                     bFilter=0,                                   # global search box on/off
                                     bInfo=0                                      # information on/off (how many records filtered, etc
                      ))

    })

    output$male_female_age_fig <- renderPlot({

        bin_plots()[[1]]
    })
    output$female_male_age_fig <- renderPlot({

        bin_plots()[[2]]
    })
    output$male_female_tt_fig <- renderPlot({

        bin_plots()[[3]]
    })
    output$female_male_tt_fig <- renderPlot({

        bin_plots()[[4]]
    })

    output$male_female_mr <- renderPlot({
        output <- mr_results()[[1]]$mr_plot
        output
    })
    output$female_male_mr <- renderPlot({
        output <- mr_results()[[2]]$mr_plot
        output
    })


}

