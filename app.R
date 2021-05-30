#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(shinythemes)
library(shinycssloaders)
library(shinyBS)
library(mathjaxr)



source("DDMsimulate.R")

# Define UI for application that draws a histogram
ui <- fluidPage(
    # Theme from shinythemes
    theme = shinytheme("flatly"),

    # Application title
    titlePanel("Diffusion Decision Model"),

    # Hide Warnings
    tags$style(type="text/css",
               ".shiny-output-error { visibility: hidden; }",
               ".shiny-output-error:before { visibility: hidden; }"
    ),

    # Sidebar with a slider input for number of bins
    sidebarLayout(
        sidebarPanel(
            sliderInput("nsims",
                        "Number of Simulated Paths:",
                        min = 30,
                        max = 500,
                        value = 100),
            sliderInput("drift",
                        "Drift Rate (v):",
                        min = -2,
                        max = 2,
                        step=0.1,
                        value = 0.7),
            sliderInput("boundary",
                        "Boundary Separation (a):",
                        min = 2,
                        max = 4,
                        step=0.1,
                        value = 3),
            sliderInput("bias",
                        "Start Point Bias (b):",
                        min = 0.1,
                        max = 0.9,
                        step=0.05,
                        value = 0.5),
            radioButtons("hist", " Response Time Plot",
                         c("Histogram" = TRUE,
                           "Density" = FALSE)),

            actionButton("btn", "Click for Reference Parameters"),
            bsTooltip(id="btn", title="Drift rate, v = 0.7</br> Boundary separation, a = 3 </br> Start Point Bias, b = 0.5 </br> Note: Reference plots are in grey",
                      trigger="click",placement="bottom")
        ),

        # Show a plot of the generated distribution
        mainPanel(
            tabsetPanel(
                tabPanel("Plot",shinycssloaders::withSpinner(plotOutput("distPlot"))),
                tabPanel("Information",shinycssloaders::withSpinner(uiOutput("info")))
            )
        )
   )
)

# Define server logic required to draw a histogram
server <- function(input, output) {

    output$distPlot <- renderPlot({
        # generate bins based on input$bins from ui.R
        res <- suppressWarnings(sim_paths(v=input$drift,a=input$boundary,b=input$bias,nsim=input$nsims))
        plot_ddmsim(res,histogram=input$hist)

    })
    output$info <- renderUI({
        withMathJax(helpText('
        INTRODUCTION
        $$ $$
        The purpose of this web app is to help you visualise how changing the drift, boundary, and start point parameters
        in the standard diffusion decision model affects choice probabilities and response time distributions.
        The simulated data that is generated with parameters that you choose are in blue, and a reference set of simulated data generated
        with set parameters (\\(v=0.7\\), \\( a=3\\), \\(b=0.5\\)) is in grey.
        $$ $$
        THE MODEL
        $$ $$
        The diffusion decision model (Ratcliff, 1978) assumes that 2-alternative forced choice decisions are made by
        a noisy process where evidence is accumulated over time towards one of two response boundaries. The time it takes for the accumulator
        to reach a boundary (first passage time) determines the response time and the boundary it reaches determines the response. In recent years, the diffusion model
        has gained traction as the model is able to disentangle the decision process and facilitate inference about
        the strength of evidence, response caution, bias, and non-decision processes such as motor response execution.
        $$ $$
        In this model, the accumulation of evidence is described as a diffusion process, defined by the stochastic differential equation:

        $$dX_{t} = vdt + sdW_{t}$$

        where \\(dX_{t}\\) is the change in the accumulated evidence that occurs in the small time interval,
        \\(dt\\). The change is the sum of a deterministic and a random component. The deterministic component is proportional
        to the drift rate, \\(v\\), which can be interpreted as the relative strength of evidence for the upper boundary response.
        The random component is defined by \\(dW_{t}\\), the differential of the Wiener process (or Brownian motion) and is proportional to
        a small standard deviation, \\(s\\).
        $$  $$
        In this parameterisation of the model, the lower boundary is set at 0 and the upper boundary is set at \\(a\\). The boundary separation, \\(a\\) can be
        interpreted as response caution, where a higher \\(a\\) indicates that more evidence has to be accumulated
        before a response is made. The start point, \\(z = a \\times b\\), is parameterised in terms of response bias, \\(b \\in (0,1)\\), where \\(b = 0.5\\) indicates
        that there is no bias. Finally, there is a non-decision time paramter, \\(T_{er}\\), which is the time it takes for processes outside the decision process, such as
        encoding and motor response. The total response time is given by \\(RT=T_D + T_{er}\\), where \\(T_D\\) is the decision time.
        $$ $$
        For more information, please refer to: Ratcliff, R. (1978). A theory of memory retrieval.\\(\\textit{ Psychological Review}\\), 85, 59-108.
        '))
    })
}

# Run the application
shinyApp(ui = ui, server = server)
