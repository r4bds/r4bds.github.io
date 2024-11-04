# Load the shiny library
library("shiny")

# Define the gene_dna function
gene_dna <- function(length, base_probs = c(0.25, 0.25, 0.25, 0.25)) {
  if (length %% 3 != 0) {
    stop("The argument to the parameter 'length' has to be divisible by 3")
  }
  if (!(abs(sum(base_probs) - 1) < 1e-10)) {
    stop("The sum of the vector-argument to the base_probs-paramter must be 1")
  }
  dna_vector <- sample(
    x = c("A", "T", "C", "G"),
    size = length,
    replace = TRUE,
    prob = base_probs
  )
  dna_string <- paste0(x = dna_vector, collapse = "")
  return(dna_string)
}

# Define the User Interface (Frontend)
ui <- fluidPage(
  titlePanel("Virtual Gene Generator"),
  sidebarLayout(
    sidebarPanel(
      sliderInput(inputId = "n_bases",
                  label = "Number of bases:",
                  min = 1,
                  max = 60,
                  value = 30),
      numericInput(inputId = "prob_A",
                   label = "Probability of A",
                   value = 0.25,
                   min = 0,
                   max = 1,
                   step = 0.1),
      numericInput(inputId = "prob_T",
                   label = "Probability of T",
                   value = 0.25,
                   min = 0,
                   max = 1,
                   step = 0.1),
      numericInput(inputId = "prob_C",
                   label = "Probability of C",
                   value = 0.25,
                   min = 0,
                   max = 1,
                   step = 0.1),
      numericInput(inputId = "prob_G",
                   label = "Probability of G",
                   value = 0.25,
                   min = 0,
                   max = 1,
                   step = 0.1)
    ),
    mainPanel(
      textOutput(outputId = "dna")
    )
  )
)

# Define the Server (Backend)
server <- function(input, output) {
  output$dna <- renderText({
    gene_dna(length = input$n_bases,
             base_probs = c(input$prob_A, input$prob_T, input$prob_C, input$prob_G))
  })
}

# Launch the shiny app
shinyApp(ui = ui, server = server)
