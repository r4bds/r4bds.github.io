# Load the shiny library
library("shiny")

# Load the library containing the `card()`-function
library("bslib")

# Define the gene_dna function
gene_dna <- function(length, base_probs = c(0.25, 0.25, 0.25, 0.25)) {
  if (length %% 3 != 0) {
    stop("The argument to the parameter 'length' has to be divisible by 3")
  }
  if (!(abs(sum(base_probs) - 1) < 1e-10)) {
    stop("The sum of the vector-argument to the base_probs-parameter must be 1")
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
ui <- page_fluid(
  layout_columns(
    col_widths = 12,
    card(
      titlePanel("Virtual Central Dogma"),
      style = "background-color: #f0f0f0; padding: 15px;"
    )),
  layout_columns(
    col_widths = 12,
    card(
      titlePanel("About"),
      helpText("Describe what your app does...")
    )),
  layout_columns(
    col_widths = 12,
    card(
      card_header("Virtual Gene Generator"),
      sliderInput(inputId = "n_bases",
                  label = "Number of bases:",
                  min = 1,
                  max = 60,
                  value = 30,
                  width = "100%"),
      layout_columns(
        col_widths = c(3, 3, 3, 3),
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
    ))),
  layout_columns(
    col_widths = 12,
    card(
      card_header("Virtual Gene output"),
      mainPanel(
        verbatimTextOutput(outputId = "dna")
      )
    ))
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
