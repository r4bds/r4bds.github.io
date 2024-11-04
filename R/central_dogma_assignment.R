# Load the shiny library --------------------------------------------------
library("shiny")



# Load the library containing the `card()`-function -----------------------
library("bslib")



# Define app functions ----------------------------------------------------
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
transcribe_dna <- function(dna){
  rna <- gsub(
    pattern = "T",
    replacement = "U",
    x = dna)
  return(rna)
}
translate_rna <- function(rna){
  if( is.null(rna) || rna == "" ){ return("") }
  l <- nchar(x = rna)
  firsts <- seq(
    from = 1,
    to = l,
    by = 3)
  lasts <- seq(
    from = 3,
    to = l,
    by = 3)
  codons <- substring(
    text = rna,
    first = firsts,
    last = lasts)
  codon_table <- c(
    "UUU" = "F", "UCU" = "S", "UAU" = "Y", "UGU" = "C",
    "UUC" = "F", "UCC" = "S", "UAC" = "Y", "UGC" = "C",
    "UUA" = "L", "UCA" = "S", "UAA" = "_", "UGA" = "_",
    "UUG" = "L", "UCG" = "S", "UAG" = "_", "UGG" = "W",
    "CUU" = "L", "CCU" = "P", "CAU" = "H", "CGU" = "R",
    "CUC" = "L", "CCC" = "P", "CAC" = "H", "CGC" = "R",
    "CUA" = "L", "CCA" = "P", "CAA" = "Q", "CGA" = "R",
    "CUG" = "L", "CCG" = "P", "CAG" = "Q", "CGG" = "R",
    "AUU" = "I", "ACU" = "T", "AAU" = "N", "AGU" = "S",
    "AUC" = "I", "ACC" = "T", "AAC" = "N", "AGC" = "S",
    "AUA" = "I", "ACA" = "T", "AAA" = "K", "AGA" = "R",
    "AUG" = "M", "ACG" = "T", "AAG" = "K", "AGG" = "R",
    "GUU" = "V", "GCU" = "A", "GAU" = "D", "GGU" = "G",
    "GUC" = "V", "GCC" = "A", "GAC" = "D", "GGC" = "G",
    "GUA" = "V", "GCA" = "A", "GAA" = "E", "GGA" = "G",
    "GUG" = "V", "GCG" = "A", "GAG" = "E", "GGG" = "G")
  protein <- paste0(
    x = codon_table[codons],
    collapse = "")
  return(protein)
}
base_freqs <- function(dna){
  if (is.null(dna) || dna == "" ){
    return( data.frame(dna_vec = factor(c("A", "C", "G", "T")),
                       Freq = c(0, 0, 0, 0)) ) }
  dna_vec <- strsplit(x = dna,
                      split = "")
  base_counts <- table(dna_vec)
  return( as.data.frame.table(base_counts) )
}



# Define the User Interface (Frontend) ------------------------------------
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
      helpText("This Shiny app emulates the central dogma of biology:",
               "DNA => RNA => Protein")
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
    )),
  layout_columns(
    col_widths = 12,
    card(
      card_header("Virtual RNA polymerase"),
      textInput(inputId = "gene_dna",
                label = "Copy/paste the virtual gene to here",
                value = "",
                width = "100%")
    )),
  layout_columns(
    col_widths = 12,
    card(
      card_header("Virtual RNA polymerase output"),
      mainPanel(
        verbatimTextOutput(outputId = "rna")
      )
    )),
  layout_columns(
    col_widths = 12,
    card(
      card_header("Virtual Ribosome"),
      textInput(inputId = "gene_rna",
                label = "Copy/paste the virtual RNA to here",
                value = "",
                width = "100%")
    )),
  layout_columns(
    col_widths = 12,
    card(
      card_header("Virtual RNA polymerase output"),
      mainPanel(
        verbatimTextOutput(outputId = "protein")
      )
    )),
  layout_columns(
    col_widths = 12,
    card(
      card_header("Analyse Nucleotide Frequencies"),
      mainPanel(
        plotOutput(outputId = "freq_analysis_plot")
      )
    ))
)



# Define the Server (Backend) ---------------------------------------------
server <- function(input, output) {
  output$dna <- renderText({
    gene_dna(length = input$n_bases,
             base_probs = c(input$prob_A, input$prob_T, input$prob_C, input$prob_G))
  })
  output$rna <- renderText({
    transcribe_dna(dna = input$gene_dna)
  })
  output$protein <- renderText({
    translate_rna(rna = input$gene_rna)
  })
  output$freq_analysis_plot <- renderPlot({
    input$gene_dna |> 
      base_freqs() |> 
    ggplot2::ggplot(
      mapping = ggplot2::aes(
        x = dna_vec,
        y = Freq / sum(Freq),
        fill = dna_vec,
        label = stringr::str_c(round(Freq / sum(Freq) * 100,
                                     digits = 2), "%"))) +
      ggplot2::geom_col(colour = "black") +
      ggplot2::geom_hline(yintercept = 0) +
      ggplot2::geom_text(nudge_y = 0.05) +
      ggplot2::scale_y_continuous(limits = c(0, 1)) +
      ggplot2::theme_minimal() +
      ggplot2::theme(legend.position = "none") +
      ggplot2::labs(x = "Nucleotide",
                    y = "Relative Frequency")
  })
}



# Launch the shiny app ----------------------------------------------------
shinyApp(ui = ui, server = server)
