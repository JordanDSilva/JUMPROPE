library(shiny)
library(magicaxis)
library(Rfits)
library(Rwcs)

# Define UI for slider demo app ----
ui <- fluidPage(
  
  # App title ----
  titlePanel("Sliders"),
  
  # Sidebar layout with input and output definitions ----
  fluidPage(
    
    textInput("VID", "type VID",
              value = c("")),
    
    # Sidebar to demonstrate various slider options ----
    sidebarPanel(
      sliderInput("red_channel", "Red channel:",
                  min = -6, max = 0,
                  value = c(-5, -1), step = 0.0001),
      sliderInput("green_channel", "Green channel:",
                  min = -6, max = 0,
                  value = c(-5, -1), step = 0.0001),
      sliderInput("blue_channel", "Blue channel:",
                  min = -6, max = 0,
                  value = c(-5, -1), step = 0.0001),
      width = 15
    ),
    # Main panel for displaying outputs ----
    fillPage(
      plotOutput("distPlot")
    )
  )
)

# Define server logic for slider examples ----
server <- function(input, output) {
  
  # Reactive expression to create data frame of all input values ----
  sliderValues <- reactive({
    
    data.frame(
      Name = c("Integer"),
      Value = as.character(c(input$integer)),
      stringsAsFactors = FALSE)
    
  })
  
  data <- reactive({
    
    rgb_dir = paste0("/Volumes/JWST")
    
    list(
      R = Rfits_read(input$Rimage, pointer = F),
      G = Rfits_read(input$Gimage, pointer = F),
      B = Rfits_read(input$Bimage, pointer = F)
    )
  })
  
  red_range <- reactive({
    c(input$red_channel[1], input$red_channel[2])
  })
  green_range <- reactive({
    c(input$green_channel[1], input$green_channel[2])
  })
  blue_range <- reactive({
    c(input$blue_channel[1], input$blue_channel[2])
  })
  
  output$distPlot <- renderPlot({
    par(mfrow = c(1,1), mar = rep(0,4), oma = rep(0,4))
    
    # plot(
    #   data()$R$image,
    #   locut = c(red_range()[1]),
    #   hicut = c(red_range()[2])
    # )
    Rwcs_imageRGB(
      data()$R$image[,],
      data()$G$image[,],
      data()$B$image[,],

      locut = 10^c(red_range()[1], green_range()[1], blue_range()[1]),
      hicut = 10^c(red_range()[2], green_range()[2], blue_range()[2]),
      
      type = "num",

      sparse = 10,

      decorate = F
    )
  })
  
}

shinyApp(ui = ui, server = server)
