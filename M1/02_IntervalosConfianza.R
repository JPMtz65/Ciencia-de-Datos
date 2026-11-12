library(shiny)
library(ggplot2)

ui <- fluidPage(
  titlePanel("Cobertura de Intervalos de Confianza (100 muestras)"),
  
  sidebarLayout(
    sidebarPanel(
      selectInput(
        "dist",
        "Distribución:",
        choices = c("Normal", "Poisson", "Bernoulli")
      ),
      
      conditionalPanel(
        condition = "input.dist == 'Normal'",
        selectInput(
          "param_normal",
          "Parámetro:",
          choices = c("Media", "Varianza")
        ),
        numericInput("mu", "Media verdadera (μ):", value = 0),
        numericInput("sigma", "Desviación estándar (σ):", value = 1, min = 0.01)
      ),
      
      conditionalPanel(
        condition = "input.dist == 'Poisson'",
        numericInput("lambda", "Lambda verdadero (λ):", value = 5, min = 0.1)
      ),
      
      conditionalPanel(
        condition = "input.dist == 'Bernoulli'",
        numericInput("p", "Probabilidad verdadera (p):", value = 0.5, min = 0.01, max = 0.99)
      ),
      
      numericInput("n", "Tamaño de muestra:", value = 30, min = 2),
      sliderInput("conf", "Nivel de confianza:",
                  min = 0.80, max = 0.99, value = 0.95, step = 0.01)
    ),
    
    mainPanel(
      plotOutput("icPlot", height = "600px"),
      verbatimTextOutput("coverageText")
    )
  )
)

server <- function(input, output) {
  
  ic_data <- reactive({
    
    B <- 100
    alpha <- 1 - input$conf
    
    lower <- upper <- numeric(B)
    true <- NA
    param_name <- ""
    
    for (i in 1:B) {
      
      if (input$dist == "Normal") {
        
        x <- rnorm(input$n, mean = input$mu, sd = input$sigma)
        
        if (input$param_normal == "Media") {
          
          xbar <- mean(x)
          s <- sd(x)
          tval <- qt(1 - alpha/2, df = input$n - 1)
          
          lower[i] <- xbar - tval * s / sqrt(input$n)
          upper[i] <- xbar + tval * s / sqrt(input$n)
          
          true <- input$mu
          param_name <- "μ"
          
        } else {
          
          s2 <- var(x)
          chi_low <- qchisq(1 - alpha/2, df = input$n - 1)
          chi_up  <- qchisq(alpha/2, df = input$n - 1)
          
          lower[i] <- (input$n - 1) * s2 / chi_low
          upper[i] <- (input$n - 1) * s2 / chi_up
          
          true <- input$sigma^2
          param_name <- "σ²"
        }
        
      } else if (input$dist == "Poisson") {
        
        x <- rpois(input$n, lambda = input$lambda)
        
        lower[i] <- qchisq(alpha/2, 2 * sum(x)) / (2 * input$n)
        upper[i] <- qchisq(1 - alpha/2, 2 * (sum(x) + 1)) / (2 * input$n)
        
        true <- input$lambda
        param_name <- "λ"
        
      } else {
        
        x <- rbinom(input$n, size = 1, prob = input$p)
        phat <- mean(x)
        
        z <- qnorm(1 - alpha/2)
        se <- sqrt(phat * (1 - phat) / input$n)
        
        lower[i] <- max(0, phat - z * se)
        upper[i] <- min(1, phat + z * se)
        
        true <- input$p
        param_name <- "p"
      }
    }
    
    data.frame(
      id = 1:B,
      lower = lower,
      upper = upper,
      cover = (lower <= true & upper >= true),
      true = true,
      param = param_name
    )
  })
  
  output$icPlot <- renderPlot({
    
    d <- ic_data()
    
    ggplot(d, aes(y = id)) +
      geom_segment(aes(x = lower, xend = upper, yend = id, color = cover),
                   linewidth = 1) +
      geom_vline(xintercept = unique(d$true),
                 linetype = "dashed",
                 linewidth = 1.2) +
      scale_color_manual(
        values = c("TRUE" = "steelblue", "FALSE" = "firebrick"),
        labels = c("FALSE" = "No cubre", "TRUE" = "Cubre")
      ) +
      labs(
        x = paste("Parámetro", unique(d$param)),
        y = "Intervalo",
        color = "Cobertura",
        title = "100 Intervalos de Confianza"
      ) +
      theme_minimal()
  })
  
  output$coverageText <- renderText({
    d <- ic_data()
    paste0(
      "Cobertura empírica: ",
      round(mean(d$cover), 3),
      "  |  Nivel nominal: ",
      input$conf
    )
  })
}

shinyApp(ui, server)
