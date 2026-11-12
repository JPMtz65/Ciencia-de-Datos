# Cargar librerías necesarias
library(shiny)
library(ggplot2)
library(dplyr)
library(tidyr)

# Interfaz de usuario (UI)
ui <- fluidPage(
  titlePanel("Estimadores Bayesianos: Familias Conjugadas"),
  
  sidebarLayout(
    sidebarPanel(
      # Barra de selección para la familia de distribuciones
      selectInput("family", "Selecciona el Modelo / Familia Conjugada:",
                  choices = c("Proporción (Beta - Binomial)" = "beta_binom",
                              "Tasa (Gamma - Poisson)" = "gamma_pois",
                              "Media (Normal - Normal)" = "norm_norm")),
      
      # Entradas dinámicas generadas desde el servidor
      uiOutput("dynamic_inputs")
    ),
    
    mainPanel(
      # Gráfico de las distribuciones
      plotOutput("dist_plot", height = "400px"),
      br(),
      
      # Tabla de resumen de estadísticos de la posterior
      h4("Estadísticos de la Distribución Posterior", align = "center"),
      tableOutput("post_stats")
    )
  )
)

# Lógica del servidor
server <- function(input, output, session) {
  
  # 1. Renderizar inputs dinámicos dependiendo de la familia elegida
  output$dynamic_inputs <- renderUI({
    req(input$family)
    
    if (input$family == "beta_binom") {
      tagList(
        h4("Parámetros Prior (Distribución Beta)"),
        numericInput("prior_a", "Alfa (\u03B1):", value = 2, min = 0.01, step = 0.5),
        numericInput("prior_b", "Beta (\u03B2):", value = 2, min = 0.01, step = 0.5),
        hr(),
        h4("Datos de la Muestra (Verosimilitud Binomial)"),
        numericInput("data_n", "Número de ensayos (n):", value = 20, min = 1, step = 1),
        numericInput("data_y", "Número de éxitos (y):", value = 14, min = 0, step = 1)
      )
    } else if (input$family == "gamma_pois") {
      tagList(
        h4("Parámetros Prior (Distribución Gamma)"),
        numericInput("prior_shape", "Forma (\u03B1):", value = 2, min = 0.01, step = 0.5),
        numericInput("prior_rate", "Tasa (\u03B2):", value = 1, min = 0.01, step = 0.5),
        hr(),
        h4("Datos de la Muestra (Verosimilitud Poisson)"),
        numericInput("data_n_pois", "Tamaño de muestra (n):", value = 5, min = 1, step = 1),
        numericInput("data_sumx", "Suma de los conteos (\u03A3x):", value = 15, min = 0, step = 1)
      )
    } else if (input$family == "norm_norm") {
      tagList(
        h4("Parámetros Prior (Distribución Normal)"),
        numericInput("prior_mu", "Media a priori (\u03BC\u2080):", value = 0, step = 1),
        numericInput("prior_sd", "Desviación Est. a priori (\u03C3\u2080):", value = 5, min = 0.1, step = 0.5),
        hr(),
        h4("Datos de la Muestra (Verosimilitud Normal)"),
        numericInput("data_n_norm", "Tamaño de muestra (n):", value = 10, min = 1, step = 1),
        numericInput("data_mean", "Media muestral (\u0304x):", value = 2, step = 0.5),
        numericInput("data_sd", "Desviación Est. conocida (\u03C3):", value = 2, min = 0.1, step = 0.5)
      )
    }
  })
  
  # 2. Cálculos reactivos de los parámetros de las distribuciones
  calc_params <- reactive({
    req(input$family)
    
    if (input$family == "beta_binom") {
      req(input$prior_a, input$prior_b, input$data_n, input$data_y)
      # Para evitar errores: y no puede ser mayor que n
      y <- min(input$data_y, input$data_n) 
      
      post_a <- input$prior_a + y
      post_b <- input$prior_b + input$data_n - y
      
      # Rango del parámetro p (0 a 1)
      x <- seq(0.001, 0.999, length.out = 500)
      y_prior <- dbeta(x, input$prior_a, input$prior_b)
      y_post <- dbeta(x, post_a, post_b)
      # Verosimilitud normalizada a densidad (Beta(y+1, n-y+1))
      y_like <- dbeta(x, y + 1, input$data_n - y + 1)
      
      # Estadísticos posteriores
      mean_post <- post_a / (post_a + post_b)
      median_post <- qbeta(0.5, post_a, post_b)
      mode_post <- ifelse(post_a > 1 & post_b > 1, 
                          (post_a - 1) / (post_a + post_b - 2), 
                          NA) # La moda puede estar en los extremos si a,b <= 1
      
      list(x = x, prior = y_prior, likelihood = y_like, posterior = y_post,
           stats = data.frame("Media" = mean_post, "Mediana" = median_post, "Moda" = mode_post))
      
    } else if (input$family == "gamma_pois") {
      req(input$prior_shape, input$prior_rate, input$data_n_pois, input$data_sumx)
      
      post_shape <- input$prior_shape + input$data_sumx
      post_rate <- input$prior_rate + input$data_n_pois
      
      # Rango (hasta el percentil 99.9 de la posterior)
      upper_bound <- qgamma(0.999, post_shape, post_rate) * 1.5
      x <- seq(0.01, upper_bound, length.out = 500)
      
      y_prior <- dgamma(x, input$prior_shape, input$prior_rate)
      y_post <- dgamma(x, post_shape, post_rate)
      # Verosimilitud normalizada a densidad (Gamma(sumx + 1, n))
      y_like <- dgamma(x, input$data_sumx + 1, input$data_n_pois)
      
      # Estadísticos posteriores
      mean_post <- post_shape / post_rate
      median_post <- qgamma(0.5, post_shape, post_rate)
      mode_post <- ifelse(post_shape >= 1, (post_shape - 1) / post_rate, 0)
      
      list(x = x, prior = y_prior, likelihood = y_like, posterior = y_post,
           stats = data.frame("Media" = mean_post, "Mediana" = median_post, "Moda" = mode_post))
      
    } else if (input$family == "norm_norm") {
      req(input$prior_mu, input$prior_sd, input$data_n_norm, input$data_mean, input$data_sd)
      
      # Precisión (1/varianza)
      prec_prior <- 1 / (input$prior_sd^2)
      prec_data <- input$data_n_norm / (input$data_sd^2)
      prec_post <- prec_prior + prec_data
      
      post_sd <- sqrt(1 / prec_post)
      post_mu <- (prec_prior * input$prior_mu + prec_data * input$data_mean) / prec_post
      
      # Rango de visualización (centrado en un rango amplio de la prior y posterior)
      min_x <- min(input$prior_mu - 4*input$prior_sd, post_mu - 4*post_sd)
      max_x <- max(input$prior_mu + 4*input$prior_sd, post_mu + 4*post_sd)
      x <- seq(min_x, max_x, length.out = 500)
      
      y_prior <- dnorm(x, input$prior_mu, input$prior_sd)
      y_post <- dnorm(x, post_mu, post_sd)
      # Verosimilitud normalizada a densidad (N(media_muestral, sd/sqrt(n)))
      y_like <- dnorm(x, input$data_mean, input$data_sd / sqrt(input$data_n_norm))
      
      # Estadísticos posteriores (En distribución normal Media = Mediana = Moda)
      mean_post <- post_mu
      median_post <- post_mu
      mode_post <- post_mu
      
      list(x = x, prior = y_prior, likelihood = y_like, posterior = y_post,
           stats = data.frame("Media" = mean_post, "Mediana" = median_post, "Moda" = mode_post))
    }
  })
  
  # 3. Generación del gráfico
  output$dist_plot <- renderPlot({
    req(calc_params())
    data <- calc_params()
    
    # Preparar los datos en formato largo para ggplot
    df <- data.frame(
      x = data$x,
      Prior = data$prior,
      Verosimilitud = data$likelihood,
      Posterior = data$posterior
    )
    
    df_long <- pivot_longer(df, cols = c(Prior, Verosimilitud, Posterior),
                            names_to = "Distribucion", values_to = "Densidad")
    
    # Ordenar los factores para que aparezcan lógicamente en la leyenda
    df_long$Distribucion <- factor(df_long$Distribucion, 
                                   levels = c("Prior", "Verosimilitud", "Posterior"))
    
    # Colores solicitados
    paleta_colores <- c("Prior" = "#006080", 
                        "Verosimilitud" = "#800060", 
                        "Posterior" = "#608000")
    
    ggplot(df_long, aes(x = x, y = Densidad, color = Distribucion, fill = Distribucion)) +
      geom_line(size = 1.2) +
      geom_area(alpha = 0.2, position = "identity") +
      scale_color_manual(values = paleta_colores) +
      scale_fill_manual(values = paleta_colores) +
      labs(title = "Actualización Bayesiana: Prior, Verosimilitud y Posterior",
           x = "Parámetro (\u03B8)", y = "Densidad",
           color = "Distribución", fill = "Distribución") +
      theme_minimal(base_size = 15) +
      theme(legend.position = "top")
  })
  
  # 4. Tabla de estadísticos posteriores
  output$post_stats <- renderTable({
    req(calc_params())
    calc_params()$stats
  }, digits = 4, align = 'c', width = "100%")
  
}

# Ejecutar la aplicación
shinyApp(ui = ui, server = server)
