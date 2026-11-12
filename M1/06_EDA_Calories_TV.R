############ Visualización estadística de variables y pruebas ###################
# Cargamos las bibliotecas necesarias
pacman::p_load(ggplot2, gridExtra, ggcorrplot, vcd, dplyr, broom, MASS)

# Paleta de colores solicitada
custom_palette <- c("#AC9FB8", "#2E364C", "#D18E54", "#B99257", 
                    "#87A6BB", "#AE7268", "#154E52", "#63331C", "#101D29")

# Función auxiliar para extraer variables
var_extract <- function(dataframe, variable){
  return(dataframe[[variable]])
}

# 1. Función de estadísticas
statcats <- function(x){
  minimo <- round(min(x, na.rm = TRUE), 4)
  maximo <- round(max(x, na.rm = TRUE), 4)
  media <- round(mean(x, na.rm = TRUE), 4)
  mediana <- round(median(x, na.rm = TRUE), 4)
  ds <- round(sd(x, na.rm = TRUE), 4)
  
  Estadisticas <- c("Minimo", "Maximo", "Media", "Mediana", "Desv. est.")
  Valores <- c(minimo, maximo, media, mediana, ds)
  
  stats <- data.frame(Estadisticas, Valores) %>% as_tibble()
  return(stats)
}

# ==============================================================================
# 2. ANÁLISIS UNIVARIADO
# ==============================================================================
# Se añade el parámetro 'type' ("numeric" o "categorical")
analyze_univar <- function(data, variable, type){
  if (!variable %in% names(data)) stop("La variable no se encuentra en el dataframe.")
  
  var <- var_extract(data, variable)
  
  # Caso 1. Variable numérica indicada por el usuario
  if(type == "numeric"){
    var <- as.numeric(var) # Coerción basada en la indicación del usuario
    
    plot_1 <- ggplot(data = NULL, aes(x = var)) + 
      theme_minimal() + ylab(NULL) + xlab(NULL) + ggtitle(variable) +
      geom_density(fill = custom_palette[2], alpha = 0.8)
    
    boxplot <- ggplot(data = NULL, aes(x = 0, y = var)) + 
      geom_boxplot(fill = custom_palette[7]) +
      theme_minimal() + xlab(NULL) + ylab(NULL) + scale_x_continuous(breaks = 0)
    
    print(grid.arrange(plot_1, boxplot, ncol = 2))
    return(statcats(var))
  }
  
  # Caso 2. Variable categórica indicada por el usuario
  else if(type == "categorical"){
    var <- as.factor(var) # Coerción basada en la indicación del usuario
    n_levels <- length(levels(var))
    colors_used <- rep(custom_palette, length.out = n_levels)
    
    bar_plot <- ggplot(data = NULL, aes(x = var, fill = var)) + 
      geom_bar() +
      scale_fill_manual(values = colors_used) +
      theme_minimal() + ylab(NULL) + xlab(NULL) + ggtitle(variable) +
      theme(legend.position = "none")
    
    new_tab <- as.data.frame(table(var)) %>% setNames(c("Category", "Freq"))
    
    pie_plot <- ggplot(new_tab, aes(x = "", y = Freq, fill = Category)) +
      geom_bar(stat = "identity", width = 1, color = "white") +
      coord_polar("y", start = 0) + 
      scale_fill_manual(values = colors_used) +
      theme_void()
    
    print(grid.arrange(bar_plot, pie_plot, ncol = 2))
    return(table(var))
  } else {
    stop("Tipo no válido. Por favor, usa 'numeric' o 'categorical'.")
  }
}

# ==============================================================================
# 3. ANÁLISIS BIVARIADO
# ==============================================================================
# Se añaden 'type1' y 'type2'
analyze_bivar <- function(data, variable1, variable2, type1, type2){
  if (!all(c(variable1, variable2) %in% names(data))) {
    stop("Alguna variable no se encuentra en el dataframe.")
  }
  
  var1 <- var_extract(data, variable1)
  var2 <- var_extract(data, variable2)
  
  # Caso 1. Ambas numéricas
  if(type1 == "numeric" & type2 == "numeric"){
    var1 <- as.numeric(var1)
    var2 <- as.numeric(var2)
    df <- data.frame(v1 = var1, v2 = var2)
    
    methods <- c("pearson", "kendall", "spearman")
    cor_tests <- data.frame(
      Pearson = cor.test(var1, var2, method = "pearson")$estimate,
      Kendall = cor.test(var1, var2, method = "kendall")$estimate,
      Spearman = cor.test(var1, var2, method = "spearman")$estimate
    )
    
    points_plot <- ggplot(df, aes(x = v1, y = v2)) + 
      geom_point(size = 3, col = custom_palette[2], alpha = 0.7) + 
      xlab(variable1) + ylab(variable2) + theme_minimal()
    
    print(points_plot)
    return(t(cor_tests))
  }
  
  # Caso 2. Ambas categóricas
  else if(type1 == "categorical" & type2 == "categorical"){
    var1 <- as.factor(var1)
    var2 <- as.factor(var2)
    tabla <- table(var1, var2)
    df_tabla <- as.data.frame(tabla)
    
    heatmap_plot <- ggplot(df_tabla, aes(x = var1, y = var2, fill = Freq)) +
      geom_tile(color = "white") +
      scale_fill_gradientn(colors = c("#1e6262", "white", "#621e1e") ) + 
      theme_minimal() + 
      labs(x = variable1, y = variable2, fill = "Frecuencia") +
      ggtitle(paste("Mapa de Calor:", variable1, "vs", variable2))
    
    print(heatmap_plot)
    
    # Pruebas estadísticas
    ji_test <- chisq.test(tabla)
    f_test <- fisher.test(tabla, simulate.p.value = TRUE)
    coef_as <- assocstats(tabla)
    
    res_stats <- list(
      Chi_Cuadrada = ji_test,
      Fisher = f_test,
      Asociacion = data.frame(
        Coeficiente = c("Phi", "Contingencia", "V de Cramer"),
        Valores = c(coef_as$phi, coef_as$contingency, coef_as$cramer)
      )
    )
    return(res_stats)
  }
  
  # Caso 3. Mixto (Numérica y Categórica)
  else if((type1 == "numeric" & type2 == "categorical") | 
          (type1 == "categorical" & type2 == "numeric")) {
    
    is_num1 <- (type1 == "numeric")
    
    v_num <- if(is_num1) as.numeric(var1) else as.numeric(var2)
    v_cat <- if(is_num1) as.factor(var2) else as.factor(var1)
    name_num <- if(is_num1) variable1 else variable2
    name_cat <- if(is_num1) variable2 else variable1
    
    kw_test <- kruskal.test(x = v_num, g = v_cat)
    anova_t <- aov(v_num ~ v_cat)
    tukey_result <- TukeyHSD(anova_t)
    tukey_df <- as.data.frame(tukey_result$v_cat)
    
    dens_plot <- ggplot(data = NULL, aes(x = v_num, fill = v_cat)) + 
      geom_density(alpha = 0.6) + 
      scale_fill_manual(values = rep(custom_palette, length.out = length(levels(v_cat)))) +
      theme_minimal() + xlab(name_num) + labs(fill = name_cat) + ylab("Densidad")
    
    int_plot <- ggplot(tukey_df, aes(y = reorder(rownames(tukey_df), diff), x = diff, xmin = lwr, xmax = upr)) + 
      geom_errorbar(linewidth = 1, width = 0.2, col = custom_palette[2]) + 
      geom_pointrange(size = 1, col = custom_palette[6]) + 
      labs(x = "Diferencia media (Tukey)", y = "Pares de grupos") + theme_minimal()
    
    print(grid.arrange(dens_plot, int_plot, ncol = 2))
    
    return(list(Kruskal_Wallis = kw_test, ANOVA = summary(anova_t)))
  } else {
    stop("Combinación de tipos no válida. Usa 'numeric' y/o 'categorical'.")
  }
}

# ==============================================================================
# 4. ANÁLISIS MULTIVARIADO
# ==============================================================================
analyze_multivar <- function(data, vars_num, type = "numeric"){
  if(type != "numeric"){
    stop("El análisis multivariado actual solo soporta variables numéricas.")
  }
  
  if (!all(vars_num %in% names(data))) {
    stop("Algunas variables especificadas no se encuentran en el dataframe.")
  }
  
  df_num <- data %>% dplyr::select(all_of(vars_num))
  
  # Forzamos a numérico ya que el usuario indicó que lo son
  df_num <- df_num %>% mutate_all(as.numeric)
  
  cor_matrix <- cor(df_num, use = "complete.obs")
  
  cor_plot <- ggcorrplot(cor_matrix, hc.order = TRUE, type = "full", 
                         lab = TRUE, method = "square",
                         colors = c(custom_palette[1], "white", custom_palette[5])) +
    ggtitle("Matriz de Correlación") + labs(x = "", y = "") +
    theme_minimal()
  
  print(cor_plot)
  
  lista_estadisticas <- lapply(df_num, statcats)
  
  return(list(Matriz_Correlacion = cor_matrix, Estadisticas = lista_estadisticas))
}

# ==============================================================================
# EJEMPLOS DE USO
# ==============================================================================
data(iris)
data("diamonds")
data("mtcars")
data(crabs)

# Univariado
analyze_univar(iris, "Sepal.Length", type = "numeric")
analyze_univar(iris, "Species", type = "categorical")

# Bivariado
analyze_bivar(diamonds, "x", "y", type1 = "numeric", type2 = "numeric")
analyze_bivar(diamonds, "cut", "color", type1 = "categorical", type2 = "categorical")

analyze_bivar(crabs, "sp", "sex", type1 = "categorical", type2 = "categorical")
analyze_bivar(iris, "Petal.Length", "Species", type1 = "numeric", type2 = "categorical")

# Multivariado
analyze_multivar(iris, c("Sepal.Width", "Sepal.Length", "Petal.Width", "Petal.Length"), type = "numeric")

