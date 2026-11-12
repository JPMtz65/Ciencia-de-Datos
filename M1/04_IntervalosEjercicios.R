#Intervalos de confianza 

#Parametros fijos 
set.seed(1989) 
alpha <- 0.05 #(1-alpha)*100% nivel de confianza
z_val <- qnorm(1 - alpha/2) #Cuantil normal estandar


#Bibliotecas ===================================================================
pacman::p_load(tidyverse)
data("iris")


# ==============================================================================
# 1. CANDIDATO PRESIDENCIAL
# ==============================================================================
n1 <- 1000
k <- 520
p_real <- k / n1

muestra_votos <- rbinom(n1, 1, p_real)

#Estimador puntual (Media muestral)
p_hat <- mean(muestra_votos)

#Suma de los 1s
exitos <- sum(muestra_votos)

# IC TLC (Asintotico)
error_est1 <- sqrt(p_hat * (1 - p_hat) / n1)

ic_tlc_1 <- c(p_hat - z_val * error_est1, p_hat + z_val * error_est1)

# IC Clopper-Pearson
ic_cp_1 <- c(qbeta(alpha/2, exitos, n1 - exitos + 1), 
             qbeta(1 - alpha/2, exitos + 1, n1 - exitos))

df_1 <- data.frame(
  Metodo = c("TLC (Asintótico)", "Clopper-Pearson (Exacto)"),
  Estimador = p_hat,
  Linf = c(ic_tlc_1[1], ic_cp_1[1]),
  Lsup = c(ic_tlc_1[2], ic_cp_1[2])
)

p1 <- ggplot(df_1, aes(y = Metodo, x = Estimador)) +
  geom_errorbar(aes(xmin = Linf, xmax = Lsup), width = 0.2, linewidth = 1,color = "#833C63") +
  geom_vline(xintercept = p_real, linetype = "dashed", color = "#B69BB6", linewidth = 1) +
  labs(title = "Proporción de Votos", x = "Proporción Estimada", y = "", subtitle = paste("p estimada = ", p_real)) +
  theme_minimal() + scale_x_continuous(n.breaks = 10)

print(p1)


# ==============================================================================
# 2. EMBOTELLADORA (Normal - Media y Desviación Estándar)
# ==============================================================================
n2 <- 30
mu_real <- 500
sigma_real <- 2.1
muestra_botellas <- rnorm(n2, mu_real, sigma_real)

x_bar <- mean(muestra_botellas)
s_dev <- sd(muestra_botellas)
t_val <- qt(1 - alpha/2, df = n2 - 1)

# MEDIA
# IC Pivotal (t-Student)
ic_media_pivotal <- c(x_bar - t_val * (s_dev/sqrt(n2)), x_bar + t_val * (s_dev/sqrt(n2)))
# IC Asintótico (Z)
ic_media_asintotico <- c(x_bar - z_val * (s_dev/sqrt(n2)), x_bar + z_val * (s_dev/sqrt(n2)))

df_2_media <- data.frame(
  Metodo = c("Pivotal (t-Student)", "Asintótico (Normal)"),
  Estimador = x_bar,
  Linf = c(ic_media_pivotal[1], ic_media_asintotico[1]),
  Lsup = c(ic_media_pivotal[2], ic_media_asintotico[2])
)

p2_media <- ggplot(df_2_media, aes(y = Metodo, x = Estimador)) +
  geom_point(size = 3, color = "#833C63") +
  geom_errorbar(aes(xmin = Linf, xmax = Lsup), width = 0.2, linewidth = 1, color = "#833C63") +
  geom_vline(xintercept = mu_real, linetype = "dashed", color = "#B69BB6") +
  labs(title = "Media Poblacional de Llenado", 
       x = "Media (ml)", y = "",
       ) +
  theme_minimal() + scale_x_continuous(n.breaks = 10)

p2_media

# DESVIACIÓN ESTÁNDAR
# IC Pivotal (Ji-cuadrada)
ji_inf <- qchisq(1 - alpha/2, df = n2 - 1)
ji_sup <- qchisq(alpha/2, df = n2 - 1)
ic_sd_pivotal <- c(sqrt((n2 - 1) * s_dev^2 / ji_inf), sqrt((n2 - 1) * s_dev^2 / ji_sup))

# IC Asintótico (Delta method)
error_est_sd <- s_dev / sqrt(2 * n2)
ic_sd_asintotico <- c(s_dev - z_val * error_est_sd, s_dev + z_val * error_est_sd)

df_2_sd <- data.frame(
  Metodo = c("Pivotal (Ji-cuadrada)", "Asintótico (Método Delta)"),
  Estimador = s_dev,
  Linf = c(ic_sd_pivotal[1], ic_sd_asintotico[1]),
  Lsup = c(ic_sd_pivotal[2], ic_sd_asintotico[2])
)

p2_sd <- ggplot(df_2_sd, aes(y = Metodo, x = Estimador)) +
  geom_point(size = 3, color = "#833C63") +
  geom_errorbar(aes(xmin = Linf, xmax = Lsup), width = 0.2, linewidth = 1, color = "#833C63") +
  geom_vline(xintercept = sigma_real, linetype = "dashed", color = "#B69BB6") +
  labs(title = "Desviación Estándar de Llenado", x = "Desviación Estándar (ml)", y = "") +
  theme_minimal() + scale_x_continuous(n.breaks = 10)

p2_media
p2_sd

# ==============================================================================
# 3. URGENCIAS HOSPITALARIAS (Poisson)
# ==============================================================================
n3 <- 30
lambda_real <- 15
muestra_poisson <- rpois(n3, lambda_real)
lambda_hat <- mean(muestra_poisson)
suma_x <- sum(muestra_poisson)

# IC Asintótico
error_est3 <- sqrt(lambda_hat / n3)
ic_pois_asintotico <- c(lambda_hat - z_val * error_est3, lambda_hat + z_val * error_est3)

# IC Aproximación Ji-cuadrada (Exacto de Garwood)
ic_pois_ji <- c(qchisq(alpha/2, 2 * suma_x) / (2 * n3),
                 qchisq(1 - alpha/2, 2 * (suma_x + 1)) / (2 * n3))

df_3 <- data.frame(
  Metodo = c("Asintótico (Wald)", "Exacto (Ji-cuadrada)"),
  Estimador = lambda_hat,
  Linf = c(ic_pois_asintotico[1], ic_pois_ji[1]),
  Lsup = c(ic_pois_asintotico[2], ic_pois_ji[2])
)

p3 <- ggplot(df_3, aes(y = Metodo, x = Estimador)) +
  geom_point(size = 3, color = "#833C63") +
  geom_errorbar(aes(xmin = Linf, xmax = Lsup), width = 0.2, linewidth = 1, color = "#833C63") +
  geom_vline(xintercept = lambda_real, linetype = "dashed", color = "#B69BB6") +
  labs(title = "Tasa Media de Ingresos a Urgencias", 
       x = "Pacientes / hora", y = "",
       subtitle = paste("Lambda estimada: ", lambda_hat)) +
  theme_minimal() + scale_x_continuous(n.breaks = 10)

p3

# ==============================================================================
# 4. LONGITUD PÉTALO IRIS (Distribución Desconocida - Bootstrap)
# ==============================================================================
data_iris <- iris$Petal.Length
n4 <- length(data_iris)
media_iris <- mean(data_iris)

B <- 5000 # Número de remuestreos
medias_boot <- numeric(B)
for (i in 1:B) {
  muestra_boot <- sample(data_iris, size = n4, replace = TRUE)
  medias_boot[i] <- mean(muestra_boot)
}

ic_bootstrap <- quantile(medias_boot, probs = c(alpha/2, 1 - alpha/2))

df_4 <- data.frame(
  Metodo = "Bootstrap (Percentil)",
  Estimador = media_iris,
  Linf = ic_bootstrap[1],
  Lsup = ic_bootstrap[2]
)

p4 <- ggplot(df_4, aes(x = Metodo, y = Estimador)) +
  geom_point(size = 3, color = "#833C63") +
  geom_errorbar(aes(ymin = Linf, ymax = Lsup), width = 0.2, linewidth = 1, color = "#833C63") +
  labs(title = "Longitud Media del Pétalo (Iris)", y = "Longitud (cm)", x = "") +
  theme_minimal() + scale_y_continuous(n.breaks=10) + coord_flip()

p4
