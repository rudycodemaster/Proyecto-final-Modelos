# Cargar librerías necesarias
library(readxl)
library(MASS)  
library(VGAM)  
library(dplyr) # Para manipulación de datos
library(fitdistrplus) 
library(ggplot2)     # Para gráficos
library(goftest)     # Para pruebas KS y Anderson-Darling
library(actuar)
library(tibble)


datos<-read_excel("C:/Users/omarp/Downloads/DBMS.xlsx")
# Filtrar para eliminar niveles con solo un valor 
sapply(lapply(datos, unique), length)
datos$Marca              <- factor(datos$Marca)
datos$Cobertura          <- factor(datos$Cobertura)
datos$Tipo_de_Perdida    <- factor(datos$Tipo_de_Perdida)
datos$Causa_del_siniestro<- factor(datos$Causa_del_siniestro)
datos$Causa_del_siniestro
######### Datos totales Frecuencia
#### Recomendación
marcas_unicas <- unique(datos$Marca)
resumen_marcas <- datos %>%
  group_by(Marca) %>%
  summarise(
    media_frec = mean(Numero_de_Vehiculos, na.rm = TRUE),
    varianza_frec = var(Numero_de_Vehiculos, na.rm = TRUE),
    .groups = "drop"
  )

# Función para recomendar distribución según media y varianza
resumen_marcas <- resumen_marcas %>%
  mutate(
    recomendacion = case_when(
      abs(media_frec - varianza_frec) < 0.1 * media_frec ~ "Poisson",
      varianza_frec < media_frec ~ "Binomial",
      varianza_frec > media_frec ~ "Binomial Negativa o Geométrica"
    )
  )

# Estimación de parámetros usando el método de momentos
resumen_marcas <- resumen_marcas %>%
  mutate(
    # ──────────────────────────────
    # POISSON: λ = media
    # ──────────────────────────────
    lambda_poisson = ifelse(recomendacion == "Poisson", media_frec, NA),
    
    # ──────────────────────────────
    # BINOMIAL:
    # p = 1 - (s² / media),    n = media / p
    # Solo válido si varianza < media
    # ──────────────────────────────
    p_binomial = ifelse(recomendacion == "Binomial" & varianza_frec < media_frec,
                        1 - (varianza_frec / media_frec), NA),
    n_binomial = ifelse(recomendacion == "Binomial" & !is.na(p_binomial) & p_binomial > 0,
                        media_frec / p_binomial, NA),
    
    # ──────────────────────────────
    # BINOMIAL NEGATIVA:
    # r = media² / (varianza - media),
    # p = media / varianza
    # Solo válido si varianza > media
    # ──────────────────────────────
    r_binomial_negativa = ifelse(recomendacion == "Binomial Negativa o Geométrica" & varianza_frec > media_frec,
                                 (media_frec^2) / (varianza_frec - media_frec), NA),
    # Asegurar que r no sea menor a 1
    r_binomial_negativa = ifelse(!is.na(r_binomial_negativa) & r_binomial_negativa < 1,
                                 1, r_binomial_negativa),
    p_binomial_negativa = ifelse(recomendacion == "Binomial Negativa o Geométrica" & varianza_frec > media_frec,
                                 media_frec / varianza_frec, NA),
    
    # ──────────────────────────────
    # GEOMÉTRICA (caso especial binomial negativa con r = 1)
    # p = 1 / media
    # ──────────────────────────────
    p_geom = ifelse(recomendacion == "Binomial Negativa o Geométrica" & !is.na(media_frec),
                    1 / media_frec, NA)
  )

# Mostrar la tabla con los resultados
print(resumen_marcas)

################################ severidad ###############################

montos<-datos$Monto_de_Siniestros
sum(montos)
head(montos)
nobs     <- length(montos)
sumlog   <- sum(log(montos))
media    <- mean(montos)
# Lista para almacenar resultados por marca
resultados_por_marca <- list()

for (marca in marcas_unicas) {
  
  cat("\n\n==============================\n")
  cat("Resultados para la Marca:", marca, "\n")
  cat("==============================\n")
  
  # Filtrar montos de la marca específica
  montos <- datos %>%
    filter(Marca == marca) %>%
    pull(Monto_de_Siniestros) %>%
    as.numeric()
  
  # Saltar si no hay suficientes datos
  if (length(montos) < 5) {
    cat("Insuficientes datos para ajustar modelos.\n")
    next
  }
  
  nobs   <- length(montos)
  sumlog <- sum(log(montos))
  media  <- mean(montos)
  
  ## AJUSTE GAMMA
  fit_gamma_mme <- fitdist(montos, "gamma", method = "mme")
  alpha1 <- fit_gamma_mme$estimate["shape"]
  
  fg <- function(a) {
    - ( nobs * ( a * (log(a) - log(media) - 1) ) -
          nobs * lgamma(a) +
          (a - 1) * sumlog )
  }
  alpha_hat <- nlm(fg, alpha1)$estimate
  theta_hat <- media / alpha_hat
  
  fit_gamma_mle <- fitdist(montos, "gamma",
                           method = "mle",
                           start = list(shape = alpha_hat, scale = theta_hat))
  
  ## AJUSTE LOG-NORMAL
  fit_logn <- fitdist(montos, "lnorm", method = "mle")
  
  ## AJUSTE PARETO II
  fit_pareto <- fitdist(montos, "pareto", start = list(shape = 1, scale = median(montos)))
  
  ## AJUSTE WEIBULL
  fit_weibull <- fitdist(montos, "weibull", method = "mle")
  
  ## AJUSTE LOG-GAMMA
  fit_lgamma <- fitdist(montos, "lgamma", method = "mle")
  
  ## AJUSTE LOG-LOGÍSTICA
  fit_llog <- fitdist(montos, "llogis", method = "mle")
  
  ## AJUSTE INVERSA WEIBULL
  fit_invw <- fitdist(montos, "invweibull", method = "mle")
  
  ## TABLA COMPARATIVA
  modelos <- list(
    gamma    = fit_gamma_mle,
    logn     = fit_logn,
    pareto   = fit_pareto,
    weibull  = fit_weibull,
    lgamma   = fit_lgamma,
    llog     = fit_llog,
    invw     = fit_invw
  )
  
  loglik <- sapply(modelos, function(x) x$loglik)
  aic    <- sapply(modelos, function(x) x$aic)
  bic    <- sapply(modelos, function(x) x$bic)
  
  salida <- rbind(LogLik = loglik,
                  AIC    = aic,
                  BIC    = bic)
  
  print(t(salida))
  
  # Guardar en la lista
  resultados_por_marca[[marca]] <- t(salida)
}
#######como se observa todos las marcas tienen como frecuencia una distribucion geometrica y de severidad una lgamma
####### esto facilita el calculo de las metricas de la variable agregada
library(tibble)
prob_geometrica <- tibble(
  Marca = c("GENERAL MOTORS", "CHRYSLER", "NISSAN", "FORD", "VOLKSWAGEN", "HONDA", "TOYOTA"),
  p_geom = c(0.02915220, 0.06180088, 0.01994933, 0.05748766, 0.03079305, 0.0272682, 0.02893794) # ejemplo
)
marcas_unicas <- unique(datos$Marca)
parametros_gamma <- list()

for (marca in marcas_unicas) {
  
  montos <- datos %>%
    filter(Marca == marca) %>%
    pull(Monto_de_Siniestros) %>%
    as.numeric()
  
  montos <- montos[!is.na(montos) & montos > 0 & is.finite(montos)]
  
  if (length(montos) < 5) {
    next
  }
  
  media_montos <- mean(montos)
  var_montos <- var(montos)
  
  if (media_montos <= 0 || var_montos <= 0) {
    next
  }
  
  # Método de momentos Gamma
  alpha_gamma <- media_montos^2 / var_montos
  theta_gamma <- var_montos / media_montos
  EX <- alpha_gamma * theta_gamma
  VarX <- alpha_gamma * theta_gamma^2
  
  # p de tabla geométrica
  p_geom <- prob_geometrica %>%
    filter(Marca == marca) %>%
    pull(p_geom)
  
  if (length(p_geom) == 0 || is.na(p_geom) || p_geom <= 0 || p_geom >= 1) {
    E_N <- NA
    Var_N <- NA
    E_S <- NA
    Var_S <- NA
    CV_S <- NA
  } else {
    E_N <- (1 - p_geom) / p_geom
    Var_N <- (1 - p_geom) / p_geom^2
    
    E_S <- E_N * EX
    Var_S <- E_N * VarX + Var_N * (EX)^2
    SD_S <- sqrt(Var_S)
    CV_S <- SD_S / E_S
  }
  
  parametros_gamma[[marca]] <- tibble(
    Marca = marca,
    Alpha = alpha_gamma,
    Theta = theta_gamma,
    `E[X]` = EX,
    `Var[X]` = VarX,
    `p_geom` = p_geom,
    `E[N]` = E_N,
    `Var[N]` = Var_N,
    `E[S]` = E_S,
    `Var[S]` = Var_S,
    `SD[S]` = SD_S,
    `CV[S]` = CV_S
  )
}

tabla_parametros_completa <- bind_rows(parametros_gamma)

# Mostrar tabla final
print(tabla_parametros_completa)

# ==============================
# GRAFICOS COMPARATIVOS POR MARCA
# ==============================

# Prima pura esperada E[S]
ggplot(tabla_parametros_completa, aes(x = reorder(Marca, `E[S]`), y = `E[S]`)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  labs(title = "Prima pura esperada (E[S]) por marca",
       x = "Marca",
       y = "E[S] (MXN)") +
  theme_minimal()

# Desviación estándar agregada SD[S]
ggplot(tabla_parametros_completa, aes(x = reorder(Marca, `SD[S]`), y = `SD[S]`)) +
  geom_col(fill = "darkgreen") +
  coord_flip() +
  labs(title = "Desviación estándar agregada (SD[S]) por marca",
       x = "Marca",
       y = "SD[S] (MXN)") +
  theme_minimal()

# Coeficiente de variación CV[S]
ggplot(tabla_parametros_completa, aes(x = reorder(Marca, `CV[S]`), y = `CV[S]`)) +
  geom_col(fill = "tomato") +
  coord_flip() +
  labs(title = "Coeficiente de variación (CV[S]) por marca",
       x = "Marca",
       y = "CV[S]") +
  theme_minimal()

# Parámetros de simulación
n_sim <- 100000  # número de simulaciones por marca

# Inicializar lista de resultados
resultados_mc <- list()

for (i in 1:nrow(tabla_parametros_completa)) {
  
  marca <- tabla_parametros_completa$Marca[i]
  p_geom <- tabla_parametros_completa$p_geom[i]
  alpha <- tabla_parametros_completa$Alpha[i]
  theta <- tabla_parametros_completa$Theta[i]
  
  # Verificar que los parámetros sean válidos
  if (is.na(p_geom) || p_geom <= 0 || p_geom >= 1 ||
      is.na(alpha) || alpha <= 0 ||
      is.na(theta) || theta <= 0) {
    next
  }
  
  # Simulación de frecuencia N usando Geométrica
  N_sim <- rgeom(n_sim, p_geom) + 1
  
  # Simulación de pérdidas agregadas
  agg_loss <- numeric(n_sim)
  
  for (j in 1:n_sim) {
    if (N_sim[j] > 0) {
      agg_loss[j] <- sum(rgamma(N_sim[j], shape = alpha, scale = theta))
    }
  }
  
  # Calcular media y percentil 99.5%
  mean_S <- mean(agg_loss)
  p995_S <- quantile(agg_loss, 0.995)
  
  # Guardar en lista
  resultados_mc[[marca]] <- tibble(
    Marca = marca,
    Sim_Mean_S = mean_S,
    Sim_Quantile_99.5 = p995_S
  )
}

# Combinar en tabla final
tabla_simulacion_mc <- bind_rows(resultados_mc)

# Mostrar tabla de simulación
print(tabla_simulacion_mc)
# Definir el factor de desviaciones estándar a sumar

k <- 1.5  # puedes cambiar según criterio actuarial de solvencia

# Calcular prima con margen
tabla_parametros_completa <- tabla_parametros_completa %>%
  mutate(
    Prima_Pura = `E[S]`,
    Margen = k * `SD[S]`,
    Prima_Comercial = Prima_Pura + Margen
  )

# Mostrar tabla con las columnas relevantes
print(
  tabla_parametros_completa %>%
    select(Marca, Prima_Pura, Margen, Prima_Comercial)
)
