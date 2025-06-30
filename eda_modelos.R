install.packages("gt")
library(gt)
library(dplyr)
library(moments)
library(ggplot2)
install.packages("tidyr")
library(tidyr)

# Asegurar que las variables numéricas estén en el formato correcto
base_modelos_para_R <- DBMS %>%
  mutate(
    `Numero de Vehiculos Siniestrados` = as.numeric(`Numero de Vehiculos`),
    `Monto de Siniestros` = as.numeric(`Monto de Siniestros`),
    `Monto Pagado` = as.numeric(`Monto Pagado`)
  )

# Calcular las medidas
desv_est <- sapply(base_modelos_para_R[, c("Numero de Vehiculos Siniestrados", "Monto de Siniestros", "Monto Pagado")], sd)
varianza <- sapply(base_modelos_para_R[, c("Numero de Vehiculos Siniestrados", "Monto de Siniestros", "Monto Pagado")], var)
rango <- sapply(base_modelos_para_R[, c("Numero de Vehiculos Siniestrados", "Monto de Siniestros", "Monto Pagado")], function(x) max(x) - min(x))
coef_var <- sapply(base_modelos_para_R[, c("Numero de Vehiculos Siniestrados", "Monto de Siniestros", "Monto Pagado")], function(x) sd(x) / mean(x))
sesgo <- sapply(base_modelos_para_R[, c("Numero de Vehiculos Siniestrados", "Monto de Siniestros", "Monto Pagado")], skewness)
curtosis <- sapply(base_modelos_para_R[, c("Numero de Vehiculos Siniestrados", "Monto de Siniestros", "Monto Pagado")], kurtosis)
suma <- sapply(base_modelos_para_R[, c("Numero de Vehiculos Siniestrados", "Monto de Siniestros", "Monto Pagado")], sum)

# hacemos el resumen para la tabla
resumen_completo <- tibble(
  Estadística = c("Mínimo", "1er Cuartil", "Mediana", "Media", "3er Cuartil", "Máximo", 
                  "Desv. Est.", "Varianza", "Rango", "Coef. Var.", "Sesgo", "Curtosis", "Suma"),
  `Vehículos Siniestrados` = c(1.00, 1.00, 3.00, 30.14, 11.00, 7837.00, 
                             desv_est[1], varianza[1], rango[1], coef_var[1], sesgo[1], curtosis[1], suma[1]),
  `Monto de Siniestros` = c(0, 8991, 36900, 356580, 166820, 51558310, 
                          desv_est[2], varianza[2], rango[2], coef_var[2], sesgo[2], curtosis[2], suma[2]),
  `Monto Pagado` = c(0, 8737, 34428, 309196, 154296, 40957671, 
                   desv_est[3], varianza[3], rango[3], coef_var[3], sesgo[3], curtosis[3], suma[3])
)


tabla_gt <- resumen_completo %>%
  gt() %>%
  tab_header(
    title = "Resumen Completo de las Variables",
    subtitle = "Estadísticas Descriptivas: Vehículos Siniestrados, Monto de Siniestros y Monto Pagado"
  ) %>%
  fmt_number(
    columns = -Estadística, # Formatear solo las columnas numéricas
    decimals = 2
  ) %>%
  cols_label(
    `Estadística` = "Estadística",
    `Vehículos Siniestrados` = "Vehículos Siniestrados",
    `Monto de Siniestros` = "Monto de Siniestros",
    `Monto Pagado` = "Monto Pagado"
  ) %>%
  tab_options(
    table.font.size = "medium",
    column_labels.font.weight = "bold",
    column_labels.background.color = "lightgray",
    table.border.top.width = px(2),
    table.border.top.color = "black",
    table.border.bottom.width = px(2),
    table.border.bottom.color = "black"
  )

# Mostrar la tabla estilizada
print(tabla_gt)

#-------------------------------------------------------------------------
#analisis mediante hisograma y boxplot

# Transformar las variables aplicando logaritmo natural
base_modelos_para_R <- base_modelos_para_R %>%
  mutate(
    `Log Numero de Vehiculos Siniestrados` = log(`Numero de Vehiculos Siniestrados` + 1),
    `Log Monto de Siniestros` = log(`Monto de Siniestros` + 1),
    `Log Monto Pagado` = log(`Monto Pagado` + 1)
  )

# Histogramas con datos suavizados
ggplot(base_modelos_para_R, aes(x = `Log Numero de Vehiculos Siniestrados`)) +
  geom_histogram(binwidth = 0.1, fill = "darkblue", color = "black", alpha = 0.7) +
  labs(title = "Distribución Logarítmica: Número de Vehículos Siniestrados",
       x = "Log Número de Vehículos Siniestrados", y = "Frecuencia") +
  theme_minimal()

ggplot(base_modelos_para_R, aes(x = `Log Monto de Siniestros`)) +
  geom_histogram(binwidth = 0.1, fill = "cadetblue4", color = "black", alpha = 0.7) +
  labs(title = "Distribución Logarítmica: Monto de Siniestros",
       x = "Log Monto de Siniestros", y = "Frecuencia") +
  theme_minimal()

ggplot(base_modelos_para_R, aes(x = `Log Monto Pagado`)) +
  geom_histogram(binwidth = 0.1, fill = "azure3", color = "black", alpha = 0.7) +
  labs(title = "Distribución Logarítmica: Monto Pagado",
       x = "Log Monto Pagado", y = "Frecuencia") +
  theme_minimal()

# Boxplot - Log Número de Vehículos Siniestrados
boxplot_vehiculos <- ggplot(base_modelos_para_R, aes(y = `Log Numero de Vehiculos Siniestrados`)) +
  geom_boxplot(outlier.color = "deepskyblue4", outlier.size = 2, fill = "dimgray", alpha = 0.7) +
  labs(
    title = "Boxplot: Log Número de Vehículos Siniestrados",
    y = "Log Número de Vehículos",
    x = NULL
  ) +
  theme_minimal(base_size = 15)

# Mostrar y guardar el boxplot del Número de Vehículos Siniestrados
print(boxplot_vehiculos)
ggsave("Boxplot_Log_Numero_de_Vehiculos_Siniestrados.png", plot = boxplot_vehiculos, width = 8, height = 6)

# Boxplot - Log Monto de Siniestros
boxplot_siniestros <- ggplot(base_modelos_para_R, aes(y = `Log Monto de Siniestros`)) +
  geom_boxplot(outlier.color = "deepskyblue4", outlier.size = 2, fill = "darkseagreen4", alpha = 0.7) +
  labs(
    title = "Boxplot: Log Monto de Siniestros",
    y = "Log Monto de Siniestros",
    x = NULL
  ) +
  theme_minimal(base_size = 15)

# Mostrar y guardar el boxplot del Monto de Siniestros
print(boxplot_siniestros)
ggsave("Boxplot_Log_Monto_de_Siniestros.png", plot = boxplot_siniestros, width = 8, height = 6)

# Boxplot - Log Monto Pagado
boxplot_pagado <- ggplot(base_modelos_para_R, aes(y = `Log Monto Pagado`)) +
  geom_boxplot(outlier.color = "deepskyblue4", outlier.size = 2, fill = "cornsilk4", alpha = 0.7) +
  labs(
    title = "Boxplot: Log Monto Pagado",
    y = "Log Monto Pagado",
    x = NULL
  ) +
  theme_minimal(base_size = 15)

# Mostrar y guardar el boxplot del Monto Pagado
print(boxplot_pagado)
ggsave("Boxplot_Log_Monto_Pagado.png", plot = boxplot_pagado, width = 8, height = 6)

# Boxplot combinado de resumen (todas las variables)
boxplot_resumen <- base_modelos_para_R %>%
  pivot_longer(
    cols = c(`Log Numero de Vehiculos Siniestrados`, `Log Monto de Siniestros`, `Log Monto Pagado`),
    names_to = "Variable",
    values_to = "Valor"
  ) %>%
  ggplot(aes(x = Variable, y = Valor, fill = Variable)) +
  geom_boxplot(outlier.color = "deepskyblue4", outlier.size = 2) +
  scale_fill_manual(values = c(
    "Log Numero de Vehiculos Siniestrados" = "dimgray",
    "Log Monto de Siniestros" = "darkseagreen4",
    "Log Monto Pagado" = "cornsilk4"
  )) +
  labs(
    title = "Boxplot Resumen: Variables Logarítmicas",
    x = "Variable",
    y = "Valor Transformado (Log)"
  ) +
  theme_minimal(base_size = 15) +
  theme(legend.position = "none")

# Mostrar el gráfico combinado
print(boxplot_resumen)

# Guardar todos los gráficos
ggsave("Boxplot_Log_Numero_de_Vehiculos_Siniestrados.png", plot = boxplot_vehiculos, width = 8, height = 6)
ggsave("Boxplot_Log_Monto_de_Siniestros.png", plot = boxplot_siniestros, width = 8, height = 6)
ggsave("Boxplot_Log_Monto_Pagado.png", plot = boxplot_pagado, width = 8, height = 6)
ggsave("Boxplot_Log_Resumen_Variables.png", plot = boxplot_resumen, width = 10, height = 8)


#---------------------------------------------------------------------
# Cargar librerías necesarias
# Cargar librerías necesarias
library(dplyr)
library(moments)  # Para calcular curtosis y sesgo
library(gt)       # Para generar tablas estilizadas

# Generar el resumen estadístico en formato vertical
resumen_estadistico <- tibble(
  `Medida Estadística` = c(
    "Media", "Mediana", "Desv. Estándar", "Varianza",
    "Sesgo", "Curtosis", "Coef. Variación"
  ),
  `Vehículos Siniestrados` = c(
    mean(base_modelos_para_R$`Log Numero de Vehiculos Siniestrados`, na.rm = TRUE),
    median(base_modelos_para_R$`Log Numero de Vehiculos Siniestrados`, na.rm = TRUE),
    sd(base_modelos_para_R$`Log Numero de Vehiculos Siniestrados`, na.rm = TRUE),
    var(base_modelos_para_R$`Log Numero de Vehiculos Siniestrados`, na.rm = TRUE),
    skewness(base_modelos_para_R$`Log Numero de Vehiculos Siniestrados`, na.rm = TRUE),
    kurtosis(base_modelos_para_R$`Log Numero de Vehiculos Siniestrados`, na.rm = TRUE),
    sd(base_modelos_para_R$`Log Numero de Vehiculos Siniestrados`, na.rm = TRUE) /
      mean(base_modelos_para_R$`Log Numero de Vehiculos Siniestrados`, na.rm = TRUE)
  ),
  `Monto de Siniestros` = c(
    mean(base_modelos_para_R$`Log Monto de Siniestros`, na.rm = TRUE),
    median(base_modelos_para_R$`Log Monto de Siniestros`, na.rm = TRUE),
    sd(base_modelos_para_R$`Log Monto de Siniestros`, na.rm = TRUE),
    var(base_modelos_para_R$`Log Monto de Siniestros`, na.rm = TRUE),
    skewness(base_modelos_para_R$`Log Monto de Siniestros`, na.rm = TRUE),
    kurtosis(base_modelos_para_R$`Log Monto de Siniestros`, na.rm = TRUE),
    sd(base_modelos_para_R$`Log Monto de Siniestros`, na.rm = TRUE) /
      mean(base_modelos_para_R$`Log Monto de Siniestros`, na.rm = TRUE)
  ),
  `Monto Pagado` = c(
    mean(base_modelos_para_R$`Log Monto Pagado`, na.rm = TRUE),
    median(base_modelos_para_R$`Log Monto Pagado`, na.rm = TRUE),
    sd(base_modelos_para_R$`Log Monto Pagado`, na.rm = TRUE),
    var(base_modelos_para_R$`Log Monto Pagado`, na.rm = TRUE),
    skewness(base_modelos_para_R$`Log Monto Pagado`, na.rm = TRUE),
    kurtosis(base_modelos_para_R$`Log Monto Pagado`, na.rm = TRUE),
    sd(base_modelos_para_R$`Log Monto Pagado`, na.rm = TRUE) /
      mean(base_modelos_para_R$`Log Monto Pagado`, na.rm = TRUE)
  )
)

# Crear la tabla estilizada sin argumentos inválidos
tabla_estadistica <- resumen_estadistico %>%
  gt() %>%
  tab_header(
    title = "Resumen Estadístico",
    subtitle = "Análisis de Variables Transformadas con Logaritmos"
  ) %>%
  fmt_number(
    columns = where(is.numeric),
    decimals = 2
  ) %>%
  cols_label(
    `Medida Estadística` = "Medida Estadística",
    `Vehículos Siniestrados` = "Vehículos Siniestrados",
    `Monto de Siniestros` = "Monto de Siniestros",
    `Monto Pagado` = "Monto Pagado"
  ) %>%
  tab_options(
    column_labels.font.weight = "bold",
    table.border.top.width = px(2),
    table.border.bottom.width = px(2),
    data_row.padding = px(5),
    table.font.size = "medium"
  )

# Mostrar la tabla estilizada
print(tabla_estadistica)

#----------------------------------------------------------------------

# Frecuencia de siniestros por Tipo de Pérdida y Causa del Siniestro
frecuencia_tipo_causa <- base_modelos_para_R %>%
  group_by(`Tipo de Perdida`, `Causa del siniestro`) %>%
  summarise(
    Frecuencia = n(),  # Conteo de siniestros
    Costo_Total = sum(`Monto de Siniestros`, na.rm = TRUE)  # Suma de montos de siniestros
  ) %>%
  arrange(desc(Frecuencia))

# Mostrar los 10 eventos más comunes
print(head(frecuencia_tipo_causa, 10))

# Visualización: Frecuencia de siniestros por causa del siniestro
ggplot(frecuencia_tipo_causa, aes(x = reorder(`Causa del siniestro`, Frecuencia), y = Frecuencia, fill = `Tipo de Perdida`)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c(
    "Parcial" = "darkseagreen4",
    "SIPAC perdida parcial" = "palegreen4",
    "Rotura de cristales" = "olivedrab",
    "RC SIPAC" = "steelblue4",
    "Total" = "gray30",
    "SIPAC perdida total" = "gray50"
  )) +
  labs(
    title = "Frecuencia de Siniestros por Tipo de Pérdida y Causa",
    x = "Causa del Siniestro",
    y = "Frecuencia",
    fill = "Tipo de Pérdida"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.y = element_text(size = 10),
    legend.position = "top"
  ) +
  coord_flip()


#-------------------------------------------------------------------

# Comparación entre marcas: Frecuencia y costo promedio
comparacion_marcas <- base_modelos_para_R %>%
  group_by(Marca) %>%
  summarise(
    Frecuencia = n(),
    Costo_Total = sum(`Monto de Siniestros`, na.rm = TRUE),
    Costo_Promedio = mean(`Monto de Siniestros`, na.rm = TRUE)
  ) %>%
  arrange(desc(Frecuencia))

# Mostrar el resumen de las 10 marcas más afectadas
print(head(comparacion_marcas, 10))

# Visualización: Frecuencia y montos de siniestros por marca
ggplot(comparacion_marcas, aes(x = reorder(Marca, -Frecuencia), y = Frecuencia, fill = Costo_Promedio)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(
    title = "Frecuencia y Montos de Siniestros por Marca",
    x = "Marca",
    y = "Frecuencia de Siniestros"
  ) +
  theme_minimal() +
  scale_fill_gradient(low = "skyblue4", high = "darkgrey") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
#---------------------------------------------------------------------
# Análisis de cobertura: Frecuencia y costos
analisis_cobertura <- base_modelos_para_R %>%
  group_by(Cobertura) %>%
  summarise(
    Frecuencia = n(),
    Costo_Total = sum(`Monto de Siniestros`, na.rm = TRUE),
    Costo_Promedio = mean(`Monto de Siniestros`, na.rm = TRUE)
  ) %>%
  arrange(desc(Frecuencia))

# Mostrar el resumen de las coberturas más frecuentes
print(head(analisis_cobertura, 10))

# Visualización: Coberturas más utilizadas y sus costos
ggplot(analisis_cobertura, aes(x = reorder(`Cobertura`, `Frecuencia`), y = `Frecuencia`, fill = `Costo_Total`)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(
    title = "Frecuencia y Costos por Tipo de Cobertura",
    x = "Tipo de Cobertura",
    y = "Frecuencia"
  ) +
  theme_minimal() +
  scale_fill_gradient(low = "azure2", high = "cyan4") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
#---------------------------------------------------------------------
# Correlación entre Monto de Siniestros y Monto Pagado
correlacion_montos <- cor(base_modelos_para_R$`Monto de Siniestros`, base_modelos_para_R$`Monto Pagado`, use = "complete.obs")

# Mostrar el coeficiente de correlación
print(paste("Coeficiente de correlación:", correlacion_montos))

# Visualización: Relación entre Monto de Siniestros y Monto Pagado
library(ggplot2)
ggplot(base_modelos_para_R, aes(x = `Monto de Siniestros`, y = `Monto Pagado`)) +
  geom_point(alpha = 0.6, color = "dodgerblue4") +
  geom_smooth(method = "lm", color = "darkgrey", se = FALSE) +
  labs(
    title = "Relación entre Monto de Siniestros y Monto Pagado",
    x = "Monto de Siniestros",
    y = "Monto Pagado"
  ) +
  theme_minimal()
#-------------------------------------------------------------------------
# Impacto de la causa del siniestro en los montos
impacto_causa <- base_modelos_para_R %>%
  group_by(`Causa del siniestro`) %>%
  summarise(
    Frecuencia = n(),
    Costo_Total = sum(`Monto de Siniestros`, na.rm = TRUE),
    Costo_Promedio = mean(`Monto de Siniestros`, na.rm = TRUE)
  ) %>%
  arrange(desc(Costo_Promedio))

# Mostrar las 10 causas con los costos promedio más altos
print(head(impacto_causa, 10))

# Visualización: Costo promedio por causa del siniestro
ggplot(impacto_causa, aes(x = reorder(`Causa del siniestro`, `Costo_Promedio`), y = `Costo_Promedio`)) +
  geom_bar(stat = "identity", fill = "darkseagreen", alpha = 0.7) +
  labs(
    title = "Costo Promedio por Causa del Siniestro",
    x = "Causa del Siniestro",
    y = "Costo Promedio (MXN)"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
#----------------------------------------------------------------
# Comparación de montos pagados según tipo de pérdida
comparacion_perdida <- base_modelos_para_R %>%
  group_by(`Tipo de Perdida`) %>%
  summarise(
    Costo_Total = sum(`Monto Pagado`, na.rm = TRUE),
    Costo_Promedio = mean(`Monto Pagado`, na.rm = TRUE),
    Frecuencia = n()
  ) %>%
  arrange(desc(Costo_Total))

# Mostrar el resumen por tipo de pérdida
print(comparacion_perdida)

# Visualización: Montos Pagados por Tipo de Pérdida
ggplot(comparacion_perdida, aes(x = reorder(`Tipo de Perdida`, `Costo_Promedio`), y = `Costo_Promedio`, fill = `Frecuencia`)) +
  geom_bar(stat = "identity", alpha = 0.8) +
  labs(
    title = "Montos Pagados por Tipo de Pérdida",
    x = "Tipo de Pérdida",
    y = "Costo Promedio (MXN)"
  ) +
  scale_fill_gradient(low = "lightblue", high = "darkblue") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
#--------------------------------------------------------------------
# Boxplot: Dispersión de montos pagados según Tipo de Pérdida
ggplot(base_modelos_para_R, aes(x = `Tipo de Perdida`, y = `Monto Pagado`, fill = `Tipo de Perdida`)) +
  geom_boxplot(outlier.color = "darkslategray4", outlier.size = 2, alpha = 0.7) +
  labs(
    title = "Dispersión de Montos Pagados según Tipo de Pérdida",
    x = "Tipo de Pérdida",
    y = "Monto Pagado (MXN)"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
#-----------------------------------------------------------------------
# Boxplot: Dispersión de montos pagados según Causa del Siniestro
ggplot(base_modelos_para_R, aes(x = `Causa del siniestro`, y = `Monto Pagado`, fill = `Causa del siniestro`)) +
  geom_boxplot(outlier.color = "darkblue", outlier.size = 2, alpha = 0.7) +
  labs(
    title = "Dispersión de Montos Pagados según Causa del Siniestro",
    x = "Causa del Siniestro",
    y = "Monto Pagado (MXN)"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
#----------------------------------------------------------------------
# Gráfico de barras: Frecuencia de siniestros según Marca
frecuencia_marca <- base_modelos_para_R %>%
  group_by(Marca) %>%
  summarise(Frecuencia = n()) %>%
  arrange(desc(Frecuencia))

ggplot(frecuencia_marca, aes(x = reorder(Marca, -Frecuencia), y = Frecuencia, fill = Frecuencia)) +
  geom_bar(stat = "identity", alpha = 0.8) +
  labs(
    title = "Frecuencia de Siniestros según Marca",
    x = "Marca",
    y = "Frecuencia de Siniestros"
  ) +
  scale_fill_gradient(low = "lightblue", high = "darkblue") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
#-----------------------------------------------------------------------
install.packages("readxl")

# Cargar librerías necesarias
library(readxl)
install.packages("dplyr")
library(MASS)  
library(VGAM)  
library(dplyr) # Para manipulación de datos
install.packages("survival")
library(fitdistrplus) 
library(ggplot2)     # Para gráficos
library(goftest)     # Para pruebas KS y Anderson-Darling




datos<-read_excel("DBMS.xlsx")
# Filtrar para eliminar niveles con solo un valor 
sapply(lapply(datos, unique), length)
######### Datos totales Frecuencia
# Modelo Poisson
modelo_poisson <- glm(`Numero de Vehiculos` ~  `Marca` + `Cobertura` + 
                        `Tipo de Perdida` + `Causa del siniestro`, 
                      data = datos, family = poisson(link = "log"))


# Modelo Binomial Negativa
modelo_bin_neg <- glm.nb(`Numero de Vehiculos` ~  `Marca` + `Cobertura` + 
                           `Tipo de Perdida` + `Causa del siniestro`, 
                         data = datos)

# Modelo Geométrico (caso especial de Binomial Negativa con theta = 1)
modelo_geom <- vglm(
  `Numero de Vehiculos` ~ `Marca` + `Cobertura` + `Tipo de Perdida` + `Causa del siniestro`,
  data = datos,
  family = geometric()  # o geometricff()
)

# Comparar modelos con AIC (menor es mejor)
aic_values <- data.frame(
  Modelo = c("Poisson", "Binomial Negativa", "Geométrica"),
  AIC = c(AIC(modelo_poisson), AIC(modelo_bin_neg), AIC(modelo_geom))
)

print(aic_values)

#### Recomendación
resumen_marcas <- datos %>%
  group_by(Marca) %>%
  summarise(
    media_frec = mean(`Numero de Vehiculos`, na.rm = TRUE),
    varianza_frec = var(`Numero de Vehiculos`, na.rm = TRUE),
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
    # Parámetros de la distribución de Poisson
    lambda_poisson = ifelse(recomendacion == "Poisson", media_frec, NA),
    
    # Parámetros de la distribución Binomial (Método de Momentos)
    n_binomial = ifelse(recomendacion == "Binomial", 
                        round((media_frec^2) / (varianza_frec - media_frec)), NA),
    p_binomial = ifelse(recomendacion == "Binomial", 
                        media_frec / n_binomial, NA),
    
    # Parámetros de la distribución Binomial Negativa o Geométrica (Método de Momentos)
    r_binomial_negativa = ifelse(recomendacion == "Binomial Negativa o Geométrica", 
                                 round((media_frec^2) / (varianza_frec - media_frec)), NA),
    # Asegurar que r_binomial_negativa sea al menos 1
    r_binomial_negativa = ifelse(r_binomial_negativa < 1, 1, r_binomial_negativa),
    
    p_binomial_negativa = ifelse(recomendacion == "Binomial Negativa o Geométrica", 
                                 media_frec / varianza_frec, NA)
  )


# Mostrar la tabla con los resultados
print(resumen_marcas)
################################


# Bucle para entrenar un modelo por cada marca
modelos_glm <- list()

for (marca in marcas_unicas) {
  # Filtrar datos de la marca específica
  datos_marca <- filter(datos, Marca == `Marca`)
  
  # Calcular media y varianza
  media_frec <- mean(datos_marca$`Numero de Vehiculos`)
  varianza_frec <- var(datos_marca$`Numero de Vehiculos`)
  
  # Elegir distribución según relación media-varianza
  if (abs(media_frec - varianza_frec) < 0.1 * media_frec) {
    distribucion <- poisson(link = "log")
    modelo <- glm(`Numero de Vehiculos` ~  `Cobertura` + `Tipo de Perdida` + `Causa del siniestro`,
                  data = datos_marca, family = distribucion)
  } else if (varianza_frec < media_frec) {
    distribucion <- binomial(link = "logit")
    modelo <- glm(`Numero de Vehiculos` ~  `Cobertura` + `Tipo de Perdida` + `Causa del siniestro`,
                  data = datos_marca, family = distribucion)
  } else {
    distribucion <- negative.binomial(theta = 1)
    modelo <- glm.nb(`Numero de Vehiculos` ~  `Cobertura` + `Tipo de Perdida` + `Causa del siniestro`,
                     data = datos_marca)
  }
  
  # Guardar el modelo en la lista con el nombre de la marca
  modelos_glm[[marca]] <- modelo
  
  # Imprimir resumen del modelo
  cat("\nResumen del modelo para la marca:", marca, "\n")
  print(summary(modelo))
}



########### Severidad #######
# Obtener la lista de marcas únicas
marcas_unicas <- unique(datos$Marca)

# Crear una lista vacía para guardar los parámetros por marca
resultados_severidad <- list()

# Función para estimar parámetros por Método de Momentos
ajustar_mm <- function(datos) {
  media <- mean(datos)
  varianza <- var(datos)
  
  # Exponencial (lambda = 1/media)
  lambda_exp <- 1 / media
  
  # Gamma (shape = media^2 / varianza, scale = varianza / media)
  shape_gamma <- media^2 / varianza
  scale_gamma <- varianza / media
  
  # Pareto (alpha = 2 + media^2 / (varianza - media^2), beta = media * (alpha - 1))
  alpha_pareto <- ifelse(varianza > media^2, 2 + media^2 / (varianza - media^2), NA)
  beta_pareto <- ifelse(!is.na(alpha_pareto), media * (alpha_pareto - 1), NA)
  
  # Burr (aproximación, ajuste simplificado)
  c_burr <- 1 + media^2 / varianza
  k_burr <- 2
  
  return(list(
    Exponencial = list(lambda = lambda_exp),
    Gamma = list(shape = shape_gamma, scale = scale_gamma),
    Pareto = list(alpha = alpha_pareto, beta = beta_pareto),
    Burr = list(c = c_burr, k = k_burr)
  ))
}

# Bucle para entrenar modelos por cada marca
for (marca in marcas_unicas) {
  cat("\n-------------------------------\n")
  cat("Modelando severidad para la marca:", marca, "\n")
  cat("-------------------------------\n")
  
  # Filtrar datos de la marca específica
  datos_marca <- filter(datos, Marca == marca)$`Monto de Siniestros`
  
  # Filtrar valores negativos o ceros
  datos_marca <- datos_marca[datos_marca > 0]
  
  # Verificar si hay suficientes datos
  if (length(datos_marca) < 10) {
    cat("⚠ Insuficientes datos para modelar severidad en", marca, "\n")
    next
  }
  
  # 📌 Método de Momentos
  parametros_mm <- ajustar_mm(datos_marca)
  
  # 📌 Máxima Verosimilitud con validaciones
  fit_exp <- tryCatch(fitdistr(datos_marca, "exponential"), error = function(e) NULL)
  fit_gamma <- tryCatch(fitdistr(datos_marca, "gamma", start = list(shape = 2, rate = 1 / mean(datos_marca))), error = function(e) NULL)
  fit_pareto <- tryCatch(vglm(datos_marca ~ 1, family = paretoII, data = data.frame(datos_marca)), error = function(e) NULL)
  fit_burr <- tryCatch(vglm(datos_marca ~ 1, family = burr, data = data.frame(datos_marca)), error = function(e) NULL)
  
  parametros_mle <- list(
    Exponencial = if (!is.null(fit_exp)) fit_exp$estimate else "Error",
    Gamma = if (!is.null(fit_gamma)) fit_gamma$estimate else "Error",
    Pareto = if (!is.null(fit_pareto)) coef(fit_pareto) else "Error",
    Burr = if (!is.null(fit_burr)) coef(fit_burr) else "Error"
  )
  
  # Guardar resultados
  resultados_severidad[[marca]] <- list(
    Metodo_Momentos = parametros_mm,
    Maxima_Verosimilitud = parametros_mle
  )
  
  print(resultados_severidad[[marca]])
}
