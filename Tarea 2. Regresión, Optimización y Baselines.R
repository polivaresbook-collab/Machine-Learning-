
#Tarea 2: Regresión, Optimización y Baselines
#Integrantes: Allen Cruz, Pedro Olivares, Joel Vásquez, Joshua Barrientos

install.packages("frenchdata")
install.packages("quantmod")
install.packages("tidyverse")
install.packages("TTR")
install.packages("timetk")

library(frenchdata)   
library(quantmod)     
library(tidyverse)    
library(TTR)          
library(timetk)       # para convertir xts a tibble
library(ggplot2)      
library(knitr)        # para tablas formateadas

# PARTE 1: Construcción del dataset de factores

# PASO 1: Descargar factores Fama-French diarios
Base_Factores_FF     <- download_french_data("Fama/French 3 Factors [Daily]")
Factores_Fama_French <- Base_Factores_FF$subsets$data[[1]]

Factores_limpios_FF <- Factores_Fama_French |>
  rename(
    mkt_rf = `Mkt-RF`,
    smb    = SMB,
    hml    = HML,
    rf     = RF
  ) |>
  mutate(
    date   = ymd(date),
    mkt_rf = mkt_rf / 100,
    smb    = smb    / 100,
    hml    = hml    / 100,
    rf     = rf     / 100
  )

#Precios de Apple
from_famafrench <- "1990-01-01"

Base_Apple <- getSymbols("AAPL",
                         src         = "yahoo",
                         from        = from_famafrench,
                         to          = Sys.Date(),
                         auto.assign = FALSE)

#Nos quedamos solo con fecha y precio ajustado
Precios_Apple <- Base_Apple |>
  timetk::tk_tbl(rename_index = "date") |>
  select(date, price = AAPL.Adjusted)

#Retornos diarios discretos de AAPL 
# Retorno discreto: (Pt - P(t-1)) / P(t-1)
Retornos_Apple <- Precios_Apple |>
  mutate(ri = TTR::ROC(price, type = "discrete")) |>
  select(date, ri) |>
  drop_na()

#Unimos bases y calculamos exceso de retorno
#inner_join conserva solo fechas presentes en ambos datasets
#Esto elimina automáticamente fines de semana y feriados
dataset <- Retornos_Apple |>
  inner_join(Factores_limpios_FF, by = "date") |>
  mutate(
    excess_ret = ri - rf    
  ) |>
  select(date, excess_ret, mkt_rf, smb, hml) |>
  drop_na()

# Verificación del dataset final
glimpse(dataset)
head(dataset)
summary(dataset)


# PARTE 2: Estimación por OLS


#  Método 1: lm() 
modelo_ols <- lm(excess_ret ~ mkt_rf + smb + hml, data = dataset)
summary(modelo_ols)

#  Método 2: Fórmula cerrada matricial θ = (XᵀX)⁻¹Xᵀy 
X <- model.matrix(~ mkt_rf + smb + hml, data = dataset)
y <- dataset$excess_ret

theta_ols <- solve(t(X) %*% X) %*% t(X) %*% y

cat("Coeficientes OLS (fórmula cerrada):\n")
print(theta_ols)

cat("\nCoeficientes OLS (lm):\n")
print(coef(modelo_ols))

cat("\nDiferencia entre métodos (debe ser ~0):\n")
print(theta_ols[,1] - coef(modelo_ols))

#  Generación de Tabla 2.3: Resultados OLS IA 

# 1. Crear un dataframe con los resultados de ambos métodos
# Nota: theta_ols es una matriz de 4x1, por eso usamos as.numeric
tabla_comparativa_ols <- data.frame(
  Coeficiente = c("(Intercepto)", "mkt_rf", "smb", "hml"),
  Metodo_lm   = as.numeric(coef(modelo_ols)),
  Metodo_Matricial = as.numeric(theta_ols),
  Diferencia  = as.numeric(coef(modelo_ols)) - as.numeric(theta_ols)
)

# 2. Renderizar la tabla con knitr::kable
cat("\n=== Comparación de Métodos OLS: lm() vs. Fórmula Cerrada ===\n")
knitr::kable(
  tabla_comparativa_ols, 
  digits = 15,          # Usamos muchos decimales para ver la diferencia de precisión de máquina
  format = "simple",    # O "pipe" si prefieres formato Markdown estándar
  col.names = c("Coeficiente", "Estimación lm()", "Fórmula Matricial", "Diferencia (E-16)"),
  align = "lccc",
  caption = "Comparación de Coeficientes OLS"
)

# ============================================================
# PARTE 3: Batch Gradient Descent desde cero
# ============================================================

#  Preparar matrices X e y 
# Se agrega columna de unos explícitamente para el intercepto
X <- dataset |>
  mutate(intercepto = 1) |>
  select(intercepto, mkt_rf, smb, hml) |>
  as.matrix()

y <- dataset$excess_ret
n <- nrow(X)

#  3.1 Función de costo J(θ) 
# J(θ) = 1/(2n) * Σ(Xθ - y)²
costo <- function(X, y, theta) {
  residuos <- X %*% theta - y
  return((1 / (2 * n)) * sum(residuos^2))
}

#  3.2 Gradiente de J(θ) 
# ∇J(θ) = (1/n) * Xᵀ(Xθ - y)
gradiente <- function(X, y, theta) {
  residuos <- X %*% theta - y
  return((1 / n) * t(X) %*% residuos)
}

#  3.3 Algoritmo Batch Gradient Descent 
batch_gd <- function(X, y,
                     alpha,
                     max_iter   = 50000,
                     tolerancia = 1e-10) {
  
  # Inicializar theta en ceros
  theta        <- matrix(0, nrow = ncol(X))
  historial    <- numeric(max_iter)
  convergencia <- FALSE
  iter_final   <- max_iter
  
  for (t in 1:max_iter) {
    
    # Calcular costo actual
    J_actual <- costo(X, y, theta)
    
    # Protección contra divergencia numérica
    if (is.nan(J_actual) || is.infinite(J_actual)) {
      cat(sprintf("  [!] Divergencia detectada en iteración %d (alpha = %g)\n", t, alpha))
      iter_final <- t
      break
    }
    
    historial[t] <- J_actual
    
    # Actualizar theta
    theta_nuevo <- theta - alpha * gradiente(X, y, theta)
    
    # Criterio de convergencia: cambio en J(θ) menor a la tolerancia
    delta <- abs(costo(X, y, theta_nuevo) - J_actual)
    if (delta < tolerancia) {
      convergencia <- TRUE
      iter_final   <- t
      theta        <- theta_nuevo
      historial    <- historial[1:t]
      break
    }
    
    theta <- theta_nuevo
  }
  
  # Recortar historial si no convergió
  if (!convergencia && iter_final == max_iter) {
    historial <- historial[1:max_iter]
  }
  
  return(list(
    theta        = theta,
    historial    = historial,
    convergencia = convergencia,
    iteraciones  = iter_final,
    J_final      = costo(X, y, theta),
    alpha        = alpha
  ))
}

# ============================================================
# PARTE 4: Experimentos con distintos learning rates
# ============================================================

learning_rates <- c(0.01, 0.1, 1, 2, 3)

# Ejecutar GD para cada learning rate y guardar resultados
resultados_gd <- lapply(learning_rates, function(lr) {
  cat(sprintf("Ejecutando GD con alpha = %.2f...\n", lr))
  batch_gd(X, y, alpha = lr, max_iter = 50000, tolerancia = 1e-10)
})

names(resultados_gd) <- paste0("alpha_", learning_rates)

# ============================================================
# PARTE 5: Resultados
# ============================================================

# -
# 5.1 Tabla de convergencia
# -

tabla_convergencia <- map_dfr(resultados_gd, function(res) {
  tibble(
    `Learning Rate`    = res$alpha,
    `Convergió`        = ifelse(res$convergencia, "Sí", "No"),
    `Iteraciones`      = res$iteraciones,
    `Costo Final J(θ)` = res$J_final
  )
})

cat("\n=== Tabla de Convergencia ===\n")
kable(tabla_convergencia, digits = 10, align = "c")

# -
# 5.2 Evolución de la función de costo (escala logarítmica)
# -

# Construir dataframe con todos los historiales de costo
historiales <- map_dfr(resultados_gd, function(res) {
  tibble(
    iteracion     = seq_along(res$historial),
    costo         = res$historial,
    learning_rate = factor(res$alpha)
  )
})

# Filtrar valores positivos y finitos para escala logarítmica
historiales_validos <- historiales |>
  filter(is.finite(costo) & costo > 0)

grafico_convergencia <- ggplot(historiales_validos,
                               aes(x = iteracion, y = costo,
                                   color = learning_rate)) +
  geom_line(linewidth = 4) +
  scale_y_log10() +
  scale_color_brewer(palette = "Set1") +
  labs(
    title  = "Evolución de la Función de Costo J(θ) por Learning Rate",
    x      = "Iteración",
    y      = "J(θ) (escala logarítmica)",
    color  = "Learning Rate (α)"
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "bottom")

print(grafico_convergencia)

# -
# 5.3 Comparación de coeficientes GD vs OLS
# -

nombres_coef <- c("α (intercepto)", "β₁ (mkt_rf)", "β₂ (smb)", "β₃ (hml)")
coef_ols     <- as.numeric(theta_ols)

# Construir tabla comparativa
tabla_coefs <- map_dfr(resultados_gd, function(res) {
  
  coef_gd <- as.numeric(res$theta)
  
  # Si hubo divergencia, theta puede ser Inf o NaN
  if (any(!is.finite(coef_gd))) {
    return(tibble(
      `Learning Rate` = res$alpha,
      Coeficiente     = nombres_coef,
      OLS             = coef_ols,
      GD              = NA_real_,
      Diferencia      = NA_real_
    ))
  }
  
  tibble(
    `Learning Rate` = res$alpha,
    Coeficiente     = nombres_coef,
    OLS             = coef_ols,
    GD              = coef_gd,
    Diferencia      = coef_gd - coef_ols
  )
})

cat("\n=== Comparación de Coeficientes: GD vs OLS ===\n")
kable(tabla_coefs, digits = 8, align = "c")

# Gráfico de diferencias absolutas (solo learning rates convergentes)
tabla_coefs_validos <- tabla_coefs |>
  filter(!is.na(Diferencia)) |>
  mutate(`Learning Rate` = factor(`Learning Rate`))

grafico_diferencias <- ggplot(tabla_coefs_validos,
                              aes(x = Coeficiente, y = abs(Diferencia),
                                  fill = `Learning Rate`)) +
  geom_col(position = "dodge") +
  scale_fill_brewer(palette = "Set1") +
  scale_y_log10() +
  labs(
    title = "Diferencia Absoluta entre GD y OLS por Coeficiente",
    x     = "Coeficiente",
    y     = "|GD - OLS| (escala logarítmica)",
    fill  = "Learning Rate (α)"
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "bottom")

print(grafico_diferencias)
