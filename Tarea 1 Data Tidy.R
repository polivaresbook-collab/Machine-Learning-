#Tarea 1
#Integrantes Allen Cruz, Pedro Olivares

install.packages("tidyverse")
library(tidyverse)


install.packages("readr")
library(readr)


## Cargamos base
df_allocation <- read.csv("C:/Users/oliva/Desktop/UNIVERSIDAD/MAGISTER FINANZAS/MACHINE LEARNING/T1/data/asset_allocation.csv")
 
## Transformar adecuadamente el tipo y formato de las variables.

df_allocation_largos <- df_allocation %>%
  pivot_longer(
    cols = -Date,
    names_to = c("Portfolio", "Activo"),
    values_to = "Peso",
    names_pattern = "^(Portfolio_\\d+)\\.(.+)$"
  ) %>%
  mutate(Date = dmy(Date))

head(df_allocation_largos)

#Lo mismo con la base de retornos
df_returns <- read.csv("C:/Users/oliva/Desktop/UNIVERSIDAD/MAGISTER FINANZAS/MACHINE LEARNING/T1/data/asset_returns.csv")

df_returns_largos <- df_returns %>%
  pivot_longer(
    cols = -Date,
    names_to = "Activo",
    values_to = "Retorno"
  ) %>%
  mutate(Date = ymd(Date), Retorno=Retorno/100)

head(df_returns_largos)

##############################

#Calculo de los retornos acumulados


# 1. Unir ambas bases por Date y Activo
df_merged <- df_allocation_largos %>%
  left_join(df_returns_largos, by = c("Date", "Activo"))

head(df_merged)

# 2. Calcular retorno ponderado por activo (peso * retorno)
df_retorno_portafolios <- df_merged %>%
  mutate(Retorno_ponderado = Peso * Retorno) %>%
  group_by(Portfolio, Date) %>%
  summarise(Retorno_periodo = sum(Retorno_ponderado, na.rm = TRUE), .groups = "drop")

head(df_retorno_portafolios)

# 3. Calcular retorno acumulado
df_retorno_acumulado <- df_retorno_portafolios %>%
  arrange(Portfolio, Date) %>%
  group_by(Portfolio) %>%
  mutate(Retorno_acumulado = cumprod(1 + Retorno_periodo / 100) - 1)

head(df_retorno_acumulado)
# Verificar

# Gráfico
ggplot(df_retorno_acumulado, aes(x = Date, y = Retorno_acumulado, color = Portfolio)) +
  geom_line(size = 0.7) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_x_date(date_breaks = "6 months", date_labels = "%b %Y") +
  labs(
    title = "Retorno Acumulado por Portafolio",
    x = NULL,
    y = "Retorno Acumulado",
    color = "Portafolio"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
