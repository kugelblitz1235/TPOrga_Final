library(tidyverse)

data <- read_csv("../results/results_03.csv") %>%
  separate(algoritmo, c("alineamiento", "lenguaje", "implementacion"), remove = F) %>%
  rename(
    longitud1_interna = longitud1,
    longitud2_interna = longitud2
  ) %>%
  mutate(
    tiempo_seg = tiempo / 10^6,
    celdas_efectivas = longitud1_interna * longitud2_interna,
    GCUPS = celdas_efectivas / (10^9 * tiempo_seg),
    longitud1 = longitud1_interna-1,
    longitud2 = longitud2_interna-1
  )

# Experimento 1: NW con secuencias random de distinto largo --------------------
NW_random <- data %>%
  filter(experimento == "NW_random")

NW_random %>%
  group_by(algoritmo, longitud1, longitud2, celdas_efectivas) %>%
  summarise(
    tiempo_promedio = mean(tiempo_seg),
    tiempo_sd = sd(tiempo_seg)
  ) %>%
  ggplot(aes(x = celdas_efectivas, y = tiempo_promedio, group = algoritmo, color = algoritmo)) +
  geom_line() +
  geom_errorbar(aes(ymin = tiempo_promedio - tiempo_sd, ymax = tiempo_promedio + tiempo_sd)) +
  ggtitle("Needleman-Wunsch: secuencias aleatorias") +
  xlab("Cantidad de celdas") +
  ylab("Tiempo (segundos)") +
  theme_bw() +
  theme(
    panel.border = element_blank(),
    axis.line = element_line()
  )

ggsave("graficos/NW_random_t_c_O3.pdf", width = 9, height = 5)

NW_random %>%
  group_by(algoritmo, longitud1, longitud2, celdas_efectivas) %>%
  summarise(
    tiempo_promedio = mean(tiempo_seg),
    tiempo_sd = sd(tiempo_seg)
  ) %>%
  ggplot(aes(x = longitud1, y = tiempo_promedio, group = algoritmo, color = algoritmo)) +
  geom_line() +
  geom_errorbar(aes(ymin = tiempo_promedio - tiempo_sd, ymax = tiempo_promedio + tiempo_sd), width = 50) +
  ggtitle("Needleman-Wunsch: secuencias aleatorias") +
  xlab("Longitud de las secuencias") +
  ylab("Tiempo (segundos)") +
  theme_bw() +
  theme(
    panel.border = element_blank(),
    axis.line = element_line()
  )

ggsave("graficos/NW_random_t_l_O3.pdf", width = 9, height = 5)

NW_random %>%
  group_by(algoritmo, longitud1, longitud2, celdas_efectivas) %>%
  summarise(
    gcups_promedio = mean(GCUPS),
    gcups_sd = sd(GCUPS)
  ) %>%
  ggplot(aes(x = celdas_efectivas, y = gcups_promedio, group = algoritmo, color = algoritmo)) +
  geom_line() +
  geom_errorbar(aes(ymin = gcups_promedio - gcups_sd, ymax = gcups_promedio + gcups_sd)) +
  ggtitle("Needleman-Wunsch: secuencias aleatorias") +
  xlab("Cantidad de celdas") +
  ylab("Giga celdas por segundo") +
  theme_bw() +
  theme(
    panel.border = element_blank(),
    axis.line = element_line()
  )

ggsave("graficos/NW_random_gcups_c_O3.pdf", width = 9, height = 5)

NW_random %>%
  group_by(algoritmo, longitud1, longitud2, celdas_efectivas) %>%
  summarise(
    gcups_promedio = mean(GCUPS),
    gcups_sd = sd(GCUPS)
  ) %>%
  ggplot(aes(x = longitud1, y = gcups_promedio, group = algoritmo, color = algoritmo)) +
  geom_line() +
  geom_errorbar(aes(ymin = gcups_promedio - gcups_sd, ymax = gcups_promedio + gcups_sd), width = 50) +
  ggtitle("Needleman-Wunsch: secuencias aleatorias") +
  xlab("Longitud de las secuencias") +
  ylab("Giga celdas por segundo") +
  theme_bw() +
  theme(
    panel.border = element_blank(),
    axis.line = element_line()
  )

ggsave("graficos/NW_random_gcups_l_O3.pdf", width = 9, height = 5)


# Experimento 2: SW con secuencias random de distinto largo --------------------
SW_random <- data %>%
  filter(experimento == "SW_random")

SW_random %>%
  group_by(algoritmo, longitud1, longitud2, celdas_efectivas) %>%
  summarise(
    tiempo_promedio = mean(tiempo_seg),
    tiempo_sd = sd(tiempo_seg)
  ) %>%
  ggplot(aes(x = celdas_efectivas, y = tiempo_promedio, group = algoritmo, color = algoritmo)) +
  geom_line() +
  geom_errorbar(aes(ymin = tiempo_promedio - tiempo_sd, ymax = tiempo_promedio + tiempo_sd)) +
  ggtitle("Smith-Waterman: secuencias aleatorias") +
  xlab("Cantidad de celdas") +
  ylab("Tiempo (segundos)") +
  theme_bw() +
  theme(
    panel.border = element_blank(),
    axis.line = element_line()
  )

ggsave("graficos/SW_random_t_c_O3.pdf", width = 9, height = 5)

SW_random %>%
  group_by(algoritmo, longitud1, longitud2, celdas_efectivas) %>%
  summarise(
    tiempo_promedio = mean(tiempo_seg),
    tiempo_sd = sd(tiempo_seg)
  ) %>%
  ggplot(aes(x = longitud1, y = tiempo_promedio, group = algoritmo, color = algoritmo)) +
  geom_line() +
  geom_errorbar(aes(ymin = tiempo_promedio - tiempo_sd, ymax = tiempo_promedio + tiempo_sd), width = 50) +
  ggtitle("Smith-Waterman: secuencias aleatorias") +
  xlab("Longitud de las secuencias") +
  ylab("Tiempo (segundos)") +
  theme_bw() +
  theme(
    panel.border = element_blank(),
    axis.line = element_line()
  )

ggsave("graficos/SW_random_t_l_O3.pdf", width = 9, height = 5)

SW_random %>%
  group_by(algoritmo, longitud1, longitud2, celdas_efectivas) %>%
  summarise(
    gcups_promedio = mean(GCUPS),
    gcups_sd = sd(GCUPS)
  ) %>%
  ggplot(aes(x = celdas_efectivas, y = gcups_promedio, group = algoritmo, color = algoritmo)) +
  geom_line() +
  geom_errorbar(aes(ymin = gcups_promedio - gcups_sd, ymax = gcups_promedio + gcups_sd)) +
  ggtitle("Smith-Waterman: secuencias aleatorias") +
  xlab("Cantidad de celdas") +
  ylab("Giga celdas por segundo") +
  theme_bw() +
  theme(
    panel.border = element_blank(),
    axis.line = element_line()
  )

ggsave("graficos/SW_random_gcups_c_O3.pdf", width = 9, height = 5)

SW_random %>%
  group_by(algoritmo, longitud1, longitud2, celdas_efectivas) %>%
  summarise(
    gcups_promedio = mean(GCUPS),
    gcups_sd = sd(GCUPS)
  ) %>%
  ggplot(aes(x = longitud1, y = gcups_promedio, group = algoritmo, color = algoritmo)) +
  geom_line() +
  geom_errorbar(aes(ymin = gcups_promedio - gcups_sd, ymax = gcups_promedio + gcups_sd), width = 50) +
  ggtitle("Smith-Waterman: secuencias aleatorias") +
  xlab("Longitud de las secuencias") +
  ylab("Giga celdas por segundo") +
  theme_bw() +
  theme(
    panel.border = element_blank(),
    axis.line = element_line()
  )

ggsave("graficos/SW_random_gcups_l_O3.pdf", width = 9, height = 5)