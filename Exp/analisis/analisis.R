library(tidyverse)

data <- read_csv("../results/results.csv") %>%
  separate(algoritmo, c("alineamiento", "lenguaje", "implementacion"), remove = F) %>%
  filter( # Esta secuencia estaba mal al momento de correr los exp el 16/1
    secuencia1 != "../sequences/genomes/H1N1_S8_JX046930.fasta",
    secuencia2 != "../sequences/genomes/H1N1_S8_JX046930.fasta"
  ) %>%
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

ggsave("graficos/NW_random_t_c.pdf", width = 9, height = 5)

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

ggsave("graficos/NW_random_t_l.pdf", width = 9, height = 5)

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

ggsave("graficos/NW_random_gcups_c.pdf", width = 9, height = 5)

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

ggsave("graficos/NW_random_gcups_l.pdf", width = 9, height = 5)


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

ggsave("graficos/SW_random_t_c.pdf", width = 9, height = 5)

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

ggsave("graficos/SW_random_t_l.pdf", width = 9, height = 5)

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

ggsave("graficos/SW_random_gcups_c.pdf", width = 9, height = 5)

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

ggsave("graficos/SW_random_gcups_l.pdf", width = 9, height = 5)

# Experimento 3: NW con genomas reales de distinto largo -----------------------
NW_genomes <- data %>%
  filter(experimento == "NW_genomes") %>%
  extract(secuencia1, "genoma1", regex = "([[:alnum:]-]+[_[:alnum:]]*)_", remove = F) %>%
  extract(secuencia2, "genoma2", regex = "([[:alnum:]-]+[_[:alnum:]]*)_", remove = F) %>%
  mutate(
  genoma1 = case_when(
    genoma1 == "Ebola_Sudan" ~ "Ebola",
    TRUE ~ genoma1
  ),
  genoma2 = case_when(
    genoma2 == "Ebola_Sudan" ~ "Ebola",
    TRUE ~ genoma2
  )) %>%
  filter(genoma1 == genoma2) %>%
  mutate(genoma = genoma1)

NW_genomes %>%
  group_by(algoritmo, genoma) %>%
  summarise(
    tiempo_promedio = mean(tiempo_seg),
    tiempo_sd = sd(tiempo_seg)
  ) %>%
  ggplot(aes(x = genoma, y = tiempo_promedio, group = algoritmo, color = algoritmo)) +
  geom_line() +
  geom_errorbar(aes(ymin = tiempo_promedio - tiempo_sd, ymax = tiempo_promedio + tiempo_sd), width = .1) +
  geom_point() +
  scale_x_discrete(
    limits = c("H1N1_S8", "H1N1_S6", "H1N1_S1", "HPV", "HIV1", "Zika", "Ebola", "Marburg", "SARS-CoV-2", "MERS"),
    labels = c("H1N1_S8\n(863 pb)", "H1N1_S6\n(1410 pb)", "H1N1_S1\n(2280 pb)", "HPV\n(7320 pb)", "HIV1\n(9181 pb)", "Zika\n(10808 pb)", "Ebola\n(18875 pb)", "Marburg\n(19114 pb)", "SARS-CoV-2\n(29903 pb)", "MERS\n(30119 pb)"),
    ) +
  ggtitle("Needleman-Wunsch: genomas virales") +
  xlab("Genoma") +
  ylab("Tiempo (segundos)") +
  theme_bw() +
  theme(
    panel.border = element_blank(),
    axis.line = element_line()
  )

ggsave("graficos/NW_genomes_t.pdf", width = 9, height = 5)

NW_genomes %>%
  group_by(algoritmo, genoma) %>%
  summarise(
    gcups_promedio = mean(GCUPS),
    gcups_sd = sd(GCUPS)
  ) %>%
  ggplot(aes(x = genoma, y = gcups_promedio, group = algoritmo, color = algoritmo)) +
  geom_line() +
  geom_errorbar(aes(ymin = gcups_promedio - gcups_sd, ymax = gcups_promedio + gcups_sd), width = .1) +
  geom_point() +
  scale_x_discrete(
    limits = c("H1N1_S8", "H1N1_S6", "H1N1_S1", "HPV", "HIV1", "Zika", "Ebola", "Marburg", "SARS-CoV-2", "MERS"),
    labels = c("H1N1_S8\n(863 pb)", "H1N1_S6\n(1410 pb)", "H1N1_S1\n(2280 pb)", "HPV\n(7320 pb)", "HIV1\n(9181 pb)", "Zika\n(10808 pb)", "Ebola\n(18875 pb)", "Marburg\n(19114 pb)", "SARS-CoV-2\n(29903 pb)", "MERS\n(30119 pb)"),
  ) +
  ggtitle("Needleman-Wunsch: genomas virales") +
  xlab("Genoma") +
  ylab("Giga celdas por segundo") +
  theme_bw() +
  theme(
    panel.border = element_blank(),
    axis.line = element_line()
  )

ggsave("graficos/NW_genomes_gcups.pdf", width = 9, height = 5)

# Experimento 4: SW reads simulados con genomas reales de distinto largo -------
SW_reads <- data %>%
  filter(experimento == "SW_reads") %>%
  extract(secuencia1, "genoma", regex = "([[:alnum:]-]+[_[:alnum:]]*)_", remove = F) %>%
  mutate(
    genoma = case_when(
      genoma == "Ebola_Sudan" ~ "Ebola",
      TRUE ~ genoma
    ))

SW_reads %>%
  group_by(algoritmo, genoma) %>%
  summarise(
    tiempo_promedio = mean(tiempo_seg),
    tiempo_sd = sd(tiempo_seg)
  ) %>%
  ggplot(aes(x = genoma, y = tiempo_promedio, group = algoritmo, color = algoritmo)) +
  geom_line() +
  geom_errorbar(aes(ymin = tiempo_promedio - tiempo_sd, ymax = tiempo_promedio + tiempo_sd), width = .1) +
  geom_point() +
  scale_x_discrete(
    limits = c("H1N1_S8", "H1N1_S6", "H1N1_S1", "HPV", "HIV1", "Zika", "Ebola", "Marburg", "SARS-CoV-2", "MERS"),
    labels = c("H1N1_S8\n(863 pb)", "H1N1_S6\n(1410 pb)", "H1N1_S1\n(2280 pb)", "HPV\n(7320 pb)", "HIV1\n(9181 pb)", "Zika\n(10808 pb)", "Ebola\n(18875 pb)", "Marburg\n(19114 pb)", "SARS-CoV-2\n(29903 pb)", "MERS\n(30119 pb)"),
  ) +
  ggtitle("Smith-Waterman: lecturas simuladas") +
  xlab("Genoma") +
  ylab("Tiempo (segundos)") +
  theme_bw() +
  theme(
    panel.border = element_blank(),
    axis.line = element_line()
  )

ggsave("graficos/SW_reads_t.pdf", width = 9, height = 5)

SW_reads %>%
  group_by(algoritmo, genoma) %>%
  summarise(
    gcups_promedio = mean(GCUPS),
    gcups_sd = sd(GCUPS)
  ) %>%
  ggplot(aes(x = genoma, y = gcups_promedio, group = algoritmo, color = algoritmo)) +
  geom_line() +
  geom_errorbar(aes(ymin = gcups_promedio - gcups_sd, ymax = gcups_promedio + gcups_sd), width = .1) +
  geom_point() +
  scale_x_discrete(
    limits = c("H1N1_S8", "H1N1_S6", "H1N1_S1", "HPV", "HIV1", "Zika", "Ebola", "Marburg", "SARS-CoV-2", "MERS"),
    labels = c("H1N1_S8\n(863 pb)", "H1N1_S6\n(1410 pb)", "H1N1_S1\n(2280 pb)", "HPV\n(7320 pb)", "HIV1\n(9181 pb)", "Zika\n(10808 pb)", "Ebola\n(18875 pb)", "Marburg\n(19114 pb)", "SARS-CoV-2\n(29903 pb)", "MERS\n(30119 pb)"),
  ) +
  ggtitle("Smith-Waterman: lecturas simuladas") +
  xlab("Genoma") +
  ylab("Giga celdas por segundo") +
  theme_bw() +
  theme(
    panel.border = element_blank(),
    axis.line = element_line()
  )

ggsave("graficos/SW_reads_gcups.pdf", width = 9, height = 5)
