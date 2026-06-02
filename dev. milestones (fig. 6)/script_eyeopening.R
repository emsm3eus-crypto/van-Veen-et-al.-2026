# =========================
# Eye opening analysis 
# =========================

library(rstudioapi)
dataFile <- selectDirectory(
  caption = "Select Directory",
  label = "Select",
  path = getActiveProject()
)
setwd(dataFile)
getwd()
list.files()

library(readxl)
library(dplyr)
library(ordinal)
library(ggplot2)

list.files()

df <- read_excel("eye_opening.xlsx")
head(df)

df <- df %>%
  filter(!is.na(animal), !is.na(genotype), !is.na(day), !is.na(score)) %>%
  mutate(
    animal = factor(animal),
    genotype = factor(genotype, levels = c("WT", "KO")),
    day = as.numeric(day),
    score = factor(score, levels = c(0, 1, 2), ordered = TRUE)
  )

# Quick check
str(df)
table(df$genotype, df$score)

model <- clmm(score ~ genotype * day + (1 | animal), data = df)
summary(model)

library(ggplot2)
library(dplyr)

plot_df <- df %>%
  mutate(score_num = as.numeric(as.character(score))) %>%
  group_by(genotype, day) %>%
  summarise(
    mean_score = mean(score_num),
    sd = sd(score_num),
    n = n(),
    se = sd / sqrt(n),
    .groups = "drop"
  )

ggplot(plot_df, aes(x = day, y = mean_score, color = genotype, group = genotype)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = mean_score - se, ymax = mean_score + se), width = 0.2) +
  scale_y_continuous(breaks = c(0, 1, 2), limits = c(0, 2)) +
  labs(
    x = "Postnatal day",
    y = "Mean eye-opening score",
    title = "Eye opening in WT and KO pups"
  ) +
  theme_classic(base_size = 13)

df %>%
  distinct(animal, genotype) %>%
  count(genotype)




