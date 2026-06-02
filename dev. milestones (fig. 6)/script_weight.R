# =========================
# Weight analysis 
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

df <- read_excel("weight.xlsx")
head(df)

df <- read_excel("weight.xlsx") %>%
  filter(!is.na(weight)) %>%
  mutate(
    animal = factor(animal),
    genotype = factor(genotype, levels = c("WT", "KO")),
    day = as.numeric(day)
  )

str(df)

library(lme4)
library(lmerTest)

model <- lmer(weight ~ genotype * day + (day| animal), data = df)
summary(model)

library(dplyr)
library(ggplot2)

plot_df <- df %>%
  filter(!is.na(weight)) %>%
  group_by(genotype, day) %>%
  summarise(
    mean_weight = mean(weight),
    sd = sd(weight),
    n = n(),
    se = sd / sqrt(n),
    .groups = "drop"
  )

p <- ggplot(plot_df, aes(x = day, y = mean_weight, color = genotype, group = genotype)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = mean_weight - se, ymax = mean_weight + se), width = 0.2) +
  labs(
    title = "Postnatal weight in WT and KO pups",
    x = "Postnatal day",
    y = "Weight (g)",
    color = "Genotype"
  ) +
  theme_classic(base_size = 13)

print(p)

df %>%
  distinct(animal, genotype) %>%
  count(genotype)

