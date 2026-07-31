# Scratch script for testing plots and loading data

library(ggplot2)
library(ggdist)
library(dplyr)
library(tidyr)
library(bayesplot)
library(readr)

df_exponential_growth <- data.frame(
  cell_number = c(15000 , 98000),
  time = c(12 , 192)
)
df_exponential_growth
df_exponential_growth |> 
  ggplot(aes(x = time, y = cell_number)) +
  geom_point(size = 3) +
    labs(x = "Elapsed Time (hours)", y = "Cell Count")

Interpolated_exponential <- read_csv("Interpolated_exponential.csv")
Interpolated_exponential

Interpolated_exponential |> 
  ggplot(aes(x = Time, y = Cell_count, color = Type)) +
  geom_point(size = 3) +
  labs(x = "Elapsed Time (hours)", y = "Cell Count")

Cell_growth <- read_csv("Cell_growth.csv")
View(Cell_growth)

ggplot(Interpolated_exponential, aes(Time, Cell_count)) + 
  geom_line(linewidth = 1.5) + 
  geom_point(data = Cell_growth, colour = "blue", size = 3) +
  labs(x = "Elapsed Time (hours)", y = "Cell Count")