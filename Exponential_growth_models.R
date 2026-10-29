# Scratch script for testing plots and loading data

library(ggplot2)
library(ggdist)
library(dplyr)
library(tidyr)
library(bayesplot)
library(readr)
library(tidybayes)
library(cowplot)
library(rstan)
library(brms)
library(ggrepel)
library(RColorBrewer)
library(modelr)

library(magrittr)
library(dplyr)
library(purrr)
library(forcats)
library(tidyr)
library(modelr)
library(ggdist)
library(tidybayes)
library(ggplot2)
library(cowplot)
library(rstan)
library(brms)
library(ggrepel)
library(RColorBrewer)
library(gganimate)
library(posterior)
library(distributional)

theme_set(theme_tidybayes() + panel_border())

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


logistic_func <- function(t, L, U, T_D, M) {
  U + ((L - U)/(1 + 2^(T_D * (t - M))))
}

# Fit custom logistic function
Logistic_cell_growth <- read_csv("Logistic_cell_growth.csv", show_col_types = FALSE)
Logistic_cell_growth

# Define the logistic function: q_2 = U + ((L - U)/(1 + 2^(T_D * (t - M))))
logistic_func <- function(t, L, U, T_D, M) {
  U + ((L - U)/(1 + 2^(T_D * (t - M))))
}

# Fit using nls with manual start values
manual_fit <- nls(Cell_count ~ logistic_func(Time, L, U, T_D, M), 
                  data = Logistic_cell_growth, 
                  start = list(L = 5000, U = 250000, T_D = 0.023, M = 195)) 

# View parameter estimates
summary(manual_fit)
coef(manual_fit)

# Create a grid of x values for the line
logistic_fit_df <- data.frame(Time = seq(min(Logistic_cell_growth$Time), max(Logistic_cell_growth$Time), length.out = 100))

# Predict y values based on the fitted model
logistic_fit_df$y_pred <- predict(manual_fit, newdata = logistic_fit_df)   

# Plotting the fit
ggplot(Logistic_cell_growth, aes(Time, Cell_count)) + 
  geom_point(color = "green", size = 3) + 
  geom_line(data = logistic_fit_df, aes(Time, y_pred), linewidth = 1.5)
  labs(x = "Elapsed Time (hours)", y = "Cell Count")

  
Linearized_cell_counts <- read_csv("Linearized_cell_counts.csv")

Cells = brm(
  Cell_count ~ Material_param_1, 
  data = Linearized_cell_counts, 
)

plot(Cells)
print(Cells)
Cells$fit


Cells_brm = brm(
  data = Linearized_cell_counts, family = gaussian,
  Cell_count ~ 1 + Material_param_1,
  prior = c(prior(normal(30, 2), class = Intercept),
            prior(normal(0, 10), class = b),
            prior(cauchy(0, 1), class = sigma)),
  iter = 2000, warmup = 1000, chains = 4, cores = 4,
  seed = 4
)

plot(Cells_brm)
print(Cells_brm)
Cells_brm$fit


Linearized_cell_counts %>%
  data_grid(Material_param_1 = seq_range(Material_param_1, n = 101)) %>%
  add_predicted_draws(Cells_brm) %>%
  ggplot(aes(x = Material_param_1, y = Cell_count, color = Approval_status)) +
  stat_lineribbon(aes(y = .prediction), .width = c(.99, .95, .8, .5), color = brewer.pal(5, "Blues")[[5]]) +
  geom_point(data = Linearized_cell_counts) +
  scale_fill_brewer() +
  scale_color_manual(values = c("green", "red"))

Logistic_cell_growth <- read_csv("Logistic_cell_growth.csv", show_col_types = FALSE)
Logistic_cell_growth

Logistic_cell_growth |> 
  ggplot(aes(x = Time, y = Cell_count)) +
  geom_point(size = 3) +
  labs(x = "Elapsed Time (hours)", y = "Cell Count")

# Define the logistic function: q_2 = U + ((L - U)/(1 + 2^(T_D * (t - M))))
logistic_func <- function(t, L, U, T_D, M) {
  U + ((L - U)/(1 + 2^(T_D * (t - M))))
}

# Fit using nls with manual start values
manual_fit <- nls(Cell_count ~ logistic_func(Time, L, U, T_D, M), 
                  data = Logistic_cell_growth, 
                  start = list(L = 5000, U = 250000, T_D = 0.023, M = 195))   

# View parameter estimates
summary(manual_fit)

# Create a grid of x values for the line
logistic_fit_df <- data.frame(
  Time = seq(min(Logistic_cell_growth$Time), max(Logistic_cell_growth$Time), length.out = 100)
)

# Predict y values based on the fitted model
logistic_fit_df$y_pred <- predict(manual_fit, newdata = logistic_fit_df)   

# Plot the non-linear fit
ggplot(Logistic_cell_growth, aes(Time, Cell_count)) + 
  geom_point(size = 3) + 
  geom_line(data = logistic_fit_df, aes(Time, y_pred), linewidth = 1.5, color = "green") +
  labs(x = "Elapsed Time (hours)", y = "Cell Count")




logistic_fit_df %>%
  mutate(
    model1 = 5171 * 2^(0.023*Time),
    model2 = 5171 * 2^(0.022*Time),
    model3 = 5171 * 2^(0.019*Time),
    model4 = 5171 * 2^(0.015*Time)
  ) %>%
  ggplot(Logistic_cell_growth, aes(Time, Cell_count)) + 
  geom_point(size = 3) + 
  geom_line(data = logistic_fit_df, aes(Time, y_pred), linewidth = 1.5, color = "green") +
  geom_line(data = logistic_fit_df, aes(Time, model1), linewidth = 1.5, color = "blue") +
  geom_line(data = logistic_fit_df, aes(Time, model2), linewidth = 1.5, color = "red") +
  geom_line(data = logistic_fit_df, aes(Time, model3), linewidth = 1.5, color = "sienna1") +
  geom_line(data = logistic_fit_df, aes(Time, model4), linewidth = 1.5, color = "hotpink4") +
  labs(x = "Elapsed Time (hours)", y = "Cell Count")

logistic_fit_df <- logistic_fit_df %>%
  mutate(
    model1 = 5171 * 2^(0.023*Time),
    model2 = 5171 * 2^(0.022*Time),
    model3 = 5171 * 2^(0.019*Time),
    model4 = 5171 * 2^(0.015*Time)
  ) 

ggplot(logistic_fit_df) + 
  geom_line(data = logistic_fit_df, aes(Time, y_pred), linewidth = 1.5, color = "green") +
  geom_line(data = logistic_fit_df, aes(Time, model1), linetype = "dashed", linewidth = 1, color = "blue") +
  geom_line(data = logistic_fit_df, aes(Time, model2), linetype = "dashed", linewidth = 1, color = "red") +
  geom_line(data = logistic_fit_df, aes(Time, model3), linetype = "dashed", linewidth = 1, color = "sienna1") +
  geom_line(data = logistic_fit_df, aes(Time, model4), linetype = "dashed", linewidth = 1, color = "hotpink4") +
  geom_point(data = Logistic_cell_growth, aes(Time, Cell_count), size = 3) +
  labs(x = "Elapsed Time (hours)", y = "Cell Count") +
  scale_y_continuous(limits = c(0, 260000))

ggplot(logistic_fit_df) + 
  geom_line(data = logistic_fit_df, aes(Time, model1), linetype = "dashed", linewidth = 1, color = "blue") +
  geom_line(data = logistic_fit_df, aes(Time, model2), linetype = "dashed", linewidth = 1, color = "red") +
  geom_line(data = logistic_fit_df, aes(Time, model3), linetype = "dashed", linewidth = 1, color = "sienna1") +
  geom_line(data = logistic_fit_df, aes(Time, model4), linetype = "dashed", linewidth = 1, color = "hotpink4") +
  geom_point(data = Logistic_cell_growth, aes(Time, Cell_count), size = 3) +
  labs(x = "Elapsed Time (hours)", y = "Cell Count") +
  scale_y_continuous(limits = c(0, 260000))

Troubleshooting_exp_growth <- read_csv("Troubleshooting_exp_growth.csv")

Cell_lin_model2 = brm(
  data = Troubleshooting_exp_growth, family = gaussian,
  Log2_ending_cells ~ 1 + Log2_starting_cells, 
  prior = c(
    prior(normal(20,2), class = Intercept),
    prior(normal(0,10), class = b),
    prior(normal(0,1), class = sigma)
  ),
  iter = 2000, warmup = 1000, chains = 4, cores = 4
)

plot(Cell_lin_model)
print(Cell_lin_model)
Cell_lin_model$fit


Troubleshooting_exp_growth %>%
  data_grid(Log2_starting_cells = seq_range(Log2_starting_cells, n = 101)) %>%
  add_predicted_draws(Cell_lin_model) %>%
  ggplot(aes(x = Log2_starting_cells, y = Log2_ending_cells, color = Status)) +
  stat_lineribbon(aes(y = .prediction), .width = c(.99, .95, .8, .5), color = brewer.pal(5, "Blues")[[5]]) +
  geom_point(data = Troubleshooting_exp_growth) +
  scale_fill_brewer() +
  scale_color_manual(values = c("green", "red"))



Cell_lin_model2 = brm(
  data = Troubleshooting_exp_growth, family = gaussian,
  Log2_ending_cells ~ Log2_starting_cells, 
  prior = c(
    prior(normal(20,2), class = Intercept),
    prior(normal(0,10), class = b),
    prior(normal(0,1), class = sigma)
  ),
  iter = 2000, warmup = 1000, chains = 4, cores = 4
)

plot(Cell_lin_model2)
print(Cell_lin_model2)
Cell_lin_model2$fit