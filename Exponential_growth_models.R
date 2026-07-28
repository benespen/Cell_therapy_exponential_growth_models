df_exponential_growth <- data.frame(
  cell_number = c(15000 , 98000),
  elapsed_time = c(12 , 192)
)
df_exponential_growth
df_exponential_growth |> 
  ggplot(aes(x = elapsed_time, y = cell_number)) +
  geom_point(size = 3) +
    labs(x = "Elapsed Time (hours)", y = "Cell Count")