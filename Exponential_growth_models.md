Exponential Growth Models for Cell Therapy
================
Benjamin I. Espen
2026-07-28

# Exponential Models

## Cell expansion processes in cell and gene therapy

Cell expansion processes are common to most cell and gene therapy
products, as you almost always start with fewer cells than you would
like to use. The biological heart of a cellular expansion process is
cellular division: one cell divides into two cells.

This then suggests a model we might use: exponential growth. And for
these purposes, an exponential growth model will work well enough. This
model multiplies a doubling time represented by $T_D$ by an elapsed time
to calculate the number of doublings. The number of doublings is then
multipled by your starting cell count $q_1$ to produce a final cell
count $q_2$

$q_2 = q_1 * 2^{T_D * t_2 - t_1}$

``` r
# Make a data frame containing two cell counts at two elapsed times
df_exponential_growth <- data.frame(
  cell_number = c(15000 , 98000),
  elapsed_time = c(12 , 192)
)
# Show the data frame
df_exponential_growth
```

    ##   cell_number elapsed_time
    ## 1       15000           12
    ## 2       98000          192

``` r
df_exponential_growth |> 
  ggplot(aes(x = elapsed_time, y = cell_number)) +
  geom_point(size = 3) +
    labs(x = "Elapsed Time (hours)", y = "Cell Count")
```

![](Exponential_growth_models_files/figure-gfm/exponential-growth-plot-1.png)<!-- -->

However, it is not very hard to find Plot of exponential growth Formula
of exponential growth

plot of logistic growth formula of logistic growth

close enough, especially if we only have a couple of data points to fit
