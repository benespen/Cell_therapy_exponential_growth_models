Exponential Growth Models for Cell Therapy
================
Benjamin I. Espen
2026-07-28

# Cell expansion processes in cell and gene therapy

This is an exploration of how I gain insight from data. The data I will
be using is derived from real data, but has been anonymized and
manipulated to isolate it from its source. I like using real(ish) data
for this because real data has more surprises in it than simulated data
does. Part of using data is deciding how to deal with all the
complications and messiness that come with data collection. Real data is
not perfect. It also does not need to be perfect before we can make use
of it.

For this I will be using cell count data from cellular expansion
processes. Cell expansion processes are common to most cell and gene
therapy products, as you almost always start with fewer cells than you
would like to use. The biological heart of a cellular expansion process
is cellular division: one cell divides into two cells.

## Exponential growth model

This then suggests a model we might use: exponential growth. This model
multiplies a doublings per unit time parameter represented by $T_D$ by
an elapsed time to calculate the number of doublings. A standard
exponential growth model uses a base of *e*, but for cell growth the
natural base is *2*. The result from raising *2* to the number of
doublings is then multiplied by your starting cell count $q_1$ to
produce a final cell count $q_2$.

$q_2 = q_1 * 2^{T_D * (t_2 - t_1)}$

## Fitting our exponential growth model

If we have a couple of cell counts with associated times, we can
calculate $T_D$. Let’s look at an example that has been taken from a
real data set, suitably anonymized.

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
# Plot data points with ggplot
df_exponential_growth |> 
  ggplot(aes(x = elapsed_time, y = cell_number)) +
  geom_point(size = 3) +
    labs(x = "Elapsed Time (hours)", y = "Cell Count")
```

![](Exponential_growth_models_files/figure-gfm/exponential-growth-plot-1.png)<!-- -->

We can interpolate intermediate cell counts between the two measurements
using our exponential growth model. First we need to calculate $T_D$
from our two data points:

$q_2 = q_1 * 2^{T_D * (t_2 - t_1)}$

$q_2 / q_1 = 2^{T_D * (t_2 - t_1)}$

$log_2(q_2/q_1) = T_D *(t_2 - t_1))$

$\frac{log_2(q_2/q_1)}{t_2 - t_1} = T_D$

$\frac{log_2(98000/15000)}{192-12} = 0.015$

Since our time was in hours, $T_D$ is then in units of doublings per
hour. Using the formula, I can then calculate a predicted cell count. I
have done so for several intermediate times and placed that information
in a file. Let’s open up that file and plot the results.

``` r
Interpolated_exponential <- read_csv("Interpolated_exponential.csv", show_col_types = FALSE)
Interpolated_exponential
```

    ## # A tibble: 17 × 3
    ##     Time Cell_count Type 
    ##    <dbl>      <dbl> <chr>
    ##  1     0     13241. Model
    ##  2    12     15000  Data 
    ##  3    24     16993. Model
    ##  4    36     19251. Model
    ##  5    48     21810. Model
    ##  6    60     24708. Model
    ##  7    72     27991. Model
    ##  8    84     31711. Model
    ##  9    96     35924. Model
    ## 10   108     40698. Model
    ## 11   120     46106. Model
    ## 12   132     52233. Model
    ## 13   144     59174. Model
    ## 14   156     67037. Model
    ## 15   168     75945. Model
    ## 16   180     86037. Model
    ## 17   192     98099. Data

This file also has a categorical column called ‘Type’ that I am using to
show whether a given value is data from a measurement or a calculation
from a model. With that category, we can color code the resulting plot.

``` r
# Plot data points with ggplot
Interpolated_exponential |> 
  ggplot(aes(x = Time, y = Cell_count, color = Type)) +
  geom_point(size = 3) +
  labs(x = "Elapsed Time (hours)", y = "Cell Count")
```

![](Exponential_growth_models_files/figure-gfm/interp-exp-plot-1.png)<!-- -->

You might observe that I also extrapolated one point, the cell count at
*t* = 0. The reason for that will become clear shortly.

One of the nice properties of an exponential growth model is that it has
a constant doublings per unit time. If our cellular expansion process
truly is exhibiting exponential growth, we should get the same $T_D$
from any two points. I have cell count data at all the same time points
as I used for interpolation, so we can see how well the model performed.

``` r
# Load the full cell growth data set
Cell_growth <- read_csv("Cell_growth.csv", show_col_types = FALSE)
Cell_growth
```

    ## # A tibble: 17 × 3
    ##    Cell_count  Time Measurement_error
    ##         <dbl> <dbl>             <dbl>
    ##  1      5078.     0              141.
    ##  2     15363.    12              912.
    ##  3     15668.    24              905.
    ##  4     15702.    36              888.
    ##  5     15759.    48              863.
    ##  6     16593.    60              904.
    ##  7     18129.    72              977.
    ##  8     20157.    84             1050.
    ##  9     23005.    96             1166.
    ## 10     26676.   108             1339.
    ## 11     31939.   120             1576.
    ## 12     38152.   132             1866.
    ## 13     46108.   144             2180.
    ## 14     55604.   156             2641.
    ## 15     67513.   168             3134.
    ## 16     81529.   180             3727.
    ## 17     98099.   192             4475.

``` r
# Overlay real data on the model
ggplot(Interpolated_exponential, aes(Time, Cell_count)) + 
  geom_line(linewidth = 1.5) + 
  geom_point(data = Cell_growth, colour = "blue", size = 3) +
  labs(x = "Elapsed Time (hours)", y = "Cell Count")
```

![](Exponential_growth_models_files/figure-gfm/exp-model-comparison-1.png)<!-- -->

Now you can see one of the complications of real data that we need to
address. The *t* = 0 data point is unexpectedly low. In the real data
set this was derived from, this pattern was seen repeatedly. I
interpreted this an error in the cell count method because I didn’t
think it was biologically plausible that the cell count would truly jump
up like that and then flatten out. In the absence of a good measurement,
I ignored the *t* = 0 data point, starting my modeling at *t* = 12, and
I then extrapolated a more plausible *t* = 0 from the model.

However, it is also clear that the exponential growth model is not quite
right for this data. From 12 hours to 48 hours, the cell count changes
very little. After that point, something that looks very similar to
exponential growth occurs. If we think biologically again for a moment,
it isn’t surprising that the cells are not dividing as normal. In
circumstances like this, cells have often just endured some insult, such
as enzymatic detachment or centrifugation. It is not uncommon to see a
delay while the cells recover.

This suggests an improvement. We could shift the beginning of our
exponential model to 48 hours, and substitute the cell count at 48 hours
for all preceding time points. Our doublings per hour value needs to
change in order to make this happen.

$\frac{log_2(98000/15759)}{192-48} = 0.018$

``` r
Shifted_exp_model <- read_csv("Shifted_exp_model.csv", show_col_types = FALSE)
Shifted_exp_model
```

    ## # A tibble: 17 × 3
    ##     Time Cell_count Type 
    ##    <dbl>      <dbl> <chr>
    ##  1     0     15759. Model
    ##  2    12     15759. Model
    ##  3    24     15759. Model
    ##  4    36     15759. Model
    ##  5    48     15759. Data 
    ##  6    60     18395. Model
    ##  7    72     21420. Model
    ##  8    84     24941. Model
    ##  9    96     29042. Model
    ## 10   108     33817. Model
    ## 11   120     39377. Model
    ## 12   132     45851. Model
    ## 13   144     53389. Model
    ## 14   156     62167. Model
    ## 15   168     72388. Model
    ## 16   180     84289. Model
    ## 17   192     98099. Data

``` r
ggplot(Shifted_exp_model, aes(Time, Cell_count)) + 
  geom_line(color = "green", linewidth = 1.5) + 
  geom_point(data = Cell_growth, colour = "blue", size = 3) +
  labs(x = "Elapsed Time (hours)", y = "Cell Count")
```

![](Exponential_growth_models_files/figure-gfm/shifted-exp-model-1.png)<!-- -->

That is much closer, although you can still see that the cell counts are
systematically lower than the model predictions up until the last data
point. The reason for this is that cell expansion processes only
approximately follow exponential growth over some time domains. The
problem is that the exponential growth model knows no limits. However,
at some point, the cells will run out of room, run out of food, run out
of oxygen, or make too much waste, and growth will slow down or stop.

This kind of growth process that slows and then stops can be modeled as
logistic growth. There are a lot of ways to parametrize logistic
functions, here I’m going to follow along with what I used for
exponential growth. *L* is the lower asymptote, or our starting cell
count, *U* is the upper asymptote, $T_D$ is again our rate parameter,
and *M* is the inflection point of the logistic curve in hours.

$q_2 = U + \frac{L-U}{1+2^{T_D*(t-M)}}$

However, it is not very hard to find Plot of exponential growth Formula
of exponential growth

plot of logistic growth formula of logistic growth

close enough, especially if we only have a couple of data points to fit
