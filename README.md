Generalized SIR and HawkesN processes
================

## Introduction

This repo hosts the code or links for reproducing the experiments
presented in \[1\]. Specifically, we summarize the necessary components
here:

  - Datasets: we used three large Twitter retweet cacade datasets in our
    experiments. These datasets can be found at the links:
    [ActiveYT](https://github.com/computationalmedia/hip-popularity/tree/master/data),
    [NEWS](https://github.com/computationalmedia/featuredriven-hawkes/tree/master/data)
    and [SEISMIC](http://snap.stanford.edu/seismic/).
  - Generalized SIR Simulation: code appeared in this repo at
    `R/sir-simulation.R`.
  - HawkesN Simulation and Fitting: code available in the package
    [evenlty](https://github.com/behavioral-ds/evently).
  - HawkesN Goodness-of-fit Tests: code appeared in this repo at
    `R/goodness-of-fit.R`.

## Tutorial 1: fitting HawkesN on simulated generalized SIR processes

In this tutorial, we are going to first simulate some stochastic SIR
processes with an exponential recovery distribution and then fit them
using EXPN (i.e. HawkesN processes with an exponential kernel).

Let’s first load all required packages and functions.

``` r
source('R/package.R')
source('R/goodness-of-fit.R')
source('R/sir-simulation.R')
set.seed(2)
```

We then simulate some stochastic SIR processes

``` r
sim_no <- 50 # simulate 100 processes
params.S <- c(beta = 5, gamma = 0.5, I.0 = 1, N = 100) # chosen SIR parameters 
sims <- lapply(seq(sim_no), function(i) {
  # let's hide the printing from generate.stochastic.sir
  sink(tempfile()) 
  on.exit(sink())
  generate.stochastic.sir(params = params.S, distribution.type = 'EXP')
})

tail(sims[[1]])
```

    ##          time S I   R   C
    ## 195  8.347471 0 5  95 100
    ## 196  8.519146 0 4  96 100
    ## 197  8.760964 0 3  97 100
    ## 198  9.783959 0 2  98 100
    ## 199  9.989227 0 1  99 100
    ## 200 10.696395 0 0 100 100

A simulated process is represented by a `data.frame` where each row is
an event, i.e. an infection or a recovery event. The columns are event
times, suscetible counts \(S_t\), infected count \(I_t\), recovered
counts \(R_t\) and cumulative infected counts \(C_t\), respectively. For
example, the row \(200\) represents an recovery event at time
\(10.696395\) and the process stopped as there is no infected individual
in the population anymore.

Before fitting an EXPN on these simulated processes, we first need to
convert SIR processes to the Hawkes-friendly format, i.e. dropping all
recovery events

``` r
sims.H <- lapply(sims, SIR2HAWKES.stochastic.event.series)

head(sims.H[[1]])
```

    ##   magnitude       time
    ## 1         1 0.00000000
    ## 2         1 0.07137156
    ## 3         1 0.25335955
    ## 4         1 0.25735855
    ## 5         1 0.26684724
    ## 6         1 0.29202057

The converted `data.frame` now only keeps the infection event
timestamps.

Let’s now fit the EXPN
model

``` r
model <- fit_series(sims.H, model_type = 'EXPN', cores = 10, observation_time = Inf)

cat(str_interp('According the param equivalence: 
beta = ${model$par[["K"]] * model$par[["theta"]]},
gamma = ${model$par[["theta"]]}'))
```

    ## According the param equivalence: 
    ## beta = 5.19875453811,
    ## gamma = 0.711039

You’ll notice the modeling bias as discussed in \[1\].

## Tutorial 2: goodness-of-fit tests on real retweet cascades

In this tutorial, we run our goodness-of-fit experiments on a sub-sample
of the NEWS dataset. Let’s first load this
dataset

``` r
index <- read_csv(file = 'https://raw.githubusercontent.com/computationalmedia/featuredriven-hawkes/master/data/index.csv', col_types = 'dd')
data_file <- read_csv(file = 'https://raw.githubusercontent.com/computationalmedia/featuredriven-hawkes/master/data/data.csv', col_types = 'dd')
head(index)
```

    ## # A tibble: 6 x 2
    ##   start_ind end_ind
    ##       <dbl>   <dbl>
    ## 1         1      63
    ## 2        64     126
    ## 3       127     242
    ## 4       243     295
    ## 5       296     365
    ## 6       366     452

``` r
head(data_file)
```

    ## # A tibble: 6 x 2
    ##    time magnitude
    ##   <dbl>     <dbl>
    ## 1     0      2086
    ## 2   136       357
    ## 3   153       340
    ## 4   168       724
    ## 5   191      5553
    ## 6   221       345

Each row in `index.csv` is a single retweet event cascade whose events
can be found in `data.csv` indexed by the `start_ind` and `end_ind`
fields. Each row in `data.csv` is a tweet/retweet where `time` is the
time relative to the initial tweet and `magnitude` is the number of
followers the corresponding user had.

There are in total 20,093 cascades in this dataset, but, in this
tutorial, we will only use the first 5 cascades as a demonstration. We
then fit the `EXPN` and `PLN` models on these complete cascades.

``` r
selected_index <- index[seq(5), ]

# get the corresponding casade events
cascades <- pmap(selected_index, function(start_ind, end_ind) data_file[start_ind:end_ind, ])
model_types <- c('mEXPN', 'mPL')

fitted_results <- map(model_types, function(model) {
  map(cascades, ~fit_series(list(.x), model_type = model, cores = 10, observation_time = Inf))
}) %>%
  set_names(model_types)
```

Now with the fitted results, let’s test how good the fits are on the
cascades

``` r
goodness_tests <- map_df(model_types, function(model) {
  map_dfr(seq_along(cascades), function(i) {
    # get the corresponding fitted model
    fitted_model <- fitted_results[[model]][[i]]
    
    # please refer to the Eq.15 in the paper
    times <- Lambda(fitted_model$data[[1]], fitted_model$par, cascades[[i]]$time[-1], kernel.type = model)

    ks_test_res <- ks.test(times, 'pexp')
    
    return(list(cascade_no = i, model_type = model, p_value = ks_test_res$p.value, distance = ks_test_res$stat[['D']]))
  })
})
head(goodness_tests)
```

    ## # A tibble: 6 x 4
    ##   cascade_no model_type   p_value distance
    ##        <int> <chr>          <dbl>    <dbl>
    ## 1          1 mEXPN      0.00505      0.216
    ## 2          2 mEXPN      0.0507       0.172
    ## 3          3 mEXPN      0.0000832    0.209
    ## 4          4 mEXPN      0.00963      0.227
    ## 5          5 mEXPN      0.000259     0.251
    ## 6          1 mPL        0.500        0.103

We can further aggregate the results in terms of the p\_values to check
the passing rate given the significance level is 0.05

``` r
goodness_tests %>%
  group_by(model_type) %>%
  summarise(passing_rate = sum(p_value > 0.05)/length(p_value))
```

    ## # A tibble: 2 x 2
    ##   model_type passing_rate
    ##   <chr>             <dbl>
    ## 1 mEXPN               0.2
    ## 2 mPL                 1

However, we note this is just a small sample results and please refer to
the paper for results on the whole dataset.

## Reference

If you find the code helpful, please consider citing the paper:

> \[1\] Kong, Quyu, Marian-Andrei Rizoiu, and Lexing Xie. “Modeling
> Information Cascades with Self-exciting Processes via Generalized
> Epidemic Models.” arXiv preprint arXiv:1910.05451 (2019).
