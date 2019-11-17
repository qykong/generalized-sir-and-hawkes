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
    [HawkesFit](https://github.com/qykong/HawkesFit).
  - HawkesN Goodness-of-fit Tests: code appeared in this repo at
    `R/goodness-of-fit.R`.

## Tutorial 1: fitting HawkesN on simulated generalized SIR processes

In this tutorial, we are going to first simulate some stochastic SIR
processes with an exponential recovery distribution and then fit them
using EXPN (i.e. HawkesN processes with an exponential kernel).

Let’s first load all required packages and functions.

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

TODO

## Reference

If you find the code helpful, please consider citing the paper:

> \[1\] Kong, Quyu, Marian-Andrei Rizoiu, and Lexing Xie. “Modeling
> Information Cascades with Self-exciting Processes via Generalized
> Epidemic Models.” arXiv preprint arXiv:1910.05451 (2019).
