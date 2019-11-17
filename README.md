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
params.S <- c(beta = 0.9, gamma = 0.5, I.0 = 1, N = 100) # chosen SIR parameters 
sims <- lapply(seq(sim_no), function(i) {
  # let's hide the printing from generate.stochastic.sir
  sink(tempfile()) 
  on.exit(sink())
  generate.stochastic.sir(params = params.S, distribution.type = 'EXP')
})

head(sims[[1]])
```

    ##        time  S I R C
    ## 1 0.0000000 99 1 0 1
    ## 2 0.3965086 98 2 0 2
    ## 3 1.4075530 97 3 0 3
    ## 4 1.4297697 96 4 0 4
    ## 5 1.4824847 95 5 0 5
    ## 6 1.5091091 95 4 1 5

A simulated process is represented by a `data.frame` where each row is
an event, i.e. an infection or a recovery event. The columns are event
times, suscetible counts \(S_t\), infected count \(I_t\), recovered
counts \(R_t\) and cumulative infected counts \(C_t\), respectively. For
example, the row \(6\) represents an recovery event at time
\(1.5091091\) and the process stopped as there is no infected individual
in the population anymore.

Before fitting an EXPN on these simulated processes, we first need to
convert SIR processes to the Hawkes-friendly format, i.e. dropping all
recovery events

``` r
sims.H <- lapply(sims, SIR2HAWKES.stochastic.event.series)

head(sims.H[[1]])
```

    ##   magnitude      time
    ## 1         1 0.0000000
    ## 2         1 0.3965086
    ## 3         1 1.4075530
    ## 4         1 1.4297697
    ## 5         1 1.4824847
    ## 6         1 1.6931406

The converted `data.frame` now only keeps the infection event
timestamps.

Let’s now fit the EXPN model

``` r
fit_series(sims.H, model_type = 'EXPN', observation_time = Inf)
```

    ## Model: EXPN 
    ## No. of cascades: 50 
    ## init_par
    ##   K 9.54e-01; theta 1.02e+00; N 4.76e+04
    ## par
    ##   K 9.54e-01; theta 1.02e+00; N 4.76e+04
    ## Neg Log Likelihood: -210.544 
    ## lower_bound
    ##   K 1.00e-100; theta 1.00e-100; N 1.00e+00
    ## upper_bound
    ##   K 1.00e+04; theta 3.00e+02; N 1.00e+07
    ## convergence: 0

## Reference

If you find the code helpful, please consider citing the paper:

> \[1\] Kong, Quyu, Marian-Andrei Rizoiu, and Lexing Xie. “Modeling
> Information Cascades with Self-exciting Processes via Generalized
> Epidemic Models.” arXiv preprint arXiv:1910.05451 (2019).
