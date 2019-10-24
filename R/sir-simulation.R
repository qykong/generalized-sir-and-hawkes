# This script contains code for simulating stochastic SIR with generalized recovery distributions.
# currently two distributions are implemented: Exponential (classical stochastic SIR) and Power-law.
# two simulation algorithms are implemented for sampling a new infection event: rejection-sampling
# and sampling by converting intensity function to single unit intensity poisson process as shown in the paper.
# NOTE: For exponential recovery distribution, there exists efficient sampling algoritm. The
# algorithm implemented here is just for verification purpose.

sample.recovery.intervals <- function(params, n, distribution.type) {
  
  if (distribution.type == 'EXP') {
    u1 <- runif(n = n, min = 0, max = 1)
    return(-log(u1)/params[['gamma']])
  }
  if (distribution.type == 'PL') {
    u <- runif(n = n, min = 0, max = 1)
    res <-  params[['c']] / ((1 - u) ^ (1/(1+params[['theta']]))) - params[['c']]
    
    return(res)
  }
}

recovery_distribution <- function(params, t, distribution.type) {
  if (distribution.type == 'EXP') {
    return(params[['gamma']] * exp(-params[['gamma']] * t))
  }
  if (distribution.type == 'PL') {
    res <- params[['c']]^(1+params[['theta']]) * (1+params[['theta']]) * (t + params[['c']])^(-2-params[['theta']])
    return(res)
  }
}

# sample infective event interval given recovery events
sample.infected.interval <- function(params, recovery.times, S) {
  infectives <- length(recovery.times)
  lambda_I_part <- params[['beta']] * S / params[['N']]
  
  t <- 0
  lambda.max <- lambda_I_part * infectives
  j <- 1
  repeat {
    u <- runif(n = 2, min = 0, max = 1)
    sampled.t <- -log(u[1])/lambda.max
    t <- sampled.t + t
    while(t >= recovery.times[j]) {
      j <- j + 1
      if (j > length(recovery.times)) return(Inf)
    }
    current.lambda <- lambda_I_part * (length(recovery.times) - j + 1)
    
    if (u[2] <= current.lambda/lambda.max) {
      return(t)
    } else {
      lambda.max <- current.lambda
    }
  }
}
sample.infected.interval.inverted <- function(params, recovery.times, S) {
  
  infectives <- length(recovery.times)
  
  lambda_I_part <- params[['beta']] * S / params[['N']]
  
  upper <- recovery.times[1] * infectives
  s <- - log(runif(1)) / lambda_I_part
  
  while (s > upper) {
    s <- s - upper
    infectives <- infectives - 1
    if (infectives == 0) return(Inf)
    upper <- infectives * (recovery.times[length(recovery.times) - infectives + 1] - recovery.times[length(recovery.times) - infectives])
  }
  
  if (length(recovery.times) - infectives == 0) return(s / infectives)
  t <- recovery.times[length(recovery.times) - infectives] + s / infectives
  t
}

#' @param params a named vector of parameters for simulation. I.0 defines the number of initial infections,
#' N defines the population size. Beta and gamma defines the exponential recovery distribution, whereas 
#' beta, c and theta determines the power law recovery distribution.
#' @param Tmax this program will simulate the SIR process until Tmax. Defaults to infinity
#' @param distribution.type a string determines which recovery distribution for simulation, PL - powerlaw, EXP - Exponential
#' @param inverted if false then rejection sampling will be used. The sampling strategy shown in the paper will be applied otherwise
#' 
#' @return a data.frame where each row is an infection or recovery event
generate.stochastic.sir.power.law <- function(params, Tmax = Inf, distribution.type = 'PL', inverted = T) {
  ## start at time 0
  state <- data.frame(t = 0, I = params["I.0"], S = params["N"] - params["I.0"], C = params["I.0"])
  rownames(state) <- NULL
  j <- 1
  infected.time <- 0
  recovery.times <- sort(sample.recovery.intervals(params, n = params[['I.0']], distribution.type = distribution.type))
  
  repeat {
    infected.interval <- ifelse(inverted, sample.infected.interval.inverted(params, recovery.times - infected.time, state$S[j]),
                                sample.infected.interval(params, recovery.times - infected.time, state$S[j]))
    infected.time <- infected.interval + infected.time
    
    # append all recovery events
    for (time in recovery.times[recovery.times < infected.time]) {
      nextState <- data.frame(t = time, I = state$I[j]-1, S = state$S[j], C = state$C[j])
      state <- rbind(state, nextState)
      j <- j + 1
      cat(sprintf("\rCurrent simulation time: %.3f / %.3f (S=%d, I=%d, R=%d, C=%d).", state$t[j], Tmax, state$S[j], state$I[j], params["N"]-state$S[j]-state$I[j], state$C[j]))
    }
    recovery.times <- recovery.times[recovery.times >= infected.time]
    
    # append infective event and update infective times for the new infective event
    if (!is.infinite(infected.time)) {
      nextState <- data.frame(t = infected.time, I = state$I[j]+1, S = state$S[j]-1, C = state$C[j]+1)
      state <- rbind(state, nextState)
      j <- j + 1
      cat(sprintf("\rCurrent simulation time: %.3f / %.3f (S=%d, I=%d, R=%d, C=%d).", state$t[j], Tmax, state$S[j], state$I[j], params["N"]-state$S[j]-state$I[j], state$C[j]))
      recovery.times <- sort(c(recovery.times, sample.recovery.intervals(params, n = 1, distribution.type = distribution.type) + infected.time))
    } else {
      break
    }
  }
  
  rownames(state) <- NULL
  state$R <- state$C - state$I
  names(state) <- c("time", "I", "S", "C", "R")
  state <- state[c("time", "S", "I", "R", "C")]
  
  cat(sprintf("\n--> Simulation done!\n"))
  return(state)
}