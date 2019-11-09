# the following functions are used for goodness-of-fit tests.
# including the functions for computing Lambda of all model types 
Lambda <- function(history, params, t, marked = T, kernel.type = 'EXPN') {
  if (any(is.na(params))) stop('params are NA!')
  params[['beta']] <- ifelse(marked, params[['beta']], 0)
  switch(kernel.type,
         PL = Lambda.PL(history, params, t),
         EXP = Lambda.EXP(history, params, t),
         EXPN = Lambda.EXPN(history, params, t),
         PLN = Lambda.PLN(history, params, t),
         QEXP = Lambda.QEXP(history, params, t),
         QEXPN = Lambda.QEXPN(history, params, t),
         SI = Lambda.SI(history, params, t))
}

Lambda.EXPN <- function(history, params, t) {
  #K * sum {i in 1..L-1} (exp(beta * log(magnitude[i] + 1e-100)) * sum{j in i..L-1}((N - j) / N * (exp(-1*theta*(time[j] - time[i])) - exp(-1*theta * (time[j+1] - time[i])))))
  max_t <- max(t)
  mid <- params[['K']] / (params[['N']])
  
  map_dbl(t, function(time) {
    if (time == 0) return(0)
    timestamps <- history$time[history$time < time]
    sum(map_dbl(seq_along(timestamps), function(j) {
      history$magnitude[j]^(params[['beta']]) *
        (exp(-params[['theta']] * (timestamps[length(timestamps)] - timestamps[j]) ) -
           exp(-params[['theta']] * (time - timestamps[j]) ) )
    })) * (params[['N']] - (length(timestamps))) * mid
  })
}

Lambda.SI <- function(history, params, t) {
  #K * sum {i in 1..L-1} (exp(beta * log(magnitude[i] + 1e-100)) * sum{j in i..L-1}((N - j) / N * (exp(-1*theta*(time[j] - time[i])) - exp(-1*theta * (time[j+1] - time[i])))))
  max_t <- max(t)
  mid <- params[['K']] / (params[['N']])
  
  map_dbl(t, function(time) {
    if (time == 0) return(0)
    timestamps <- history$time[history$time < time]
    sum(map_dbl(seq_along(timestamps), function(j) {
      history$magnitude[j]^(params[['beta']]) *
        ( time - timestamps[length(timestamps)] )
    })) * (params[['N']] - (length(timestamps))) * mid
  })
}

Lambda.QEXPN <- function(history, params, t) {
  #K / ((2 - theta) * N) * sum {i in 1..L-1} (magnitude[i]^beta * sum{j in i..L-1}( (N - j) * ( (1 + (theta-1)*(time[j] - time[i]))^((2 - theta) / (1-theta)) - (1 + (theta-1)*(time[j+1] - time[i]))^((2-theta)/(1-theta) ) )))
  max_t <- max(t)
  mid <- params[['K']] / ((2 - params[['theta']]) * params[['N']])
  theta_tmp <- (2 - params[['theta']]) / (1 - params[['theta']])
  map_dbl(t, function(time) {
    if (time == 0) return(0)
    timestamps <- history$time[history$time < time]
    sum(map_dbl(seq_along(timestamps), function(j) {
      history$magnitude[j]^(params[['beta']]) *
        ((1 + (params[['theta']] - 1) * (timestamps[length(timestamps)] - timestamps[j]))^theta_tmp -
           (1 + (params[['theta']] - 1) * (time - timestamps[j]))^theta_tmp)
    })) * (params[['N']] - (length(timestamps))) * mid
  })
}

Lambda.EXP <- function(history, params, t) {
  # K * sum {i in 1..L} (magnitude[i]^beta * (1 - exp(-1 * theta * (time[L] - time[i]))))
  max_t <- max(t)
  magnitude <- history$magnitude[history$time < max_t]^(params[['beta']])
  
  map_dbl(t, function(time) {
    if (time == 0) return(0)
    timestamps <- history$time[history$time < time]
    sum(map_dbl(seq_along(timestamps), function(j) {
      magnitude[j] * (exp(-params[['theta']] * (timestamps[length(timestamps)] - timestamps[j])) - exp(-params[['theta']] * (time - timestamps[j])))
    })) * params[['K']]
  })
}

Lambda.QEXP <- function(history, params, t) {
  # K / ((2 - theta) ) * sum {i in 1..L-1} (magnitude[i]^beta * sum{j in i..L-1}(  ( (1 + (theta-1)*(time[j] - time[i]))^((2 - theta) / (1-theta)) - (1 + (theta-1)*(time[j+1] - time[i]))^((2-theta)/(1-theta) ) ))
  max_t <- max(t)
  magnitude <- history$magnitude[history$time < max_t]^(params[['beta']])
  theta_tmp <- (2 - params[['theta']]) / (1 - params[['theta']])
  mid <- params[['K']] / ((2 - params[['theta']]))
  
  map_dbl(t, function(time) {
    if (time == 0) return(0)
    timestamps <- history$time[history$time < time]
    sum(map_dbl(seq_along(timestamps), function(j) {
      magnitude[j] * ((1 + (params[['theta']] - 1) * (timestamps[length(timestamps)] - timestamps[j]))^theta_tmp -
                        (1 + (params[['theta']] - 1) * (time - timestamps[j]))^theta_tmp)
    })) * mid
  })
}

Lambda.PLN <- function(history, params, t) {
  #K * sum {i in 1..L-1} ( magnitude[i]^beta * sum{j in i..L-1}( (N - j) * ( (time[j] - time[i] + c)^(-1*theta) - (time[j+1] - time[i] + c)^(-1*theta) ))) / (theta * N)
  max_t <- max(t)
  magnitude <- history$magnitude[history$time < max_t]^(params[['beta']])
  mid <- params[['K']] / (params[['N']] * params[['theta']])
  
  map_dbl(t, function(time) {
    if (time == 0) return(0)
    timestamps <- history$time[history$time < time]
    sum(map_dbl(seq_along(timestamps), function(j) {
      magnitude[j] * 
        ((timestamps[length(timestamps)] - timestamps[j] + params[['c']])^(-params[['theta']]) -
           (time - timestamps[j] + params[['c']])^(-params[['theta']]) )
    })) *  (params[['N']] - (length(timestamps))) * mid
  })
}

Lambda.PL <- function(history, params, t) {
  # K * sum {i in 1..L} (magnitude[i]^beta * ( (1 / c)^theta - ( 1 / (time[L] + c - time[i]))^theta )) / theta
  max_t <- max(t)
  magnitude <- history$magnitude[history$time < max_t]^(params[['beta']])
  mid <- params[['K']] / params[['theta']]
  map_dbl(t, function(time) {
    if (time == 0) return(0)
    timestamps <- history$time[history$time < time]
    sum(map_dbl(seq_along(timestamps), function(j) {
      magnitude[j] * ((timestamps[length(timestamps)] - timestamps[j] + params[['c']])^(-params[['theta']]) - (time - timestamps[j] + params[['c']])^(-params[['theta']]))
    })) * mid
  })
}

ed.test <- function(samples) {
  (1 - pnorm(sqrt(length(samples)) * (var(samples) - 1) / sqrt(8)))
}

goodness_of_fit_tests <- function(history, params, marked = F, kernel.type) {
  times <- Lambda(history, params, history$time, marked, kernel.type = kernel.type)
  ks <- ks.test(unique(times), 'pexp')
  lb <- Box.test(unique(times), type = 'Ljung-Box')
  ed <- ed.test(times)
  return(list(ks = list(p.value = ks$p.value,
                        statistic = ks$statistic[['D']]),
              ed = list(p.value = ed),
              lb = list(p.value = lb$p.value)))
}
