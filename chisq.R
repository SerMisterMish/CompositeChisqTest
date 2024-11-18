library(ggplot2)
library(tidyr)
library(dplyr)

group.states <- function(freq.pract, prob.theor, M) {
  state.names <- names(freq.pract)
  K <- length(freq.pract)
  N <- sum(freq.pract)
  K.new <- K
  
  for (i in 1:(K - 1)) {
    if (N * prob.theor[i] < M) {
      prob.theor[i + 1] <- prob.theor[i] + prob.theor[i + 1]
      freq.pract[i + 1] <- freq.pract[i] + freq.pract[i + 1]
      prob.theor[i] <- -1
      freq.pract[i] <- -1
      
      if (grepl(">", state.names[i + 1]))
        state.names[i + 1] <- paste0(">", strsplit(state.names[i], ",")[[1]][1])
      else
        state.names[i + 1] <- paste0(state.names[i], ", ", state.names[i + 1])
      state.names[i] <- NA
      
      K.new <- K.new - 1
    }
  }
  
  prob.theor <- prob.theor[prob.theor >= 0]
  freq.pract <- freq.pract[freq.pract >= 0]
  state.names <- state.names[!is.na(state.names)]
  
  if (N * prob.theor[K.new] < M) {
    prob.theor[K.new - 1] <- prob.theor[K.new] + prob.theor[K.new - 1]
    freq.pract[K.new - 1] <- freq.pract[K.new] + freq.pract[K.new - 1]
    prob.theor[K.new] <- -1
    freq.pract[K.new] <- -1
    
    if (grepl(">", state.names[K.new]))
      state.names[K.new - 1] <- paste0(">", strsplit(state.names[K.new - 1], ",")[[1]][1])
    else
      state.names[K.new - 1] <- paste0(state.names[K.new - 1], ", ", state.names[K.new])
    state.names[K.new] <- NA
    
    K.new <- K.new - 1
  }
  
  names(prob.theor) <- state.names[!is.na(state.names)]
  names(freq.pract) <- state.names[!is.na(state.names)]
  
  list(freq.pract = freq.pract[freq.pract >= 0],
       prob.theor = prob.theor[prob.theor >= 0],
       K = K.new)
}

chisq.pval <- function(x,
                       states,
                       mle.parameters,
                       model.parameters,
                       dist.fun,
                       M = 5,
                       group.fun = "group.states") {
  freq.pract <- numeric(length(states))
  names(freq.pract) <- states
  tmp <- table(x)
  freq.pract[names(tmp)] <- tmp
  
  if (tail(states, 1) == ">") {
    prob.theor <- numeric(length(states))
    numeric.states <- as.numeric(states[-length(states)])
    args <- append(list(x = numeric.states), mle.parameters) |> append(model.parameters)
    prob.theor[-length(states)] <- do.call(dist.fun, args)
    prob.theor[length(states)] <- 1 - sum(prob.theor)
  }
  else
  {
    args <- append(list(x = states), mle.parameters) |> append(model.parameters)
    prob.theor <- do.call(dist.fun, args)
  }
  names(prob.theor) <- states
  
  group.states <- get(group.fun)
  grouped <- group.states(freq.pract, prob.theor, M)
  freq.pract <- grouped$freq.pract
  prob.theor <- grouped$prob.theor
  K <- grouped$K
  
  t <- chisq.test(freq.pract, p = prob.theor)$statistic
  names(t) <- NULL
  df <- K - 1 - length(mle.parameters)
  
  if (df <= 0){
    p <- NA
    warning("Non-positive degrees of freedom")
  }
  else
    p <-  1 - pchisq(t, df = df)
  
  list(
    t = t,
    p = p,
    grouped.sample = grouped,
    df = df
  )
}

plot.pval <- function(pvals,
                      Ms,
                      main,
                      sub,
                      xlab = expression(alpha),
                      ylab = expression(F[n](alpha)),
                      ...) {
  pvals |> lapply(t) |> as.data.frame() |>
    `colnames<-`(paste0(c("p", "states"), rep(Ms, each = 2))) |>
    select(!starts_with("states")) |>
    pivot_longer(
      cols = starts_with("p"),
      names_to = "M",
      names_prefix = "p",
      values_to = "p"
    ) |>
    ggplot() + stat_ecdf(aes(x = p, colour = M)) +
    geom_abline(slope = 1,
                intercept = 0,
                col = "grey") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5)) + labs(title = main, subtitle = sub) +
    xlab(xlab) + ylab(ylab) + xlim(0, 1) + ylim(0, 1)
}

model.pval <- function(rgen,
                       gen.parameters,
                       h0.dist,
                       mle.fun,
                       model.parameters,
                       right.border,
                       M) {
  x <- do.call(rgen, gen.parameters)
  if (!is.null(mle.fun))
    mle.parameters <- mle.fun(x, model.parameters)
  else 
    mle.parameters <- NULL
  
  states <- 0:max(x)
  if (is.infinite(right.border))
    states <- c(states, ">")
  else if (tail(states, 1) < right.border)
    states <- 0:right.border
  
  chisq.pval(x, states, mle.parameters, model.parameters, h0.dist, M)
}
