source("chisq.R")

# Binomial distribution
modeling.times <- 500
sample.size <- 100

gen.parameters <- list(n = sample.size, size = 10, prob = 0.1)
model.parameters <- list(size = 10)
mle.p <- function(x, model.parameters) {
  mean(x) / model.parameters$size
}

Ms <- c(5, 0, 1, 10, 20)

set.seed(5)


pvals <- lapply(Ms, (\(M)
                     replicate(modeling.times, {
                       tmp <- model.pval(rbinom,
                                         gen.parameters,
                                         dbinom,
                                         mle.p,
                                         model.parameters,
                                         model.parameters$size,
                                         M)
                       c(tmp$p, tmp$grouped.sample$K)
                     })))

plot.pval(
  pvals,
  Ms = Ms,
  main = "Выборочная функция распределения p-v при верной H0",
  sub = paste0(
    "Bin(n = ",
    gen.parameters$size,
    ", p = ",
    gen.parameters$prob,
    "), sample size = ",
    sample.size,
    ", p-v size = ",
    modeling.times
  )
)
for (m in seq_along(Ms))
  print(paste("M = ", Ms[m], ", med = ", median(pvals[[m]][2, ])))


# Geometric distribution
modeling.times <- 500
sample.size <- 100

gen.parameters <- list(n = sample.size, prob = 0.3)
model.parameters <- NULL
mle.q <- function(x, model.parameters) {
  1 / (mean(x) + 1)
}

Ms <- c(5, 0, 1, 10, 20)

set.seed(5)


pvals <- lapply(Ms, (\(M)
                     replicate(modeling.times, {
                       tmp <- model.pval(rgeom, gen.parameters, dgeom, mle.q, model.parameters, Inf, M)
                       c(tmp$p, tmp$grouped.sample$K)
                     })))

plot.pval(
  pvals,
  Ms = Ms,
  main = "Выборочная функция распределения p-v при верной H0",
  sub = paste0(
    "Geom(q = ",
    gen.parameters$prob,
    "), sample size = ",
    sample.size,
    ", p-v size = ",
    modeling.times
  )
)
for (m in seq_along(Ms))
  print(paste("M = ", Ms[m], ", med = ", median(pvals[[m]][2, ])))


# H0: Geometric(0.3), exact test
modeling.times <- 500
sample.size <- 100

gen.parameters <- list(n = sample.size, prob = 0.3)
model.parameters <- list(prob = 0.3)


Ms <- c(0, 1, 5, 10, 20)

set.seed(5)


pvals <- lapply(Ms, (\(M)
                     replicate(modeling.times, {
                       tmp <- model.pval(rgeom, gen.parameters, dgeom, NULL, model.parameters, Inf, M)
                       c(tmp$p, tmp$grouped.sample$K)
                     })))

plot.pval(
  pvals,
  Ms = Ms,
  main = "Выборочная функция распределения p-v при H0: x~Geom(0.3)",
  sub = paste0(
    "Geom(q = ",
    gen.parameters$prob,
    "), sample size = ",
    sample.size,
    ", p-v size = ",
    modeling.times
  )
)
for (m in seq_along(Ms))
  print(paste("M = ", Ms[m], ", med = ", median(pvals[[m]][2, ])))

# H0: Geom; reality: Bin; sample.size = 1000
modeling.times <- 500
sample.size <- 1000

gen.parameters <- list(n = sample.size, size = 10, prob = 0.01)
model.parameters <- NULL
mle.q <- function(x, model.parameters) {
  1 / (mean(x) + 1)
}

Ms <- c(1, 5, 10)

set.seed(5)

pvals <- lapply(Ms, (\(M)
                     replicate(modeling.times, {
                       tmp <- model.pval(rbinom, gen.parameters, dgeom, mle.q, model.parameters, Inf, M)
                       c(tmp$p, tmp$grouped.sample$K)
                     })))

plot.pval(
  pvals,
  Ms = Ms,
  main = "Выборочная функция распределения p-v при H0: x ~ Geom(q)",
  sub = paste0(
    "Binom(n = ",
    gen.parameters$size,
    ", p = ",
    gen.parameters$prob,
    "), sample size = ",
    sample.size,
    ", p-v size = ",
    modeling.times
  )
)
for (m in seq_along(Ms))
  print(paste("M = ", Ms[m], ", med = ", median(pvals[[m]][2, ])))


# H0: Geom; reality: Bin; sample.size = 5000
modeling.times <- 500
sample.size <- 5000

gen.parameters <- list(n = sample.size, size = 10, prob = 0.01)
model.parameters <- NULL
mle.q <- function(x, model.parameters) {
  1 / (mean(x) + 1)
}

Ms <- c(1, 5, 10)

set.seed(5)

pvals <- lapply(Ms, (\(M)
                     replicate(modeling.times, {
                       tmp <- model.pval(rbinom, gen.parameters, dgeom, mle.q, model.parameters, Inf, M)
                       c(tmp$p, tmp$grouped.sample$K)
                     })))

plot.pval(
  pvals,
  Ms = Ms,
  main = "Выборочная функция распределения p-v при H0: x ~ Geom(q)",
  sub = paste0(
    "Binom(n = ",
    gen.parameters$size,
    ", p = ",
    gen.parameters$prob,
    "), sample size = ",
    sample.size,
    ", p-v size = ",
    modeling.times
  )
)
for (m in seq_along(Ms))
  print(paste("M = ", Ms[m], ", med = ", median(pvals[[m]][2, ])))
