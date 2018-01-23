


data_fun <- function(mean.vec = apply(X, 2, mean),
                     prec.mat = diag(1 / apply(X, 2, var))) {
  res <- as.data.frame(rmvnorm(500, mean = mean.vec, sigma = solve(prec.mat)))
  colnames(res) <- c("x1", "x2")
  res
}

ggplot(faithful, aes(eruptions, waiting)) +
  geom_point() +
  theme_bw() -> p
# p + geom_contour(data = plot.df, aes(x1, x2, z = z), bins = 2)
for (k in 1:K) {
  p <- p + stat_ellipse(data = data_fun(m0[, k], Lambda[[k]]), aes(x1, x2),
                        geom = "path")
}
p

# Data
meanx <- apply(X, 2, mean)
varx <- apply(X, 2, var)
sdx <- apply(X, 2, sd)

K <- 6  # number of clusters
N <- nrow(X)
D <- ncol(X)

# Hyperparameters
alpha <- rep(0.5, K)  # Dirichlet
m <- rbind(
  # runif(K, min = min(X[, 1]), max = max(X[, 1])),
  # runif(K, min = min(X[, 2]), max = max(X[, 2]))
  sample(X[, 1], size = K, replace = FALSE),
  sample(X[, 2], size = K, replace = FALSE)
)
# m0 <- t(X[sample(seq_len(N), size = K, replace = FALSE), ])
beta <- rep(1, K)
nu <- rep(ncol(X), K)
W <- list()
for (k in seq_len(K)) {
  v1 <- runif(1, min = 0.005, max = 0.2)
  v2 <- runif(1, min = varx[2] / 10, max = varx[2] / 10)
  W[[k]] <- diag(1 / c(v1, v2)) / nu[k]
}

# Initialise
pi <- rep(1, K) / K
Lambda <- mapply("*", W, nu, SIMPLIFY = FALSE)
log.rnk <- matrix()
for (k in 1:K) {
  x.minus.m <- sweep(X, MARGIN = 2, STATS = m[, k], FUN = "-")
  res <- log(pi[k]) + 0.5 * log(Lambda[[k]]) - D / (2 * beta[k]) -
    nu[k] * diag(x.minus.m %*% W[[k]] %*% t(x.minus.m)) / 2
}

