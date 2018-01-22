source("01-prelim.R")

## ---- points ----
set.seed(123)
N <- 150
f <- function(x, truth = FALSE) {
  35 * dnorm(x, mean = 1, sd = 0.8) +
    65 * dnorm(x, mean = 4, sd = 1.5) +
    (x > 4.5) * (exp((1.25 * (x - 4.5))) - 1) +
    3 * dnorm(x, mean = 2.5, sd = 0.3)
}
x <- c(seq(0.2, 1.9, length = N * 5 / 8), seq(3.7, 4.6, length = N * 3 / 8))
x <- sample(x, size = N)
x <- x + rnorm(N, sd = 0.65)  # adding random fluctuation to the x
x <- sort(x)
y.err <- rt(N, df = 1)
y <- f(x) + sign(y.err) * pmin(abs(y.err), rnorm(N, mean = 4.1))  # adding random terms to the y

# True values
x.true <- seq(-2.1, 7, length = 1000)
y.true <- f(x.true, TRUE)

# Data for plot
dat <- data.frame(x, y)
dat.truth <- data.frame(x.true, y.true)

p1 <- ggplot() +
  geom_point(data = dat, aes(x = x, y = y)) +
  scale_x_continuous(
    limits = c(min(x.true), max(x.true)),
    breaks = NULL, name = expression(italic(x))
  ) +
  scale_y_continuous(
    limits = c(min(y) - 5, max(y) + 5),
    breaks = NULL, name = expression(italic(y))
  ) +
  theme_bw()

## ---- plot.function1 ----
fnH4 <- function(x, y = NULL, l = 1) {
x <- scale(x, scale = FALSE)
if (is.vector(x))
  x <- matrix(x, ncol = 1)
n <- nrow(x)
A <- matrix(0, n, n)
index.mat <- upper.tri(A)
index <- which(index.mat, arr.ind = TRUE)
xcrossprod <- tcrossprod(x)
if (is.null(y)) {
  tmp1 <- diag(xcrossprod)[index[, 1]]
  tmp2 <- diag(xcrossprod)[index[, 2]]
  tmp3 <- xcrossprod[index]
  A[index.mat] <- tmp1 + tmp2 - 2 * tmp3
  A <- A + t(A)
  tmp <- exp(-A / (2 * l ^ 2))
} else {
  if (is.vector(y))
    y <- matrix(y, ncol = 1)
  else y <- as.matrix(y)
  y <- sweep(y, 2, attr(x, "scaled:center"), "-")
  m <- nrow(y)
  B <- matrix(0, m, n)
  indexy <- expand.grid(1:m, 1:n)
  ynorm <- apply(y, 1, function(z) sum(z ^ 2))
  xycrossprod <- tcrossprod(y, x)
  tmp1 <- ynorm[indexy[, 1]]
  tmp2 <- diag(xcrossprod)[indexy[, 2]]
  tmp3 <- as.numeric(xycrossprod)
  B[, ] <- tmp1 + tmp2 - 2 * tmp3
  tmp <- exp(-B / (2 * l ^ 2))
}
tmp
}

dev.SEkern <- function(theta, y = y) {
  alpha <- mean(y)
  lambda <- theta[1]
  psi <- theta[2]
  n <- length(y)
  H <- fnH4(x, l = theta[3])
  H2 <- H %*% H
  Vy <- psi * lambda ^ 2 * H2 + diag(1 / psi, n)
  tmp <- -(n / 2) * log(2 * pi) - (1 / 2) * determinant(Vy)$mod -
    (1 / 2) * (y - alpha) %*% solve(Vy, y - alpha)
  as.numeric(-2 * tmp)
}

plot1 <- function(kernel, no.of.draws = 100) {
  # Fit an I-prior model -------------------------------------------------------
  if (kernel == "SE") {
    mod.se <- optim(c(1, 1, 1), dev.SEkern, method = "L-BFGS", y = y,
                    lower = c(-Inf, 1e-9, 1e-9))
    n <- length(y)
    H <- fnH4(x, l = mod.se$par[3])
    H2 <- H %*% H
    alpha <- mean(y)
    lambda <- mod.se$par[1]
    psi <- mod.se$par[2]
    Vy <- psi * lambda ^ 2 * H2 + diag(1 / psi, n)
    w.hat <- psi * lambda * H %*% solve(Vy, y - alpha)
    VarY.inv <- solve(Vy)
    H.star <- fnH4(x = x, y = x.true, l = mod.se$par[3])
    y.fitted <- as.numeric(mean(y) + lambda * H.star %*% w.hat)
    y.fitted2 <- as.numeric(mean(y) + lambda * H %*% w.hat)
  } else {
    if (kernel == "Canonical") {
      mod <- iprior(kernL(y, x, model = list(kernel = "Canonical")))
      H.star <- fnH2(x = x, y = x.true)
    }
    if (kernel == "FBM") {
      mod <- fbmOptim(kernL(y, x, model = list(kernel = "FBM")))
      H.star <- fnH3(x = x, y = x.true, gamma = mod$ipriorKernel$model$Hurst[1])
    }
    # Estimated values  ----------------------------------------------------------
    y.fitted <- predict(mod, list(matrix(x.true, ncol = 1)))
    y.fitted2 <- fitted(mod)
    lambda <- mod$lambda
    psi <- mod$psi
    w.hat <- mod$w.hat
    H <- mod$ipriorKernel$Hl[[1]]
    H2 <- H %*% H
    Vy <- vary(mod)
    VarY.inv <- solve(Vy)
  }

  # Prepare random draws from prior and posterior ------------------------------
  draw.pri <- draw.pos <- matrix(NA, ncol = no.of.draws, nrow = nrow(H.star))
  L <- chol(VarY.inv)
  for (i in 1:no.of.draws) {
    w.draw <- w.hat + crossprod(L, rnorm(length(y)))
    draw.pos[, i] <- mean(y) + lambda * as.numeric(H.star %*% w.draw)
    draw.pri[, i] <- mean(y) +
      lambda * as.numeric(H.star %*% rnorm(length(y), sd = sqrt(psi)))
  }
  dat.f <- rbind(data.frame(x = x.true, y = y.fitted, type = "Posterior"),
                 data.frame(x = x.true, y = mean(y), type = "Prior"))
  melted.pos <- melt(data.frame(f = draw.pos, x = x.true), id.vars = "x")
  melted.pri <- melt(data.frame(f = draw.pri, x = x.true), id.vars = "x")
  melted <- rbind(cbind(melted.pri, type = "Prior"),
                  cbind(melted.pos, type = "Posterior"))

  # Posterior predictive covariance matrix -------------------------------------
  varystar <- psi * (lambda ^ 2) * H.star %*% t(H.star) +
    diag(1 / psi, nrow(H.star))
  covystary <- psi * (lambda ^ 2) * H.star %*% H
  VarY.stary <- varystar - covystary %*% VarY.inv %*% t(covystary)
  dat.fit <- data.frame(x.true, y.fitted, sdev = sqrt(diag(VarY.stary)),
                        type = "95% credible interval")

  # Prepare random draws for posterior predictive checks -----------------------
  VarY.hat <- Vy - (psi ^ 2) * (lambda ^ 4) * H2 %*% VarY.inv %*% H2
  ppc <- matrix(NA, ncol = no.of.draws, nrow = nrow(H))
  L <- chol(VarY.hat)
  for (i in 1:no.of.draws) {
    ppc[, i] <- y.fitted2 + crossprod(L, rnorm(n))
    # ppc[, i] <- y.fitted2 + rnorm(n, sd = sqrt(diag(VarY.hat)))
  }
  melted.ppc <- melt(data.frame(x = x, ppc = ppc), id.vars = "x")
  melted.ppc <- cbind(melted.ppc, type = "Posterior predictive check")

  # Random draws from prior and posterior function -----------------------------
  p2.tmp <- ggplot() +
    geom_point(data = dat, aes(x = x, y = y), col = "grey55", alpha = 0.5) +
    scale_x_continuous(
      limits = c(min(x.true), max(x.true)),
      breaks = NULL, name = expression(italic(x))
    ) +
    scale_y_continuous(
      limits = c(min(y) - 5, max(y) + 5),
      breaks = NULL, name = expression(italic(y))
    ) +
    theme_bw()
  p2 <- p2.tmp +
    geom_line(data = melted, aes(x = x, y = value, group = variable),
              col = "steelblue3", size = 0.19, alpha = 0.5) +
    facet_grid(type ~ .) +
    geom_line(data = dat.f, aes(x = x, y = y), size = 1, linetype = 2, col = "grey10")
  p2.prior <- p2.tmp +
    geom_line(data = subset(melted, type == "Prior"),
              aes(x = x, y = value, group = variable),
              col = "steelblue3", size = 0.19, alpha = 0.5) +
    facet_grid(type ~ .)
  p2.prior.line <- p2.prior +
    geom_line(data = subset(dat.f, type == "Prior"), aes(x = x, y = y),
              size = 1, linetype = 2, col = "grey10")
  p2.posterior <- p2.tmp +
    geom_line(data = subset(melted, type == "Posterior"),
              aes(x = x, y = value, group = variable),
              col = "steelblue3", size = 0.19, alpha = 0.5) +
    facet_grid(type ~ .)
  p2.posterior.line <- p2.posterior +
    geom_line(data = subset(dat.f, type == "Posterior"), aes(x = x, y = y),
              size = 1, linetype = 2, col = "grey10")

  # Confidence band for predicted values  --------------------------------------
  p3 <- p1 +
    geom_line(data = dat.fit, aes(x = x.true, y = y.fitted), col = "grey50",
              size = 0.9, linetype = 2) +
    geom_ribbon(data = dat.fit, fill = "grey70", alpha = 0.5,
                aes(x = x.true, ymin = y.fitted - 1.96 * sdev,
                    ymax = y.fitted + 1.96 * sdev)) +
    facet_grid(type ~ .)

  p4 <- p2 +
    geom_line(data = dat.truth, aes(x = x.true, y = y.true, col = "Fitted"),
              size = 1, alpha = 0.75) + theme(legend.position = "none")

  # Posterior predictive checks ------------------------------------------------
  p5 <- ggplot() +
  scale_x_continuous(breaks = NULL, name = expression(italic(y))) +
  scale_y_continuous(breaks = NULL) +
  geom_line(data = melted.ppc,
            aes(x = value, group = variable, col = "yrep", size = "yrep"),
            stat = "density", alpha = 0.5) +
  geom_line(data = dat, aes(x = y, col = "y", size = "y"), stat = "density") +
    theme(legend.position = "bottom") +
  scale_colour_manual(
    name = NULL, labels = c("Observed", "Replications"),
    values = c("grey10", "steelblue3")
  ) +
  scale_size_manual(
    name = NULL, labels = c("Observed", "Replications"),
    values = c(1.1, 0.19)
  ) +
  facet_grid(type ~ .) +
  theme_bw() +
  theme(legend.position = c(0.9, 0.5))

  list(p2 = p2, p2.prior = p2.prior, p2.posterior = p2.posterior,
       p2.prior.line = p2.prior.line,
       p2.posterior.line = p2.posterior.line, p3 = p3, p4 = p4, p5 = p5)
}

## ---- canonical.kernel ----
plot.can <- plot1("Canonical")

## ---- fbm.kernel ----
plot.fbm <- plot1("FBM")

## ---- se.kernel.mle ----
plot.se <- plot1("SE")

## ---- save.plots.for.presentation ----
ggsave("../figure/points.pdf", p1, width = 6.5, height = 6.5 / 2.25)
# ggsave("figure/can-prior.pdf", plot.can$p2.prior.line,
#        width = 6.5, height = 6.5 / 1.5)
# ggsave("figure/can-posterior.pdf", plot.can$p2.posterior.line,
#        width = 6.5, height = 6.5 / 1.5)
ggsave("../figure/fbm-prior.pdf", plot.fbm$p2.prior.line,
       width = 6.5, height = 6.5 / 1.5)
ggsave("../figure/fbm-posterior.pdf", plot.fbm$p2.posterior.line,
       width = 6.5, height = 6.5 / 1.5)
ggsave("../figure/fbm-posterior-truth.pdf", {
  plot.fbm$p2.posterior +
    geom_line(data = dat.truth, aes(x = x.true, y = y.true), size = 1,
              alpha = 0.75, col = "red3") +
    annotate("text", label = "Truth", col = "red3", x = max(x.true),
             y = max(y.true) + 1)
  }, width = 6.5, height = 6.5 / 1.5)
p1 <- p1 +
  geom_line(data = dat.truth, aes(x = x.true, y = y.true), size = 1,
            alpha = 0.75, col = "red3") +
  scale_x_continuous(
    breaks = NULL, name = NULL
  ) +
  scale_y_continuous(
    breaks = NULL, name = NULL
  ) + theme_classic()
ggsave("../figure/plot-line.pdf", p1, width = 4, height = 4 / 2.25)
ggsave("../figure/credible-interval.pdf", plot.fbm$p3, width = 6.5, height = 6.5 / 1.5)
ggsave("../figure/ppc.pdf", plot.fbm$p5, width = 6.5, height = 6.5 / 1.5)
