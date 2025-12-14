# -------------------------------
# Hybrid tree simulation: multiplicative Beta + LogNormal for non-overlapping nodes
# -------------------------------
library(nnls)
library(nleqslv)
library(ggplot2)
library(ggpubr)
library(latex2exp)
library(truncdist)

# -------------------------------
# Example tree
# -------------------------------
parent <- c(NA, 1, 1, 3, 2, 2, 6, 6, 8, 9)
t_lo <- c(10.0, 9.4, 4.8, 4.0, 8.0, 6.8, 1.8, 2.0, 1.7, 1.6)
t_hi <- c(10.5, 9.6, 5.0, 5.0, 9.5, 8.0, 2.5, 3.2, 2.5, 2.5)
N <- length(parent)
p_cov <- 0.90
n_sim <- 50000

# -------------------------------
# Helper: log-moment targets
# -------------------------------
compute_log_targets <- function(t_l, t_u, p) {
  z <- qnorm((1 + p)/2)
  sigma <- (log(t_u) - log(t_l)) / (2*z)
  mu <- log(t_l) + z*sigma
  list(mu = mu, sigma2 = sigma^2)
}

# Invert log-moments to Beta parameters
invert_logmoments_to_beta <- function(m, v){
  Ey <- exp(m + v/2)
  Vy <- exp(2*m + v)*(exp(v)-1)
  if(Vy >= Ey*(1-Ey)) Vy <- 0.99*Ey*(1-Ey)
  ab_sum <- Ey*(1-Ey)/Vy -1
  a0 <- max(1e-3, Ey*ab_sum)
  b0 <- max(1e-3, (1-Ey)*ab_sum)
  F <- function(par){
    a <- par[1]; b <- par[2]
    c(digamma(a)-digamma(a+b)-m,
      trigamma(a)-trigamma(a+b)-v)
  }
  sol <- nleqslv(c(a0,b0), F)
  list(alpha = sol$x[1], beta = sol$x[2])
}

safe_rlnorm_trunc <- function(meanlog, sdlog, upper_vec) {
  n <- length(upper_vec)
  out <- numeric(n)
  
  for (k in seq_len(n)) {
    u <- upper_vec[k]
    if (is.na(u) || u <= 0) {
      out[k] <- NA_real_
      next
    }
    
    # check if truncation range is numerically valid
    p_upper <- plnorm(u, meanlog, sdlog)
    if (p_upper < 1e-8) {
      # parent too small relative to distribution scale
      out[k] <- runif(1, 0.5, 0.9) * u
    } else {
      out[k] <- rtrunc(1, "lnorm", a = 0, b = u, meanlog = meanlog, sdlog = sdlog)
    }
  }
  out
}


# -------------------------------
# Step 1: classify nodes
# -------------------------------
is_beta_node <- rep(TRUE, N)
is_beta_node[1] <- FALSE  # root always LogNormal

for(i in 2:N){
  par_i <- parent[i]
  if(!is.na(par_i)){
    # If child interval does not overlap parent, treat as new root
    if(t_hi[i] <= t_lo[par_i]){
      is_beta_node[i] <- FALSE
    }
  }
}

# -------------------------------
# Step 2: compute log-targets
# -------------------------------
log_targets <- lapply(1:N, function(i) compute_log_targets(t_lo[i], t_hi[i], p_cov))
mu_vec <- sapply(log_targets, `[[`, "mu")
sigma2_vec <- sapply(log_targets, `[[`, "sigma2")

# -------------------------------
# Step 3: compute edges for beta nodes
# -------------------------------
edges <- which(is_beta_node & !is.na(parent))
n_edges <- length(edges)
edge_names <- paste(parent[edges], edges, sep="_")

A_mean <- matrix(0, nrow=n_edges, ncol=n_edges)
rownames(A_mean) <- as.character(edges)
colnames(A_mean) <- edge_names
b_mean <- numeric(n_edges)
b_var <- numeric(n_edges)

for(idx in seq_along(edges)){
  i <- edges[idx]
  par_i <- parent[i]
  m_edge <- mu_vec[i] - mu_vec[par_i]
  v_edge <- sigma2_vec[i]  # approx, relative variance
  b_mean[idx] <- m_edge
  b_var[idx] <- v_edge
  A_mean[idx, idx] <- 1
}

# -------------------------------
# Step 4: solve for beta parameters
# -------------------------------
m_hat <- solve(A_mean, b_mean)
v_hat <- nnls(A_mean, b_var)$x
v_hat <- pmax(v_hat, 1e-6)

beta_params <- lapply(1:n_edges, function(j) invert_logmoments_to_beta(m_hat[j], v_hat[j]))
edge_alpha <- sapply(beta_params, `[[`, "alpha")
edge_beta <- sapply(beta_params, `[[`, "beta")
names(edge_alpha) <- edge_names
names(edge_beta) <- edge_names

# -------------------------------
# Step 5: simulation
# -------------------------------
Wsim <- matrix(NA, nrow = n_sim, ncol = N)
Wsim[, 1] <- rlnorm(n_sim, mu_vec[1], sqrt(sigma2_vec[1]))  # root

for (i in 2:N) {
  par_i <- parent[i]
  if (is.na(par_i)) next
  
  if (is_beta_node[i]) {
    key <- paste(par_i, i, sep = "_")
    Wsim[, i] <- Wsim[, par_i] * rbeta(n_sim, edge_alpha[key], edge_beta[key])
  } else {
    Wsim[, i] <- safe_rlnorm_trunc(mu_vec[i], sqrt(sigma2_vec[i]), Wsim[, par_i])
  }
}


# -------------------------------
# Step 6: coverage
# -------------------------------
coverage <- sapply(1:N, function(i) mean(Wsim[,i] >= t_lo[i] & Wsim[,i] <= t_hi[i]))
print(round(coverage,3))

# -------------------------------
# Step 7: plots
# -------------------------------
plots <- list()
for(i in 1:N){
  df <- data.frame(value=Wsim[, i])
  cov_i <- round(coverage[i],3)
  plots[[i]] <- ggplot(df, aes(x=value)) +
    geom_density(fill="grey", alpha=0.5) +
    geom_vline(xintercept = t_lo[i], color="red", linetype="dashed") +
    geom_vline(xintercept = t_hi[i], color="red", linetype="dashed") +
    labs(title=paste0("Node ", i, " (Coverage=", cov_i, ")"),
         x = TeX(paste0("$T_{", i, "}$")),
         y="Density") +
    theme_minimal()  +
    theme(
      plot.title = element_text(size = 10, face = "bold")
    )
}

panel <- ggarrange(plotlist=plots, ncol=2, nrow=ceiling(N/2))
print(panel)

# -------------------------------
# Outputs
# -------------------------------
list(
  beta_nodes = edges,
  edge_alpha = edge_alpha,
  edge_beta = edge_beta,
  Wsim = Wsim,
  coverage = coverage
)
