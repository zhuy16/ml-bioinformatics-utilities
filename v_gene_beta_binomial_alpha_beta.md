# V-Gene–Specific Beta–Binomial Model with Direct Alpha/Beta Parameterization (Numerically Stable)

> Goal: Model sequencing **error counts** `y_i` out of `n_i` trials for each V-gene group in the 260 bp FR1–3 region using a Beta–Binomial likelihood. Parameterize directly with \(\alpha_g\) and \(\beta_g\) for each V-gene group to avoid transformations (like logit) that can cause numerical instability. Fit via Stan with ~2,000 posterior draws.

---

## 1) Data layout
For each observation \(i = 1, \dots, N\):
- `y_i`: number of erroneous bases
- `n_i`: total bases assessed
- `g_i`: V-gene group index (1,…,G)

---

## 2) Model formulation (direct alpha/beta)
**Likelihood:**
\[
 y_i \sim \text{BetaBinomial}(n_i, \alpha_{g_i}, \beta_{g_i})
\]

**Priors:**
\[
 \alpha_g \sim \mathrm{Gamma}(a_\alpha, b_\alpha), \quad \beta_g \sim \mathrm{Gamma}(a_\beta, b_\beta)
\]
Gamma priors ensure positivity and avoid compound transformations.

---

## 3) Stan implementation
```stan
data {
  int<lower=1> N;            // observations
  int<lower=0> y[N];         // error counts
  int<lower=0> n[N];         // coverage
  int<lower=1> G;            // number of V-gene groups
  int<lower=1,upper=G> g[N]; // group index for each observation
  real<lower=0> a_alpha;
  real<lower=0> b_alpha;
  real<lower=0> a_beta;
  real<lower=0> b_beta;
}

parameters {
  vector<lower=0>[G] alpha_g;  // alpha per V-gene group
  vector<lower=0>[G] beta_g;   // beta per V-gene group
}

model {
  // Priors
  alpha_g ~ gamma(a_alpha, b_alpha);
  beta_g  ~ gamma(a_beta,  b_beta);

  // Likelihood
  for (i in 1:N)
    y[i] ~ beta_binomial(n[i], alpha_g[g[i]], beta_g[g[i]]);
}

generated quantities {
  vector[G] mu_g;
  vector[G] phi_g;
  vector[G] rho_g;
  for (k in 1:G) {
    mu_g[k]  = alpha_g[k] / (alpha_g[k] + beta_g[k]);
    phi_g[k] = alpha_g[k] + beta_g[k];
    rho_g[k] = 1 / (1 + phi_g[k]);
  }
}
```

---

## 4) Running the model in R (CmdStanR)
```r
library(cmdstanr)

stan_data <- list(
  N = length(y),
  y = y,
  n = n,
  G = max(g),
  g = g,
  a_alpha = 1,
  b_alpha = 1,
  a_beta  = 1,
  b_beta  = 1
)

mod <- cmdstan_model("v_gene_beta_binomial_alpha_beta.stan")
fit <- mod$sample(
  data = stan_data,
  seed = 1234,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 1000,
  iter_sampling = 500
)

print(fit$summary(c("mu_g", "phi_g", "rho_g")))
```

---

## 5) Interpretation
- **\(\alpha_g, \beta_g\)**: Beta distribution parameters for the error probability in V-gene \(g\).
- **\(\mu_g\)**: mean error rate \(\alpha_g / (\alpha_g + \beta_g)\)
- **\(\phi_g\)**: concentration parameter (\(\alpha_g + \beta_g\))
- **\(\rho_g\)**: overdispersion metric \(1 / (1 + \phi_g)\)

This approach avoids logit or mean–concentration reparameterizations, improving numerical stability and interpretability.

