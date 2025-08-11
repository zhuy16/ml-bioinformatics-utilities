# V-gene–wise Error Rate Modeling Pipeline (Allele Clustering → Beta–Binomial in Stan)

## 1) Inputs & Definitions
- **Assemblies**: per-sample reconstructed sequences annotated with a **V-gene** (e.g., TRBV5-1).
- **Region of interest**: concatenated **FR1–FR3** segment within the V-gene region (constant length L once aligned; see below).
- **Goal**: for each **(V-gene, allele group)**, estimate an **error rate** θ via a Beta–Binomial model. We nest allele groups under each V-gene but keep the modeling V-gene–specific (allele groups share the same V-gene prior).

> Practical alignment note: If assemblies are already aligned to the same V-gene reference coordinates, FR1–FR3 substrings should have equal length (Hamming distance is then OK). If not, do a global alignment to the group centroid and count mismatches + gaps as errors.

---

## 2) Per–V-gene Preprocessing
For each V-gene v:

### 2.1 Extract FR1–FR3 sequences
1. Collect all assemblies annotated with V-gene v.
2. Extract FR1–FR3 subsequences (use consistent coordinate definitions; if varying, realign to a V-reference first).

### 2.2 Distance matrix
Compute pairwise **edit distances** D(i,j) between sequences:
- If all sequences are same length and aligned → **Hamming distance**.
- Otherwise → **Levenshtein** (edit distance) or global alignment–based mismatch+gap count.

### 2.3 Two-allele clustering (sequence-specific)
We want at most **two** clusters (two putative alleles). Use **k-medoids (PAM)** with k=2 on the distance matrix because:
- It handles arbitrary distance matrices (like edit distance).
- The **medoid** sequence (centroid by medoid) is an actual sequence—useful for reference.

**Procedure**
1. Run k-medoids (k=2) on D.
2. Compute cluster sizes n₁, n₂. If the smaller cluster has **n < 2**, **discard** that cluster and keep only the larger one (single-allele case for this V-gene instance).
3. For each retained cluster c, take its **medoid** m_c as the **centroid**.

> Alternative: run agglomerative clustering with average linkage on D and cut the tree at 2 groups; still discard any cluster with n<2. k-medoids gives you a natural centroid.

### 2.4 Error counting per allele group
For each allele group c with centroid (medoid) m_c and sequences {s_i}:
- If lengths are equal and aligned: **errors_i = Hamming(s_i, m_c)**
- Else: **errors_i = edit_distance(s_i, m_c)** or mismatch+gap count from a global alignment
- **trials_i = length(m_c)** (or aligned length considered for mismatches)

Aggregate per allele group:
- **y_c = Σ_i errors_i**  (total mismatches across sequences in group c)
- **n_c = Σ_i trials_i**  (total bases compared in group c)
- Store record: (v_gene=v, allele_id=c, errors=y_c, trials=n_c)

**Output table schema**
```
v_gene, allele_id, n_sequences, centroid, errors, trials
TRBV5-1, A1, 23, ACTG...,  187,  23*L
TRBV5-1, A2,  7,  ACTA...,   51,   7*L
TRBV6-2, A1, 12,  ...
```

---

## 3) Hierarchical Beta–Binomial Model (no logit)
We model **V-gene–specific** error-rate distributions; each allele group’s error rate θ_{v,c} is drawn from that V-gene’s Beta prior.

**Likelihood** (per allele group c under V-gene v):
- y_{v,c} ~ Binomial(n_{v,c}, θ_{v,c})

**Prior**: parameterize via mean–concentration
- θ_{v,c} ~ Beta( μ_v * κ_v, (1−μ_v) * κ_v )

**Hyperpriors** (weakly informative)
- μ_v ~ Beta( a_μ, b_μ )
- κ_v ~ Gamma( a_κ, b_κ )

---

## 4) Stan Model
```stan
data {
  int<lower=1> V;                  // number of V-genes
  int<lower=1> C;                  // number of allele groups total
  int<lower=1, upper=V> v_idx[C];  // map allele group -> V-gene index
  int<lower=0> y[C];               // errors per allele group
  int<lower=0> n[C];               // trials per allele group
  real<lower=0> a_mu;
  real<lower=0> b_mu;
  real<lower=0> a_kappa;
  real<lower=0> b_kappa;
}
parameters {
  vector<lower=0, upper=1>[V] mu;
  vector<lower=0>[V] kappa;
  vector<lower=0, upper=1>[C] theta;
}
model {
  mu ~ beta(a_mu, b_mu);
  kappa ~ gamma(a_kappa, b_kappa);
  for (c in 1:C) {
    int v = v_idx[c];
    theta[c] ~ beta(mu[v] * kappa[v], (1 - mu[v]) * kappa[v]);
  }
  y ~ binomial(n, theta);
}
```

---

## 5) Practical Notes
- **BLASTn vs edit distance?** Edit distance is enough for FR1–FR3; BLASTn is overkill unless you expect large indels or off-targets.
- **Centroid choice:** Use medoid to avoid synthetic bias.
- **Discarding tiny clusters:** Use n<2 rule to reduce noise.
- **Two alleles vs more:** Cap at two unless you specifically need more.

---

## 6) Outputs to Save
- Aggregated table: `v_gene, allele_id, n_sequences, centroid, errors, trials`
- Optional QC table: per-sequence errors vs centroid
- Stan posterior summaries for μ_v, κ_v, θ_{v,c}
