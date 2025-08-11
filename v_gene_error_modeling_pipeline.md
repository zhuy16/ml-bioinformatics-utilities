# V-gene–wise Error Rate Modeling Pipeline

## Overview
This document outlines a method for modeling error rates in V-gene–specific assemblies by clustering sequences into allelic groups, identifying a centroid for each group, counting errors relative to the centroid, and using Bayesian inference with a Beta–Binomial model.

## Steps

### 1. Data Extraction
- Identify all assemblies annotated with a specific V-gene.
- Extract the FR1–FR3 region sequences from each assembly.

### 2. Clustering into Alleles
- Compute a pairwise sequence distance matrix using edit distance (Hamming if aligned; Levenshtein otherwise).
- Perform clustering using k-medoids (PAM) with k=2.
- If one cluster has fewer than 2 assemblies, discard it.
- The medoid of each remaining cluster serves as its centroid.

### 3. Error Counting
- For each allele group, compare each sequence to the centroid.
- Count mismatches (edit distance) as errors.
- Record:
  - `errors`: total mismatches across all sequences in the allele group.
  - `trials`: total bases compared (aligned length × number of sequences).

### 4. Output Table
Example schema:
```
v_gene, allele_id, n_sequences, centroid, errors, trials
TRBV5-1, A1, 23, ACTG..., 187, 23*L
TRBV5-1, A2, 7, ACTA..., 51, 7*L
```

### 5. Modeling: Beta–Binomial Hierarchical Model
- **Likelihood**:
  \[ y_{v,c} \sim \text{Binomial}(n_{v,c}, \theta_{v,c}) \]
- **Prior**:
  \[ \theta_{v,c} \sim \text{Beta}(\mu_v \kappa_v, (1 - \mu_v) \kappa_v) \]
- **Hyperpriors**:
  - \( \mu_v \sim \text{Beta}(a_\mu, b_\mu) \)
  - \( \kappa_v \sim \text{Gamma}(a_\kappa, b_\kappa) \)

### 6. Bayesian Inference
- Fit using Stan.
- Use MCMC to obtain posterior distributions for \(\mu_v, \kappa_v, \theta_{v,c}\).
- Perform posterior predictive checks.

### 7. Notes
- Discard small clusters to avoid noise.
- Edit distance calculation is sufficient for FR1–FR3; BLASTn is not needed.
- This method can be extended to allow >2 alleles if necessary.

---

**End goal:** A set of posterior distributions for V-gene–specific error rates, ready for interpretation and further modeling.
