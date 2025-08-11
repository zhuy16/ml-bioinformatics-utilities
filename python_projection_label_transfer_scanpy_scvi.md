# Python Single-Cell Projection & Label Transfer
**Scanpy `ingest` (PCA/UMAP-based) vs. scVI/SCANVI + scArches (VAE-based)**

Goal: replicate the Seurat projection workflow in Python by (1) building a reference, (2) **projecting** a new query into the **same embedding**, and (3) **transferring labels**. Two viable paths:

- **Path A — Scanpy `ingest`**: fast PCA/UMAP projection + label transfer (no neural networks).  
- **Path B — scVI/SCANVI + scArches**: VAE latent space with robust batch handling, UMAP projection via `.transform`, and probabilistic label transfer.

---

## Prereqs & assumptions
- `adata_ref`: reference AnnData with curated labels in `adata_ref.obs['t_subtype']` (or your label key).  
- `adata_q`: query AnnData (new patient TILs).  
- Genes **must** be aligned (intersection, same order).  
- Keep preprocessing **consistent** across reference and query for the chosen path.

---

## Path A — Scanpy `ingest`: PCA/UMAP projection + label transfer
Closest to Seurat’s PCA/UMAP-based `MapQuery`: fit PCA/UMAP on reference, **project** query into the same spaces, and transfer labels via nearest neighbors.

```python
import scanpy as sc

# ---------- Reference ----------
# 1) normalize/log
sc.pp.normalize_total(adata_ref, target_sum=1e4)
sc.pp.log1p(adata_ref)

# 2) choose HVGs (reference defines the gene set)
sc.pp.highly_variable_genes(adata_ref, n_top_genes=3000, flavor="seurat_v3")
adata_ref = adata_ref[:, adata_ref.var['highly_variable']].copy()

# 3) scale, PCA, neighbors, UMAP
sc.pp.scale(adata_ref, max_value=10)
sc.tl.pca(adata_ref, n_comps=50)
sc.pp.neighbors(adata_ref, n_neighbors=30, n_pcs=50)
sc.tl.umap(adata_ref)

# ---------- Query ----------
# Subset to the reference HVGs and apply the same preprocessing
adata_q = adata_q[:, adata_ref.var_names].copy()
sc.pp.normalize_total(adata_q, target_sum=1e4)
sc.pp.log1p(adata_q)
sc.pp.scale(adata_q, max_value=10)

# ---------- Projection + label transfer ----------
# `obs='t_subtype'` tells ingest which label to transfer
sc.tl.ingest(adata_q, adata_ref, obs='t_subtype')

# Results:
# - adata_q.obsm['X_pca']   : projected PCs (into reference PCA)
# - adata_q.obsm['X_umap']  : query on the **reference UMAP**
# - adata_q.obs['t_subtype']: transferred labels (nearest-neighbor vote)
```

**When to use:** you want a **simple, transparent** projection equivalent to Seurat’s PCA/UMAP mapping without training a model.

---

## Path B — scVI/SCANVI + scArches: VAE latent + UMAP transform + probabilistic labels
Train a **reference SCVI** model on raw counts → obtain a stable latent; optionally train **SCANVI** for label transfer. Use **scArches** to map a new query into the **same latent**. Fit UMAP on the reference latent; `.transform()` query latent into the same 2D coordinates.

```python
import scanpy as sc
import scvi
import umap
import numpy as np

# --------------------
# A) Reference (SCVI)
# --------------------
# If counts are in a layer, use: scvi.model.SCVI.setup_anndata(adata_ref, layer="counts", batch_key="batch")
scvi.model.SCVI.setup_anndata(adata_ref, batch_key="batch")  # donor/experiment as batch
scvi_ref = scvi.model.SCVI(adata_ref, n_latent=30)
scvi_ref.train(max_epochs=400)

# Reference latent
adata_ref.obsm['X_scvi'] = scvi_ref.get_latent_representation()

# Fit UMAP ON THE LATENT and keep the fitted model
umap_ref = umap.UMAP(n_neighbors=30, min_dist=0.3, random_state=1)
adata_ref.obsm['X_umap_scvi'] = umap_ref.fit_transform(adata_ref.obsm['X_scvi'])

# --------------------
# B) Optional label transfer model (SCANVI)
# --------------------
# Ensure labels exist; set unlabeled to "Unknown" if needed.
# adata_ref.obs['t_subtype'] = adata_ref.obs['t_subtype'].cat.add_categories(["Unknown"]).fillna("Unknown")

scanvi_ref = scvi.model.SCANVI.from_scvi_model(
    scvi_ref, unlabeled_category="Unknown", labels_key="t_subtype"
)
scanvi_ref.train(max_epochs=20, n_samples_per_label=100)

# --------------------
# C) Map query into the same latent (scArches)
# --------------------
# 1) Align genes
adata_q = adata_q[:, adata_ref.var_names].copy()

# 2) Attach query to reference model (scArches); light fine-tuning
scvi_q = scvi.model.SCVI.load_query_data(adata_q, scvi_ref)
scvi_q.train(max_epochs=100)  # small

# 3) Get query latent & project to REFERENCE UMAP
adata_q.obsm['X_scvi'] = scvi_q.get_latent_representation(adata_q)
adata_q.obsm['X_umap_scvi'] = umap_ref.transform(adata_q.obsm['X_scvi'])

# --------------------
# D) Predict labels (SCANVI)
# --------------------
scanvi_q = scvi.model.SCANVI.load_query_data(adata_q, scanvi_ref)
scanvi_q.train(max_epochs=20)  # short tune
pred_labels = scanvi_q.predict(adata_q)            # categorical labels
pred_probs  = scanvi_q.predict(adata_q, soft=True) # per-class probabilities

adata_q.obs['t_subtype_pred'] = pred_labels
adata_q.obs['t_subtype_conf'] = np.max(pred_probs, axis=1)
```

**When to use:** you need **robust generalization** across batches/platforms, want a consistent **latent space** for many studies, and prefer **probabilistic label transfer** with confidence scores.

---

## Comparison (quick)
| Aspect | Scanpy `ingest` | scVI/SCANVI + scArches |
|---|---|---|
| Embedding transfer | Yes (PCA/UMAP) | Yes (SCVI latent; UMAP `.transform`) |
| Model training | No | Yes (VAE) |
| Batch handling | Limited (depends on preprocessing) | Strong (learned latent, batch correction) |
| Label transfer | kNN vote | SCANVI classifier (+ probabilities) |
| Speed/complexity | Fast & simple | Slower; more infra but more robust |

---

## Tips & gotchas
- **Gene alignment:** `adata_q = adata_q[:, adata_ref.var_names]` *before* projection/mapping.  
- **Consistency:** Use the **same** normalization and scaling choices in Path A; for Path B, use **raw counts** with `setup_anndata` (`layer='counts'` if needed).  
- **UMAP reproducibility:** Fit UMAP **once** on the reference (PCA or SCVI latent), then call `.transform` for queries.  
- **Confidence thresholds:** With SCANVI, filter or flag cells with low `t_subtype_conf`.  
- **Versioning:** Save fitted objects: `umap_ref`, SCVI/SCANVI state dicts, and reference gene list.

---

## Which to pick?
- Want a **Seurat-like projection** with minimal overhead? → **Scanpy `ingest`**.  
- Need a **durable reference** across many cohorts/platforms with **probabilistic labels**? → **scVI/SCANVI + scArches**.
