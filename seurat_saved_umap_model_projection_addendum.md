# Seurat Projection With Saved UMAP Model (No Harmony)

**Context:** For projection-based annotation you asked about saving the **UMAP model** from the reference, then projecting a new query onto that *same* embedding. In Seurat, this is done by running `RunUMAP(..., return.model = TRUE)` on the **reference**. The saved model is then used automatically by `MapQuery()` when you pass `reduction.model = "umap"`.

---

## Key points

- Use `return.model = TRUE` (not `save=TRUE`) in `RunUMAP` **on the reference** to store the fitted **uwot** model inside the UMAP reduction.
- The saved model lives in the UMAP reduction's `misc` slot (e.g., `ref.t[["umap"]]@misc`).
- `MapQuery()` will **detect and use** this saved model when you set `reduction.model = "umap"`, so your query gets **projected into the reference UMAP** (no refitting).

---

## Minimal code snippet (drop-in for your pipeline)

```r
# 1) Build the integrated reference as before (SCT anchors, no Harmony)
# ... ref.t <- IntegrateData(...) ; RunPCA(ref.t, npcs = 50)

# 2) Fit UMAP on the reference **and keep the model**
ref.t <- RunUMAP(ref.t, dims = 1:50, return.model = TRUE)

# (Optional) sanity check: confirm a uwot model is stored
# This prints the uwot model summary; structure may vary by Seurat version
ref_umap_model <- ref.t[["umap"]]@misc$model
print(ref_umap_model)

# 3) Find transfer anchors and project the query onto the reference UMAP
anchors.map <- FindTransferAnchors(
  reference = ref.t,
  query = til.query,
  normalization.method = "SCT",
  reference.reduction = "pca",
  dims = 1:50
)

til.mapped <- MapQuery(
  anchorset = anchors.map,
  query = til.query,
  refdata = list(t_subtype = ref.t$t_subtype),
  reference = ref.t,
  refdr = "pca",
  reduction.model = "umap"  # <-- uses the saved uwot model from ref.t
)

# Results:
# - til.mapped$predicted.t_subtype : transferred labels
# - Embeddings(til.mapped, "ref.umap") : query coordinates in the **reference UMAP**
```

> If `return.model = TRUE` was not used on the reference, `MapQuery()` cannot project onto the original UMAP; it will still transfer labels, but the query will not share the exact 2D embedding. Re-run `RunUMAP` on the reference with `return.model = TRUE` and try again.

---

## Where the model is stored

- After `RunUMAP(..., return.model = TRUE)`, Seurat places the fitted **uwot** model inside the UMAP reduction:
  - Example access: `ref.t[["umap"]]@misc$model`
- You generally **don’t need** to pull this object manually — `MapQuery()` reads it when `reduction.model = "umap"`.

---

## Version notes

- The `return.model` behavior is available in Seurat v3+ (uwot backend).  
- Exact internal slot names can vary slightly across versions, but accessing via `ref.t[["umap"]]@misc` is stable enough for a quick check.

```

