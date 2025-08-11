# Seurat PCA Anchor-Based Projection (No Harmony)
**Goal:** Harmonize multiple scRNA-seq datasets with SCT integration, build a **T-cell reference** with subtype labels, and **project** a new TIL dataset onto the reference to obtain (i) predicted T-cell subtypes and (ii) the **UMAP projection in the reference embedding** — **without using Harmony** (keep PCA as the reference reduction).

---

## Inputs
- Multiple Seurat objects with raw counts: `seu1, seu2, seu3, ...` (donors/cohorts)
- One new query Seurat object with raw counts: `til.query` (new patient TILs)

> Replace placeholders with your own objects or loading code.

---

## 0) Setup
```r
set.seed(1)
library(Seurat)
library(dplyr)
# Optional parallelization:
# library(future)
# plan("multisession", workers = 4)
```

---

## 1) Light QC per dataset
Adjust thresholds to your biology/platform.
```r
qc_seurat <- function(obj, mt.pattern = "^MT-|^mt-",
                      min.features = 200, max.features = 8000, max.mt = 15) {
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = mt.pattern)
  subset(obj, subset = nFeature_RNA >= min.features &
                      nFeature_RNA <= max.features &
                      percent.mt <= max.mt)
}

seu.list <- list(seu1, seu2, seu3) |> lapply(qc_seurat)
```

---

## 2) SCTransform per dataset
Use the same regression variables across datasets for stability.
```r
seu.list <- lapply(seu.list, function(x) {
  SCTransform(x, vars.to.regress = "percent.mt", verbose = FALSE)
})
```

---

## 3) Select features & prep for SCT integration
```r
features <- SelectIntegrationFeatures(object.list = seu.list, nfeatures = 3000)
seu.list <- PrepSCTIntegration(object.list = seu.list, anchor.features = features)
```

---

## 4) Find anchors & integrate (SCT)
**No Harmony here** — we keep PCA as the reference reduction.
```r
anchors <- FindIntegrationAnchors(object.list = seu.list,
                                  normalization.method = "SCT",
                                  anchor.features = features)

ref.all <- IntegrateData(anchorset = anchors, normalization.method = "SCT")
```

---

## 5) Dimensional reduction & clustering on the integrated reference
Compute PCA/UMAP on the integrated object. Keep PCA as the reference reduction.
```r
ref.all <- RunPCA(ref.all, npcs = 50, verbose = FALSE)
ref.all <- RunUMAP(ref.all, dims = 1:50, verbose = FALSE)
ref.all <- FindNeighbors(ref.all, dims = 1:50)
ref.all <- FindClusters(ref.all, resolution = 0.5)
```

---

## 6) Build a T-cell reference (subsetting & labeling)
Curate subtype labels **now**, then we’ll use this as the reference for projection.

**A. Identify T cells** (use your own markers/annotation procedures and refine manually):
```r
DefaultAssay(ref.all) <- "SCT"
ref.all <- AddModuleScore(ref.all, features = list(c("CD3D","TRAC")), name = "TScore")
# Inspect distribution and pick a sensible threshold for your data
ref.t <- subset(ref.all, subset = TScore1 > 0)
```

**B. Create T-cell subtype labels** (replace with your curated labels):
```r
ref.t <- AddModuleScore(ref.t, features = list(c("IL7R","CCR7")), name = "NaiveLike")
ref.t <- AddModuleScore(ref.t, features = list(c("GZMB","NKG7","PRF1")), name = "CytoLike")

ref.t$t_subtype <- dplyr::case_when(
  ref.t$NaiveLike1 - ref.t$CytoLike1 > 0.1 ~ "T_NaiveLike",
  ref.t$CytoLike1 - ref.t$NaiveLike1 > 0.1 ~ "T_CytotoxicLike",
  TRUE ~ "T_Mixed"
)
Idents(ref.t) <- "t_subtype"

# Optionally rebuild PCA/UMAP on T cells alone for a tighter embedding
ref.t <- RunPCA(ref.t, npcs = 50, verbose = FALSE)
ref.t <- RunUMAP(ref.t, dims = 1:50, verbose = FALSE)

# saveRDS(ref.t, "reference_Tcells_integrated.rds")
```

---

## 7) Prepare the new TIL query dataset
QC + SCTransform with the same approach.
```r
til.query <- qc_seurat(til.query)
til.query <- SCTransform(til.query, vars.to.regress = "percent.mt", verbose = FALSE)
```

---

## 8) Find transfer anchors (reference = T-cell atlas, reduction = **PCA**)
**Key point:** Use PCA as `reference.reduction` so we can do **projection** into the fixed reference space.
```r
anchors.map <- FindTransferAnchors(
  reference = ref.t,
  query = til.query,
  normalization.method = "SCT",
  reference.reduction = "pca",
  dims = 1:50
)
```

---

## 9) Projection-based annotation **and** UMAP projection
`MapQuery()` transfers labels *and* projects the query into the **reference UMAP** (keeps reference fixed).
```r
til.mapped <- MapQuery(
  anchorset = anchors.map,
  query = til.query,
  refdata = list(t_subtype = ref.t$t_subtype),  # transfer subtype labels
  reference = ref.t,
  refdr = "pca",                 # project via the reference PCA
  reduction.model = "umap"       # and into the reference UMAP
)

# Results:
# - til.mapped$predicted.t_subtype : predicted T-cell subtype labels
# - Embeddings(til.mapped, "ref.umap") : projected UMAP coordinates
```

---

## 10) Quick visualization
```r
p.ref <- DimPlot(ref.t, reduction = "umap", group.by = "t_subtype", label = TRUE) +
  ggtitle("Reference T cells (UMAP)")

p.query <- DimPlot(til.mapped, reduction = "ref.umap", group.by = "predicted.t_subtype") +
  ggtitle("New TILs projected into reference UMAP")

# print(p.ref); print(p.query)
```

---

## 11) Optional: retrieve prediction scores explicitly
`MapQuery()` runs label transfer internally. If you also want the score matrix:
```r
pred <- TransferData(
  anchorset = anchors.map,
  refdata = ref.t$t_subtype,
  query = til.query,
  dims = 1:50
)
head(pred)  # per-class scores
til.mapped@meta.data <- cbind(til.mapped@meta.data, pred)
```

---

## Why **no Harmony** here?
- Harmony is excellent for **within-dataset** batch correction but does **not** yield a reusable mapping function for *new* data.  
- Using PCA as the reference reduction enables **anchor-based projection** (`MapQuery`) that keeps your reference fixed and places new cells into the same PCA/UMAP space.

---

## Tips / Pitfalls
- Keep `nfeatures`, `dims`, and QC thresholds consistent across runs.  
- Curate T-cell subtypes carefully (markers + manual review); the better the reference labels, the better the transfer.  
- If platform drift is high, increase `nfeatures`, check anchor counts, and consider rebalancing donors in the reference.  
- Save `ref.t` to reuse the exact same reference atlas across projects.

---

## Versions
- Tested pattern with Seurat v4/v5 API. Check `MapQuery()` docs for any version-specific args.  
- Set seeds to improve reproducibility for PCA/UMAP.

```

