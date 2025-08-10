# Fine-tuning scGPT for T-cell Subtype Prediction

This guide documents the likely workflow used to fine-tune a pretrained **scGPT** model on **TIL single-cell transcriptomes** to predict **T-cell subtypes**.  
It includes conceptual steps, pseudo-code, GPU training setup, and AWS deployment tips.

---

## 1. Data Preparation & Tokenization

```python
# inputs
counts_matrix = load_h5ad("til_scRNA.h5ad")        # cells x genes
labels = load_labels("tcell_subtypes.csv")         # e.g., CD4, CD8, Treg...

# tokenize: genes → tokens (scGPT uses gene vocab)
gene_vocab = load_vocab("scgpt_gene_vocab.json")
tokens, expr_bins = tokenize_cells(counts_matrix, gene_vocab, topk=2000, binning="log1p-quantile")

# split
train, val, test = split_dataset(tokens, expr_bins, labels, stratify=labels)
```

---

## 2. Load Pretrained scGPT & Replace Classifier Head

```python
device = "cuda" if torch.cuda.is_available() else "cpu"

model = load_pretrained_scgpt("scgpt_pbmc_base.pt")   # encoder weights
freeze_some_layers(model.encoder, ratio=0.0)          # optional: partial freeze

num_subtypes = len(unique(labels))
model.classifier_head = nn.Linear(model.hidden_dim, num_subtypes)  # new head
model.to(device)
```

---

## 3. Fine-Tuning Loop (GPU, Mixed Precision)

```python
optimizer = torch.optim.AdamW(model.parameters(), lr=2e-5, weight_decay=0.01)
scheduler = get_cosine_annealing(optimizer, T_max=10)
scaler = torch.cuda.amp.GradScaler()

accum_steps = 4
best_val = -float("inf")

for epoch in range(epochs):
    model.train()
    optimizer.zero_grad()
    for step, (tok, bins, y) in enumerate(train_loader):
        tok, bins, y = tok.to(device), bins.to(device), y.to(device)
        with torch.cuda.amp.autocast():
            logits = model(tok, expr_bins=bins)
            loss = F.cross_entropy(logits, y)

        scaler.scale(loss / accum_steps).backward()
        if (step + 1) % accum_steps == 0:
            scaler.step(optimizer); scaler.update()
            optimizer.zero_grad()
            scheduler.step()

    # validation
    model.eval()
    with torch.no_grad():
        val_logits, val_y = [], []
        for tok, bins, y in val_loader:
            logits = model(tok.to(device), expr_bins=bins.to(device))
            val_logits.append(logits.cpu()); val_y.append(y)
    val_acc = accuracy(cat(val_logits), cat(val_y))
    if val_acc > best_val:
        best_val = val_acc
        save_checkpoint("scgpt_tcell_head_best.pt", model, optimizer, epoch, val_acc)
```

---

## 4. Inference After Fine-Tuning

```python
model = load_checkpoint("scgpt_tcell_head_best.pt").to(device)
model.eval()

with torch.no_grad(), torch.cuda.amp.autocast():
    for tok, bins in test_loader:
        logits = model(tok.to(device), expr_bins=bins.to(device))
        preds = logits.argmax(dim=1)
        write_predictions(preds, "tcell_subtype_predictions.csv")
```

---

## 5. Artifact Management
- **Model registry**: push `scgpt_tcell_head_best.pt` + a `model_card.md`
- **Metrics**: accuracy/F1/confusion matrix
- **Data lineage**: S3 path to input h5ad & label version
- **Environment**: `environment.yml` (conda) or `requirements.txt`

---

## 6. AWS GPU Setup (EC2)

### Launch Instance
1. Choose GPU instance:
   - `g5.2xlarge` (A10G 24GB) – good balance
   - `p3.2xlarge` (V100 16GB) – older, cheaper
   - `p4d.24xlarge` (A100 40GB x8) – high end
2. AMI: **Deep Learning AMI (Ubuntu)** (pre-installed CUDA + PyTorch)
3. Storage: 200–500 GB gp3 EBS
4. Security group: allow SSH (port 22) from your IP
5. IAM role: S3 read/write if pulling data

### Connect & Verify
```bash
ssh -i ~/.ssh/your_key.pem ubuntu@ec2-XX-XX-XX-XX.compute-1.amazonaws.com
nvidia-smi
```

### Conda Env & Install
```bash
conda create -n scgpt python=3.10 -y
conda activate scgpt

pip install torch torchvision --index-url https://download.pytorch.org/whl/cu121
pip install scgpt anndata scanpy pandas numpy scikit-learn transformers
```

### Run
```bash
aws s3 sync s3://your-bucket/til_data ./til_data
python train_scgpt_tcell.py --config configs/tcell.yaml
```

---

## 7. AWS SageMaker (Optional)

```python
import sagemaker
from sagemaker.pytorch import PyTorch

sess = sagemaker.Session()
role = sagemaker.get_execution_role()

est = PyTorch(
    entry_point="train_scgpt_tcell.py",
    source_dir="src/",
    role=role,
    framework_version="2.2",
    py_version="py310",
    instance_type="ml.g5.2xlarge",
    instance_count=1,
    hyperparameters={"lr":2e-5, "epochs":10, "accum_steps":4},
    use_spot_instances=True,
    max_wait=36000,
)

est.fit({"train": "s3://your-bucket/til/train/", "val": "s3://your-bucket/til/val/"})
est.deploy(instance_type="ml.g5.xlarge", initial_instance_count=1)
```

---

## 8. Interview Blurb

> *"I started from a PBMC-pretrained scGPT encoder and swapped in a fresh classification head for our **T-cell subtype** labels. I fine-tuned on our TIL single-cell data with **mixed-precision** on GPU, logged metrics, and exported a model artifact that predicts **only T-cell subtypes**—by design, because the output layer was trained on that restricted label set."*
