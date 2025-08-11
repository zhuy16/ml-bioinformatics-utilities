# RAG & Hugging Face Utilities Overview

This document combines:
1. **How to Implement a Context-Aware RAG Pipeline**
2. **Hugging Face Utilities Overview**

---

## 1. How to Implement a Context-Aware RAG Pipeline

### Concept Overview
- **RAG** = Retrieve relevant chunks from your own data → feed them into an LLM → generate grounded answers.
- **Context-aware retrieval** = Use metadata filters, hybrid search, rerankers, and query rewriting to retrieve the *most relevant* data for a query.

---

### Typical Workflow

1. **Ingest & normalize**  
   - Sources: PDFs (OCR), HTML, CSV/Parquet, SQL tables, imaging metadata.  
   - Extract text + structured fields (e.g., patient cohort, assay type). Handle PHI appropriately.

2. **Chunk & enrich**  
   - Split docs into ~500–1,000 token chunks with small overlaps.
   - Attach metadata (e.g., project_id, tumor type, assay version).

3. **Embed & index**  
   - Compute embeddings (OpenAI, Hugging Face, etc.).
   - Store in a vector database: FAISS, pgvector (Postgres), OpenSearch, or Pinecone.

4. **Retrieve (context-aware)**  
   - **Hybrid search** = combine keyword (BM25) + semantic vector retrieval.
   - Apply **metadata filters** (e.g., `tumor=NSCLC`, `assay=WGS`).
   - Optional: rewrite query using conversation history.

5. **Rerank (optional)**  
   - Use a cross-encoder reranker to reorder top results by relevance.

6. **Generate with citations**  
   - Feed top chunks into the LLM prompt.
   - Return answer + source citations for auditability.

7. **Evaluate & trace**  
   - Log queries, retrieved docs, and outputs.
   - Maintain golden test sets for groundedness checks.

---

### Tech Stack Options

#### Open-source / self-hosted
- **Frameworks**: LangChain, LlamaIndex
- **Vector DB**: FAISS, pgvector (Postgres), OpenSearch
- **Hybrid Search**: BM25 (OpenSearch) + vectors
- **LLMs**: OpenAI API, Azure OpenAI, Anthropic, local models

#### Cloud-native
- **Azure AI Search** (hybrid search + vector) + **Azure OpenAI**
- Managed options: Pinecone, Weaviate, Milvus

---

### Minimal Pseudocode

```python
# 1) Ingest + chunk
docs = load_docs(["/reports/*.pdf", "/sql/assays", "/notes/*.txt"])  # OCR where needed
chunks = chunk(docs, size=800, overlap=100, metadata=extract_metadata(docs))

# 2) Embed + index
emb = EmbeddingModel()
vectors = [emb.encode(c.text) for c in chunks]
index = VectorIndex(store="pgvector", dim=emb.dim)   # or FAISS/OpenSearch
index.upsert(vectors, metadata=[c.meta for c in chunks])

# 3) Retrieve (context-aware)
query = history_aware_rewrite(user_query, chat_history)   # optional
hits = hybrid_search(query, index, bm25_backend="opensearch", filters={"tumor":"NSCLC"})

# 4) Rerank (optional)
reranked = cross_encoder_rerank(query, hits)[:k]

# 5) Generate with citations
answer = llm.generate(prompt_with_citations(query, reranked))
return answer
```

---

### Diagram: Context-Aware RAG Pipeline

```mermaid
flowchart LR
    A[Data Sources] -->|Ingest & OCR| B[Document Loader]
    B --> C[Chunking & Metadata Enrichment]
    C --> D[Embedding Model]
    D --> E[Vector DB]
    C --> F[Keyword Index (BM25)]
    E --> H[Retriever]
    F --> H
    H -->|Hybrid Search + Filters| I[Reranker (optional)]
    I --> J[LLM with Retrieved Context]
    J --> K[Answer + Citations]
    
    subgraph Data Sources
        A1[PDFs] 
        A2[SQL Tables] 
        A3[NoSQL Data]
        A4[Scanned Docs]
    end
```

---

## 2. Hugging Face Utilities Overview

Hugging Face is **not just a model registry** — it’s a full **ecosystem** for loading, fine-tuning, deploying, and sharing machine learning models.

---

### Core Components

#### A. **Model Hub**
- Repository of 500k+ open-source models (NLP, vision, audio, multimodal, biology).
- Pretrained weights, model cards, usage examples.
- Models are versioned and can be shared publicly or privately.

**Example:**
```python
from transformers import pipeline
classifier = pipeline("sentiment-analysis")
print(classifier("I love precision medicine!"))
```

---

#### B. **Framework & APIs**
- **`transformers`**: Unified API for inference and fine-tuning across PyTorch, TensorFlow, JAX.
- Handles tokenization, training loops, distributed training, mixed precision, and model export.
- Supports dozens of architectures: BERT, GPT, T5, ViT, Whisper, etc.

**Example:**
```python
from transformers import AutoTokenizer, AutoModelForSequenceClassification

tokenizer = AutoTokenizer.from_pretrained("distilbert-base-uncased")
model = AutoModelForSequenceClassification.from_pretrained("distilbert-base-uncased")

inputs = tokenizer("AI in healthcare is transformative", return_tensors="pt")
outputs = model(**inputs)
```

---

#### C. **Deployment & MLOps Tools**
- **Datasets**: Load and process thousands of public datasets with a common API.
- **Evaluate**: Common metrics (accuracy, BLEU, ROUGE, etc.).
- **Optimum**: Model optimization for speed (quantization, ONNX export, hardware acceleration).
- **Inference Endpoints**: Managed model hosting with auto-scaling APIs.
- **Spaces**: Deploy demos (Gradio/Streamlit) easily.
- **Model Cards**: Standardized documentation for ML models.
- **Hub integrations**: Works with MLflow, LangChain, Weights & Biases.

---

### Why It Matters for AI/GenAI Engineering
Fluency in Hugging Face typically means:
- Loading and running inference from the Model Hub.
- Fine-tuning models on domain-specific data.
- Correctly tokenizing inputs for the chosen architecture.
- Pushing trained models to the Hub for sharing/versioning.
- Integrating with production pipelines (self-hosted or managed).
- Using optimization and deployment tools for real-time applications.

---

### Useful Links
- [Transformers library](https://huggingface.co/docs/transformers/index)
- [Model Hub](https://huggingface.co/models)
- [Datasets library](https://huggingface.co/docs/datasets/index)
- [Evaluate library](https://huggingface.co/docs/evaluate/index)
- [Optimum](https://huggingface.co/docs/optimum/index)
- [Inference Endpoints](https://huggingface.co/inference-endpoints)
- [Spaces](https://huggingface.co/spaces)
