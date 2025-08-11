# Nextflow DSL2 — Minimal HPC Skeleton (Cheat Sheet)

## Core ideas
- **Process** = one *task template* (shell script or containerized step).
  - Declares `input:`, `output:`, optional resources (`cpus`, `memory`, `time`), and a `script:` section (bash).
  - Each unique set of inputs → one job on your executor (local, SLURM, LSF, SGE, PBS, etc.).
- **Channel** = a *typed stream* of values that connects processes.
  - Create from files/values (e.g., `Channel.fromPath(...)`, `Channel.of(...)`), then transform with operators (`map`, `groupTuple`, `filter`, etc.).
  - Processes consume from input channels and emit to output channels.
- **Parameters (`params`)**
  - Defaults defined in `nextflow.config` (or in `main.nf`), override via CLI `--name value` or `-params-file`.
  - Precedence: **CLI** > **params-file** > **config** > **script defaults**.

---

## Minimal project layout
```
myflow/
├─ main.nf
└─ nextflow.config
```

---

## `main.nf` (self-contained minimal pipeline)
```groovy
// main.nf
// Minimal Nextflow DSL2 example: count lines in input files and publish results
nextflow.enable.dsl=2

// --------- defaults (optional; can also live in nextflow.config) ----------
params.reads  = params.reads  ?: 'data/*.txt'   // glob for demo text files
params.outdir = params.outdir ?: 'results'      // output folder

// ----------------------------- PROCESS ------------------------------------
process COUNT_LINES {
  tag "$sample_id"
  cpus 1
  memory '1 GB'
  time '30m'

  // publish outputs to a clean location
  publishDir "${params.outdir}/count_lines", mode: 'copy'

  input:
    tuple val(sample_id), path(file)

  output:
    tuple val(sample_id), path("${sample_id}.lines.txt"), emit: out

  script:
    """
    # count lines and save just the integer
    wc -l ${file} | awk '{print \$1}' > ${sample_id}.lines.txt
    """
}

// ----------------------------- WORKFLOW -----------------------------------
workflow {
  // 1) materialize files into tuples: (sample_id, path)
  reads_ch = Channel
              .fromPath(params.reads, checkIfExists: true)
              .map { f -> tuple(f.baseName, f) }

  // 2) run the process like a function
  COUNT_LINES(reads_ch)

  // 3) for demo: print (sample_id, line_count) to stdout
  COUNT_LINES.out.view { sid, result_file =>
    def n = result_file.text.trim()
    "sample=${sid}\tlines=${n}"
  }
}
```

---

## `nextflow.config` (local + HPC profiles, containers)
```groovy
// nextflow.config
nextflow.enable.dsl = 2

params {
  // defaults (can be overridden via CLI or -params-file)
  reads  = 'data/*.txt'
  outdir = 'results'
}

profiles {
  // run: nextflow run . -profile local
  local {
    process.executor = 'local'
    trace { enabled = true }
  }

  // run on SLURM cluster: nextflow run . -profile hpc --reads '/path/*.txt'
  hpc {
    process.executor = 'slurm'
    process.queue    = 'norm'       // change to your partition
    process.maxForks = 100

    // Apptainer/Singularity preferred on HPC (optional in this toy example)
    singularity.enabled    = true
    singularity.autoMounts = true
    singularity.cacheDir   = "$HOME/.singularity_cache"

    // keep heavy I/O off shared home when possible
    workDir = "/scratch/$USER/nextflow_work"
  }
}
```

---

## How to run
```bash
# local (uses defaults in nextflow.config)
nextflow run . -profile local

# override parameters
nextflow run . -profile local --reads "data/*.txt" --outdir "nf_results" -resume

# SLURM HPC example
nextflow run . -profile hpc --reads "/data/project/*.txt" --outdir "results" -resume
```

### Using a params file
`params.yaml`
```yaml
reads: "/data/project/*.txt"
outdir: "results"
```
Run:
```bash
nextflow run . -profile hpc -params-file params.yaml -resume
```

---

## Quick reference (interview-ready lines)
- **Process = task** (resources + script); **Channel = data stream**.
- **Immutability + caching**: completed tasks are reused with `-resume`.
- **Per-process resources** (`cpus/memory/time`) → efficient scheduling on HPC.
- **Profiles** in `nextflow.config` switch executors/envs instantly (`local`, `slurm`, etc.).
- **Containers**: prefer **Singularity/Apptainer** on HPC; one flag to enable.
- **Params**: define in config; override via `--name` or `-params-file` (YAML/JSON).

---

## Common channel patterns
```groovy
// single values
Channel.of('A','B','C')

// files by glob (error if none found)
Channel.fromPath('data/*.fastq.gz', checkIfExists: true)

// pair reads into tuples (sample_id, R1, R2)
Channel
  .fromFilePairs('data/*_{R1,R2}.fastq.gz', flat: true)
  .map { sid, pair -> tuple(sid, pair[0], pair[1]) }
```
