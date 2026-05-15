# HiFiMAP
High resolution fast identity-by-descent (IBD) mapping
<br />
Current version: 2.0.0 (Optimized Streaming Pipeline)

## Overview
The new HiFiMAP pipeline has been completely rebuilt for massive scalability and memory efficiency. It uses a highly optimized Python preparser to convert IBD segments into chunked sparse matrices (`.mtx`) and differential updates (`.diff`). The main R package utilizes C++ (`Rcpp`) to stream these updates, allowing for fast, parallelized genome-wide association testing without frequently loading massive IBD matrices into RAM.

## Quick Installation 

`HiFiMAP` can be downloaded via:
```bash
git clone [https://github.com/baihongguo/HiFiMAP](https://github.com/baihongguo/HiFiMAP)
```

## Preparations
The software is developed and tested in Linux HPC environments. 

**Python Dependencies (Python >= 3.6):**
* `numpy`

**R Dependencies (R >= 4.0.0):**
* `data.table`
* `Matrix`
* `fastmatch`
* `GMMAT`
* `Rcpp`
* `CompQuadForm`

---

## Usage & Example Pipeline

Suppose you have an `example` folder containing the following files for Chromosome 21:
* `chr21_p_smoother_hap_ibd_res_3cm.ibd`: IBD segments generated from hapIBD.
* `chr21_reference.vcf.gz`: The corresponding VCF file used to generate the IBDs.
* `phenotype.txt`: Contains individual IDs, phenotype, and covariates (e.g., age, sex).
* `global_kinship.RData`: The global IBD (kinship) matrix.

### Step 0: Pre-parse the IBD Segments (Python)
This step transforms the raw IBD segments into high-performance sparse matrices and delta differences. It automatically chunks the chromosome to allow for parallel processing in Step 2.

```bash
python3 reformat_ibds_v1.3.py \
    --ibd chr21_p_smoother_hap_ibd_res_3cm.ibd \
    --vcf chr21_reference.vcf.gz \
    --output ./ibd_prep/chr21 \
    --n-checkpoints 20
```

| Argument | Description |
| ------------- | ------------- |
| `--ibd` | The IBD segment file output from hapIBD. |
| `--vcf` | The reference VCF file to extract exact SNP positions and sample IDs. |
| `--output` | The directory where the chunked matrices and diffs will be saved. |
| `--n-checkpoints` | Number of chunks to divide the chromosome into (e.g., `20`). This controls parallelization scaling in Step 2. Set this number to be less than the number of CPU you have for better performance, recommended is 20-30. |

*(This will generate a folder `./ibd_prep/chr21` containing `sites.txt`, `samples.txt`, and the chunked `.mtx` / `.diff` files).*

<br />

### Step 1: Fit the Null GLMM Model (R)
Next, fit the null model using your phenotypes and global kinship matrix. **This step only needs to be run once per phenotype**, regardless of how many chromosomes you analyze.

```bash
Rscript Step1_Preprocess_HiFiMAP.R \
    phenotype.txt \
    global_kinship.RData \
    toy_pheno \
    12345
```

| Argument | Description |
| ------------- | ------------- |
| `phenotype.txt` | The file containing ID, Phenotype, and covariates. |
| `global_kinship.RData` | The global kinship matrix saved in `.RData` format. |
| `toy_pheno` | Prefix for the output file. |
| `12345` | Random seed for reproducibility. |

*(This generates the null model object: `toy_pheno_glmmkin2randomvec.rds`, the subjects in the phenotype.txt can be a subset of the individuals included in your vcf files used for step0 and hapIBD).*

<br />

### Step 2: Run Parallelized HiFiMAP Scan (Bash wrapper)
To run the actual association scan efficiently across all chunks, create a bash wrapper script (e.g., `Run_HiFiMAP_Parallel.sh`) with the following contents:

```bash
#!/bin/bash

# ==========================================
# CONFIGURATION
# ==========================================
n_analysis_chunks=20   # MUST MATCH the --n-checkpoints used in Step 0
min_jobs=4             # Throttle: launch next chunks when active jobs < min_jobs
seed=12345
threads=4
target_chromosomes="21"

OUT_BASE="./results"
IBD_PREP_BASE="./ibd_prep"
GLMM_RDS="toy_pheno_glmmkin2randomvec.rds"
script_name="Step2_Run_HiFiMAP_Direct_Parallel_flexible.R"

mkdir -p "$OUT_BASE"

# ==========================================
# EXECUTION
# ==========================================
for chr in $target_chromosomes; do
    echo "Launching chunks for Chromosome $chr..."
    
    for chunk in $(seq 1 $n_analysis_chunks); do
        
        # Queue Management
        running_jobs=$(pgrep -f "$script_name" | wc -l)
        if [ "$running_jobs" -ge "$min_jobs" ]; then
            while [ "$running_jobs" -ge "$min_jobs" ]; do
                sleep 5
                running_jobs=$(pgrep -f "$script_name" | wc -l)
            done
        fi
        
        # Launch Chunk in background
        out_file="${OUT_BASE}/HiFiMAP_toy_pheno_chr${chr}_chunk${chunk}.txt"
        
        Rscript "$script_name" \
            "$chr" "$seed" "toy_pheno" "$threads" \
            "$GLMM_RDS" "${IBD_PREP_BASE}/chr${chr}" \
            "$chunk" "$n_analysis_chunks" > "${out_file}.log" 2>&1 &
            
    done
    
    # Wait for all chunks on this chromosome to finish
    wait 
    
    # Merge chunks
    echo "Merging results for Chromosome $chr..."
    final_out="${OUT_BASE}/HiFiMAP_toy_pheno_chr${chr}.txt"
    chunk1="${OUT_BASE}/HiFiMAP_toy_pheno_chr${chr}_chunk1.txt"
    
    if [ -f "$chunk1" ]; then
        head -n 1 "$chunk1" > "$final_out"
        for chunk in $(seq 1 $n_analysis_chunks); do
            tail -n +2 -q "${OUT_BASE}/HiFiMAP_toy_pheno_chr${chr}_chunk${chunk}.txt" >> "$final_out"
            rm "${OUT_BASE}/HiFiMAP_toy_pheno_chr${chr}_chunk${chunk}.txt"* # Cleanup
        done
    fi
done

echo "Pipeline Complete!"
```

Execute the wrapper script:
```bash
bash Run_HiFiMAP_Parallel.sh
```

### Final Output
The final merged result file (`results/HiFiMAP_toy_pheno_chr21.txt`) will contain the association statistics for every tested site:

```text 
chr     pos       n.ibd.segs    p.value
21      9411239   104523        0.83120
21      9411255   104526        0.82901
21      9411261   104550        0.65103
21      9411300   104811        0.00312
```
